"""
Initialization module for BHMC cluster sampling.

This module handles:
1. Loading molecular structures from XYZ files
2. Fragmenting into submolecules based on connectivity
3. Optimizing individual submolecules
4. Creating simulation boxes
5. Generating initial random configurations
"""

import numpy as np
from typing import List, Tuple, Optional
from dataclasses import dataclass
import matplotlib.pyplot as plt


from molecule_class import Molecule
from box import SimulationBox
from psi4_interface import Psi4Calculator, Config
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# Import sobol sequence generator
from scipy.stats import qmc


# Logger import
from logger import Logger

@dataclass
class InitializationConfig:
    """Configuration for cluster initialization."""
    method: str = "hf"              # QM method for optimization
    basis: str = "sto-3g"           # Basis set for optimization
    box_type: str = "sphere"        # "sphere" or "cube"
    box_scale_factor: float = 2.5   # Scaling factor for box size
    min_distance: float = 2.0       # Minimum distance between submolecules (Angstrom)
    max_placement_attempts: int = 1000  # Max attempts to place molecules
    optimize_submolecules: bool = True  # Whether to optimize submolecules
    verbose: bool = True

@dataclass
class InitEnergyConfig:
    backend: str = "xtb" # TODO Psi4 implementation
    method: str = "hf"
    basis: str = "sto-3g"
    xtb_method: str = "GFN2-xTB"
    verbose: bool = False 

class InitEnergyEvaluator:
    """ 
    Energy evaluator for initialization pre-screening
    """
    EV_TO_HARTREE = 0.0367493  # 1 eV in Hartree

    def __init__(self, config: InitEnergyConfig):
        """ 
        Initialize the energy evaluator with the specified configuration
        """
        self.config = config
        self.backend = config.backend.lower()
        self._calculator = None
        if self.backend in {"xtb", "ase_emt"}:
            self._init_ase_backend()
    
    def _log(self, msg: str) -> None:
        if self.config.verbose:
            print(msg)
        
    def _init_ase_backend(self) -> None:
        try:
            from ase.calculators.emt import EMT
            if self.backend == "xtb":
                from xtb.ase.calculator import XTB
                self._calculator = XTB(method=self.config.xtb_method)
            else:
                self._calculator = EMT()
        except ImportError as e:
            raise RuntimeError(
                f"Failed to initialize backend {self.backend}: {e}. "
            )
    def evaluate_energy(self, molecule: Molecule) -> Optional[float]:
        """  
        Returns the energy 
        """
        try: 
            if self.backend in {"xtb", "ase_emt"}:
                from ase import Atoms
                symbols = molecule.get_symbols()
                if symbols is None:
                    return None
                atoms = Atoms(symbols=symbols, positions=molecule.coordinates)
                atoms.calc = self._calculator
                energy_ev = atoms.get_potential_energy()
                return float(energy_ev) * self.EV_TO_HARTREE
            
            self._log(f"Unknown backend for energy evaluation: {self.backend}")
            return None
        except Exception as e:
            self._log(f"Energy evaluation failed for backend {self.backend}: {e}")
            return None
        
    def rank_candidates(
            self,
            molecules: List[Molecule],
            n_keep: int
        ) -> List[Tuple[float, Molecule]]:
        """  
        Evaluate all candidates and return lowest n_keep as (energy, molecule)
        """
        scored: List[Tuple[Optional[float], Molecule]] = []
        for mol in molecules:
            e = self.evaluate_energy(mol)
            if e is not None:
                scored.append((e, mol))
            
        # Sort according to lowest energy and return top n_keep
        scored.sort(key=lambda x: x[0])
        return scored[:n_keep]

    


class ClusterInitializer:
    """Handles initialization of molecular clusters for BHMC sampling."""
    
    def __init__(self, config: InitializationConfig, logger: Optional[Logger] = None):
        """Initialize the cluster initializer.
        
        Args:
            config: Initialization configuration
        """
        self.config = config
        self.logger = logger 
        self.calculator = None
        if self.config.optimize_submolecules:
            qm_config = Config(method=config.method, basis=config.basis)
            self.calculator = Psi4Calculator(config=qm_config, verbose=config.verbose)

    def _log(self, msg: str, level: str = "info") -> None:
        """   
        Helper Method to log messages
        """
        if self.config.verbose:
            print(msg)
        if self.logger:
           getattr(self.logger, level)(msg) 

    def initialize_from_xyz(
        self,
        xyz_file: str,
        n_workers: int = 1,
        n_configurations: int = 10000,
        placing_method: str = "random",
        energy_backend: str= "xtb",
        energy_xtb_method: str = "GFN2-xTB",
    ) -> Tuple[Molecule, List[List[int]], SimulationBox]:
        """
        Complete initialization workflow from XYZ file.
        First fragments the molecule, then optimizes the submolecules using psi4.
        Then we create a simulation box and generate n_configurations random initial configurations
        Next we evaluate the energy for each configuration using the specified energy backend
        and return the lowest n_workers configurations as the initial molecules for the BHMC sampling         

        Args:
            xyz_file: Path to XYZ file containing molecular structure
            n_workers: Number of configurations to return (e.g, one per BHMC worker)
            n_configurations: Number of random configurations to generate for pre-screening
            placing_method: Placement method for initial coordinates ("random" or "sobol")
                energy_backend: Backend to use for energy evaluation ("xtb", "ase_emt", or "psi4")
                    energy_xtb_method: Method to use for xTB energy evaluation (e.g., "GFN2-xTB")

        Returns:
            Tuple of (initial_molecules, submolecule_indices, simulation_box)
        """
        if self.logger:
            self.logger.header("Cluster Initialization")
            msg = (
                f"Configuration:\n"
                f"  QM Method (Optimization): {self.config.method}\n"
                f"  Basis Set (Optimization): {self.config.basis}\n"
                f"  Simulation Box Type: {self.config.box_type}\n"
                f"  Simulation Box Scale Factor: {self.config.box_scale_factor}\n"
                f"  Minimum Distance Between Submolecules: {self.config.min_distance} Angstrom\n"
                f"  Candidate configurations generated: {n_configurations}\n"
                f"  Workers (n_keep) for energy ranking: {n_workers}\n"
                f"  Energy Evaluation Backend: {energy_backend}\n"
                f"  xTB Method (if applicable): {energy_xtb_method}\n"
            )
            self.logger.info(msg)

        
        # Step 1: Load molecule
        molecule = self._load_molecule(xyz_file)
    
        # Step 2: Fragment into submolecules
        submolecules = self._fragment_molecule(molecule)
    
        # Step 2b: Store submolecule indices BEFORE optimization
        submol_indices = [submol.get_index_in_parent() for submol in submolecules]
    
        # Step 3: Optimize submolecules (optional)
        if self.config.optimize_submolecules:
            submolecules = self._optimize_submolecules(submolecules)
    
        # Step 4: Create simulation box
        simulation_box = self._create_simulation_box(submolecules)
    
        # Step 5: Generate random initial configuration
        initial_molecules = self._generate_random_configurations(
            submolecules = submolecules,
            simulation_box = simulation_box,
            n_configurations = n_configurations,
            placing_method = placing_method,
        )

        # Step 5b: Energy prescreening and ranking
        energy_cfg = InitEnergyConfig(
            backend=energy_backend,
            method=self.config.method,
            basis=self.config.basis,
            xtb_method=energy_xtb_method,
            verbose=self.config.verbose
        )
        evaluator = InitEnergyEvaluator(config=energy_cfg)
        scored: List[Tuple[float, Molecule]] = []
        failed = 0
        for mol in initial_molecules:
            e = evaluator.evaluate_energy(mol)
            if e is not None:
                scored.append((e, mol))
            else:
                failed += 1
        self._log(f"\nEnergy evaluation complete: {len(scored)} successful, {failed} failed")
        scored.sort(key=lambda x: x[0])
        selected_molecules = [mol for _, mol in scored[:n_workers]]
        e_min = scored[0][0]
        e_max = scored[n_workers - 1][0]
        self._log(
            f"Pre-screening complete: valid = {len(scored)}, failed = {failed}, "
            f"Selected energies in range {e_min:.6f} to {e_max:.6f} Hartree ?? check unit :))"
        )

        summary_lines = [
            f"{'='*60}",
            f"Initialization Complete!",
            f"Total atoms: {len(initial_molecules[0].coordinates)}",
            f"Submolecules: {len(submolecules)}",
            f"Simulation box type: {simulation_box.box_type}",
            f"Candidates generated: {n_configurations}, Valid: {len(scored)}, Failed: {failed}",
            f"Sekected for workers: {n_workers}, Energy range: {e_min:.6f} to {e_max:.6f} Hartree",
        ]
        if simulation_box.box_type == "sphere":
            summary_lines.append(f"Box radius: {simulation_box.radius:.2f} Angstrom")
            summary_lines.append(f"Box volume: {simulation_box.get_volume():.2f} Angstrom^3")
        else:
            summary_lines.append(f"Box dimensions: {simulation_box.box_dimensions}")
            summary_lines.append(f"Box volume: {simulation_box.get_volume():.2f} Angstrom^3")
        
        summary = "\n".join(summary_lines)
        self._log(summary)
    
        return selected_molecules, submol_indices, simulation_box 
    
    def _generate_random_configurations(
            self,
            submolecules: List[Molecule],
            simulation_box: SimulationBox,
            n_configurations: int = 1,
            placing_method: str = "random",
        ) -> List[Molecule]:
        """   
        Helper function to generate multiple random initial configurations inside the simulation box

        Args:
            submolecules: List of submolecules to place
            simulation_box: Simulation box to place molecules in
            n_configurations: Number of random configurations to generate
            placing_method: Placement method for initial coordinates ("random" or "sobol")
        """

        configurations = []
        for config_idx in range(n_configurations):
            mol = self._generate_configuration(
                submolecules,
                simulation_box,
                placing_method=placing_method,
            )
            configurations.append(mol)
        return configurations 



    def _load_molecule(self, xyz_file: str) -> Molecule:
        """Load molecule from XYZ file."""
        self._log(f"\n[1/5] Loading molecule from {xyz_file}") 

        
        with open(xyz_file, 'r') as f:
            xyz_content = f.read()
        
        molecule = Molecule.from_xyz(xyz_content)

        self._log(f" Loaded: {len(molecule.coordinates)} atoms")
        
        return molecule
    
    def _fragment_molecule(self, molecule: Molecule) -> List[Molecule]:
        """Fragment molecule into connected components."""
        self._log(f"\n[2/5] Fragmenting molecule into submolecules")
        

        molecule.compute_bonds()
        submolecules = molecule.fragment_by_connectivity()
        
        self._log(f" Found {len(submolecules)} submolecules:")
        for i, submol in enumerate(submolecules):
            self._log(f"  Submolecule {i+1}: {len(submol.coordinates)} atoms")
        

        
        return submolecules
    
    def _optimize_submolecules(self, submolecules: List[Molecule]) -> List[Molecule]:
        """Optimize each submolecule individually (parallelized)."""
        self._log(f"\n[3/5] Optimizing Submolecules (parallelized)")

    
        n_workers = min(len(submolecules), multiprocessing.cpu_count())
        self._log(f"Workers: {n_workers}")

        optimized = [None] * len(submolecules)
    
        # Extract serializable data from submolecules
        submol_data = [
            {
                'atom_labels': submol.atom_labels,
                'coordinates': submol.coordinates,
                'index': i
            }
            for i, submol in enumerate(submolecules)
        ]
    
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            # Submit jobs with serializable data
            future_to_idx = {
                executor.submit(
                    self._optimize_single_static,
                    data['atom_labels'],
                    data['coordinates'],
                    self.config.method,
                    self.config.basis
                ): data['index']
                for data in submol_data
            }
            
            # Collect results
            for future in as_completed(future_to_idx):
                idx = future_to_idx[future]
                try:
                    result_mol = future.result()
                    optimized[idx] = result_mol
                    self._log(f"  Submolecule {idx+1} optimized successfully")
                except Exception as e:
                    self._log(f"  Submolecule {idx+1} optimization failed: {e}", level="error")
                    # Fall back to original submolecule
                    optimized[idx] = submolecules[idx]
    
        return optimized

    @staticmethod
    def _optimize_single_static(
        atom_labels: List[str],
        coordinates: np.ndarray,
        method: str,
        basis: str
    ) -> Molecule:
        """
        Optimize a single submolecule (static method for parallel execution).
        
        Args:
            atom_labels: List of atom symbols
            coordinates: Atomic coordinates
            method: QM method
            basis: Basis set
        
        Returns:
            Optimized molecule or original if optimization fails
        """
        from psi4_interface import Psi4Calculator, Config
        from molecule_class import Molecule
        
        # Create molecule from data
        submol = Molecule.from_labels_and_coords(
            atom_labels=atom_labels,
            coordinates=coordinates
        )
        
        # Create calculator inside worker process
        qm_config = Config(method=method, basis=basis)
        calc = Psi4Calculator(config=qm_config, verbose=False)
        
        try:
            result = calc.optimize(submol)
            return result.molecule if result.success else submol
        except Exception:
            return submol
        
    
    
    def _create_simulation_box(self, submolecules: List[Molecule]) -> SimulationBox:
        """Create simulation box based on submolecules."""
        self._log(f"\n[4/5] Creating simulation box")
        
        # Collect all covalent radii
        all_cov_radii = []
        total_atoms = 0
        
        for submol in submolecules:
            all_cov_radii.extend(submol.covalent_radii)
            total_atoms += len(submol.coordinates)
        
        # Create box
        simulation_box = SimulationBox.from_covalent_radii(
            covalent_radii=all_cov_radii,
            n_atoms=total_atoms,
            box_type=self.config.box_type,
            scale_factor=self.config.box_scale_factor
        )
        self._log(f" Box type:  {simulation_box.box_type}")
        if simulation_box.box_type == "sphere":
            self._log(f" Radius: {simulation_box.radius:.2f} Angstrom")
            self._log(f" Volume: {simulation_box.get_volume():.2f} Angstrom^3")
        else:
            self._log(f" Dimensions: {simulation_box.box_dimensions}")
            self._log(f" Volume: {simulation_box.get_volume():.2f} Angstrom^3")
        
        return simulation_box
    
    def _generate_configuration(
        self, 
        submolecules: List[Molecule], 
        simulation_box: SimulationBox,
        placing_method: str = "random",
    ) -> Molecule:
        """
        Generate random initial configuration inside simulation box.
        
        Args:
            submolecules: List of submolecules to place
            simulation_box: Simulation box to place molecules in
            placing_method: Method to place molecules:
                - "random" for random placement with distance checks
                - "grid" for systematic grid placement (TODO: Implement)
                - "sobol" for low-discrepancy sequence placement (TODO: Implement)

        """
        
        if len(submolecules) == 0:
            raise ValueError("No submolecules to place")
        

        method = placing_method.lower().strip()
        
        if method not in {"random", "grid", "sobol"}:
            raise ValueError(f"Invalid placing method: {placing_method}. "
                             "Choose from 'random', 'grid', or 'sobol'.")

        sobol_engine = None
        if method == "sobol":
            sobol_engine = qmc.Sobol(d=3, scramble=True)
            sobol_engine_rotation = qmc.Sobol(d=3, scramble=True)
        

        placed_coords = np.empty((0, 3))
        placed_atoms = []

        for i, submol in enumerate(submolecules):
            if self.config.verbose:
                print(f"  Placing submolecule {i+1}...", end=" ", flush=True)
            
            # Get submolecule COM
            submol_com = np.average(submol.coordinates, axis=0, weights=submol.masses)
            centered_coords = submol.coordinates - submol_com
            
            placed = False
            for attempt in range(self.config.max_placement_attempts):
                
                if method == "random":
                    position = simulation_box.get_random_position(n_points=1)
                elif method == "sobol":
                    u = sobol_engine.random(n=1)[0]  # Get one point

                    if simulation_box.box_type == "cube":
                        half_dims = simulation_box.box_dimensions / 2.0
                        position = simulation_box.center + (2.0 * u - 1.0) * half_dims

                    elif simulation_box.box_type == "sphere":
                        # Uniform in volume:
                        # r = R * u1^(1/3), cos(theta) = 2*u2 - 1, phi = 2*pi*u3
                        r = simulation_box.radius * (u[0] ** (1.0/3.0))
                        cos_theta = 2.0 * u[1] - 1.0
                        sin_theta = np.sqrt(max(0.0, 1.0 - cos_theta * cos_theta)) 
                        phi = 2.0 * np.pi * u[2]

                        direction = np.array([
                            sin_theta * np.cos(phi),
                            sin_theta * np.sin(phi),
                            cos_theta
                        ])  
                        position = simulation_box.center + r * direction
                    else:
                        raise ValueError(f"Unsupported box type for Sobol placement: {simulation_box.box_type}")
                
                # Sample Random Rotation
                if method == "sobol":
                    rotation_matrix = self._sobol_rotation_matrix(sobol_engine_rotation)
                else:
                    rotation_matrix = self._random_rotation_matrix()
            
                rotated_coords = (rotation_matrix @ centered_coords.T).T
                
                # Translate to new position
                new_coords = rotated_coords + position
                
                # Check distance constraints
                if i == 0 or self._check_min_distance(new_coords, placed_coords):
                    placed_coords = np.vstack([placed_coords, new_coords])                     
                    placed_atoms.extend(submol.atom_labels.tolist())
                    placed = True
                    break
            
            if not placed:
                raise RuntimeError(
                    f"Could not place submolecule {i+1} after "
                    f"{self.config.max_placement_attempts} attempts. "
                    "Try increasing box size or reducing min_distance."
                )
        
        initial_molecule = Molecule.from_labels_and_coords(
            atom_labels=placed_atoms,
            coordinates=placed_coords.copy(),  
        )
        
        return initial_molecule
    
    def _check_min_distance(
        self, 
        new_coords: np.ndarray, 
        existing_coords: np.ndarray 
    ) -> bool:
        """Check if new coordinates satisfy minimum distance constraint."""
        if len(existing_coords) == 0:
            return True
        
        # existing_coords is already a numpy array, no conversion needed
        diff = new_coords[:, np.newaxis, :] - existing_coords[np.newaxis, :, :]
        distances = np.linalg.norm(diff, axis=2)
        return np.all(distances >= self.config.min_distance)



    @staticmethod
    def _random_rotation_matrix() -> np.ndarray:
        """Generate a random rotation matrix using quaternions (vecotized)"""
        u = np.random.rand(3)
        q = np.array([
            np.sqrt(1 - u[0]) * np.sin(2 * np.pi * u[1]),
            np.sqrt(1 - u[0]) * np.cos(2 * np.pi * u[1]),
            np.sqrt(u[0]) * np.sin(2 * np.pi * u[2]),
            np.sqrt(u[0]) * np.cos(2 * np.pi * u[2])
        ])
        q0, q1, q2, q3 = q 
        return np.array([
            [1 - 2*(q2**2 + q3**2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
            [2*(q1*q2 + q0*q3), 1 - 2*(q1**2 + q3**2), 2*(q2*q3 - q0*q1)],
            [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1**2 + q2**2)]
        ])
    
    def _sobol_rotation_matrix(self, sobol_engine) -> np.ndarray:
        """ 
        Generate a quasi-random rotation matrix usign the sobol sequencse
        """
        u = sobol_engine.random(n=1)[0]

        # Shoemake's algorithm like above
        q = np.array([
            np.sqrt(1 - u[0]) * np.sin(2 * np.pi * u[1]),
            np.sqrt(1 - u[0]) * np.cos(2 * np.pi * u[1]),
            np.sqrt(u[0]) * np.sin(2 * np.pi * u[2]),
            np.sqrt(u[0]) * np.cos(2 * np.pi * u[2])
        ])
        q0, q1, q2, q3 = q
        return np.array([
            [1 - 2*(q2**2 + q3**2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
            [2*(q1*q2 + q0*q3), 1 - 2*(q1**2 + q3**2), 2*(q2*q3 - q0*q1)],
            [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1**2 + q2**2)]
        ])

        

    
# =============================== Debugging and Testing ===============================

def test_initializer(xyz_file: str,
                     method: str = "hf",
                     basis: str = "cc-pvdz",
                     optimize: bool = True,
                     box_type: str = "sphere",
                     box_scale: float = 2.5,
                     min_distance: float = 2.0,
                     placing_method: str = "random",
                     n_configurations: int = 1,
                     save_output: bool = True):
    """   
    Test the ClusterInitializer with a given XYZ file

    This function performs a complete initialization workflow and
    outputs some diagnostics
    
    Args:
        xyz_file: Path to XYZ file
        method: QM method for optimization
        basis: Basis set for optimization
        optimize: Whether to optimize submolecules
        box_type: Type of simulation box ("sphere" or "cube")
        box_scale: Scaling factor for box size
        min_distance: Minimum distance between submolecules (Angstrom)
        placing_method: Placement method ("random" or "sobol")
        save_output: Whether to save output files (initial configuration and representatives)
    """
    print(f"\n{'='*80}")
    print(f"Testing ClusterInitializer with {xyz_file}")
    print(f"{'='*80}\n")
    print(f"Method: {method}, Basis: {basis}, Optimize: {optimize}")
    print(f"Box type: {box_type}, Box scale: {box_scale}, Min distance: {min_distance} Angstrom")
    print(f"Placement method: {placing_method}")

    config = InitializationConfig(
        method=method,
        basis=basis,
        optimize_submolecules=optimize,
        box_type=box_type,
        box_scale_factor=box_scale,
        min_distance=min_distance,
        verbose=True
    )

    initializer = ClusterInitializer(config=config)
    try:
        initial_molecules, submol_indices, simulation_box = initializer.initialize_from_xyz(
            xyz_file,
            n_configurations=n_configurations,
            placing_method=placing_method,
        )
        
        # Use first configuration for diagnostics
        initial_molecule = initial_molecules[0]
        
        print(f"\nInitialization successful!")
        print(f"Generated {len(initial_molecules)} initial configuration(s)")
        print(f"Molecule Information:")
        print(f"  Total atoms: {len(initial_molecule.coordinates)}")
        print(f" Total mass: {np.sum(initial_molecule.masses):.2f} amu")
        print(f" Center of mass: {np.average(initial_molecule.coordinates, axis=0, weights=initial_molecule.masses)}")
        print(f"\n Submolecule Information:")
        print(f"  Number of submolecules: {len(submol_indices)}")
        for i, indices in enumerate(submol_indices):
            print(f"    Submolecule {i+1}: {len(indices)} atoms")
        print(f"\n Simulation Box Information:")
        print(f"  Box type: {simulation_box.box_type}")
        if simulation_box.box_type == "sphere":
            print(f"  Box radius: {simulation_box.radius:.2f} Angstrom")
            print(f"  Box volume: {simulation_box.get_volume():.2f} Angstrom^3")
        else:
            print(f"  Box dimensions: {simulation_box.box_dimensions}")
            print(f"  Box volume: {simulation_box.get_volume():.2f} Angstrom^3")
        all_inside = simulation_box.is_inside(initial_molecule.coordinates)
        print(f"  All atoms inside box: {'Yes' if all_inside else 'No'}")

        # Calculate pairwise distances 
        print(f"\n Distance Check:")
        min_dist = np.inf 
        max_dist = 0.0 
        coords = initial_molecule.coordinates
        for i in range(len(coords)):
            for j in range(i+1, len(coords)):
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < min_dist:
                    min_dist = dist
                if dist > max_dist:
                    max_dist = dist
        print(f"  Minimum interatomic distance: {min_dist:.2f} Angstrom")
        print(f"  Maximum interatomic distance: {max_dist:.2f} Angstrom")

        # Calculate center of mass distance between submolecules
        print(f"\n Submolecule Distance Check:")
        if len(submol_indices) > 1:
            coms = []
            for indices in submol_indices:
                sub_coords = initial_molecule.coordinates[indices]
                sub_masses = initial_molecule.masses[indices]
                com = np.average(sub_coords, axis=0, weights=sub_masses)
                coms.append(com)
            for i in range(len(coms)):
                for j in range(i+1, len(coms)):
                    dist = np.linalg.norm(coms[i] - coms[j])
                    print(f"  Distance between submolecule {i+1} and {j+1}: {dist:.2f} Angstrom")
        else:
            print(f"  Only one submolecule, skipping distance check.")
        # Save output files
        if save_output:
            # Write all initial configuration into 1 xyz file
            with open("initial_configurations.xyz", "w") as f:
                for idx, mol in enumerate(initial_molecules):
                    f.write(f"{len(mol.coordinates)}\n")
                    f.write(f"Configuration {idx+1}\n")
                    for label, coord in zip(mol.atom_labels, mol.coordinates):
                        f.write(f"{label} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")
            print(f"\nInitial configurations saved to initial_configurations.xyz") 
        
    except Exception as e:
        print(f"\nInitialization failed with error: {e}")

if __name__ == "__main__":
    """  
    Run test when module is executed directly
    """
    xyz_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/nh3_2.xyz"
    #xyz_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/h2o_n2.xyz"
    #xyz_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/h2o_2.xyz"
    
    H2O_2_XYZ = """
    6
    Coordinates from ORCA-job h2o_2 E -152.102897751726
      O           0.20131422818946     -0.13419863189991     -0.38118207664628
      H           1.10826926774619      0.03054527004232     -0.16898516572928
      H          -0.26386315635905     -0.06924940339370      0.43390806815981
      O           3.06513444516272      0.52224610599392      0.05902763481169
      H           3.25817767296947      1.29105665797100     -0.44996761253291
      H           3.66908154229119     -0.14039999871364     -0.23001884806303    
    """



    method = "hf"
    basis = "cc-pvdz"
    optimize = True
    box_type = "sphere"
    box_scale = 2.5
    min_distance = 2.0
    save_output = True
    n_configurations = 28

    # Numerical test sample 2^10 points in [0,1] x [0,1] using Sobol sequence
    sampler = qmc.Sobol(d=2, scramble=True)
    sample = sampler.random(n=2**10)
    # Sample 1000 random points in [0,1] x [0,1] using uniform distribution
    uniform_sample = np.random.rand(2**10, 2)
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.scatter(sample[:, 0], sample[:, 1], s=5, color='blue')
    plt.title('Sobol Sequence Sample (1000 points)')
    plt.subplot(1, 2, 2)
    plt.scatter(uniform_sample[:, 0], uniform_sample[:, 1], s=5, color='red')
    plt.title('Uniform Random Sample (1000 points)')
    plt.tight_layout()
    plt.savefig("figures/sobol_vs_uniform.png") 

    # Try the Energy Evaluator
    energy_config = InitEnergyConfig(
        backend="xtb",
        method="GFN2-xTB",
        verbose=True
    )
    energy_evaluator = InitEnergyEvaluator(config=energy_config)
    # Create a simple molecule for testing
    molecule = Molecule.from_xyz(H2O_2_XYZ)
    energy = energy_evaluator.evaluate_energy(molecule) 
    print(f"\nEnergy evaluation test:")
    print(f"  Evaluated energy: {energy:.6f} Hartree")

#    test_initializer(
#        xyz_file=xyz_file,
#        method=method,
#        basis=basis,
#        optimize=optimize,
#        box_type=box_type,
#        box_scale=box_scale,
#        min_distance=min_distance,
#        save_output=save_output,
#        n_configurations=n_configurations
#    )


