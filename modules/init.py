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

from molecule_class import Molecule
from box import SimulationBox
from psi4_interface import Psi4Calculator, Config
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

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


class ClusterInitializer:
    """Handles initialization of molecular clusters for BHMC sampling."""
    
    def __init__(self, config: InitializationConfig):
        """Initialize the cluster initializer.
        
        Args:
            config: Initialization configuration
        """
        self.config = config
        self.calculator = None
        if self.config.optimize_submolecules:
            qm_config = Config(method=config.method, basis=config.basis)
            self.calculator = Psi4Calculator(config=qm_config, verbose=config.verbose)
    
    def initialize_from_xyz(self, xyz_file: str) -> Tuple[Molecule, List[List[int]], SimulationBox]:
        """Complete initialization workflow from XYZ file.
        
        Args:
            xyz_file: Path to XYZ file containing molecular structure
            
        Returns:
            Tuple of (initial_molecule, submolecule_indices, simulation_box)
        """
        if self.config.verbose:
            print(f"\n{'='*70}")
            print(f"Cluster Initialization")
            print(f"{'='*70}")
    
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
        initial_molecule = self._generate_random_configuration(
            submolecules, simulation_box
        )
    
        # Note: submol_indices already captured at step 2b
    
        if self.config.verbose:
            print(f"{'='*70}")
            print(f"Initialization Complete!")
            print(f"  Total atoms: {len(initial_molecule.coordinates)}")
            print(f"  Submolecules: {len(submolecules)}")
            print(f"  Simulation box: {simulation_box.box_type}")
            if simulation_box.box_type == "sphere":
                print(f"  Box radius: {simulation_box.radius:.2f} Angstrom")
            else:
                print(f"  Box dimensions: {simulation_box.box_dimensions}")
            print(f"{'='*70}\n")
    
        return initial_molecule, submol_indices, simulation_box
    
    def _load_molecule(self, xyz_file: str) -> Molecule:
        """Load molecule from XYZ file."""
        if self.config.verbose:
            print(f"\n[1/5] Loading molecule from {xyz_file}")
        
        with open(xyz_file, 'r') as f:
            xyz_content = f.read()
        
        molecule = Molecule.from_xyz(xyz_content)
        
        if self.config.verbose:
            print(f"  Loaded: {len(molecule.coordinates)} atoms")
        
        return molecule
    
    def _fragment_molecule(self, molecule: Molecule) -> List[Molecule]:
        """Fragment molecule into connected components."""
        if self.config.verbose:
            print(f"\n[2/5] Fragmenting molecule by connectivity")
        
        molecule.compute_bonds()
        submolecules = molecule.fragment_by_connectivity()
        
        if self.config.verbose:
            print(f"  Found {len(submolecules)} submolecules:")
            for i, submol in enumerate(submolecules):
                atoms = len(submol.coordinates)
                print(f"    Submolecule {i+1}: {atoms} atoms") 
        
        return submolecules
    
    def _optimize_submolecules(self, submolecules: List[Molecule]) -> List[Molecule]:
        """Optimize each submolecule individually (parallelized)."""
        if self.config.verbose:
            print(f"\n[3/5] Optimizing individual submolecules (parallel)")
    
        n_workers = min(len(submolecules), multiprocessing.cpu_count())
    
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
                    
                    if self.config.verbose:
                        print(f"  Submolecule {idx+1} complete ✓")
                except Exception as e:
                    if self.config.verbose:
                        print(f"  Submolecule {idx+1} failed: {e}")
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
        if self.config.verbose:
            print(f"\n[4/5] Creating simulation box")
        
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
        
        if self.config.verbose:
            print(f"  Box type: {self.config.box_type}")
            if simulation_box.box_type == "sphere":
                print(f"  Radius: {simulation_box.radius:.2f} Angstrom")
                print(f"  Volume: {simulation_box.get_volume():.2f} Angstrom^3")
            else:
                print(f"  Dimensions: {simulation_box.box_dimensions}")
                print(f"  Volume: {simulation_box.get_volume():.2f} Angstrom^3")
        
        return simulation_box
    
    def _generate_random_configuration(
        self, 
        submolecules: List[Molecule], 
        simulation_box: SimulationBox
    ) -> Molecule:
        """Generate random initial configuration inside simulation box."""
        if self.config.verbose:
            print(f"\n[5/5] Generating random initial configuration")
        
        if len(submolecules) == 0:
            raise ValueError("No submolecules to place")
        
        # Use numpy arrays from the start instead of lists
        placed_coords = np.empty((0, 3))  # <-- CHANGED: Start with empty 2D array
        placed_atoms = []
        
        for i, submol in enumerate(submolecules):
            if self.config.verbose:
                print(f"  Placing submolecule {i+1}...", end=" ", flush=True)
            
            # Get submolecule COM
            submol_com = np.average(submol.coordinates, axis=0, weights=submol.masses)
            centered_coords = submol.coordinates - submol_com
            
            # Try to place submolecule
            placed = False
            for attempt in range(self.config.max_placement_attempts):
                # Random position
                position = simulation_box.get_random_position(n_points=1)
                # Random rotation
                rotation_matrix = self._random_rotation_matrix()
                rotated_coords = (rotation_matrix @ centered_coords.T).T
                new_coords = rotated_coords + position
                
                # Check distance constraints
                if i == 0 or self._check_min_distance(new_coords, placed_coords):
                    placed_coords = np.vstack([placed_coords, new_coords])  # <-- CHANGED: vstack instead of extend
                    placed_atoms.extend(submol.atom_labels.tolist())
                    placed = True
                    if self.config.verbose:
                        print(f"✓ (attempt {attempt+1})")
                    break
            
            if not placed:
                if self.config.verbose:
                    print(f"✗ Failed after {self.config.max_placement_attempts} attempts")
                raise RuntimeError(
                    f"Could not place submolecule {i+1} after "
                    f"{self.config.max_placement_attempts} attempts. "
                    "Try increasing box size or reducing min_distance."
                )
        
        # Create new molecule with explicit copy of coordinates
        initial_molecule = Molecule.from_labels_and_coords(
            atom_labels=placed_atoms,
            coordinates=placed_coords.copy(),  # <-- CHANGED: Explicit copy
        )
        
        return initial_molecule
    
    def _check_min_distance(
        self, 
        new_coords: np.ndarray, 
        existing_coords: np.ndarray  # <-- CHANGED: Type hint to numpy array
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

    
# =============================== Debugging and Testing ===============================

def test_initializer(xyz_file: str,
                     method: str = "hf",
                     basis: str = "cc-pvdz",
                     optimize: bool = True,
                     box_type: str = "sphere",
                     box_scale: float = 2.5,
                     min_distance: float = 2.0,
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
        save_output: Whether to save output files (initial configuration and representatives)
    """
    print(f"\n{'='*80}")
    print(f"Testing ClusterInitializer with {xyz_file}")
    print(f"{'='*80}\n")
    print(f"Method: {method}, Basis: {basis}, Optimize: {optimize}")
    print(f"Box type: {box_type}, Box scale: {box_scale}, Min distance: {min_distance} Angstrom")

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
        initial_molecule, submol_indices, simulation_box = initializer.initialize_from_xyz(xyz_file)
        
        print(f"\nInitialization successful!")
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
            # Write xyz file for initial configuration
            initial_xyz = initial_molecule.to_xyz()
            with open("initial_configuration.xyz", "w") as f:
                f.write(initial_xyz)
            print(f"\nInitial configuration saved to initial_configuration.xyz")
        
    except Exception as e:
        print(f"\nInitialization failed with error: {e}")

if __name__ == "__main__":
    """  
    Run test when module is executed directly
    """
    xyz_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/nh3_2.xyz"
    #xyz_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/h2o_n2.xyz"
    #xyz_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/h2o_2.xyz"
    
    method = "hf"
    basis = "cc-pvdz"
    optimize = True
    box_type = "sphere"
    box_scale = 2.5
    min_distance = 2.0
    save_output = True

    test_initializer(
        xyz_file=xyz_file,
        method=method,
        basis=basis,
        optimize=optimize,
        box_type=box_type,
        box_scale=box_scale,
        min_distance=min_distance,
        save_output=save_output
    )


