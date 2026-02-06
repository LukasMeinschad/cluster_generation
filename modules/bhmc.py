import numpy as np  
from typing import List, Optional, Tuple, Dict, Union, Callable
from dataclasses import dataclass
from molecule_class import Molecule
from transformations import Transformation, GeometryOps, Quaternion
from psi4_interface import Psi4Calculator, Psi4Config
import random
from multiprocessing import Pool, cpu_count
from functools import partial
import matplotlib.pyplot as plt
@dataclass
class BHMCConfig:
    """   
    Configuration for the Basin Hopping Monte Carlo
    """
    temperature: float = 400.0 # Temperature in Kelvin
    max_steps: int = 10
    step_size: float = 0.2 # Step size for random perturbations

    # Psi4 settings 
    method: str = "hf"
    basis: str = "sto-3g"


    # Steps for different Phases
    phase_a_steps: int = 50 # Exploration Phase with Nonlocal Operators


class SimulationBox:
    """   
    Defines a Simulation box to constrain the molecular structures during the BHMC algorithm.
    Prevents unrealistic structures and ensure the algorithm explores a reasonable region of the PES.
    """
    def __init__(self,
                 box_type: str = "sphere",
                 radius: Optional[float] = None,
                 box_dimensions: Optional[np.ndarray] = None,
                 center: Optional[np.ndarray] = None):
        """    
        Initializes the Simulation Box

        Args:
            box_type: Type of box ("sphere" or "cube")
            radius: Radius for spherical box (required if box_type is "sphere")
            box_dimensions: Dimensions for cubic box (required if box_type is "cube")
            center: Center of the box (default: origin)
        """
        self.box_type = box_type
        self.radius = radius
        self.box_dimensions = box_dimensions
        self.center = center if center is not None else np.zeros(3)

        if box_type == "sphere" and radius is None:
            raise ValueError("Radius must be provided for spherical box.")
        if box_type == "cube" and box_dimensions is None:
            raise ValueError("Box dimensions must be provided for cubic box.")
        
    @staticmethod
    def from_cov_radii_molecule(molecule: Molecule, padding: float = 2.0) -> 'SimulationBox':
        """   
        Generates a spherical simulation box based on the covalent radii of the atoms

        R = 2 R_C * [1/2 + (3N/4π)^(1/3)] + padding

        here R_C is the average covalent radius of the M chemical species R_c = (1/M) sum_alpha R_c^alpha
        """
        average_radii_list = molecule.get_average_covalent_radii()
        average_radius = np.mean(average_radii_list)
        n_atoms = len(molecule.coordinates)
        radius = 2 * average_radius * (0.5 + ((3* n_atoms) / (4 * np.pi)) ** (1/3)) + padding
        return SimulationBox(box_type="sphere", radius=radius, center=np.zeros(3))

    def is_inside(self, coordinates: np.ndarray) -> bool:
        """   
        Checks if a given set of coordinates is inside the simulation box.
        """
        if self.box_type == "sphere":
            distances = np.linalg.norm(coordinates - self.center, axis=1)
            return np.all(distances <= self.radius)
        elif self.box_type == "cube":
            lower_bound = self.center - self.box_dimensions / 2
            upper_bound = self.center + self.box_dimensions / 2
            return np.all((coordinates >= lower_bound) & (coordinates <= upper_bound))
        else:
            raise ValueError(f"Invalid box type: {self.box_type}")


    # Testing function
    def plot_simulation_box(self):
        """   
        Plots the Simulation box and prints its parameters for verification.
        """
        if self.box_type == "sphere":
            print(f"Simulation Box: Sphere with radius {self.radius:.2f} Å centered at {self.center}")
            # Plotting the sphere
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            x = self.radius * np.outer(np.cos(u), np.sin(v)) + self.center[0]
            y = self.radius * np.outer(np.sin(u), np.sin(v)) + self.center[1]
            z = self.radius * np.outer(np.ones(np.size(u)), np.cos(v)) + self.center[2]

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(x, y, z, color='b', alpha=0.3)
            ax.set_title("Spherical Simulation Box")
            plt.show()
            plt.close()

    



class LocalOptimizer:
    """   
    Local geometry optimizer using Psi4
    """
    def __init__(self, method: str = "hf", basis: str = "sto-3g", verbose: bool = False):
        """   
        Initializes the local optimizer with specified quantum chemistry method and basis set.
        """
        self.method = method
        self.basis_set = basis

        # Setup Psi4 calculator
        config = Psi4Config(method=self.method, basis_set=self.basis_set)
        self.calculator = Psi4Calculator(config=config, verbose=verbose)

    def optimize(self, molecule: Molecule) -> Tuple[Optional[Molecule], Optional[float]]:
        """   
        Performs a local optimization on the given molecule.
        """
        mol_copy = molecule.copy()
        
        # Run optimization
        energy = self.calculator.geometry_optimization(mol_copy)
        
        if energy is None:
            return None, None

        return mol_copy, energy
    

class EnergyEvaluator:
    """  
    Energy Evaluator for BHMC using Psi4
    """
    def __init__(self, method: str = "hf", basis: str = "sto-3g"):
        """   
        Initializes the energy evaluator with specified quantum chemistry method and basis set.
        """
        self.method = method
        self.basis_set = basis

        # Setup Psi4 calculator
        config = Psi4Config(method=self.method, basis_set=self.basis_set)
        self.calculator = Psi4Calculator(config)

    def evaluate(self, molecule: Molecule) -> float:
        """   
        Evaluates the energy of the given molecule.
        """
        energy = self.calculator.single_point_energy(molecule)
        return energy


class LocalOperators:
    """   
    Defines a set of Local Operators for Perturbation of Molecular Geometries
    """
    def __init__(self):
        """   
        Initializes the local operators.
        """
        pass

    def random_displacement_submolecule(self, molecule: Molecule, submolecule_indices: List[List[int]], 
                                       delta_range: Tuple[float, float] = (0.3, 0.8)) -> Molecule:
        """   
        Randomly selects a submolecule and displaces it by a random vector.
        
        Args:
            molecule: Input molecule
            submolecule_indices: List of submolecule indices [[int], [int], ...]
            delta_range: Range for displacement magnitude in Angstrom
        """
        mol_copy = molecule.copy()
        n_atoms = len(mol_copy.coordinates)

        submol_idx = random.choice(submolecule_indices)
        # Generate random displacement magnitude
        delta_value = random.uniform(delta_range[0], delta_range[1])


        # Check if indices are valid
        if any(idx >= n_atoms or idx < 0 for idx in submol_idx):
            print(f"Invalid submolecule indices: {submol_idx} for molecule with {n_atoms} atoms.")

            raise ValueError("Submolecule indices are out of bounds.")
        
        # Generate random direction (unit vector)
        direction = np.random.uniform(-1, 1, size=3)
        direction = direction / np.linalg.norm(direction)
        
        # Apply same displacement to all atoms in submolecule
        displacement = direction * delta_value
        
        for atom_idx in submol_idx:
            mol_copy.coordinates[atom_idx] += displacement

        return mol_copy
    
    def random_atomic_displacement(self, molecule: Molecule, delta_range: float = 0.2) -> Molecule:
        """
        Apply random displacement to each atom independently.
        
        Args:
            molecule: Input molecule
            delta_range: Maximum displacement per atom in Angstrom
        """
        mol_copy = molecule.copy()
        
        n_atoms = len(mol_copy.coordinates)
        displacements = np.random.uniform(-delta_range, delta_range, size=(n_atoms, 3))
        
        mol_copy.coordinates += displacements
        
        return mol_copy


class NonLocalOperators:
    """   
    Defines a nonlocal set of operators for Basin Hopping and larger structural perturbations
    """
    def __init__(self):
        self.transformer = Transformation()

    def twist_operator(self, molecule: Molecule, submolecule_indices: List[List[int]],
                      rotation_angle_range: Tuple[float, float] = (0, 360)) -> Molecule:
        """  
        Twist Operator - rotates one submolecule around the center of mass of another.
        
        Args:
            molecule: Input molecule
            submolecule_indices: List of submolecule indices [[int], [int], ...]
            rotation_angle_range: Range of rotation angles in degrees
        """
        if len(submolecule_indices) < 2:
            raise ValueError("At least two submolecules are required for the twist operator.")

        mol_copy = molecule.copy()

        # Select two different submolecules
        ref_idx, rot_idx = random.sample(range(len(submolecule_indices)), 2)
        ref_submol = submolecule_indices[ref_idx]
        rot_submol = submolecule_indices[rot_idx]

        # Calculate the center of mass of reference submolecule
        ref_coords = mol_copy.coordinates[ref_submol]
        ref_masses = mol_copy.masses[ref_submol]
        ref_center = GeometryOps.center_of_mass(ref_coords, ref_masses)

        # Generate random rotation axis
        axis = self._random_unit_vector()
        
        # Generate random rotation angle
        angle_deg = np.random.uniform(*rotation_angle_range)
        angle_rad = np.radians(angle_deg)
        
        # Get rotation matrix
        R = GeometryOps.rotation_matrix_rodrigues(axis, angle_rad)
        
        # Rotate the rotating submolecule around reference center
        for atom_idx in rot_submol:
            coord = mol_copy.coordinates[atom_idx] - ref_center
            coord_rotated = coord @ R.T
            mol_copy.coordinates[atom_idx] = coord_rotated + ref_center
        
        return mol_copy

    def large_displacement(self,molecule,submolecule_indices: List[List[int]]) -> Molecule:
        """   
        Large Displacement Operator for Global Exploration of PES
        """
        mol_copy = molecule.copy()
        submol_idx = random.choice(submolecule_indices)

        # Large Displacement 1-3 Angstroms
        delta_value = random.uniform(1.0, 3.0)
        direction = self._random_unit_vector()
        displacement = direction * delta_value

        for atom_idx in submol_idx:
            mol_copy.coordinates[atom_idx] += displacement
        return mol_copy

    def _random_unit_vector(self) -> np.ndarray:
        """   
        Generates a random unit vector in 3D space.
        """
        phi = np.random.uniform(0, 2 * np.pi)
        costheta = np.random.uniform(-1, 1)
        theta = np.arccos(costheta) 

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        return np.array([x, y, z])
    
    def mirror_operator(self, molecule: Molecule, submolecule_indices: List[List[int]]) -> Molecule:
        """   
        Mirror Operator that mirrors a submolecule across a random plane.

        Algorithm:

        1. Select fixed submolecule
        2. Construct reference frame for fixed submolecule
        3. Select one of the coordinate planes (XY, YZ, XZ)
        4. Build Reflection Matrix H =  I - 2*n*n^T
        5. Apply reflection to the other submolecule
        """
        mol_copy = molecule.copy()
        if len(submolecule_indices) < 2:
            raise ValueError("At least two submolecules are required for the mirror operator.")
        # Select two different submolecules
        ref_idx, mir_idx = random.sample(range(len(submolecule_indices)), 2)
        ref_submol = submolecule_indices[ref_idx]
        mir_submol = submolecule_indices[mir_idx]
        # Construct reference frame for fixed submolecule
        ref_frame = self.transformer.set_reference_frame_from_indices(mol_copy, atom_indices=ref_submol)

        # Select random plane
        plane = random.choice(["xy", "yz", "xz"])
        if plane == "xy":
            # Select xy plane from reference frame
            n = ref_frame.z_axis
        elif plane == "yz":
            # Select yz plane from reference frame
            n = ref_frame.x_axis
        elif plane == "xz":
            # Select xz plane from reference frame
            n = ref_frame.y_axis
        else:
            raise ValueError(f"Invalid plane selection: {plane}")
        n = n / np.linalg.norm(n)  # Ensure normal vector is unit length
        # Build Reflection Matrix H = I - 2*n*n^T
        H = np.eye(3) - 2 * np.outer(n, n)
        # Apply reflection to the other submolecule
        for atom_idx in mir_submol:
            coord = mol_copy.coordinates[atom_idx] - ref_frame.origin
            coord_reflected = coord @ H.T
            mol_copy.coordinates[atom_idx] = coord_reflected + ref_frame.origin
        return mol_copy
    
    def roto_reflection_operator(self, molecule: Molecule, submolecule_indices: List[List[int]]) -> Molecule:
        """   
        Roto-Reflection Operator that combines rotation and reflection for more complex nonlocal transformations.
        """
        mol_copy = molecule.copy()
        if len(submolecule_indices) < 2:
            raise ValueError("At least two submolecules are required for the roto-reflection operator.")
        # Select two different submolecules
        ref_idx, rot_ref_idx = random.sample(range(len(submolecule_indices)), 2)
        ref_submol = submolecule_indices[ref_idx]
        rot_ref_submol = submolecule_indices[rot_ref_idx]
        # Construct reference frame for fixed submolecule
        ref_frame = self.transformer.set_reference_frame_from_indices(mol_copy, atom_indices=ref_submol)
        # Select random plane for reflection
        plane = random.choice(["xy", "yz", "xz"])
        if plane == "xy":
            n = ref_frame.z_axis
        elif plane == "yz":
            n = ref_frame.x_axis
        elif plane == "xz":
            n = ref_frame.y_axis
        else:
            raise ValueError(f"Invalid plane selection: {plane}")
        n = n / np.linalg.norm(n)
        # Build Reflection Matrix H = I - 2*n*n^T
        H = np.eye(3) - 2 * np.outer(n, n)
        # Generate random rotation angle
        angle_deg = np.random.uniform(0, 360)
        angle_rad = np.radians(angle_deg)
        # Get rotation matrix
        R = GeometryOps.rotation_matrix_rodrigues(n, angle_rad)
        # Combined Roto-Reflection Matrix
        M = R @ H
        # Apply combined transformation to the other submolecule
        for atom_idx in rot_ref_submol:
            coord = mol_copy.coordinates[atom_idx] - ref_frame.origin
            coord_transformed = coord @ M.T
            mol_copy.coordinates[atom_idx] = coord_transformed + ref_frame.origin
        return mol_copy

@dataclass
class OperatorConfig:
    """   
    Configuration of a Single Operator
    """
    name: str 
    operator_func: Callable
    operator_type: str  # "local" or "nonlocal"
    weight: float = 1.0 # Relative weight for selection probability

    # Stats
    n_applied: int = 0
    n_accepted: int = 0

    def acceptance_rate(self) -> float:
        """  
        Computes acceptance rate of the operator.
        """
        return self.n_accepted / self.n_applied if self.n_applied > 0 else 0.0
    
# ================================ Multiprocessing Worker for Independent BHMC Chains
def phase_a_worker(args: Tuple) -> List[Tuple]:
    """   
    Worker function that runs an independent BHMC chain for Phase A.
    Each worker generates n_structures by applying operators and evaluating energies.
    
    Args:
        args: Tuple containing (initial_molecule, submolecule_indices, n_structures, 
                                operator_configs, config_dict, worker_id)
    
    Returns:
        List of (Molecule, energy) tuples for accepted structures
    """
    initial_molecule, submolecule_indices, n_structures, operator_configs, config_dict, worker_id = args
    
    # Initialize operators for this worker
    nonlocal_ops = NonLocalOperators()
    k_B = 8.617333262145e-5  # Boltzmann constant in eV/K
    
    # Recreate operator functions
    operators = []
    weights = []
    for op_config in operator_configs:
        if op_config['name'] == 'twist':
            func = lambda m, s: nonlocal_ops.twist_operator(m, s)
        elif op_config['name'] == 'large_displacement':
            func = lambda m, s: nonlocal_ops.large_displacement(m, s)
        elif op_config['name'] == 'mirror':
            func = lambda m, s: nonlocal_ops.mirror_operator(m, s)
        elif op_config['name'] == 'roto_reflection':
            func = lambda m, s: nonlocal_ops.roto_reflection_operator(m, s)
        else:
            continue
        operators.append({'name': op_config['name'], 'func': func})
        weights.append(op_config['weight'])
    
    weights = np.array(weights)
    weights /= weights.sum()
    
    # Initialize energy evaluator
    evaluator = EnergyEvaluator(method=config_dict['method'], basis=config_dict['basis'])
    
    # Start BHMC chain
    current_structure = initial_molecule.copy()
    current_energy = evaluator.evaluate(current_structure)
    
    accepted_structures = []
    n_accepted = 0
    n_rejected = 0
    
    print(f"Worker {worker_id}: Starting with energy {current_energy:.6f} Hartree")
    
    for step in range(n_structures):
        # Select operator
        op_idx = np.random.choice(len(operators), p=weights)
        operator = operators[op_idx]
        
        try:
            # Apply operator
            new_structure = operator['func'](current_structure, submolecule_indices)
            
            # Evaluate energy
            new_energy = evaluator.evaluate(new_structure)
            
            # Metropolis acceptance
            if new_energy < current_energy:
                accept = True
            else:
                hartree_to_ev = 27.2114
                delta_e_ev = (new_energy - current_energy) * hartree_to_ev
                beta = 1.0 / (k_B * config_dict['temperature'])
                prob = np.exp(-beta * delta_e_ev)
                accept = random.random() < prob
            
            if accept:
                current_structure = new_structure
                current_energy = new_energy
                accepted_structures.append((current_structure.copy(), current_energy))
                n_accepted += 1
            else:
                n_rejected += 1
                
        except Exception as e:
            print(f"Worker {worker_id} step {step}: Error with operator {operator['name']}: {e}")
            n_rejected += 1
            continue
        
        # Progress update every 50 steps
        if (step + 1) % 50 == 0:
            print(f"Worker {worker_id}: {step + 1}/{n_structures} steps, "
                  f"Accepted: {n_accepted}, Rejected: {n_rejected}, "
                  f"Current E: {current_energy:.6f} Ha")
    
    print(f"Worker {worker_id}: Completed. Total accepted: {n_accepted}/{n_structures}")
    
    return accepted_structures


class MultiPhaseBHMC:
    """    
    Multi-Phase Basin Hopping Monte Carlo:

    Phase A (Exploration): Nonlocal Optimization, No Local Optimization Generation of Grid for Parallel Local Optimizations
    """
    def __init__(self, config: BHMCConfig):
        """   
        Initializes the Multi-Phase BHMC with the given configuration.
        """
        self.config = config
        self.optimizer = LocalOptimizer(method=config.method, basis=config.basis, verbose=False)
        self.energy_evaluator = EnergyEvaluator(method=config.method, basis=config.basis)
        self.local_ops = LocalOperators()
        self.nonlocal_ops = NonLocalOperators()

        self.k_B = 8.617333262145e-5  # Boltzmann constant in eV/K

        # Statistics
        self.n_accepted = 0
        self.n_rejected = 0
        self.energy_history = []
        self.best_energy = np.inf 
        self.best_structure = None
        self.current_structure = None
        self.current_energy = None
        self.trajectory = []

        # Result of Phase A: Collected Candidate Structures
        self.phase_a_structures = []

    def _initialize_operators(self) -> List[OperatorConfig]:
        """    
        Initialize the Operator Registry.
        """
        return [
            # Nonlocal Operators (A)
            OperatorConfig(
                name="twist",
                operator_func=lambda m,s: self.nonlocal_ops.twist_operator(m, s),
                operator_type="nonlocal",
                weight=1.5
            ),
            OperatorConfig(
                name="large_displacement",
                operator_func=lambda m,s: self.nonlocal_ops.large_displacement(m, s),
                operator_type="nonlocal",
                weight=1.0 
            ), 
            OperatorConfig(
                name="mirror",
                operator_func=lambda m,s: self.nonlocal_ops.mirror_operator(m, s),
                operator_type="nonlocal",
                weight=1.0
            ),
            OperatorConfig(
                name="roto_reflection",
                operator_func=lambda m,s: self.nonlocal_ops.roto_reflection_operator(m, s),
                operator_type="nonlocal",
                weight=0.5
            )
        ]


    def metropolis_accept(self, energy_new: float, energy_old: float) -> bool:
        """   
        Metropolis Acceptance Criterion
        """
        if energy_new < energy_old:
            return True 
        
        hartree_to_ev = 27.2114
        delta_e_ev = (energy_new - energy_old) * hartree_to_ev
        beta = 1.0 / (self.k_B * self.config.temperature)
        prob = np.exp(-beta * delta_e_ev)
        accept = random.random() < prob
        return accept
    

    def run_phase_a(self, 
                    initial_molecule: Molecule,
                    submolecule_indices: List[List[int]],
                    n_structures_per_worker: int = 300,
                    n_processes: int = 10):
        """   
        Runs Phase A with parallel independent BHMC chains.
        Each worker runs n_structures_per_worker steps independently.
        
        Args:
            initial_molecule: Starting molecular structure
            submolecule_indices: List of submolecule atom indices
            n_structures_per_worker: Number of BHMC steps per worker (default: 300)
            n_processes: Number of parallel workers (default: 10)
            
        Returns:
            List of all accepted (Molecule, energy) tuples from all workers
        """
        print(f"\n{'='*60}")
        print(f"Phase A: Global Exploration with Parallel BHMC Chains")
        print(f"{'='*60}")
        print(f"Number of workers: {n_processes}")
        print(f"Steps per worker: {n_structures_per_worker}")
        print(f"Total structures to generate: {n_processes * n_structures_per_worker}")
        print(f"Method: {self.config.method}/{self.config.basis}")
        print(f"Temperature: {self.config.temperature} K")
        print(f"{'='*60}\n")
        
        # Get operator configurations
        nonlocal_ops = [op for op in self._initialize_operators() if op.operator_type == "nonlocal"]
        operator_configs = [
            {'name': op.name, 'weight': op.weight}
            for op in nonlocal_ops
        ]
        
        # Configuration dictionary
        config_dict = {
            'method': self.config.method,
            'basis': self.config.basis,
            'temperature': self.config.temperature
        }
        
        # Prepare arguments for each worker
        args_list = [
            (
                initial_molecule,
                submolecule_indices,
                n_structures_per_worker,
                operator_configs,
                config_dict,
                worker_id
            )
            for worker_id in range(n_processes)
        ]
        
        # Run parallel workers
        print("Starting parallel BHMC chains...\n")
        all_accepted_structures = []
        
        with Pool(processes=n_processes) as pool:
            results = pool.map(phase_a_worker, args_list)
        
        # Collect all results
        for worker_id, worker_results in enumerate(results):
            print(f"Worker {worker_id}: Collected {len(worker_results)} accepted structures")
            all_accepted_structures.extend(worker_results)
        
        self.phase_a_structures = all_accepted_structures
        
        # Statistics
        total_generated = n_processes * n_structures_per_worker
        total_accepted = len(all_accepted_structures)
        
        print(f"\n{'='*60}")
        print(f"Phase A Complete")
        print(f"{'='*60}")
        print(f"Total structures generated: {total_generated}")
        print(f"Total structures accepted: {total_accepted}")
        print(f"Overall acceptance rate: {total_accepted/total_generated*100:.2f}%")
        
        if all_accepted_structures:
            energies = [e for _, e in all_accepted_structures]
            print(f"Energy range: [{min(energies):.6f}, {max(energies):.6f}] Hartree")
            print(f"Mean energy: {np.mean(energies):.6f} Hartree")
            print(f"Std energy: {np.std(energies):.6f} Hartree")
        
        print(f"{'='*60}\n")
        
        return self.phase_a_structures
