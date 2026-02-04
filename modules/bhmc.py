import numpy as np  
from typing import List, Optional, Tuple, Dict, Union, Callable
from dataclasses import dataclass
from molecule_class import Molecule
from transformations import Transformation, GeometryOps, Quaternion
from psi4_interface import Psi4Calculator, Psi4Config
import random

@dataclass
class BHMCConfig:
    """   
    Configuration for the Basin Hopping Monte Carlo
    """
    temperature: float = 300.0 # Temperature in Kelvin
    max_steps: int = 10
    step_size: float = 0.2 # Step size for random perturbations

    # Psi4 settings 
    method: str = "hf"
    basis: str = "sto-3g"


    # Steps for different Phases
    phase_a_steps: int = 200 # Exploration Phase with Nonlocal Operators


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
    

    def run_phase_a(self, initial_molecule: Molecule, submolecule_indices: List[List[int]]):
        """   
        Runs Phase A the Global Exploration Phase using Nonlocal Operators without Local Optimization 
        """
        nonlocal_ops = [op for op in self._initialize_operators() if op.operator_type == "nonlocal"]
        weights = np.array([op.weight for op in nonlocal_ops])
        weights /= weights.sum()  # Normalize to probabilities
        current_sturcture = initial_molecule.copy()
        current_energy = self.energy_evaluator.evaluate(current_sturcture)
        
        for step in range(self.config.phase_a_steps):
            print(f"Phase A - Step {step + 1}/{self.config.phase_a_steps}")

            # Select operator based on weights
            operator = np.random.choice(nonlocal_ops, p=weights)
            print(f"  Applying Operator: {operator.name}")
            try:
                new_structure = operator.operator_func(current_sturcture, submolecule_indices)
            except Exception as e:
                print(f"  Operator {operator.name} failed with error: {e}")
                self.n_rejected += 1
                continue
            new_energy = self.energy_evaluator.evaluate(new_structure)
            print(f"  New Energy: {new_energy:.6f} Hartree, Current Energy: {current_energy:.6f} Hartree")
            
            operator.n_applied += 1

            accept = self.metropolis_accept(new_energy, current_energy)
            if accept:
                print("  Move accepted.")
                current_sturcture = new_structure
                current_energy = new_energy
                operator.n_accepted += 1
                self.phase_a_structures.append((current_sturcture.copy(), current_energy))
            else:
                print("  Move rejected.")
                self.n_rejected += 1
            
        return self.phase_a_structures




class BHMC:
    """"   
    Basin Hopping Monte Carlo for Global Optimization of Molecular Geometries
    """
    def __init__(self, config: BHMCConfig):
        """   
        Initializes the BHMC with the given configuration.
        """
        self.config = config
        self.optimizer = LocalOptimizer(method=config.method, basis=config.basis, verbose=False)
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

        # Initialize Psi4 Calculator for energy evaluations
        psi4_config = Psi4Config(method=config.method, basis_set=config.basis)
        self.calculator = Psi4Calculator(psi4_config, verbose=False)

    def metropolis_accept(self, energy_new: float, energy_old: float) -> bool:
        """   
        Metropolis acceptance criterion

        Args:
            energy_new: Energy of the new structure (Hartree)
            energy_old: Energy of the old structure (Hartree)

        Returns:
            True if accepted, False otherwise
        """
        if energy_new < energy_old:
            return True
        
        hartree_to_ev = 27.2114
        delta_e_ev = (energy_new - energy_old) * hartree_to_ev

        beta = 1.0 / (self.k_B * self.config.temperature)
        prob = np.exp(-beta * delta_e_ev)
        
        accept = random.random() < prob
        
        print(f"  Metropolis: ΔE = {delta_e_ev:.4f} eV, prob = {prob:.4f}, accept = {accept}")
        
        return accept
    
    def run(self, initial_molecule: Molecule, submolecule_indices: List[List[int]], 
            operator_type: str = "local") -> Dict:
        """   
        Runs the BHMC algorithm starting from the initial molecule.

        Following Algorithm:

        1. Local Optimization of the Initial Structure
        2. Iterative BHMC Steps
            a. Apply Random Perturbation
            b. Local Optimization
            c. Metropolis Acceptance Criterion
        """
        self.current_structure, self.current_energy = self.optimizer.optimize(initial_molecule)
        print(self.current_structure.atom_labels)
        if self.current_structure is None or self.current_energy is None:
            raise RuntimeError("Initial local optimization failed.")
        self.best_energy = self.current_energy
        self.best_structure = self.current_structure.copy()
        self.energy_history.append(self.current_energy)
        self.trajectory.append(self.current_structure.copy())
        print(f"Initial Energy: {self.current_energy:.6f} Hartree")

        for step in range(self.config.max_steps):
            print(f"\nBHMC Step {step + 1}/{self.config.max_steps}")

            # TODO implement different operators
            if operator_type == "local":
                new_structure = self.local_ops.random_displacement_submolecule(
                    self.current_structure, submolecule_indices, delta_range=(0.1, 0.3))
            elif operator_type == "nonlocal":
                new_structure = self.nonlocal_ops.twist_operator(
                    self.current_structure, submolecule_indices, rotation_angle_range=(0, 360))
            else:
                raise ValueError(f"Unknown operator type: {operator_type}")
            # Local Optimization
            optimized_structure, optimized_energy = self.optimizer.optimize(new_structure)
            if optimized_structure is None or optimized_energy is None:
                print("  Local optimization failed, rejecting move.")
                self.n_rejected += 1
                self.energy_history.append(self.current_energy)
                self.trajectory.append(self.current_structure.copy())
                continue
            print(f"  Optimized Energy: {optimized_energy:.6f} Hartree")
            # Metropolis Acceptance Criterion
            if self.metropolis_accept(optimized_energy, self.current_energy):
                self.current_structure = optimized_structure
                self.current_energy = optimized_energy
                self.n_accepted += 1
                print("  Move accepted.")
                # Update best structure
                if optimized_energy < self.best_energy:
                    self.best_energy = optimized_energy
                    self.best_structure = optimized_structure.copy()
            else:
                self.n_rejected += 1
                print("  Move rejected.")
            self.energy_history.append(self.current_energy)
            self.trajectory.append(self.current_structure.copy())
        print(f"\nBHMC Completed: Accepted {self.n_accepted}, Rejected {self.n_rejected}")
        return self.trajectory, self.energy_history, self.trajectory
    

class BHMCAnalyzer:
    """  
    Helper Class for Analyzing the BHMC Results
    """
    def __init__(self):
        pass 

    def plot_energy_distribution_phase_a(self, phase_a_structures: List[Tuple[Molecule, float]], filename: str = "figures/phase_a_energy_distribution.png"):
        """
        Plots the energy distribution of structures collected in Phase A.
        """
        import matplotlib.pyplot as plt

        energies = [energy for _, energy in phase_a_structures]
        
        plt.figure(figsize=(10, 6))
        plt.hist(energies, bins=20, color='skyblue', edgecolor='black')
        plt.xlabel("Energy (Hartree)", fontsize=12)
        plt.ylabel("Frequency", fontsize=12)
        plt.title("Energy Distribution of Phase A Structures", fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()


# =============================== TEST Functions ==================================

def test_twist_operator():
    """   
    Test function for the twist operator makes a simple molecule and applies is 50 times and 
    writes a trajectory file.
    """
    h2o_molecule=""" 
    6
    H2O Nonopt
    O     -4.651077    1.267697   -0.000000
    H     -5.057317    1.928988   -0.547872
    H     -4.969101    1.355490    0.890872
    O     -6.998837    4.503183    0.000001
    H     -7.632732    4.320290   -0.683540
    H     -7.009572    5.431464    0.201689
    """
    molecule = Molecule.from_xyz(h2o_molecule, name="h2o_dimer")
    molecule.compute_bonds()
    submolecules = molecule.fragment_by_connectivity()
    submol_indices = [submol.get_index_in_parent() for submol in submolecules]
    NonlocalOperators = NonLocalOperators()
    traj = []
    for i in range(50):
        twisted_mol = NonlocalOperators.twist_operator(molecule=molecule, submolecule_indices=submol_indices)
        traj.append(twisted_mol)
    with open("twist_operator_test.xyz", "w") as file:
        for mol in traj:
            file.write(mol.to_xyz_string())
            file.write("\n")

def test_random_displacement_submolecule():
    """   
    Test function for random displacement of submolecule.
    """
    h2o_molecule=""" 
    6
    H2O Nonopt
    O     -4.651077    1.267697   -0.000000
    H     -5.057317    1.928988   -0.547872
    H     -4.969101    1.355490    0.890872
    O     -6.998837    4.503183    0.000001
    H     -7.632732    4.320290   -0.683540
    H     -7.009572    5.431464    0.201689
    """
    molecule = Molecule.from_xyz(h2o_molecule, name="h2o_dimer")
    molecule.compute_bonds()
    submolecules = molecule.fragment_by_connectivity()
    submol_indices = [submol.get_index_in_parent() for submol in submolecules]
    local_ops = LocalOperators()
    traj = []
    for i in range(50):
        displaced_mol = local_ops.random_displacement_submolecule(molecule=molecule, submolecule_indices=submol_indices)
        traj.append(displaced_mol)
    with open("random_displacement_submolecule_test.xyz", "w") as file:
        for mol in traj:
            file.write(mol.to_xyz_string())
            file.write("\n")