import numpy as np  
from typing import List, Optional, Tuple, Dict, Union, Callable
from dataclasses import dataclass
from molecule_class import Molecule
from transformations import LocalOperators, NonLocalOperators
from cluster import BHMCAnalyzer


from box import SimulationBox
from psi4_interface import Psi4Calculator, Config, direct_energy 
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
        config = Config(method=self.method, basis=basis)
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
        config = Config(method=method, basis=basis)
        self.calculator = Psi4Calculator(config)

    def evaluate(self, molecule: Molecule) -> float:
        """   
        Evaluates the energy of the given molecule.
        """
        energy = direct_energy(molecule, method=self.method, basis=self.basis_set) 
        return energy

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
    """Worker function for Phase A"""
    (initial_molecule, submolecule_indices, n_structures,
     operator_configs, config_dict, worker_id, sim_box_dict) = args
    
    # Setup
    sim_box = SimulationBox.from_dict(sim_box_dict) if sim_box_dict else None
    nonlocal_ops = NonLocalOperators(simulation_box=sim_box)
    k_B = 8.617333262145e-5
    
    # Setup evaluator
    evaluator = EnergyEvaluator(
        method=config_dict.get('method', 'hf'),
        basis=config_dict.get('basis', 'sto-3g'),
    )
    
    # Recreate operator functions
    operators = []
    weights = []
    op_map = {
        'twist': nonlocal_ops.twist_operator,
        'large_displacement': nonlocal_ops.large_displacement,
        'mirror': nonlocal_ops.mirror_operator,
        'roto_reflection': nonlocal_ops.roto_reflection_operator,
        'exchange': nonlocal_ops.exchange_operator,
       # 'inversion': nonlocal_ops.inversion_operator,
       # 'screw_axis': nonlocal_ops.screw_axis_operator,
    }
    
    for op_config in operator_configs:
        if op_config['name'] in op_map:
            operators.append({
                'name': op_config['name'],
                'func': op_map[op_config['name']]
            })
            weights.append(op_config['weight'])
    
    weights = np.array(weights)
    weights /= weights.sum()
    
    # Start BHMC chain
    current_structure = initial_molecule.copy()
    current_energy = evaluator.evaluate(current_structure)
    
    # WICHTIG: Prüfen ob Energie berechnet werden konnte
    if current_energy is None:
        print(f"Worker {worker_id}: ERROR - Initial energy calculation failed!")
        print(f"Worker {worker_id}: Molecule has {len(current_structure.coordinates)} atoms")
        print(f"Worker {worker_id}: Config: {config_dict}")
        return []  # Leere Liste zurückgeben wenn erste Berechnung fehlschlägt
    
    accepted_structures = []
    n_accepted = 0
    n_rejected = 0
    
    print(f"Worker {worker_id}: Starting with energy {current_energy:.6f} Hartree")
    
    for step in range(n_structures):
        op_idx = np.random.choice(len(operators), p=weights)
        operator = operators[op_idx]
        
        try:
            new_structure = operator['func'](current_structure, submolecule_indices)
            new_energy = evaluator.evaluate(new_structure)
            
            # Skip wenn Energie-Berechnung fehlgeschlagen
            if new_energy is None:
                n_rejected += 1
                continue
            
            # Metropolis acceptance
            if new_energy < current_energy:
                accept = True
            else:
                hartree_to_ev = 27.2114
                delta_e_ev = (new_energy - current_energy) * hartree_to_ev
                beta = 1.0 / (k_B * config_dict.get('temperature', 300.0))
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
            print(f"Worker {worker_id} step {step}: Error: {e}")
            n_rejected += 1
            continue
        
        if (step + 1) % 10 == 0:
            print(f"Worker {worker_id}: {step + 1}/{n_structures}, Accepted: {n_accepted}")
    
    print(f"Worker {worker_id}: Completed. Accepted: {n_accepted}/{n_structures}")
    
    return accepted_structures


class MultiPhaseBHMC:
    """    
    Multi-Phase Basin Hopping Monte Carlo:

    Phase A (Exploration): Nonlocal Optimization, No Local Optimization Generation of Grid for Parallel Local Optimizations
    """
    def __init__(self, config: BHMCConfig, simulation_box: Optional[SimulationBox] = None):
        """   
        Initializes the Multi-Phase BHMC with the given configuration.
        
        Args:
            config: BHMC configuration
            simulation_box: Optional simulation box for constraining structures
        """
        self.config = config
        self.simulation_box = simulation_box
        self.optimizer = LocalOptimizer(method=config.method, basis=config.basis, verbose=False)
        self.energy_evaluator = EnergyEvaluator(method=config.method, basis=config.basis)
        
        # Initialize operators with simulation box
        self.local_ops = LocalOperators(simulation_box=simulation_box)
        self.nonlocal_ops = NonLocalOperators(simulation_box=simulation_box)

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
        Initialize the Operator Registry with all available operators.
        """
        return [
            # Nonlocal Operators
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
                weight=1.2 
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
                weight=0.8
            ),
            OperatorConfig(
                name="exchange",
                operator_func=lambda m,s: self.nonlocal_ops.exchange_operator(m, s),
                operator_type="nonlocal",
                weight=0.7
            ),
            OperatorConfig(
                name="inversion",
                operator_func=lambda m,s: self.nonlocal_ops.inversion_operator(m, s),
                operator_type="nonlocal",
                weight=0.5
            ),
            OperatorConfig(
                name="screw_axis",
                operator_func=lambda m,s: self.nonlocal_ops.screw_axis_operator(m, s),
                operator_type="nonlocal",
                weight=0.5
            ),
            # Local Operators
            OperatorConfig(
                name="local_displacement",
                operator_func=lambda m,s: self.local_ops.random_displacement_submolecule(m, s),
                operator_type="local",
                weight=1.0
            ),
            OperatorConfig(
                name="local_rotation",
                operator_func=lambda m,s: self.local_ops.random_rotation_submolecule(m, s),
                operator_type="local",
                weight=0.8
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
        if self.simulation_box:
            print(f"Simulation Box: {self.simulation_box}")
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
        
        # Serialize simulation box
        sim_box_dict = self.simulation_box.to_dict() if self.simulation_box else None
        
        # Prepare arguments for each worker
        args_list = [
            (
                initial_molecule,
                submolecule_indices,
                n_structures_per_worker,
                operator_configs,
                config_dict,
                worker_id,
                sim_box_dict
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
    
    @staticmethod
    def analyse_phase_a_results(phase_a_structures: List[Tuple[Molecule, float]], submolecule_indices: List[List[int]] = None) -> List[Molecule]:
        """   
        Analyses the results from the exploration phase
        Args:
            phase_a_structures: List of (Molecule, energy) tuples from Phase A
            submolecule_indices: Optional list of submolecule indices for clustering analysis
        
        """
        analyzer = BHMCAnalyzer(submolecule_indices=submolecule_indices)
        analyzer.add_structures_batch(phase_a_structures)
        analyzer.rmsd_filtering() # Filter out duplicates based on RMSD 
        analyzer.plot_energy_distribution()
        analyzer.plot_energy_vs_rmsd()
        analyzer.plot_rg_vs_energy()
        analyzer.plot_int_d_vs_e()
        analyzer.plot_pca_agglomerative(n_clusters=10)

        # Run Clustering 
        analyzer.AgglomerativeClustering(n_clusters=10)
        representatives = analyzer.get_cluster_representatives()

        return representatives

