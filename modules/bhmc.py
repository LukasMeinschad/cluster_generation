"""
Basin Hopping Monte Carlo (BHMC) implementation for molecular cluster exploration.

This module provides a parallel BHMC implementation for exploring the potential
energy surface of molecular clusters using various nonlocal transformation operators.
"""

import numpy as np  
import random
from typing import List, Optional, Tuple, Dict
from dataclasses import dataclass
from multiprocessing import Pool

# Module imports
from molecule_class import Molecule
from transformations import NonLocalOperators
from cluster import BHMCAnalyzer
from box import SimulationBox
from psi4_interface import direct_energy


@dataclass
class BHMCConfig:
    """Configuration for Basin Hopping Monte Carlo sampling."""
    temperature: float = 400.0  # Temperature in Kelvin
    method: str = "hf"          # Quantum chemistry method
    basis: str = "sto-3g"       # Basis set


class EnergyEvaluator:
    """Energy evaluator for molecular structures using Psi4."""
    
    def __init__(self, method: str = "hf", basis: str = "sto-3g"):
        """Initialize the energy evaluator.
        
        Args:
            method: Quantum chemistry method
            basis: Basis set
        """
        self.method = method
        self.basis_set = basis

    def evaluate(self, molecule: Molecule) -> Optional[float]:
        """Evaluate the energy of a molecule.
        
        Args:
            molecule: Molecule to evaluate
            
        Returns:
            Energy in Hartree, or None if calculation failed
        """
        return direct_energy(molecule, method=self.method, basis=self.basis_set)


# ================================ Multiprocessing Worker ================================

def _phase_a_worker(args: Tuple) -> List[Tuple[Molecule, float]]:
    """Worker function for parallel BHMC Phase A exploration.
    
    Each worker runs an independent BHMC chain using nonlocal operators.
    
    Args:
        args: Tuple containing (initial_molecule, submolecule_indices, n_structures,
              operator_list, config_dict, worker_id, sim_box_dict)
              
    Returns:
        List of (Molecule, energy) tuples for accepted structures
    """
    (initial_molecule, submolecule_indices, n_structures,
     operator_list, config_dict, worker_id, sim_box_dict) = args
    
    # Constants
    k_B = 8.617333262145e-5  # Boltzmann constant in eV/K
    HARTREE_TO_EV = 27.2114
    
    # Setup simulation box and operators
    sim_box = SimulationBox.from_dict(sim_box_dict) if sim_box_dict else None
    nonlocal_ops = NonLocalOperators(simulation_box=sim_box)
    
    # Map operator names to functions
    op_map = {
        'twist': nonlocal_ops.twist_operator,
        'large_displacement': nonlocal_ops.large_displacement,
        'mirror': nonlocal_ops.mirror_operator,
        'roto_reflection': nonlocal_ops.roto_reflection_operator,
        'exchange': nonlocal_ops.exchange_operator,
        'random_so3': nonlocal_ops.random_so3_operator,
    }
    
    # Build operator list with weights
    operators = []
    weights = []
    for op_name, weight in operator_list:
        if op_name in op_map:
            operators.append({'name': op_name, 'func': op_map[op_name]})
            weights.append(weight)
    
    if not operators:
        print(f"Worker {worker_id}: ERROR - No valid operators configured!")
        return []
    
    weights = np.array(weights)
    weights /= weights.sum()
    
    # Setup energy evaluator
    evaluator = EnergyEvaluator(
        method=config_dict['method'],
        basis=config_dict['basis']
    )
    
    # Initialize BHMC chain
    current_structure = initial_molecule.copy()
    current_energy = evaluator.evaluate(current_structure)
    
    if current_energy is None:
        print(f"Worker {worker_id}: ERROR - Initial energy calculation failed!")
        return []
    
    print(f"Worker {worker_id}: Starting with energy {current_energy:.6f} Hartree")
    
    # BHMC loop
    accepted_structures = []
    n_accepted = 0
    n_rejected = 0
    temperature = config_dict['temperature']
    beta = 1.0 / (k_B * temperature)
    
    for step in range(n_structures):
        # Select and apply operator
        op_idx = np.random.choice(len(operators), p=weights)
        operator = operators[op_idx]
        
        try:
            new_structure = operator['func'](current_structure, submolecule_indices)
            new_energy = evaluator.evaluate(new_structure)
            
            if new_energy is None:
                n_rejected += 1
                continue
            
            # Metropolis acceptance criterion
            if new_energy < current_energy:
                accept = True
            else:
                delta_e_ev = (new_energy - current_energy) * HARTREE_TO_EV
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
            print(f"Worker {worker_id} step {step}: Error applying {operator['name']}: {e}")
            n_rejected += 1
            continue
        
        # Progress update
        if (step + 1) % 10 == 0:
            print(f"Worker {worker_id}: {step + 1}/{n_structures}, Accepted: {n_accepted}, "
                  f"Rate: {n_accepted/(n_accepted+n_rejected)*100:.1f}%")
    
    print(f"Worker {worker_id}: Completed. Accepted: {n_accepted}/{n_structures} "
          f"({n_accepted/n_structures*100:.1f}%)")
    
    return accepted_structures


# ================================ Main BHMC Class ================================

class MultiPhaseBHMC:
    """Multi-Phase Basin Hopping Monte Carlo for molecular cluster exploration.
    
    This implementation uses parallel BHMC chains in Phase A to explore the potential
    energy surface using various nonlocal transformation operators.
    """
    
    # Default operator configuration (name, weight)
    DEFAULT_OPERATORS = [
        ('twist', 1.5),
        ('large_displacement', 1.2),
        ('mirror', 1.0),
        ('roto_reflection', 0.8),
        ('exchange', 0.7),
        ('random_so3', 0.5),
    ]
    
    def __init__(self, 
                 config: BHMCConfig,
                 simulation_box: Optional[SimulationBox] = None,
                 operators: Optional[List[Tuple[str, float]]] = None):
        """Initialize the BHMC sampler.
        
        Args:
            config: BHMC configuration
            simulation_box: Optional simulation box for constraining structures
            operators: Optional list of (operator_name, weight) tuples. 
                      Uses DEFAULT_OPERATORS if not provided.
        """
        self.config = config
        self.simulation_box = simulation_box
        self.operators = operators or self.DEFAULT_OPERATORS
        
        # Storage for Phase A results
        self.phase_a_structures: List[Tuple[Molecule, float]] = []
    
    def run_phase_a(self, 
                    initial_molecule: Molecule,
                    submolecule_indices: List[List[int]],
                    n_structures_per_worker: int = 300,
                    n_processes: int = 10) -> List[Tuple[Molecule, float]]:
        """Run Phase A: Parallel exploration with independent BHMC chains.
        
        Each worker runs an independent BHMC chain using nonlocal operators to
        explore the potential energy surface.
        
        Args:
            initial_molecule: Starting molecular structure
            submolecule_indices: List of submolecule atom indices for operators
            n_structures_per_worker: Number of BHMC steps per worker
            n_processes: Number of parallel workers
            
        Returns:
            List of (Molecule, energy) tuples for all accepted structures
        """
        print(f"\n{'='*70}")
        print(f"Phase A: Global Exploration with Parallel BHMC Chains")
        print(f"{'='*70}")
        print(f"Configuration:")
        print(f"  Workers: {n_processes}")
        print(f"  Steps per worker: {n_structures_per_worker}")
        print(f"  Total structures: {n_processes * n_structures_per_worker}")
        print(f"  Method/Basis: {self.config.method}/{self.config.basis}")
        print(f"  Temperature: {self.config.temperature} K")
        print(f"  Operators: {', '.join(op[0] for op in self.operators)}")
        if self.simulation_box:
            print(f"  Simulation Box: {self.simulation_box}")
        print(f"{'='*70}\n")
        
        # Configuration dictionary for workers
        config_dict = {
            'method': self.config.method,
            'basis': self.config.basis,
            'temperature': self.config.temperature
        }
        
        # Serialize simulation box if present
        sim_box_dict = self.simulation_box.to_dict() if self.simulation_box else None
        
        # Prepare arguments for each worker
        args_list = [
            (
                initial_molecule,
                submolecule_indices,
                n_structures_per_worker,
                self.operators,
                config_dict,
                worker_id,
                sim_box_dict
            )
            for worker_id in range(n_processes)
        ]
        
        # Run parallel workers
        print("Starting parallel BHMC chains...\n")
        
        with Pool(processes=n_processes) as pool:
            results = pool.map(_phase_a_worker, args_list)
        
        # Collect results
        all_accepted_structures = []
        for worker_id, worker_results in enumerate(results):
            n_accepted = len(worker_results)
            all_accepted_structures.extend(worker_results)
            print(f"Worker {worker_id}: Collected {n_accepted} structures")
        
        self.phase_a_structures = all_accepted_structures
        
        # Print statistics
        total_generated = n_processes * n_structures_per_worker
        total_accepted = len(all_accepted_structures)
        
        print(f"\n{'='*70}")
        print(f"Phase A Complete")
        print(f"{'='*70}")
        print(f"Total structures generated: {total_generated}")
        print(f"Total structures accepted: {total_accepted}")
        print(f"Overall acceptance rate: {total_accepted/total_generated*100:.2f}%")
        
        if all_accepted_structures:
            energies = [e for _, e in all_accepted_structures]
            print(f"\nEnergy Statistics (Hartree):")
            print(f"  Range: [{min(energies):.6f}, {max(energies):.6f}]")
            print(f"  Mean: {np.mean(energies):.6f}")
            print(f"  Std: {np.std(energies):.6f}")
        
        print(f"{'='*70}\n")
        
        return self.phase_a_structures
    
    @staticmethod
    def analyse_phase_a_results(
            phase_a_structures: List[Tuple[Molecule, float]], 
            submolecule_indices: Optional[List[List[int]]] = None,
            n_clusters: int = 10) -> List[Molecule]:
        """Analyze Phase A results and extract representative structures.
        
        Performs RMSD filtering, clustering, and various analyses on the
        Phase A structures to identify unique cluster representatives.
        
        Args:
            phase_a_structures: List of (Molecule, energy) tuples from Phase A
            submolecule_indices: Optional submolecule indices for clustering
            n_clusters: Number of clusters to identify
        
        Returns:
            List of representative Molecule structures from each cluster
        """
        print(f"\n{'='*70}")
        print(f"Analyzing Phase A Results")
        print(f"{'='*70}")
        print(f"Total structures: {len(phase_a_structures)}")
        print(f"Target clusters: {n_clusters}")
        print(f"{'='*70}\n")
        
        analyzer = BHMCAnalyzer(submolecule_indices=submolecule_indices)
        analyzer.add_structures_batch(phase_a_structures)
        
        # Filter duplicates
        analyzer.rmsd_filtering()
        
        # Generate analysis plots
        analyzer.plot_energy_distribution()
        analyzer.plot_energy_vs_rmsd()
        analyzer.plot_rg_vs_energy()
        analyzer.plot_int_d_vs_e()
        analyzer.plot_pca_agglomerative(n_clusters=n_clusters)
        
        # Perform clustering and extract representatives
        analyzer.AgglomerativeClustering(n_clusters=n_clusters)
        representatives = analyzer.get_cluster_representatives()
        
        print(f"\nExtracted {len(representatives)} cluster representatives\n")
        
        return representatives


