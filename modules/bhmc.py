"""
Basin Hopping Monte Carlo (BHMC) implementation for molecular cluster exploration.

This module provides a parallel BHMC implementation for exploring the potential
energy surface of molecular clusters using various nonlocal transformation operators.
"""

import numpy as np  
import random
import logging
from typing import List, Optional, Tuple, Dict
from dataclasses import dataclass
from multiprocessing import Pool
import copy

# Module imports
from molecule_class import Molecule
from transformations import NonLocalOperators
from cluster import BHMCAnalyzer
from box import SimulationBox
from psi4_interface import direct_energy

# Logger import 
from logger import Logger

@dataclass
class BHMCConfig:
    """Configuration for Basin Hopping Monte Carlo sampling."""
    temperature: float = 400.0  # Temperature in Kelvin
    method: str = "hf"          # Quantum chemistry method
    basis: str = "sto-3g"       # Basis set
    verbose: bool = False       # Enable debug output in workers
    adaptive_operators: bool = True  # Use adaptive scaling for operators


class EnergyEvaluator:
    """Energy evaluator for molecular structures using Psi4."""
    
    def __init__(self, method: str = "hf", basis: str = "sto-3g"):
        self.method = method
        self.basis_set = basis

    def evaluate(self, molecule: Molecule) -> Tuple[Optional[float], Optional[List[float]], Optional[float]]:
        """ 
        Evaluate energy and the dipole moment using Psi4

        Returns:
            Tuple of (energy in Hartree, dipole vector [dx,dy,dz], dipole magnitude) or (None, None, None) if evaluation fails.
        """
        return direct_energy(molecule, method=self.method, basis=self.basis_set)


# =============================== Worker Logger ================================

def _get_worker_logger(log_file: Optional[str] = None, worker_id: int = 0) -> Optional[Logger]:
    """Create a logger instance inside a worker process.
    
    Args:
        log_file: Path to the worker log file, or None to disable.
        worker_id: Unique identifier for the worker.
        
    Returns:
        Logger instance or None if logging is disabled.
    """
    if log_file is None:
        return None
    return Logger(
        name=f"bhmc_worker_{worker_id}",
        log_file=log_file,
        file_mode="a"
    )


# ================================ Multiprocessing Worker ================================

def _phase_a_worker(args: Tuple) -> List[Tuple[Molecule, float, List[float], float]]:
    """Worker function for parallel BHMC Phase A exploration.
    
    Args:
        args: Tuple containing (initial_molecule, submolecule_indices, n_structures,
              operator_list, config_dict, worker_id, sim_box_dict, worker_log_file)
              
    Returns:
       Dictionary with keys:
              'accepted_structures': List of (Molecule, energy, dipole_vec, dipole_magnitude) tuples
              'energy_trajectory': List of (step, current_energy) for every step
              'worker_id': Unique identifier for the worker

    """
    (initial_molecule, submolecule_indices, n_structures,
     operator_list, config_dict, worker_id, sim_box_dict, worker_log_file) = args
    
    # Create worker-local logger — writes to bhmc_workers.log, NOT cluster_gen.out
    logger = _get_worker_logger(worker_log_file, worker_id)

    def _log(msg: str, level: str = "info") -> None:
        if config_dict.get('verbose', False):
            print(msg)
        if logger:
            getattr(logger, level)(msg)

    # Constants
    k_B = 8.617333262145e-5  # Boltzmann constant in eV/K
    HARTREE_TO_EV = 27.2114
    
    # Setup simulation box and operators
    sim_box = SimulationBox.from_dict(sim_box_dict) if sim_box_dict else None
    nonlocal_ops = NonLocalOperators(simulation_box=sim_box)
    
    adaptive = config_dict.get('adaptive_operators', True)
    
    adaptive_operators = {'twist', 'large_displacement', 'roto_reflection'}
    
    op_map = {
        'twist': nonlocal_ops.twist_operator,
        'large_displacement': nonlocal_ops.large_displacement,
        'mirror': nonlocal_ops.mirror_operator,
        'roto_reflection': nonlocal_ops.roto_reflection_operator,
        'exchange': nonlocal_ops.exchange_operator,
        'random_so3': nonlocal_ops.random_so3_operator,
    }
    
    operators = []
    weights = []
    for op_name, weight in operator_list:
        if op_name in op_map:
            operators.append({
                'name': op_name, 
                'func': op_map[op_name],
                'supports_adaptive': op_name in adaptive_operators
            })
            weights.append(weight)
    
    if not operators:
        _log(f"Worker {worker_id}: No valid operators found!", level="error")
        return {'accepted_structures': [], 'energy_trajectory': [], 'worker_id': worker_id}
    
    
    weights = np.array(weights)
    weights /= weights.sum()
    
    evaluator = EnergyEvaluator(
        method=config_dict['method'],
        basis=config_dict['basis']
    )
    
    current_structure = copy.deepcopy(initial_molecule)
    eval_result = evaluator.evaluate(current_structure)
    
    if eval_result is None:
        _log(f"Worker {worker_id}: Failed to evaluate initial energy!", level="error")
        return {'accepted_structures': [], 'energy_trajectory': [], 'worker_id': worker_id}
    
    current_energy, current_dipole_vec, current_dipole_mag = eval_result
    
    if current_energy is None:
        _log(f"Worker {worker_id}: Failed to evaluate initial energy!", level="error")
        return {'accepted_structures': [], 'energy_trajectory': [], 'worker_id': worker_id}
    
    _log(f"Worker {worker_id}: Starting energy = {current_energy:.6f} Hartree")

    accepted_structures = []
    energy_trajectory = [(0, current_energy)] # Track energy at each step
    n_accepted = 0
    n_rejected = 0
    temperature = config_dict['temperature']
    beta = 1.0 / (k_B * temperature)
    
    for step in range(n_structures):
        op_idx = np.random.choice(len(operators), p=weights)
        operator = operators[op_idx]
        
        try:
            if operator['supports_adaptive']:
                new_structure = operator['func'](
                    current_structure, submolecule_indices, adaptive=adaptive
                )
            else:
                new_structure = operator['func'](
                    current_structure, submolecule_indices
                )
            
            if np.allclose(new_structure.coordinates, current_structure.coordinates):
                _log(f"Worker {worker_id} step {step}: {operator['name']} — no coordinate change (box rejection?)", level="warning")
                n_rejected += 1
                energy_trajectory.append((step+1, current_energy))
                continue
            
            new_eval = evaluator.evaluate(new_structure)
            
            if new_eval is None:
                _log(f"Worker {worker_id} step {step}: Energy eval failed after {operator['name']}", level="warning")
                n_rejected += 1
                energy_trajectory.append((step+1, current_energy))
                continue
            
            new_energy, new_dipole_vec, new_dipole_mag = new_eval
            
            if new_energy is None:
                _log(f"Worker {worker_id} step {step}: Energy eval returned None after {operator['name']}", level="warning")
                n_rejected += 1
                energy_trajectory.append((step+1, current_energy))
                continue
            
            if new_energy < current_energy:
                accept = True
            else:
                delta_e_ev = (new_energy - current_energy) * HARTREE_TO_EV
                prob = np.exp(-beta * delta_e_ev)
                accept = random.random() < prob
            
            if accept:
                current_structure = copy.deepcopy(new_structure)
                current_energy = new_energy
                current_dipole_vec = new_dipole_vec
                current_dipole_mag = new_dipole_mag
                accepted_structures.append((
                    copy.deepcopy(current_structure), 
                    current_energy,
                    current_dipole_vec,
                    current_dipole_mag
                ))
                n_accepted += 1
                _log(f"Worker {worker_id} step {step}: ACCEPTED {operator['name']} E={current_energy:.6f} Ha dipole={current_dipole_mag:.4f} D")
            else:
                n_rejected += 1
            
            # Record energy trajectory for every step, accepted or not
            energy_trajectory.append((step+1, current_energy))
                
        except Exception as e:
            _log(f"Worker {worker_id} step {step}: Exception in {operator['name']}: {e}", level="error")
            n_rejected += 1
            energy_trajectory.append((step+1, current_energy))
            continue
        
        if (step + 1) % 10 == 0:
            rate = n_accepted / (n_accepted + n_rejected) * 100 if (n_accepted + n_rejected) > 0 else 0.0
            _log(f"Worker {worker_id}: {step+1}/{n_structures} — accepted={n_accepted}, rate={rate:.1f}%")
    
    total = n_accepted + n_rejected
    rate = n_accepted / total * 100 if total > 0 else 0.0
    _log(f"Worker {worker_id}: DONE — {n_accepted}/{total} accepted ({rate:.1f}%)")
    
    return {
        'accepted_structures': accepted_structures,
        'energy_trajectory': energy_trajectory,
        'worker_id': worker_id}


# ================================ Main BHMC Class ================================

class MultiPhaseBHMC:
    """Multi-Phase Basin Hopping Monte Carlo for molecular cluster exploration."""
    
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
                 operators: Optional[List[Tuple[str, float]]] = None,
                 logger: Optional[Logger] = None,
                 worker_log_file: str = "bhmc_workers.log"):
        """Initialize the BHMC sampler.
        
        Args:
            config: BHMC configuration
            simulation_box: Optional simulation box for constraining structures
            operators: Optional list of (operator_name, weight) tuples.
            logger: Optional Logger for summary output (cluster_gen.out)
            worker_log_file: Path to log file for worker-level detail.
                            Set to None to disable worker logging.
        """
        self.config = config
        self.simulation_box = simulation_box
        self.operators = operators or self.DEFAULT_OPERATORS
        self.logger = logger
        self.worker_log_file = worker_log_file
        self.worker_trajectories = {}  # Store energy trajectories for each worker
        
        # Storage for Phase A results
        self.phase_a_structures: List[Tuple[Molecule, float, List[float], float]] = []

    def _log(self, msg: str, level: str = "info") -> None:
        """Log a summary message to the main log file."""
        if self.config.verbose:
            print(msg)
        if self.logger:
            getattr(self.logger, level)(msg)

    def run_phase_a(self, 
                    initial_molecules: List[Molecule],
                    submolecule_indices: List[List[int]],
                    n_structures_per_worker: int = 300,
                    n_processes: int = 10) -> List[Tuple[Molecule, float, List[float], float]]:
        """Run Phase A: Parallel exploration with independent BHMC chains.
        
        Each worker starts from a different initial configuration for better PES coverage.
        """
        if n_processes is None:
            n_processes = len(initial_molecules)
        
        if len(initial_molecules) != n_processes:
            raise ValueError(f"Number of initial molecules ({len(initial_molecules)}) must match n_processes ({n_processes})")
        
        if self.logger:
            self.logger.header("Phase A: Global PES Exploration")
        
        self._log(f"Configuration:")
        self._log(f"  Workers: {n_processes}")
        self._log(f"  Steps per worker: {n_structures_per_worker}")
        self._log(f"  Total steps: {n_processes * n_structures_per_worker}")
        self._log(f"  Method/Basis: {self.config.method}/{self.config.basis}")
        self._log(f"  Temperature: {self.config.temperature} K")
        self._log(f"  Adaptive operators: {self.config.adaptive_operators}")
        self._log(f"  Operators: {', '.join(op[0] for op in self.operators)}")
        self._log(f"  Independent initial configs: {n_processes}")
        if self.simulation_box:
            self._log(f"  Simulation Box: {self.simulation_box}")
        if self.worker_log_file:
            self._log(f"  Worker log: {self.worker_log_file}")
        
        config_dict = {
            'method': self.config.method,
            'basis': self.config.basis,
            'temperature': self.config.temperature,
            'verbose': self.config.verbose,
            'adaptive_operators': self.config.adaptive_operators,
        }
        
        sim_box_dict = self.simulation_box.to_dict() if self.simulation_box else None
        
        args_list = [
            (
                initial_molecules[worker_id],
                submolecule_indices,
                n_structures_per_worker,
                self.operators,
                config_dict,
                worker_id,
                sim_box_dict,
                self.worker_log_file,
            )
            for worker_id in range(n_processes)
        ]
        
        self._log("Starting parallel BHMC chains...")
        
        with Pool(processes=n_processes) as pool:
            results = pool.map(_phase_a_worker, args_list)
        
        # Collect results and trajectories
        all_accepted_structures = []
        self.worker_trajectories = {}  # Store for plotting
        
        for worker_result in results:
            wid = worker_result['worker_id']
            accepted = worker_result['accepted_structures']
            trajectory = worker_result['energy_trajectory']
            
            all_accepted_structures.extend(accepted)
            self.worker_trajectories[wid] = trajectory
            self._log(f"  Worker {wid}: {len(accepted)} structures accepted, {len(trajectory)} trajectory points")
        
        self.phase_a_structures = all_accepted_structures
        
        # Summary statistics
        total_generated = n_processes * n_structures_per_worker
        total_accepted = len(all_accepted_structures)
        
        if self.logger:
            self.logger.section("Phase A Results")
        
        self._log(f"Total steps: {total_generated}")
        self._log(f"Total accepted: {total_accepted}")
        if total_generated > 0:
            self._log(f"Overall acceptance rate: {total_accepted/total_generated*100:.2f}%")
        
        if all_accepted_structures:
            energies = [e for _, e, _, _ in all_accepted_structures]
            dipoles = [d for _, _, _, d in all_accepted_structures]
            self._log(f"Energy Statistics (Hartree):")
            self._log(f"  Min:  {min(energies):.6f}")
            self._log(f"  Max:  {max(energies):.6f}")
            self._log(f"  Mean: {np.mean(energies):.6f}")
            self._log(f"  Std:  {np.std(energies):.6f}")
            self._log(f"Dipole Statistics (Debye):")
            self._log(f"  Min:  {min(dipoles):.4f}")
            self._log(f"  Max:  {max(dipoles):.4f}")
            self._log(f"  Mean: {np.mean(dipoles):.4f}")
        
        # Plot worker trajectories
        self.plot_worker_trajectories()
        
        return self.phase_a_structures
    
    def analyse_phase_a_results(
            self,
            phase_a_structures: List[Tuple[Molecule, float]], 
            submolecule_indices: Optional[List[List[int]]] = None,
            n_clusters: int = 10,
            simulation_box: Optional[SimulationBox] = None,
            logger: Optional[Logger] = None
            ) -> List[Molecule]:
        """Analyze Phase A results and extract representative structures.
        
        Args:
            phase_a_structures: List of (Molecule, energy) tuples from Phase A
            submolecule_indices: Optional submolecule indices for clustering
            n_clusters: Number of clusters to identify
            simulation_box: Optional simulation box for trajectory plots
        
        Returns:
            List of representative Molecule structures from each cluster
        """
        if self.logger:
            self.logger.header("Phase A Analysis")
        
        self._log(f"Total structures: {len(phase_a_structures)}")
        self._log(f"Target clusters: {n_clusters}")
        
        analyzer = BHMCAnalyzer(submolecule_indices=submolecule_indices, logger=logger)
        analyzer.add_structures_batch(phase_a_structures)
        
        # Get energy statistics for analysis log
        energy_stats = analyzer.get_energy_statistics()
        
        # Log dipole statistics
        dipoles = np.array([s.dipole_magnitude for s in analyzer.structures])
        if np.any(dipoles > 0):
            self._log(f"Dipole Magnitude Statistics (Debye):")
            self._log(f"  Min:  {np.min(dipoles):.4f}")
            self._log(f"  Max:  {np.max(dipoles):.4f}")
            self._log(f"  Mean: {np.mean(dipoles):.4f}")
            self._log(f"  Std:  {np.std(dipoles):.4f}")
        
     


        analyzer.plot_com_trajectory_2d_projection(
            simulation_box=simulation_box,
            save_path="figures/phase_a_com_trajectory.png",
            separate_submolecules=True
        )

        analyzer.rmsd_filtering()
        analyzer.plot_rmsd_heatmap()
        
        # Precompute the feature matrix
        analyzer.feature_matrix(normalize=True)  # Normalized features for clustering


        analyzer.plot_energy_distribution()
        analyzer.plot_energy_vs_rmsd()
        analyzer.plot_rg_vs_energy()
        analyzer.plot_int_d_vs_e()
        analyzer.plot_pca_explained_variance()
        analyzer.plot_pca_agglomerative(n_clusters=n_clusters)
        analyzer.plot_tsne_agglomerative(n_clusters=n_clusters)
        analyzer.plot_umap_agglomerative(n_clusters=n_clusters)
        

        analyzer.agglomerative_clustering(n_clusters=n_clusters)
        representatives = analyzer.get_cluster_representatives()


        analyzer.plot_cluster_populations()
        analyzer.plot_cluster_graph(
                representatives=representatives,
                
        ) 

        
        return representatives
    
    # Additional plotting utilities for worker trajectories
    def plot_worker_trajectories(self,
                                 save_path: str = "figures/worker_energy_trajectories.png") -> None:
        """  
        Plot energy vs step for each BHMC worker
        """
        import matplotlib.pyplot as plt
        from pathlib import Path

        if not hasattr(self, 'worker_trajectories') or not self.worker_trajectories:
            self._log("No worker trajectories to plot.", level="warning")
            return
        
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        n_workers = len(self.worker_trajectories)

        cmap = plt.get_cmap('tab20', n_workers)

        plt.figure(figsize=(12, 6))
        for worker_id, trajectory in self.worker_trajectories.items():
            steps, energies = zip(*trajectory)
            plt.plot(steps, energies, label=f"Worker {worker_id}", color=cmap(worker_id))

        plt.xlabel("Step")
        plt.ylabel("Energy (Hartree)")
        plt.title("Energy Trajectories of BHMC Workers")
        plt.legend(
            ncols = 4,
            loc = "upper right",
            columnspacing = 0.5,
        )
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()



# =================================== Benchmarking =========================================

def benchmark_temperature_acceptance(
    initial_molecule: Molecule,
    submolecule_indices: List[List[int]],
    simulation_box: Optional[SimulationBox] = None,
    temperatures: Optional[List[float]] = None,
    method: str = "hf",
    basis: str = "sto-3g",
    n_steps: int = 100,
    n_trials: int = 5,
    operators: Optional[List[Tuple[str, float]]] = None,
    save_plot: bool = True,
    plot_filename: str = "acceptance_vs_temperature.png",
    logger: Optional[Logger] = None,
    verbose: bool = True
) -> Dict[float, Dict[str, float]]:
    """Benchmark acceptance rates across different temperatures.
    
    Args:
        initial_molecule: Starting molecular structure
        submolecule_indices: List of submolecule atom indices
        simulation_box: Optional simulation box
        temperatures: List of temperatures to test (K).
        method: QM method for energy evaluation
        basis: Basis set
        n_steps: Number of BHMC steps per trial
        n_trials: Number of independent trials per temperature
        operators: Optional operator list.
        save_plot: Whether to save the plot
        plot_filename: Output filename for plot
        logger: Optional Logger for summary output.
        verbose: Whether to print to stdout.
        
    Returns:
        Dictionary mapping temperature to acceptance/energy statistics.
    """
    import matplotlib.pyplot as plt
    
    def _log(msg: str, level: str = "info") -> None:
        if verbose:
            print(msg)
        if logger:
            getattr(logger, level)(msg)
    
    if temperatures is None:
        temperatures = [50, 100, 200, 300, 400, 500, 750, 1000, 1500]
    
    if operators is None:
        operators = MultiPhaseBHMC.DEFAULT_OPERATORS
    
    if logger:
        logger.header("Temperature Benchmark")
    
    _log(f"Testing {len(temperatures)} temperatures with {n_trials} trials each")
    _log(f"Steps per trial: {n_steps}")
    _log(f"Method/Basis: {method}/{basis}")
    
    k_B = 8.617333262145e-5
    HARTREE_TO_EV = 27.2114
    
    sim_box = simulation_box
    nonlocal_ops = NonLocalOperators(simulation_box=sim_box)
    
    op_map = {
        'twist': nonlocal_ops.twist_operator,
        'large_displacement': nonlocal_ops.large_displacement,
        'mirror': nonlocal_ops.mirror_operator,
        'roto_reflection': nonlocal_ops.roto_reflection_operator,
        'exchange': nonlocal_ops.exchange_operator,
        'random_so3': nonlocal_ops.random_so3_operator,
    }
    
    operator_funcs = []
    weights = []
    for op_name, weight in operators:
        if op_name in op_map:
            operator_funcs.append({'name': op_name, 'func': op_map[op_name]})
            weights.append(weight)
    
    weights = np.array(weights)
    weights /= weights.sum()
    
    evaluator = EnergyEvaluator(method=method, basis=basis)
    results = {}
    
    for temp in temperatures:
        _log(f"\nT = {temp} K:")
        beta = 1.0 / (k_B * temp)
        
        trial_acceptance_rates = []
        trial_energy_changes = []
        
        for trial in range(n_trials):
            current_structure = initial_molecule.copy()
            current_energy = evaluator.evaluate(current_structure)
            
            if current_energy is None:
                _log(f"  Trial {trial+1}: Initial energy failed", level="warning")
                continue
            
            n_accepted = 0
            energy_changes = []
            
            for step in range(n_steps):
                op_idx = np.random.choice(len(operator_funcs), p=weights)
                operator = operator_funcs[op_idx]
                
                try:
                    new_structure = operator['func'](current_structure, submolecule_indices)
                    new_energy = evaluator.evaluate(new_structure)
                    
                    if new_energy is None:
                        continue
                    
                    delta_e = new_energy - current_energy
                    delta_e_ev = delta_e * HARTREE_TO_EV
                    
                    if delta_e < 0:
                        accept = True
                    else:
                        prob = np.exp(-beta * delta_e_ev)
                        accept = random.random() < prob
                    
                    if accept:
                        energy_changes.append(abs(delta_e))
                        current_structure = new_structure
                        current_energy = new_energy
                        n_accepted += 1
                        
                except Exception:
                    continue
            
            acceptance_rate = n_accepted / n_steps
            avg_energy_change = np.mean(energy_changes) if energy_changes else 0.0
            
            trial_acceptance_rates.append(acceptance_rate)
            trial_energy_changes.append(avg_energy_change)
            
            _log(f"  Trial {trial+1}/{n_trials}: {acceptance_rate*100:.1f}%")
        
        if trial_acceptance_rates:
            results[temp] = {
                'mean_acceptance': np.mean(trial_acceptance_rates),
                'std_acceptance': np.std(trial_acceptance_rates),
                'mean_energy_change': np.mean(trial_energy_changes),
                'std_energy_change': np.std(trial_energy_changes)
            }
            _log(f"  Average: {results[temp]['mean_acceptance']*100:.1f}% "
                 f"± {results[temp]['std_acceptance']*100:.1f}%")
        else:
            _log(f"  All trials failed for T = {temp} K", level="warning")
    
    # Plot
    if save_plot and results:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        temps = sorted(results.keys())
        mean_acc = [results[t]['mean_acceptance'] * 100 for t in temps]
        std_acc = [results[t]['std_acceptance'] * 100 for t in temps]
        mean_energy = [results[t]['mean_energy_change'] * HARTREE_TO_EV * 1000 for t in temps]
        std_energy = [results[t]['std_energy_change'] * HARTREE_TO_EV * 1000 for t in temps]
        
        ax1.errorbar(temps, mean_acc, yerr=std_acc, marker='o', capsize=5, 
                    linewidth=2, markersize=8, color='steelblue')
        ax1.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='50% target')
        ax1.axhline(y=30, color='orange', linestyle='--', alpha=0.5, label='30% minimum')
        ax1.set_xlabel('Temperature (K)', fontsize=12)
        ax1.set_ylabel('Acceptance Rate (%)', fontsize=12)
        ax1.set_title('BHMC Acceptance Rate vs Temperature', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        ax1.set_ylim(0, 100)
        
        ax2.errorbar(temps, mean_energy, yerr=std_energy, marker='s', capsize=5,
                    linewidth=2, markersize=8, color='darkorange')
        ax2.set_xlabel('Temperature (K)', fontsize=12)
        ax2.set_ylabel('Average Energy Change (meV)', fontsize=12)
        ax2.set_title('Average Accepted Energy Change', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        _log(f"Plot saved to {plot_filename}")
        plt.close()
    
    # Summary
    if logger:
        logger.section("Benchmark Summary")
    
    _log(f"{'Temp (K)':<12} {'Acceptance (%)':<22} {'Avg ΔE (meV)':<20}")
    _log(f"{'-'*54}")
    for temp in sorted(results.keys()):
        ma = results[temp]['mean_acceptance'] * 100
        sa = results[temp]['std_acceptance'] * 100
        me = results[temp]['mean_energy_change'] * HARTREE_TO_EV * 1000
        se = results[temp]['std_energy_change'] * HARTREE_TO_EV * 1000
        _log(f"{temp:<12.0f} {ma:>6.1f} ± {sa:<5.1f}        {me:>6.2f} ± {se:<6.2f}")
    
    if results:
        target = 0.45
        best_temp = min(sorted(results.keys()), 
                       key=lambda t: abs(results[t]['mean_acceptance'] - target))
        best_acc = results[best_temp]['mean_acceptance'] * 100
        _log(f"Recommended temperature: {best_temp} K (acceptance: {best_acc:.1f}%)")
    
    return results


