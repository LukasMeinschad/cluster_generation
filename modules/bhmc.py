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
from transformations import NonLocalOperators, LocalOperators
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


    # In-worker adaptive box control (Phase) A
    adaptive_box: bool = True
    box_update_interval: int = 10  # Update every N steps
    box_target_acceptance: float = 0.6 # Desired acceptance rate 
    box_acceptance_window: float = 0.05 # no update of box rate inside [box_target_acceptance - window, box_target_acceptance + window]
    box_growth_kp: float = 0.6 # propotional gain
    box_growth_max: float = 1.15 # cap per update
    box_max_scale: float = 4.0 # Relative to initial box size
    box_stable_windows: int = 3 # stop growing after this many stable windows

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
    
    
    # Adaptive box_settings (worker-local)
    adaptive_box = bool(config_dict.get('adaptive_box', False)) and (sim_box is not None)
    box_update_interval = int(config_dict.get('box_update_interval', 10))
    box_target = float(config_dict.get('box_target_acceptance', 0.6))
    box_acceptance_window = float(config_dict.get('box_acceptance_window', 0.05))
    box_kp = float(config_dict.get('box_growth_kp', 0.6))  
    box_growth_max = float(config_dict.get('box_growth_max', 1.15)) 
    box_max_scale = float(config_dict.get('box_max_scale', 4.0))
    box_stable_windows = int(config_dict.get('box_stable_windows', 3))
    
    if adaptive_box:
        if sim_box.box_type == "sphere":
            initial_box_size = sim_box.radius
        elif sim_box.box_type == "cube":
            initial_box_size = float(np.max(sim_box.box_dimensions))
        else:
            initial_box_size = 1.0
    else:
        initial_box_size = 1.0
    
    # Iterators for adaptive box_control
    stable_windows = 0
    last_update_step = 0
    accepted_since_update = 0
    rejected_since_update = 0

    # Operator setup
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
        'principal_axis_rotation': nonlocal_ops.principal_axis_rotation_operator,
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
                accepted_since_update += 1
                _log(f"Worker {worker_id} step {step}: ACCEPTED {operator['name']} E={current_energy:.6f} Ha dipole={current_dipole_mag:.4f} D")
            else:
                n_rejected += 1
                rejected_since_update += 1
            


            # Record energy trajectory for every step, accepted or not
            energy_trajectory.append((step+1, current_energy))
                
        except Exception as e:
            _log(f"Worker {worker_id} step {step}: Exception in {operator['name']}: {e}", level="error")
            n_rejected += 1
            energy_trajectory.append((step+1, current_energy))
            continue

        # Worker-local adaptive box control
        if adaptive_box and (step + 1) % box_update_interval == 0:
            window_total = accepted_since_update + rejected_since_update
            if window_total > 0:
                window_acc = accepted_since_update / window_total
            else:
                window_acc = 0.0
            
            lower = box_target - box_acceptance_window
            upper = box_target + box_acceptance_window
            
            grew = False
            if window_acc < lower:
                # proportional growth with cap
                error = lower - window_acc
                growth = min(box_growth_max, 1.0 + box_kp * error)

                if sim_box.box_type == "sphere":
                    current_scale = sim_box.radius / max(initial_box_size, 1e-12)
                    if current_scale < box_max_scale:
                        new_radius = min(sim_box.radius * growth, initial_box_size * box_max_scale)
                        if new_radius > sim_box.radius:
                            sim_box.radius = new_radius
                            grew = True
                elif sim_box.box_type == "cube":
                    current_scale = max(sim_box.box_dimensions) / max(initial_box_size, 1e-12)
                    if current_scale < box_max_scale:
                        new_dimensions = np.minimum(sim_box.box_dimensions * growth, initial_box_size * box_max_scale)
                        if np.any(new_dimensions > sim_box.box_dimensions):
                            sim_box.box_dimensions = new_dimensions
                            grew = True
            if grew:
                stable_windows = 0
                if sim_box.box_type == "sphere":
                    _log(f"Worker {worker_id} step {step}: Adaptive box growth — new radius = {sim_box.radius:.2f} Å (acc={window_acc:.2%}) ")
                else:
                    _log(f"Worker {worker_id} step {step}: Adaptive box growth — new dimensions = {sim_box.box_dimensions} Å (acc={window_acc:.2%}) ")
            else:
                if lower <= window_acc <= upper:
                    stable_windows += 1
                    _log(f"Worker {worker_id} step {step}: Box stable window {stable_windows}/{box_stable_windows} (acc={window_acc:.2%})")
                else:
                    stable_windows = 0

            # Reset the window counters
            accepted_since_update = 0
            rejected_since_update = 0
            last_update_step = step + 1

            if stable_windows >= box_stable_windows:
                _log(f"Worker {worker_id} step {step}: Box growth stabilized after {stable_windows} windows. Stopping growth.", level="info")
                adaptive_box = False  # Stop further growth
        
        
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


def _phase_b_worker(args: Tuple) -> Dict:
    """ 
    Worker function for parallel BHMC Phase B local refinement.

    Each worker starts from a representative structure found in Phase A and applies
    local operators to refine within that basin of attraction

    Args:
        args: Tuple containing (initial_molecule, initial_energy, submolecule_indices,
        n_steps, operator_list, config_dict, worker_id, sim_box_dict, worker_log_file)
    Returns:
       Dictionary with keys:
                'best_strucutre': (Molecule, energy, dipole_vec, dipole_magnitude) tuple of the best structure found
                'accepted_structures': List of (Molecule, energy, dipole_vec, dipole_mag) tuples
                'energy_trajectory': List of (step, current_energy) for every step
                'worker_id': Unique identifier for the worker
    """
    (initial_molecule, initial_energy, submolecule_indices,
     n_steps, operator_list, config_dict, worker_id, sim_box_dict, worker_log_file) = args
    
    logger = _get_worker_logger(worker_log_file, worker_id)

    def _log(msg: str, level: str = "info") -> None:
        if config_dict.get('verbose', False):
            print(msg)
        if logger:
            getattr(logger, level)(msg)

    # Constants
    k_B = 8.617333262145e-5  # Boltzmann constant in eV/K
    HARTREE_TO_EV = 27.2114

    # Setup
    sim_box = SimulationBox.from_dict(sim_box_dict) if sim_box_dict else None
    local_ops = LocalOperators(simulation_box=sim_box)
    adaptive = config_dict.get('adaptive_operators', True)


    op_map = {
        "random_displacement": local_ops.random_displacement_submolecule,
        "random_rotation": local_ops.random_rotation_submolecule,
        "local_twist": local_ops.local_twist_submolecule,
        "correlated_displacement": local_ops.correlated_displacement,
        "small_principal_axis_rotation": local_ops.small_principal_axis_rotation
    }

    operators = []
    weights = []
    for op_name, weight in operator_list:
        if op_name in op_map:
            operators.append({"name": op_name, "func": op_map[op_name]})
            weights.append(weight)

    if not operators:
        _log(f"Worker {worker_id}: No valid local operator found!", level="error")
        return {
            "best_structure": None,
            "accepted_structures": [],
            "energy_trajectory": [],
            "worker_id": worker_id
        }
    weights = np.array(weights)
    weights /= weights.sum()
    evaluator = EnergyEvaluator(
        method=config_dict['method'],
        basis=config_dict['basis']
    )

    current_structure = copy.deepcopy(initial_molecule)

    # TODO this can be improved
    eval_result = evaluator.evaluate(current_structure)
    if eval_result is None or eval_result[0] is None:
        _log(f"Worker {worker_id}: Failed to evaluate initial energy!", level="error")
        return {
            "best_structure": None,
            "accepted_structures": [],
            "energy_trajectory": [],
            "worker_id": worker_id
        }
    current_energy, current_dipole_vec, current_dipole_mag = eval_result
    _log(f"Worker {worker_id}: Starting energy = {current_energy:.6f} Hartree")

    # Track best structure found
    best_energy = current_energy
    best_structure = copy.deepcopy(current_structure)
    best_dipole_vec = current_dipole_vec
    best_dipole_mag = current_dipole_mag
    
    accepted_structures = []
    energy_trajectory = [(0, current_energy)] # Track energy at each step
    n_accepted = 0
    n_rejected = 0
    temperature = config_dict['temperature']
    beta = 1.0 / (k_B * temperature)

    for step in range(n_steps):
        op_idx = np.random.choice(len(operators), p=weights)
        operator = operators[op_idx]

        try:
            new_structure = operator['func'](
                current_structure, submolecule_indices, adaptive=adaptive
            )
            if np.allclose(new_structure.coordinates, current_structure.coordinates):
                _log(f"Worker {worker_id} step {step}: {operator['name']} — no coordinate change (box rejection?)", level="warning")
                n_rejected += 1
                energy_trajectory.append((step+1, current_energy))
                continue
            new_eval = evaluator.evaluate(new_structure)
            if new_eval is None or new_eval[0] is None:
                _log(f"Worker {worker_id} step {step}: Energy eval failed after {operator['name']}", level="warning")
                n_rejected += 1
                energy_trajectory.append((step+1, current_energy))
                continue
            new_energy, new_dipole_vec, new_dipole_mag = new_eval

            # Metropolis acceptance
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
                
                if current_energy < best_energy:
                    best_energy = current_energy
                    best_structure = copy.deepcopy(current_structure)
                    best_dipole_vec = current_dipole_vec
                    best_dipole_mag = current_dipole_mag
            else:
                n_rejected += 1
            energy_trajectory.append((step+1, current_energy))

        except Exception as e:
            _log(f"Worker {worker_id} step {step}: Exception in {operator['name']}: {e}", level="error")
            n_rejected += 1
            energy_trajectory.append((step+1, current_energy))
            continue

        if (step + 1) % 10 == 0:
            rate = n_accepted / (n_accepted + n_rejected) * 100 if (n_accepted + n_rejected) > 0 else 0.0
            _log(f"Worker {worker_id}: {step+1}/{n_steps} — accepted={n_accepted}, rate={rate:.1f}%")

    total = n_accepted + n_rejected
    rate = n_accepted / total * 100 if total > 0 else 0.0
    _log(f"Worker {worker_id}: DONE — {n_accepted}/{total} accepted ({rate:.1f}%)")
    
    return {
        "best_structure": (best_structure, best_energy, best_dipole_vec, best_dipole_mag),
        "accepted_structures": accepted_structures,
        "energy_trajectory": energy_trajectory,
        "worker_id": worker_id
    }



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
        ('principal_axis_rotation', 0.8),
    ]

    DEFAULT_LOCAL_OPERATORS = [
        ("random_displacement", 1.8),
        ("random_rotation", 1.0),
        ("local_twist", 0.8),
        ("correlated_displacement", 0.8),
        ("small_principal_axis_rotation", 0.8)   
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

            # Adaptive Box scaling
            'adaptive_box': self.config.adaptive_box,
            'box_update_interval': self.config.box_update_interval,
            'box_target_acceptance': self.config.box_target_acceptance,
            'box_acceptance_window': self.config.box_acceptance_window,
            'box_growth_kp': self.config.box_growth_kp,
            'box_growth_max': self.config.box_growth_max,
            'box_max_scale': self.config.box_max_scale,
            'box_stable_windows': self.config.box_stable_windows
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
    
    def run_phase_b(self,
                    representatives: List["StructureData"],
                    submolecule_indices: List[List[int]],
                    n_steps_per_worker = 200,
                    local_operators: Optional[List[Tuple[str, float]]] = None,
                    n_processes: Optional[int] = None
                    ) -> List[Dict]:
        """ 
        Run Phase B: Local Refinement around representative structures from Phase A

        Each representative is assigned to one worker that applies local operators
        to refine the structure within one basin of attraction.

        Args:
            representatives: List of representative structures from Phase A
            submolecule_indices: Submolecule indices for local operators
            n_steps_per_worker: Number of BHMC steps for each local refinement worker
            local_operators: Optional list of (operator_name, weight) tuples for local refinement.
                            If None, defaults to MultiPhaseBHMC.DEFAULT_LOCAL_OPERATORS
            n_processes: Number of parallel processes to use. If None, defaults to len(representatives)

        Returns:
            List of worker results dicts with keys:
                'best_structure': (Molecule, energy, dipole_vec, dipole_mag) tuple of the best structure found
                'accepted_structures': List of (Molecule, energy, dipole_vec, dipole_mag) tuples
                'energy_trajectory': List of (step, current_energy) for every step
        """
        if local_operators is None:
            local_operators = self.DEFAULT_LOCAL_OPERATORS
        
        n_workers = len(representatives)
        if n_processes is None:
            n_processes = n_workers
        if self.logger:
            self.logger.header("Phase B: Local Refinement")
        
        self._log(f"Configuration:")
        self._log(f" Representatives: {n_workers}")
        self._log(f" Steps per worker: {n_steps_per_worker}")
        self._log(f" Total Steps: {n_workers * n_steps_per_worker}")
        self._log(f" Method/Basis: {self.config.method}/{self.config.basis}")
        self._log(f" Temperature: {self.config.temperature} K")
        self._log(f" Adaptive operators: {self.config.adaptive_operators}")
        self._log(f" Local Operators: {', '.join(op[0] for op in local_operators)}")
        if self.simulation_box:
            self._log(f" Simulation Box: {self.simulation_box}")
        if self.worker_log_file:
            self._log(f" Worker log: {self.worker_log_file}")
        config_dict = {
            "method": self.config.method,
            "basis": self.config.basis,
            "temperature": self.config.temperature,
            "verbose": self.config.verbose,
            "adaptive_operators": self.config.adaptive_operators,
        }
        
        sim_box_dict = self.simulation_box.to_dict() if self.simulation_box else None

        args_list = [
            (
                rep.molecule,
                rep.energy,
                submolecule_indices,
                n_steps_per_worker,
                local_operators,
                config_dict,
                worker_id,
                sim_box_dict,
                self.worker_log_file
            )
            for worker_id, rep in enumerate(representatives)
        ]
        self._log("Starting parallel local refinement workers...")
        with Pool(processes=n_processes) as pool:
            results = pool.map(_phase_b_worker, args_list)

        # Collect results and trajectories
        self.phase_b_trajectories = {}
        self.phase_b_results = results

        for worker_result in results:
            wid = worker_result["worker_id"]
            accepted = worker_result["accepted_structures"]
            trajectory = worker_result["energy_trajectory"]
            best = worker_result["best_structure"]

            self.phase_b_trajectories[wid] = trajectory
            best_e_str = f"{best[1]:.6f} Ha" if best else "N/A"
            self._log(f" Worker {wid}: Best energy = {best_e_str}, accepted structures = {len(accepted)}, trajectory points = {len(trajectory)}")

        if self.logger:
            self.logger.section("Phase B Results")
        
        best_structures = [r["best_structure"] for r in results if r["best_structure"] is not None]
        total_accepted = sum(len(r["accepted_structures"]) for r in results)
        total_steps = n_workers * n_steps_per_worker


        self._log(f" Total steps: {total_steps}")
        self._log(f" Total accepted: {total_accepted}")
        if total_steps > 0:
            self._log(f" Overall acceptance rate: {total_accepted/total_steps*100:.2f}%")
        
        if best_structures:
            energies = [bs[1] for bs in best_structures]
            self._log(f"Best energies across workers (Hartree):")
            self._log(f" Min:  {min(energies):.6f}")
            self._log(f" Max:  {max(energies):.6f}")
            self._log(f" Mean: {np.mean(energies):.6f}")
            self._log(f" Std:  {np.std(energies):.6f}")


        self.plot_worker_trajectories(
            save_path="figures/phase_b_worker_energy_trajectories.png",
            trajectories=self.phase_b_trajectories,
            title="Phase B Worker Energy Trajectories"
        )

        return results

    def analyse_phase_a_results(
            self,
            phase_a_structures: List[Tuple[Molecule, float]], 
            submolecule_indices: Optional[List[List[int]]] = None,
            cluster_method: str = "kmeans",
            n_clusters = None,
            simulation_box: Optional[SimulationBox] = None,
            logger: Optional[Logger] = None,
            **cluster_kwargs
            ) -> List[Molecule]:
        """Analyze Phase A results and extract representative structures.
        
        Args:
            phase_a_structures: List of (Molecule, energy) tuples from Phase A
            submolecule_indices: Optional submolecule indices for clustering
            cluster_method: Clustering method to use ("kmeans", "dbscan", etc.)
            n_clusters: Number of clusters to identify
            simulation_box: Optional simulation box for trajectory plots
            cluster_kwargs: Additional kwargs for clustering method (e.g. eps for DBSCAN) 
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
        #analyzer.plot_rmsd_heatmap()
        
        # Precompute the feature matrix
        analyzer.feature_matrix(normalize=True)  # Normalized features for clustering


        analyzer.plot_energy_distribution()
        analyzer.plot_energy_vs_rmsd()
        analyzer.plot_rg_vs_energy()
        analyzer.plot_int_d_vs_e()
        analyzer.plot_pca_explained_variance()
         
        

        if cluster_method in {"kmeans", "agglomerative"}:
            if n_clusters is None:
                raise ValueError(f"{cluster_method} clustering requires n_clusters to be specified")
            labels = analyzer.cluster(
                method=cluster_method,
                n_clusters=n_clusters,
                **cluster_kwargs
            )
        else:
            labels = analyzer.cluster(
                method=cluster_method,
                **cluster_kwargs,
            )




        analyzer.plot_pca_clustered(n_clusters=n_clusters,
                                    cluster_method=cluster_method,**cluster_kwargs)
        analyzer.plot_tsne_clustered(n_clusters=n_clusters,
                                     cluster_method=cluster_method,**cluster_kwargs)
        analyzer.plot_umap_clustered(n_clusters=n_clusters,
                                     cluster_method=cluster_method,**cluster_kwargs)

        representatives = analyzer.get_cluster_representatives()
        analyzer.log_representative_features(representatives,labels=labels) 
    
        
        return representatives
    
    # Additional plotting utilities for worker trajectories
    def plot_worker_trajectories(self,
                                 save_path: str = "figures/worker_energy_trajectories.png",
                                 trajectories: Optional[Dict] = None,
                                 title: str = "Worker Energy Trajectories"):
        """  
        Plot energy vs step for each BHMC worker
        """
        import matplotlib.pyplot as plt
        from pathlib import Path

        if trajectories is None:
            trajectories = self.worker_trajectories

        if not trajectories:
            self._log("No worker trajectories to plot.", level="warning")
            return
        
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        n_workers = len(trajectories)
        cmap = plt.get_cmap('tab20', max(n_workers, 1))

        plt.figure(figsize=(12, 6))
        for worker_id, trajectory in trajectories.items():
            steps, energies = zip(*trajectory)
            plt.plot(steps, energies, label=f"Worker {worker_id}", color=cmap(worker_id))

        plt.xlabel("Step")
        plt.ylabel("Energy (Hartree)")
        plt.title(title)
        plt.legend(
            ncols = 4,
            loc = "upper right",
            columnspacing = 0.5,
        )
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()

    # ========================= Interpolation ====

    def generate_interpolated_candidates(
            self,
            representatives: List["StructureData"],
            submolecule_indices: Optional[List[List[int]]] = None,
            n_interpolations: int = 2,
            lambdas: Optional[List[float]] = None
        ) -> List["StructureData"]:
        """
        Generate interpolated candidate structures between all pairs of representative structures.
        This uses the Kabsch algorithm
        
        Args:
            representatives: List of representative structures to interpolate between
            n_interpolations: Number of interpolated structures to generate between each pair
            lambdas: Optional list of interpolation parameters (0 to 1). If None, uses equally spaced.
        
        Returns:
        List of StructureData objects for the interpolated candidates, with energy and dipole set to None.
        """
        from transformations import Transformation
        from cluster import StructureData

        transformer = Transformation()
        interpolated_candidates = []
        n_reps = len(representatives)

        self._log(f"Generating interpolated candidates between {n_reps} representatives")
        self._log(f" Pairs: {n_reps * (n_reps - 1) // 2} ")
        self._log(f" Interpolations per pair: {n_interpolations}")
        for i in range(n_reps):
            for j in range(i + 1, n_reps):
                mol_a = representatives[i].molecule
                mol_b = representatives[j].molecule

                interpolated_mols = transformer.kabsch_interpolate(
                    mol_a, mol_b,
                    submolecule_indices=submolecule_indices,
                    lambdas = lambdas,
                    n_points= n_interpolations
                )
                for mol in interpolated_mols:
                    candidate = StructureData(
                        molecule = mol,
                        energy = 0.0,
                        phase= "interpolated",
                        metadata={"parent i": i, "parent j": j}
                    )
                    interpolated_candidates.append(candidate)
        self._log(f"Generated {len(interpolated_candidates)} interpolated candidates")
        return interpolated_candidates



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


if __name__ == "__main__":
    
    # Module Import just for testing
    import multiprocessing as mp
    import sys
    from pathlib import Path 
    import matplotlib.pyplot as plt
    import numpy as np 

    module_dir = Path(__file__).parent / "modules"
    sys.path.insert(0, str(module_dir))

    from molecule_class import Molecule
    from cluster import StructureData
    from transformations import LocalOperators
    from box import SimulationBox
    from bhmc import MultiPhaseBHMC, BHMCConfig
    from logger import Logger 



    mp.set_start_method('spawn', force=True)
    

    xyz = """
    6
    H2O Nonopt
    O     -4.651077    1.267697   -0.000000
    H     -5.057317    1.928988   -0.547872
    H     -4.969101    1.355490    0.890872
    O     -6.998837    4.503183    0.000001
    H     -7.632732    4.320290   -0.683540
    H     -7.009572    5.431464    0.201689
    """
    molecule = Molecule.from_xyz(xyz)
    submolecule_indices = [[0, 1, 2], [3, 4, 5]]
    covalent_radii = [r for r in molecule.covalent_radii]
    
    # Center molecule
    centroid = np.mean(molecule.coordinates, axis=0)
    molecule.coordinates -= centroid

    
    
    sim_box = SimulationBox.from_covalent_radii(
        covalent_radii=covalent_radii,
        n_atoms=len(molecule.atom_labels),
        box_type="sphere",
        scale_factor=5.0
    )

    # Benchmark Displacements of local operators
    local_ops = LocalOperators(simulation_box=sim_box)
    diag = local_ops.diagnose_perturbations(
        molecule, submolecule_indices, n_samples=200, adaptive=True
    )
    diag_noadapt = local_ops.diagnose_perturbations(
        molecule,submolecule_indices, n_samples=200, adaptive=False
    )


    # We need to wrap the molecule as structure data for phase b
    representative = StructureData(
        molecule=molecule,
        energy=0.0,
        phase="test"
    )
    config = BHMCConfig(
        temperature=150,
        method="hf",
        basis="6-311g(d,p)",
        adaptive_operators=False
    )

    sampler = MultiPhaseBHMC(
        config=config,
        simulation_box=sim_box,
    )
    results = sampler.run_phase_b(
        representatives=[representative],
        submolecule_indices=submolecule_indices,
        n_steps_per_worker=300,
        n_processes=1
    )
    worker_result = results[0]
    trajectory = worker_result["energy_trajectory"]
    accepted = worker_result["accepted_structures"]
    best = worker_result["best_structure"]
    steps, energies = zip(*trajectory)
    energies_ha = np.array(energies)

    print(f"Best energy: {best[1]:.6f} Ha")
    print(f"Accepted structures: {len(accepted)}")
    print(f"Energy trajectory points: {len(trajectory)}")
    print(f"Acceptance rate: {len(accepted)/len(trajectory)*100:.1f}%")
    print(f"Initial energy: {trajectory[0][1]:.6f} Ha")

    print(f"Energy change: {(best[1] - trajectory[0][1])} Ha")


    # Plot convergence of energy
    plt.figure(figsize=(8, 5))
    plt.plot(steps, energies_ha, marker='o', color='steelblue')
    plt.xlabel("Step")
    plt.ylabel("Energy (Hartree)")
    plt.title("Phase B Energy Trajectory")
    plt.grid(True)
    plt.tight_layout()
    plt.show() 
    # Save the energy trajectory
    # Save the best structure as XYZ
    best_structure = best[0]
    best_energy = best[1]
    best_dipole_vec = best[2]
    best_dipole_mag = best[3]
    best_structure.write_xyz("best_structure_phase_b.xyz")


