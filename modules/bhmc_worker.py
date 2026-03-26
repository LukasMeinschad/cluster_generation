import copy
import random
import numpy as np
from typing import List, Tuple, Dict, Optional

from molecule_class import Molecule
from operators import NonLocalOperators, LocalOperators
from box import SimulationBox
from bhmc_energy import EnergyEvaluator
from logger import Logger

# Constants
K_B = 8.617333262145e-5  # Boltzmann constant in eV/K
HARTREE_TO_EV = 27.2114  # 1 Hartree in eV

# ====================== Helper Classes =================================

class AdaptiveBoxController:
    """   
    Manages the worker-local adaptive scaling of the simulation box during Phase A.
    """
    def __init__(self, config_dict: dict, sim_box: Optional[SimulationBox] = None):
        self.sim_box = sim_box
        self.active = bool(config_dict.get('adaptive_box', False) and (sim_box is not None))
        self.interval = int(config_dict.get('box_update_interval', 10))
        self.target = float(config_dict.get('box_target_acceptance', 0.6))
        self.window = float(config_dict.get('box_acceptance_window', 0.05))
        self.kp = float(config_dict.get('box_growth_kp', 0.6))
        self.growth_max = float(config_dict.get('box_growth_max', 1.15))
        self.max_scale = float(config_dict.get('box_max_scale', 4.0))
        self.max_stable_windows = int(config_dict.get('box_stable_windows', 3))

        self.stable_windows = 0
        self.accepted_since_update = 0
        self.rejected_since_update = 0

        if self.active:
            if sim_box.box_type == "sphere":
                self.initial_size = float(sim_box.radius)
            elif sim_box.box_type == "cube":
                self.initial_size = float(np.max(sim_box.box_dimensions))
            else:
                self.initial_size = 1.0

        else:
            self.initial_size = 1.0  # Default to 1.0 if not active

    def get_current_size(self) -> Optional[float]:
        """  
        Returns the current box size (radius for sphere, max dimension for cube) or None if adaptive box is not active.
        """ 
        if not self.sim_box:
            return None 
        if self.sim_box.box_type == "sphere":
            return float(self.sim_box.radius)
        elif self.sim_box.box_type == "cube":
            return float(np.max(self.sim_box.box_dimensions))
        else:
            return None
        
    def record_step(self, accepted: bool):
        """ 
        Helper function to record the outcome of a step
        """
        if accepted:
            self.accepted_since_update += 1
        else:
            self.rejected_since_update += 1
    
    def update(self) -> Dict: 
        """  
        Updates the box size based on recent acceptance rates and returns update info
        """
        if not self.active:
            return {}
        
        total = self.accepted_since_update + self.rejected_since_update
        window_acc = self.accepted_since_update / total if total > 0 else 0.0

        lower = self.target - self.window
        upper = self.target + self.window

        grew = False 
        if window_acc < lower:
            # Too low acceptance, growth is needed
            error = lower - window_acc 
            growth = min(self.growth_max, 1.0 + self.kp * error)
            current_size = self.get_current_size() or 1.0
            current_scale = current_size / max(self.initial_size, 1e-12)

            if current_scale < self.max_scale:
                # Apply growth
                if self.sim_box.box_type == "sphere":
                    self.sim_box.radius *= growth
                    grew = True
                elif self.sim_box.box_type == "cube":
                    self.sim_box.box_dimensions = [d * growth for d in self.sim_box.box_dimensions]
                    grew = True
        
            if grew:
                self.stable_windows = 0 # Reset the stable window counter
            else:
                if lower <= window_acc <= upper:
                    self.stable_windows += 1
                    if self.stable_windows >= self.max_stable_windows:
                        # Stop growing if we've been in the target window for too long
                        self.active = False
                else:
                    # If acceptence is to high, no shrinking just reset the stable window counter
                    self.stable_windows = 0
            
        info = {
            "window_acceptance": window_acc,
            "grew": grew,
            "stable_windows": self.stable_windows,
            "current_size": self.get_current_size()
        }

        # Reset counters for next interval
        self.accepted_since_update = 0
        self.rejected_since_update = 0
        return info

class AnnealingScheduler:
    """  
    Manages temperature reduction (Simulated Annealing) for local refinement
    """
    def __init__(self, start_temp: float, end_temp: float, n_steps: int, p: float = 2.0):
        """  
        Args:
            start_temp: Initial temperature in Kelvin
            end_temp: Final temperature in Kelvin
            n_steps: Total number of steps over which to anneal
            p: Exponent for controlling the annealing curve (p=1 linear, p>1 faster decay)
        """
        self.start_temp = start_temp
        self.end_temp = end_temp
        self.n_steps = n_steps
        self.p = p

    def get_temperature(self, step: int) -> float:
        fraction = step / max(1, self.n_steps - 1)
        # Power law interpolation between start and end temperature
        return self.start_temp + (fraction ** self.p) * (self.end_temp - self.start_temp)
    
# ======================== Worker Logger =================================

def _get_worker_logger(log_file: Optional[str] = None, worker_id: int = 0) -> Optional[Logger]:
    
    if log_file is None:
        return None
    return Logger(name=f"bhmc_worker_{worker_id}", log_file=log_file,file_mode="a")

# ========================= Multiprocessing Workers ===============================

def _phase_a_worker(args: Tuple) -> Dict:
    """  
    Worker function for parallel BHMC Phase A exploration 
    """
    
    (initial_molecule, submolecule_indices, n_structures,
     operator_list, config_dict, worker_id, sim_box_dict, worker_log_file) = args
    
    logger = _get_worker_logger(worker_log_file, worker_id)

    def _log(msg: str, level: str = "info") -> None:
        if config_dict.get('verbose', False):
            print(msg)
        if logger:
            getattr(logger, level)(msg)
 
    sim_box = SimulationBox.from_dict(sim_box_dict) if sim_box_dict else None
    box_controller = AdaptiveBoxController(config_dict, sim_box)

    # For tracking and plotting
    box_updates = []

    if box_controller.active:
        box_updates.append({
            "step": 0, "box_size": box_controller.get_current_size(),
            "window_acceptance": None, "grew": False, "stable_windows": 0
        })


    # Operator Configuration
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
    
    operators, weights = [], []

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
        return {'accepted_structures': [], 'energy_trajectory': [], 'box_updates': box_updates, 'worker_id': worker_id}
    
    
    weights = np.array(weights) / sum(weights)  # Normalize weights


    
    evaluator = EnergyEvaluator(
        method=config_dict.get('method', 'hf'),
        basis=config_dict.get('basis', 'sto-3g'),
        backend=config_dict.get('backend', 'psi4'),
        xtb_method=config_dict.get('xtb_method', 'GFN2-xTB')
    )

    
    current_structure = copy.deepcopy(initial_molecule)
    current_energy, current_dipole_vec, current_dipole_mag = evaluator.evaluate(current_structure)

    if current_energy is None:
        _log(f"Worker {worker_id}: Failed to evaluate initial energy!", level="error")
        return {'accepted_structures': [], 'energy_trajectory': [], 'box_updates': box_updates, 'worker_id': worker_id}
    
    _log(f"Worker {worker_id}: Starting energy = {current_energy:.6f} Hartree")

    accepted_structures = []
    energy_trajectory = [(0, current_energy)]
    n_accepted = n_rejected = 0
    temperature = config_dict.get('temperature', 400.0)
    beta = 1.0 / (K_B * temperature)

    for step in range(n_structures):
        op_info = np.random.choice(operators, p=weights)


        # Try to apply operator
        try:
            if op_info['supports_adaptive']:
                new_structure = op_info["func"](current_structure,submolecule_indices=submolecule_indices, adaptive=True)
            else:
                new_structure = op_info["func"](current_structure,submolecule_indices=submolecule_indices)

            if np.allclose(new_structure.coordinates, current_structure.coordinates):
                _log(f"Worker {worker_id}: Operator {op_info['name']} did not change the structure. Skipping.", level="warning")
                n_rejected += 1; box_controller.record_step(False)
                energy_trajectory.append((step+1, current_energy))
                continue
            
            new_energy, new_dipole_vec, new_dipole_mag = evaluator.evaluate(new_structure)
            
            if new_energy is None:
                _log(f"Worker {worker_id}: Energy evaluation failed for proposed structure. Rejecting.", level="warning")
                n_rejected += 1; box_controller.record_step(False)
                energy_trajectory.append((step+1, current_energy))
                continue
            
            # Metropolis Criterion
            accept = False
            if new_energy < current_energy:
                accept = True
            else:
                delta_e_ev = (new_energy - current_energy) * HARTREE_TO_EV
                prob = np.exp(-beta * delta_e_ev)
                rand_val = random.random()
                accept = rand_val < prob
                if (step +1) % 10 == 0:
                    _log(f"Worker {worker_id}: Step {step+1}, ΔE = {delta_e_ev:.4f} eV, Prob = {prob:.4f}, Rand = {rand_val:.4f}")
            
            if accept:
                # Set new values
                current_structure, current_energy = copy.deepcopy(new_structure), new_energy
                current_dipole_vec, current_dipole_mag = new_dipole_vec, new_dipole_mag
                accepted_structures.append((copy.deepcopy(current_structure), current_energy, current_dipole_vec, current_dipole_mag))
                n_accepted += 1
                _log(f"Worker {worker_id}: Accepted new structure at step {step+1} with energy {current_energy:.6f} Hartree")
            else:
                n_rejected += 1
            
            box_controller.record_step(accept)
            energy_trajectory.append((step+1, current_energy))

        except Exception as e:
            _log(f"Worker {worker_id}: Exception during operator application or evaluation: {e}", level="error")
            n_rejected += 1; box_controller.record_step(False)
            energy_trajectory.append((step+1, current_energy))
            continue

        # Adaptive Box Update
        if box_controller.active and (step + 1) % box_controller.interval == 0:
            update_info = box_controller.update()

            if update_info['grew']:
                _log(f"Worker {worker_id}: Box grew at step {step+1}. New size: {update_info['current_size']:.3f}, Acceptance in window: {update_info['window_acceptance']:.3f}")
            elif update_info['stable_windows'] > 0:
                _log(f"Worker {worker_id}: Box stable for {update_info['stable_windows']} windows at step {step+1}. Current size: {update_info['current_size']:.3f}, Acceptance in window: {update_info['window_acceptance']:.3f}")

            box_updates.append({
                "step": step + 1, "box_size": update_info["current_size"],
                "window_acceptance": update_info.get("window_acceptance"),
                "grew": update_info.get("grew", False),
                "stable_windows": update_info.get("stable_windows", 0)
            })

            if not box_controller.active:
                _log(f"Worker {worker_id}: Box adaptation stopped after step {step+1}.")
            
    mode_rate = n_accepted / (n_accepted + n_rejected) * 100 if (n_accepted + n_rejected) > 0 else 0.0
    _log(f"Worker {worker_id}: Finished Phase A. Accepted {n_accepted} structures, Rejected {n_rejected} structures. Acceptance Rate: {mode_rate:.2f}%")

    return {
        'accepted_structures': accepted_structures,
        'energy_trajectory': energy_trajectory,
        'box_updates': box_updates,
        'worker_id': worker_id
    }


def _phase_b_worker(args: Tuple) -> Dict:
    """Worker function for parallel BHMC Phase B local refinement."""
    (initial_molecule, initial_energy, submolecule_indices, n_steps, 
     operator_list, config_dict, worker_id, sim_box_dict, worker_log_file) = args
    
    logger = _get_worker_logger(worker_log_file, worker_id)

    def _log(msg: str, level: str = "info") -> None:
        if config_dict.get('verbose', False):
            print(msg)
        if logger:
            getattr(logger, level)(msg)

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

    operators, weights = [], []
    for op_name, weight in operator_list:
        if op_name in op_map:
            operators.append({"name": op_name, "func": op_map[op_name]})
            weights.append(weight)

    weights = np.array(weights) / sum(weights)
    
    evaluator = EnergyEvaluator(
        method=config_dict.get('method', 'hf'),
        basis=config_dict.get('basis', 'sto-3g'),
        backend=config_dict.get('backend', 'psi4'),
        xtb_method=config_dict.get('xtb_method', 'GFN2-xTB')
    )

    current_structure = copy.deepcopy(initial_molecule)
    current_energy, current_dipole_vec, current_dipole_mag = evaluator.evaluate(current_structure)
    
    if current_energy is None:
        return {"best_structure": None, "accepted_structures": [], "energy_trajectory": [], "worker_id": worker_id}

    best_energy = current_energy
    best_structure = copy.deepcopy(current_structure)
    best_dipole_vec, best_dipole_mag = current_dipole_vec, current_dipole_mag
    
    accepted_structures = []
    energy_trajectory = [(0, current_energy)]
    n_accepted = n_rejected = 0

    # Temperature Control
    scheduler = AnnealingScheduler(
        start_temp=config_dict['temperature'], 
        end_temp=config_dict.get('max_temperature', config_dict['temperature']), 
        n_steps=n_steps
    )

    for step in range(n_steps):
        op_info = np.random.choice(operators, p=weights)
        current_temp = scheduler.get_temperature(step)
        beta = 1.0 / (K_B * current_temp)

        try:
            new_structure = op_info['func'](current_structure, submolecule_indices, adaptive=adaptive)
            if np.allclose(new_structure.coordinates, current_structure.coordinates):
                raise ValueError("No coordinate change")
                
            new_energy, new_dipole_vec, new_dipole_mag = evaluator.evaluate(new_structure)
            if new_energy is None:
                raise ValueError("Energy invalid")

            # Metropolis-Hastings Acceptance
            accept = False
            if new_energy < current_energy:
                accept = True
            else:
                delta_e_ev = (new_energy - current_energy) * HARTREE_TO_EV
                if random.random() < np.exp(-beta * delta_e_ev):
                    accept = True

            if accept:
                current_structure, current_energy = copy.deepcopy(new_structure), new_energy
                accepted_structures.append((copy.deepcopy(current_structure), current_energy, new_dipole_vec, new_dipole_mag))
                n_accepted += 1
                
                # Check for absolute best
                if current_energy < best_energy:
                    best_energy = current_energy
                    best_structure = copy.deepcopy(current_structure)
                    best_dipole_vec, best_dipole_mag = new_dipole_vec, new_dipole_mag
            else:
                n_rejected += 1
                
            energy_trajectory.append((step+1, current_energy))

        except Exception as e:
            n_rejected += 1
            energy_trajectory.append((step+1, current_energy))

    return {
        "best_structure": (best_structure, best_energy, best_dipole_vec, best_dipole_mag),
        "accepted_structures": accepted_structures,
        "energy_trajectory": energy_trajectory,
        "worker_id": worker_id
    }
