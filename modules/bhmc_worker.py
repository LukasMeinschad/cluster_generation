import copy
import random
import numpy as np
from typing import List, Tuple, Dict, Optional

import cProfile
import pstats
import io
import tracemalloc

from modules.molecule_class import Molecule
from modules.operators import NonLocalOperators, LocalOperators
from modules.box import SimulationBox
from modules.calculator import EnergyEvaluator
from modules.logger import Logger

K_B = 8.617333262145e-5   # eV / K
HARTREE_TO_EV = 27.2114


# ======================== Adaptive Box Controller ========================

class AdaptiveBoxController:
    """Manages worker-local adaptive scaling of the simulation box."""

    def __init__(self, config_dict: dict, sim_box: Optional[SimulationBox] = None):
        self.sim_box = sim_box
        self.active = bool(config_dict.get('adaptive_box', False) and sim_box is not None)
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
            self.initial_size = 1.0

    def get_current_size(self) -> Optional[float]:
        if not self.sim_box:
            return None
        if self.sim_box.box_type == "sphere":
            return float(self.sim_box.radius)
        elif self.sim_box.box_type == "cube":
            return float(np.max(self.sim_box.box_dimensions))
        return None

    def record_step(self, accepted: bool):
        if accepted:
            self.accepted_since_update += 1
        else:
            self.rejected_since_update += 1

    def update(self) -> Dict:
        if not self.active:
            return {}

        total = self.accepted_since_update + self.rejected_since_update
        window_acc = self.accepted_since_update / total if total > 0 else 0.0

        lower = self.target - self.window
        upper = self.target + self.window
        grew = False

        if window_acc < lower:
            error = lower - window_acc
            growth = min(self.growth_max, 1.0 + self.kp * error)
            current_size = self.get_current_size() or 1.0
            current_scale = current_size / max(self.initial_size, 1e-12)

            if current_scale < self.max_scale:
                if self.sim_box.box_type == "sphere":
                    self.sim_box.radius *= growth
                    grew = True
                elif self.sim_box.box_type == "cube":
                    self.sim_box.box_dimensions = [d * growth for d in self.sim_box.box_dimensions]
                    grew = True

        if grew:
            self.stable_windows = 0
        elif lower <= window_acc <= upper:
            self.stable_windows += 1
            if self.stable_windows >= self.max_stable_windows:
                self.active = False
        else:
            self.stable_windows = 0

        info = {
            "window_acceptance": window_acc,
            "grew": grew,
            "stable_windows": self.stable_windows,
            "current_size": self.get_current_size(),
        }
        self.accepted_since_update = 0
        self.rejected_since_update = 0
        return info


# ======================== Worker Logger ==================================

def _get_worker_logger(log_file: Optional[str], worker_id: int) -> Optional[Logger]:
    if log_file is None:
        return None
    return Logger(name=f"bhmc_worker_{worker_id}", log_file=log_file, file_mode="a")


# ======================== Unified BHMC Worker ============================

def _bhmc_worker(args: Tuple) -> Dict:
    """
    Single BHMC worker that draws from a unified local + non-local operator pool.
    Non-local operators provide basin hops; local operators explore within a basin.
    The balance between the two is controlled by the weights in operator_list.
    """
    import os
    import tempfile
    from pathlib import Path

    pr = cProfile.Profile()
    pr.enable()
    tracemalloc.start()

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    (initial_molecule, submolecule_indices, n_steps,
     operator_list, config_dict, worker_id, sim_box_dict, worker_log_file) = args

    # Unique psi4 scratch per worker
    base = Path.cwd() / ".psi_scratch"
    base.mkdir(parents=True, exist_ok=True)
    task_scratch = Path(tempfile.mkdtemp(prefix=f"worker_{worker_id}_", dir=base))
    os.environ["PSI_SCRATCH"] = str(task_scratch)

    result = {
        'accepted_structures': [],
        'energy_trajectory': [],
        'box_updates': [],
        'worker_id': worker_id,
        'operator_acceptance_rates': {},
    }

    try:
        logger = _get_worker_logger(worker_log_file, worker_id)

        def _log(msg: str, level: str = "info") -> None:
            if config_dict.get('verbose', False):
                print(msg)
            if logger:
                getattr(logger, level)(msg)

        # ---- Simulation box & adaptive controller ----
        sim_box = SimulationBox.from_dict(sim_box_dict) if sim_box_dict else None
        box_controller = AdaptiveBoxController(config_dict, sim_box)
        box_updates = []
        if box_controller.active:
            box_updates.append({
                "step": 0, "box_size": box_controller.get_current_size(),
                "window_acceptance": None, "grew": False, "stable_windows": 0,
            })

        # ---- Build unified operator pool ----
        nonlocal_ops = NonLocalOperators(simulation_box=sim_box)
        local_ops = LocalOperators(simulation_box=sim_box)
        adaptive = config_dict.get('adaptive_operators', True)

        nonlocal_adaptive = {'twist', 'large_displacement', 'roto_reflection'}

        nonlocal_op_map = {
            'twist':                   nonlocal_ops.twist_operator,
            'large_displacement':      nonlocal_ops.large_displacement,
            'mirror':                  nonlocal_ops.mirror_operator,
            'roto_reflection':         nonlocal_ops.roto_reflection_operator,
            'exchange':                nonlocal_ops.exchange_operator,
            'random_so3':              nonlocal_ops.random_so3_operator,
            'principal_axis_rotation': nonlocal_ops.principal_axis_rotation_operator,
            'com_com_approach':        nonlocal_ops.com_com_approach_operator,
            'com_com_separation':      nonlocal_ops.com_com_separation_operator,
        }
        local_op_map = {
            'random_displacement':          local_ops.random_displacement_submolecule,
            'random_rotation':              local_ops.random_rotation_submolecule,
            'local_twist':                  local_ops.local_twist_submolecule,
            'correlated_displacement':      local_ops.correlated_displacement,
            'small_principal_axis_rotation': local_ops.small_principal_axis_rotation,
        }

        operators, weights = [], []
        for op_name, weight in operator_list:
            if op_name in nonlocal_op_map:
                operators.append({
                    'name': op_name,
                    'func': nonlocal_op_map[op_name],
                    'local': False,
                    'supports_adaptive': op_name in nonlocal_adaptive,
                })
                weights.append(weight)
            elif op_name in local_op_map:
                operators.append({
                    'name': op_name,
                    'func': local_op_map[op_name],
                    'local': True,
                    'supports_adaptive': True,
                })
                weights.append(weight)
            else:
                _log(f"Worker {worker_id}: Unknown operator '{op_name}', skipping.", level="warning")

        if not operators:
            _log(f"Worker {worker_id}: No valid operators found!", level="error")
            return result

        weights = np.array(weights, dtype=float)
        weights /= weights.sum()

        op_accepted = {op['name']: 0 for op in operators}
        op_attempted = {op['name']: 0 for op in operators}

        # ---- Energy evaluator ----
        evaluator = EnergyEvaluator(
            backend=config_dict.get('backend', 'psi4'),
            qm_method=config_dict.get('qm_method', 'hf'),
            qm_basis=config_dict.get('qm_basis', 'sto-3g'),
            xtb_method=config_dict.get('xtb_method', 'GFN2-xTB'),
            gpaw_mode=config_dict.get('gpaw_mode', 'lcao'),
            gpaw_basis=config_dict.get('gpaw_basis', 'dzp'),
            gpaw_xc=config_dict.get('gpaw_xc', 'B3LYP'),
        )

        # ---- Initial energy ----
        current_structure = copy.deepcopy(initial_molecule)
        try:
            current_energy = evaluator.evaluate(current_structure)
        except Exception as exc:
            _log(f"Worker {worker_id}: Failed to evaluate initial energy: {exc}", level="error")
            return result

        _log(f"Worker {worker_id}: Starting energy = {current_energy:.6f} Hartree")

        accepted_structures = []
        energy_trajectory = [(0, current_energy)]
        n_accepted = n_rejected = 0
        temperature = config_dict.get('temperature', 400.0)
        beta = 1.0 / (K_B * temperature)

        # ---- Main BHMC loop ----
        for step in range(n_steps):
            op_info = operators[np.random.choice(len(operators), p=weights)]
            op_attempted[op_info['name']] += 1

            try:
                # Apply operator
                if op_info['local']:
                    new_structure = op_info['func'](
                        current_structure, submolecule_indices, adaptive=adaptive
                    )
                elif op_info['supports_adaptive'] and adaptive:
                    new_structure = op_info['func'](
                        current_structure, submolecule_indices=submolecule_indices, adaptive=True
                    )
                else:
                    new_structure = op_info['func'](
                        current_structure, submolecule_indices=submolecule_indices
                    )

                # Clash check
                if nonlocal_ops._has_inter_submolecule_clashes(
                    new_structure, submolecule_indices=submolecule_indices, scale=0.5
                ):
                    _log(f"Worker {worker_id}: Clash after '{op_info['name']}', rejecting.", level="warning")
                    n_rejected += 1
                    box_controller.record_step(False)
                    energy_trajectory.append((step + 1, current_energy))
                    continue

                # No-move guard
                if np.allclose(new_structure.coordinates, current_structure.coordinates):
                    _log(f"Worker {worker_id}: No change from '{op_info['name']}', skipping.", level="warning")
                    n_rejected += 1
                    box_controller.record_step(False)
                    energy_trajectory.append((step + 1, current_energy))
                    continue

                # Energy evaluation
                try:
                    new_energy = evaluator.evaluate(new_structure)
                except Exception:
                    new_energy = None

                if new_energy is None:
                    _log(f"Worker {worker_id}: Energy evaluation failed, rejecting.", level="warning")
                    n_rejected += 1
                    box_controller.record_step(False)
                    energy_trajectory.append((step + 1, current_energy))
                    continue

                # Metropolis acceptance
                if new_energy < current_energy:
                    accept = True
                else:
                    delta_e_ev = (new_energy - current_energy) * HARTREE_TO_EV
                    accept = random.random() < np.exp(-beta * delta_e_ev)
                    if (step + 1) % 10 == 0:
                        _log(
                            f"Worker {worker_id}: Step {step+1}, op={op_info['name']}, "
                            f"ΔE={delta_e_ev:.4f} eV, accept={accept}"
                        )

                if accept:
                    current_structure = copy.deepcopy(new_structure)
                    current_energy = new_energy
                    accepted_structures.append((copy.deepcopy(current_structure), current_energy))
                    op_accepted[op_info['name']] += 1
                    n_accepted += 1
                    _log(f"Worker {worker_id}: Step {step+1} accepted, E={current_energy:.6f} Ha")
                else:
                    n_rejected += 1

                box_controller.record_step(accept)
                energy_trajectory.append((step + 1, current_energy))

            except Exception as exc:
                _log(f"Worker {worker_id}: Exception at step {step+1}: {exc}", level="error")
                n_rejected += 1
                box_controller.record_step(False)
                energy_trajectory.append((step + 1, current_energy))
                continue

            # Adaptive box update
            if box_controller.active and (step + 1) % box_controller.interval == 0:
                update_info = box_controller.update()
                box_updates.append({
                    "step": step + 1,
                    "box_size": update_info["current_size"],
                    "window_acceptance": update_info.get("window_acceptance"),
                    "grew": update_info.get("grew", False),
                    "stable_windows": update_info.get("stable_windows", 0),
                })
                if update_info.get('grew'):
                    _log(
                        f"Worker {worker_id}: Box grew at step {step+1}. "
                        f"Size={update_info['current_size']:.3f}, "
                        f"acc={update_info['window_acceptance']:.3f}"
                    )

        acc_rate = n_accepted / (n_accepted + n_rejected) * 100 if (n_accepted + n_rejected) > 0 else 0.0
        _log(f"Worker {worker_id}: Done. Accepted={n_accepted}, Rejected={n_rejected}, Rate={acc_rate:.1f}%")

        operator_rates = {
            op: (op_accepted[op] / op_attempted[op] * 100 if op_attempted[op] > 0 else 0.0)
            for op in op_accepted
        }

        result = {
            'accepted_structures': accepted_structures,
            'energy_trajectory': energy_trajectory,
            'box_updates': box_updates,
            'worker_id': worker_id,
            'operator_acceptance_rates': operator_rates,
        }

    finally:
        pr.disable()
        s = io.StringIO()
        pstats.Stats(pr, stream=s).sort_stats('cumulative').print_stats(10)
        print(f"\n--- Profiling: Worker {worker_id} ---\n{s.getvalue()}")
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        print(f"Memory — current: {current/1e6:.2f} MB, peak: {peak/1e6:.2f} MB")

    return result


# ======================== Quick self-test ================================

if __name__ == "__main__":
    h2o_xyz = """6
H2O dimer
O     -2.429551    2.071454   -0.000000
H     -2.835791    2.732744   -0.547872
H     -2.747575    2.159247    0.890872
O      0.269950    0.440735    0.000000
H      1.078961    0.704535    0.422388
H     -0.176739   -0.198484    0.542557
"""
    mol = Molecule.from_xyz(h2o_xyz)
    submolecule_indices = [[0, 1, 2], [3, 4, 5]]

    operator_list = [
        # non-local
        ("twist",               0.6),
        ("large_displacement",  0.5),
        # local
        ("random_displacement", 3.0),
        ("random_rotation",     2.5),
    ]

    config_dict = {
        'backend': 'xtb',
        'xtb_method': 'GFN2-xTB',
        'temperature': 300.0,
        'adaptive_operators': False,
        'verbose': True,
        'adaptive_box': False,
    }

    args = (mol, submolecule_indices, 20, operator_list, config_dict, 0, None, None)
    result = _bhmc_worker(args)
    print(f"\nAccepted structures : {len(result['accepted_structures'])}")
    print(f"Trajectory points  : {len(result['energy_trajectory'])}")
    print(f"Operator rates     : {result['operator_acceptance_rates']}")
