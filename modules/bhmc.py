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
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
from tqdm import tqdm
import time

# Matplotlib
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for plotting


# Profiling and Memory Tracking
import cProfile
import pstats
import io 
import tracemalloc
import sys

# Module imports
from modules.molecule_class import Molecule
from modules.operators import NonLocalOperators, LocalOperators

# Clustering and Feature Matrix Generation
from modules.featurizer import Featurizer, FeaturizerConfig
from modules.cluster import Clustering

from modules.box import SimulationBox
from modules.bhmc_config import BHMCConfig
from modules.calculator import EnergyEvaluator
from modules.logger import Logger
from modules.bhmc_worker import _phase_a_worker, _phase_b_worker





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
        ('com_com_approach', 0.8),
        ('com_com_separation', 0.3)
    
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
        self.worker_box_updates = {}   # Store adaptive box updates for each worker
        
        # Storage for Phase A results
        self.phase_a_structures: List[Tuple[Molecule, float, List[float], float]] = []

        # Storage for Phase B results
        self.phase_b_structures: List[Tuple[Molecule, float, List[float], float]] = []


    def _log(self, msg: str, level: str = "info") -> None:
        """Log a summary message to the main log file."""
        if self.config.verbose:
            print(msg)
        if self.logger:
            getattr(self.logger, level)(msg)

    def estimate_single_point_time(self, molecule: Molecule) -> Optional[float]:
        """ 
        Helper function that estimates the time of a SP calculation
        This can then be used to estimate the total time of a Phase
        """
        evaluator = EnergyEvaluator(
            backend=self.config.backend,
            qm_method=self.config.qm_method,
            qm_basis=self.config.qm_basis,
            xtb_method=self.config.xtb_method,
            gpaw_mode=self.config.gpaw_mode,
            gpaw_basis=self.config.gpaw_basis,
            gpaw_xc=self.config.gpaw_xc,
        )
        start_time = time.time()
        try:
            energy = evaluator.evaluate(molecule)
        except Exception as exc:
            self._log(f"Warning: Single point energy evaluation failed during timing estimation: {exc}", level="warning")
            return None
        end_time = time.time()
        elapsed = end_time - start_time
        self._log(f"Estimated single point time: {elapsed:.2f} seconds (Energy: {energy:.6f} Hartree)")
        return elapsed
    

    def run_phase_a(self, 
                    initial_molecules: List[Molecule],
                    submolecule_indices: List[List[int]],
                    n_structures_per_worker: int = 300,
                    n_processes: int = 10) -> List[Tuple[Molecule, float, List[float], float]]:
        """Run Phase A: Parallel exploration with independent BHMC chains.
        
        Each worker starts from a different initial configuration for better PES coverage.
        """
        # Memory Tracking and CPU Profiling
        pr = cProfile.Profile()
        pr.enable()
        tracemalloc.start()
        try: 
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
            self._log(f"  Method/Basis: {self.config.qm_method}/{self.config.qm_basis}")
            self._log(f"  Temperature: {self.config.temperature} K")
            self._log(f"  Adaptive operators: {self.config.adaptive_operators}")
            self._log(f"  Operators: {', '.join(op[0] for op in self.operators)}")
            self._log(f"  Independent initial configs: {n_processes}")
            if self.simulation_box:
                self._log(f"  Simulation Box: {self.simulation_box}")
            if self.worker_log_file:
                self._log(f"  Worker log: {self.worker_log_file}")
                # Clear worker log file at the start of Phase A
                with open(self.worker_log_file, 'w') as f:
                    f.write(f"BHMC Phase A Worker Log\n")
                    f.write(f"Configuration: {n_processes} workers, {n_structures_per_worker} steps each\n")
                    f.write(f"Method: {self.config.qm_method}, Basis: {self.config.qm_basis}, Temperature: {self.config.temperature} K\n")
                    f.write(f"Operators: {', '.join(op[0] for op in self.operators)}\n")
                    if self.simulation_box:
                        f.write(f"Simulation Box: {self.simulation_box}\n")
                    f.write("\n")
            
            config_dict = {
                'qm_method': self.config.qm_method,
                'qm_basis': self.config.qm_basis,
                'backend': self.config.backend,
                'xtb_method': self.config.xtb_method,
                'gpaw_mode': self.config.gpaw_mode,
                'gpaw_basis': self.config.gpaw_basis,
                'gpaw_xc': self.config.gpaw_xc,
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

            # Run initial sp calculation to estimate time
            sp_molecule = initial_molecules[0]
            sp_time = self.estimate_single_point_time(sp_molecule)
            if sp_time is not None:
                est_time = sp_time * n_structures_per_worker * n_processes
                self._log(f"Estimated Phase A time: {est_time/60:.2f} minutes")
                print(f"Estimated Phase A time: {est_time/60:.2f} minutes")


            
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
            ctx = multiprocessing.get_context("spawn")
            results = [None] * n_processes
            with ProcessPoolExecutor(max_workers=n_processes, mp_context=ctx) as executor:
                future_to_wid = {
                    executor.submit(_phase_a_worker, args): args[5]
                    for args in args_list
                }
                for future in tqdm(as_completed(future_to_wid), total=n_processes, desc="Phase A Workers"):
                    wid = future_to_wid[future]
                    results[wid] = future.result()

            
            # Collect results and trajectories
            all_accepted_structures = []
            self.worker_trajectories = {}  # Store for plotting
            self.worker_box_updates = {}   # Store adaptive box updates for plotting
            self.worker_operator_acceptance = {}  # Store operator acceptance rates for analysis 


            for worker_result in results:
                wid = worker_result['worker_id']
                accepted = worker_result['accepted_structures'] 
                trajectory = worker_result['energy_trajectory']
                box_updates = worker_result.get('box_updates', [])
                operator_rates = worker_result.get('operator_acceptance_rates', {})
                
                all_accepted_structures.extend(accepted)
                self.worker_trajectories[wid] = trajectory
                self.worker_box_updates[wid] = box_updates
                self.worker_operator_acceptance[wid] = operator_rates
                self._log(
                    f"  Worker {wid}: {len(accepted)} structures accepted, "
                    f"{len(trajectory)} trajectory points, {len(box_updates)} box updates"
                    f"  Operator acceptance rates: {', '.join(f'{op}: {rate:.1f}%' for op, rate in operator_rates.items())}"
                )

            # Log size of accepted structures for Phase A
            self._log(f"Size of accepted structures (Phase A): {sys.getsizeof(all_accepted_structures)/10**6:.2f} MB")
            
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
                energies = [e for _, e in all_accepted_structures]
                self._log(f"Energy Statistics (Hartree):")
                self._log(f"  Min:  {min(energies):.6f}")
                self._log(f"  Max:  {max(energies):.6f}")
                self._log(f"  Mean: {np.mean(energies):.6f}")
                self._log(f"  Std:  {np.std(energies):.6f}")

            
            # Plot worker trajectories
            self.plot_worker_trajectories()
            self.plot_adaptive_box_updates()
            self.plot_worker_operator_acceptance()

        except Exception as exc:
            self._log(f"Phase A failed with exception: {exc}", level="error")
            import traceback
            self._log(traceback.format_exc(), level="error")
            raise
        finally:
            pr.disable()
            current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            s = io.StringIO()
            ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
            ps.print_stats(10)
            # Log memory usage
            self._log(f"Profiling Results for Phase A:")
            self._log(f"Current memory usage: {current / 10**6:.2f} MB; Peak memory usage: {peak / 10**6:.2f} MB")

            # Write profiling results to an out file
            from pathlib import Path
            Path("profiling").mkdir(parents=True, exist_ok=True)
            with open("profiling/phase_a_profile.txt", 'w') as f:
                f.write(f"--- Profiling Results for Phase A ---\n")
                f.write(s.getvalue())
                f.write(f"\n--- Memory Usage for Phase A ---\n")
                f.write(f"Current memory usage: {current / 10**6:.2f} MB; Peak memory usage: {peak / 10**6:.2f} MB\n")

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
                'best_structure': (Molecule, energy) tuple of the best structure found
                'accepted_structures': List of (Molecule, energy) tuples
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
        self._log(f" Method/Basis: {self.config.qm_method}/{self.config.qm_basis}")
        self._log(f" Temperature: {self.config.temperature} K")
        self._log(f" Adaptive operators: {self.config.adaptive_operators}")
        self._log(f" Local Operators: {', '.join(op[0] for op in local_operators)}")
        if self.simulation_box:
            self._log(f" Simulation Box: {self.simulation_box}")
        if self.worker_log_file:
            self._log(f" Worker log: {self.worker_log_file}")
            # Append Header to worker log file for Phase B
            with open(self.worker_log_file, 'a') as f:
                f.write(f"\n\nBHMC Phase B Worker Log\n")
                f.write(f"Configuration: {n_processes} workers, {n_steps_per_worker} steps each\n")
                f.write(f"Method: {self.config.qm_method}, Basis: {self.config.qm_basis}, Temperature: {self.config.temperature} K\n")
                f.write(f"Local Operators: {', '.join(op[0] for op in local_operators)}\n")
                if self.simulation_box:
                    f.write(f"Simulation Box: {self.simulation_box}\n")
                f.write("\n")
            
        config_dict = {
            "qm_method": self.config.qm_method,
            "qm_basis": self.config.qm_basis,
            "backend": self.config.backend,
            "xtb_method": self.config.xtb_method,
            "gpaw_mode": self.config.gpaw_mode,
            "gpaw_basis": self.config.gpaw_basis,
            "gpaw_xc": self.config.gpaw_xc,
            "temperature": self.config.temperature,
            "max_temperature": self.config.max_temperature,
            "verbose": self.config.verbose,
            "adaptive_operators": self.config.adaptive_operators,
        }

        # Estimate time for a single local refinement step
        sp_molecule = representatives[0].molecule
        sp_time = self.estimate_single_point_time(sp_molecule)
        if sp_time is not None:
            est_time = sp_time * n_steps_per_worker * n_workers
            self._log(f"Estimated Phase B time: {est_time/60:.2f} minutes")
            print(f"Estimated Phase B time: {est_time/60:.2f} minutes")
        

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
        ctx = multiprocessing.get_context("spawn")
        results = [None] * n_workers
        with ProcessPoolExecutor(max_workers=n_processes, mp_context=ctx) as executor:
            future_to_wid = {
                executor.submit(_phase_b_worker, args): args[6]
                for args in args_list
            }
            for future in tqdm(as_completed(future_to_wid), total=n_workers, desc="Phase B Workers"):
                wid = future_to_wid[future]
                results[wid] = future.result()

        # Collect results and trajectories
        self.phase_b_trajectories = {}
        self.phase_b_results = results
        all_accepted_structures = []


        for worker_result in results:
            wid = worker_result["worker_id"]
            accepted = worker_result["accepted_structures"]
            trajectory = worker_result["energy_trajectory"]
            best = worker_result["best_structure"]
            all_accepted_structures.extend(accepted)

            self.phase_b_trajectories[wid] = trajectory
            best_e_str = f"{best[1]:.6f} Ha" if best else "N/A"  # best = (mol, energy)
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
            energies = [bs[1] for bs in best_structures]  # (mol, energy)
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

        

        self.phase_b_structures = all_accepted_structures

        # Plot local basin exploration for each worker
        self.plot_phase_b_basin_exploration(results)

        return results

    def analyse_phase_a_results(
            self,
            phase_a_structures: Optional[List[Tuple[Molecule, float]]] = None):
        """
        Applies Clustering to the Structures obtained from Phase A to select representative structures for Phase B local refinement  
        """
        if self.logger:
            self.logger.header("Phase A Analysis")

        if not phase_a_structures:
            # Use self.phase_a_structures if no argument provided
            self._log("No Phase A structures provided for analysis, using internal storage", level="warning")
            phase_a_structures = self.phase_a_structures


        # Configure the Featurizer and Clustering
        featurizer_config = FeaturizerConfig(
            descriptor_type="SOAP",
            soap_rcut = 6,
            soap_nmax = 12,
            soap_lmax = 6,
            soap_sigma = 0.5,
        )
        featurizer = Featurizer(config=featurizer_config)
        mols, energies = zip(*phase_a_structures)
        feature_mat = featurizer.build_feature_matrix(mols, energies=energies, include_hbonds=True)
        self._log(f"Computed feature matrix for {len(mols)} structures with shape {feature_mat.shape}")
        clustering = Clustering(feature_matrix=feature_mat, normalize=True, logger=self.logger)
        # Benchmark the clustering operations per sample for algorithm selection
        _, best_algorithm_complexity = clustering.benchmark_operations_per_sample()
        print(f"Best clustering algorithm based on complexity: {best_algorithm_complexity}")  

        if best_algorithm_complexity == "dbscan_avg" or best_algorithm_complexity == "dbscan_worst":
            clustering.cluster(method="dbscan", eps=0.5, min_samples=20)
        
        embedding = clustering.umap(n_components=2)
        clustering.plot_embedding(embedding, title="Phase A Clustering Embedding", save_path="figures/phase_a_clustering_embedding.png")


        
        
        


    def analyse_phase_b_results(self, 
                                phase_b_results: List[Dict],
                                submolecule_indices: Optional[List[List[int]]] = None,
                                logger: Optional[Logger] = None
                                ) -> List[Tuple[Molecule, float, List[float], float]]:
        """  
        Analyzes the Results of the Phase B local refinement and extracts the best structures found by each worker, along with their energies and dipole moments.
        
        Args:
            phase_b_results: List of dictionaries returned by each Phase B worker, containing keys 'best_structure', 'accepted_structures', 'energy_trajectory', and 'worker_id'.
        """
        if self.logger:
            self.logger.header("Phase B Analysis")

        # Initialize an analyzer

        
       
        
    # Plot Exploration of local funnels of phase b
    def plot_phase_b_basin_exploration(self,
                                       phase_b_results: List[Dict],):
        """  
        Plots Energy vs RMSD to the best structure for each worker in Phase B to visualize the local exploration of a given funnel
        """
        import matplotlib.pyplot as plt
        from pathlib import Path

        Path("figures").mkdir(parents=True, exist_ok=True)
        
        # Compute RMSD to best structure for each accepted structure
        rmsd_data = []
        for result in phase_b_results:
            worker_id = result["worker_id"]
            accepted = result["accepted_structures"]
            best_structure = result["best_structure"][0] if result["best_structure"] else None
            if best_structure is None:
                continue
            for mol, energy, _ , _, _ in accepted:
                rmsd = best_structure.compute_rmsd(mol)
                rmsd_data.append((worker_id, energy, rmsd))
        
        # Make suplots with 4 cols depending on number of workers
        n_workers = len(phase_b_results)
        n_cols = min(4, n_workers)
        n_rows = (n_workers + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows), sharex=True, sharey=True)
        axes = axes.flatten() if n_workers > 1 else [axes]
        for worker_id, energy, rmsd in rmsd_data:
            ax = axes[worker_id]
            ax.scatter(rmsd, energy, alpha=0.6)
            ax.set_title(f"Worker {worker_id}")
            ax.set_xlabel("RMSD to best structure (A)")
            ax.set_ylabel("Energy (Hartree)")
            ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig("figures/phase_b_basin_exploration.png", dpi=300)
        plt.close()

            
            



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

    def plot_adaptive_box_updates(
            self,
            save_path: str = "figures/adaptive_box_trajectories.png",
            box_updates_by_worker: Optional[Dict[int, List[Dict]]] = None,
            title: str = "Adaptive Box Control"):
        """Plot adaptive box size and window acceptance per worker."""
        import matplotlib.pyplot as plt
        from pathlib import Path

        if box_updates_by_worker is None:
            box_updates_by_worker = self.worker_box_updates

        if not box_updates_by_worker:
            self._log("No adaptive box updates to plot.", level="warning")
            return

        valid_updates = {wid: updates for wid, updates in box_updates_by_worker.items() if updates}
        if not valid_updates:
            self._log("Adaptive box updates are empty for all workers.", level="warning")
            return

        Path(save_path).parent.mkdir(parents=True, exist_ok=True)

        fig, (ax_size, ax_acc) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
        cmap = plt.get_cmap('tab20', max(len(valid_updates), 1))

        for i, (worker_id, updates) in enumerate(sorted(valid_updates.items(), key=lambda x: x[0])):
            color = cmap(i)

            steps = [u["step"] for u in updates]
            sizes = [u["box_size"] for u in updates]
            ax_size.plot(steps, sizes, marker='o', linewidth=1.8, markersize=4,
                         label=f"Worker {worker_id}", color=color)

            acc_steps = [u["step"] for u in updates if u["window_acceptance"] is not None]
            acc_vals = [u["window_acceptance"] for u in updates if u["window_acceptance"] is not None]
            if acc_vals:
                ax_acc.plot(acc_steps, acc_vals, marker='o', linewidth=1.5, markersize=4,
                            label=f"Worker {worker_id}", color=color)

        target = self.config.box_target_acceptance
        window = self.config.box_acceptance_window
        ax_acc.axhline(target, linestyle='--', linewidth=1.2, color='black', alpha=0.8, label='Target')
        ax_acc.axhline(target - window, linestyle=':', linewidth=1.0, color='gray', alpha=0.8, label='Target window')
        ax_acc.axhline(target + window, linestyle=':', linewidth=1.0, color='gray', alpha=0.8)

        ax_size.set_ylabel("Box size (A)")
        ax_size.set_title(title)
        ax_size.grid(True, alpha=0.3)
        ax_size.legend(ncols=4, loc='best')

        ax_acc.set_xlabel("Step")
        ax_acc.set_ylabel("Window acceptance")
        ax_acc.set_ylim(0.0, 1.0)
        ax_acc.grid(True, alpha=0.3)
        ax_acc.legend(ncols=4, loc='best')

        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()

    # plot operator acceptance rates for each worker
    def plot_worker_operator_acceptance(self,
                                      save_path: str = "figures/worker_operator_acceptance.png"):
        """  
        Plots the Operator Acceptance Rates for each worker in Phase A as a bar chart
        """ 
        import matplotlib.pyplot as plt
        from pathlib import Path
        ncols = max(1, min(4, len(self.worker_operator_acceptance)))
        nrows = (len(self.worker_operator_acceptance) + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows), sharey=True)
        axes = axes.flatten() if len(self.worker_operator_acceptance) > 1 else [axes]
        for i, (worker_id, op_rates) in enumerate(sorted(self.worker_operator_acceptance.items(), key=lambda x: x[0])):
            ax = axes[i]
            operators = list(op_rates.keys())
            rates = [op_rates[op] for op in operators]
            ax.bar(operators, rates, color='skyblue')
            # Adjust tick labels for better readability
            ax.set_xticklabels(operators, rotation=45, ha='right')
            ax.set_title(f"Worker {worker_id}")
            ax.set_xlabel("Operator")
            ax.set_ylabel("Acceptance Rate (%)")
            ax.set_ylim(0, 100)
            ax.grid(True, alpha=0.3)
        # Hide any unused subplots
        for j in range(i + 1, len(axes)):
            axes[j].axis('off')
        plt.tight_layout()
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
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



# Testing

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
    from init import ClusterInitializer, InitializationConfig
    from bhmc import MultiPhaseBHMC, BHMCConfig
    from logger import Logger 



    mp.set_start_method('spawn', force=True)
    

    xyz_path = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/h2o_2_nonopt.xyz"

    init_config = InitializationConfig(
        method="hf",
        basis="cc-pvdz",
        box_type="sphere",
        box_scale_factor=3,
        min_distance=1.5,
        optimize_submolecules=True,
        verbose=False
    )
    initializer = ClusterInitializer(config=init_config)
    initial_molecules, submol_indices, simulation_box = initializer.initialize_from_xyz(
        xyz_path, n_configurations=1
    )
    bhmc_config = BHMCConfig(
        temperature=1,
        qm_method="hf",
        qm_basis="6-311G(d,p)",
        verbose=False,
        adaptive_operators=True,
        adaptive_box = True,
        box_update_interval=5, # Update box every 5 steps
        box_target_acceptance =0.6,
        box_acceptance_window = 0.05,
        box_growth_max = 1.15,
        box_max_scale = 4.0,
        box_stable_windows = 3
    )

    bhmc_sampler = MultiPhaseBHMC(
        config=bhmc_config,
        simulation_box=simulation_box,
        worker_log_file="bhmc_workers.log"
    )
    # Initialize StuctureData with one molecule
    initial_structure = StructureData(
        molecule=initial_molecules[0],
        energy=0.0,
        dipole_magnitude=0.0,
        phase="initial" 
    )
    representatives = [initial_structure]
    submolecule_indices = submol_indices
    phase_b_results = bhmc_sampler.run_phase_b(
        representatives=representatives,
        submolecule_indices=submolecule_indices,
        n_steps_per_worker=300,
        local_operators=None,
        n_processes=1
    )
    # Get best structures from results
    best_structures = [r["best_structure"] for r in phase_b_results if r["best_structure"] is not None]
    # write best structures to xyz

    # Get the trajectory
    trajectory = bhmc_sampler.phase_b_structures

    # 
    
    for i, best in enumerate(best_structures):
        mol, energy = best
        mol.write_xyz(f"best_structure_worker_{i}.xyz")
        print(f"Worker {i}: Best energy = {energy:.6f} Ha")
    
    


