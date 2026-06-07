"""
Basin Hopping Monte Carlo (BHMC) implementation for molecular cluster exploration.

Single-loop parallel BHMC: local and non-local operators are drawn from one unified
probability distribution. Non-local moves provide basin hops; local moves explore
within a basin. The balance is controlled by the weights in BHMCConfig.operators.
"""

import numpy as np
from typing import List, Optional, Tuple, Dict
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
import cProfile
import pstats
import io
import tracemalloc
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC

try:
    from modules.molecule_class import Molecule
except ImportError:
    from molecule_class import Molecule

from modules.featurizer import Featurizer, FeaturizerConfig
from modules.cluster import Clustering
from modules.box import SimulationBox
from modules.bhmc_config import BHMCConfig
from modules.calculator import EnergyEvaluator
from modules.logger import Logger
from modules.bhmc_worker import _bhmc_worker
from tqdm import tqdm


class BHMC:
    """Parallel Basin Hopping Monte Carlo for molecular cluster exploration."""

    def __init__(
        self,
        config: BHMCConfig,
        simulation_box: Optional[SimulationBox] = None,
        logger: Optional[Logger] = None,
        worker_log_file: str = "bhmc_workers.log",
    ):
        self.config = config
        self.simulation_box = simulation_box
        self.logger = logger
        self.worker_log_file = worker_log_file

        # Populated after run()
        self.accepted_structures: List[Tuple[Molecule, float]] = []
        self.worker_trajectories: Dict = {}
        self.worker_box_updates: Dict = {}
        self.worker_operator_acceptance: Dict = {}

    # ------------------------------------------------------------------ logging

    def _log(self, msg: str, level: str = "info") -> None:
        if self.config.verbose:
            print(msg)
        if self.logger:
            getattr(self.logger, level)(msg)

    # ---------------------------------------------------------- time estimation

    def estimate_single_point_time(self, molecule: Molecule) -> Optional[float]:
        """Run one SP calculation and return its wall-clock time in seconds."""
        evaluator = EnergyEvaluator(
            backend=self.config.backend,
            qm_method=self.config.qm_method,
            qm_basis=self.config.qm_basis,
            xtb_method=self.config.xtb_method,
            gpaw_mode=self.config.gpaw_mode,
            gpaw_basis=self.config.gpaw_basis,
            gpaw_xc=self.config.gpaw_xc,
        )
        t0 = time.time()
        try:
            energy = evaluator.evaluate(molecule)
        except Exception as exc:
            self._log(f"SP timing failed: {exc}", level="warning")
            return None
        elapsed = time.time() - t0
        self._log(f"SP time estimate: {elapsed:.2f} s  (E = {energy:.6f} Ha)")
        return elapsed

    # --------------------------------------------------------------- main entry

    def run(
        self,
        initial_molecules: List[Molecule],
        submolecule_indices: List[List[int]],
        n_steps_per_worker: int = 300,
        n_processes: Optional[int] = None,
    ) -> List[Tuple[Molecule, float]]:
        """
        Launch parallel BHMC chains and return all accepted (molecule, energy) pairs.

        Each worker starts from its own initial configuration and draws operators
        from the unified pool defined in config.operators.
        """
        if n_processes is None:
            n_processes = len(initial_molecules)

        if len(initial_molecules) != n_processes:
            raise ValueError(
                f"Number of initial molecules ({len(initial_molecules)}) "
                f"must match n_processes ({n_processes})"
            )

        if self.logger:
            self.logger.header("BHMC Run")

        self._log(f"Workers            : {n_processes}")
        self._log(f"Steps / worker     : {n_steps_per_worker}")
        self._log(f"Total steps        : {n_processes * n_steps_per_worker}")
        if self.config.backend == "psi4":
            self._log(f"QM method         : {self.config.qm_method} / {self.config.qm_basis}")
        if self.config.backend == "xtb":
            self._log(f"xTB method        : {self.config.xtb_method}")
        if self.config.backend == "gpaw":
            self._log(f"GPAW mode         : {self.config.gpaw_mode}")
            self._log(f"GPAW basis        : {self.config.gpaw_basis}")
            self._log(f"GPAW xc           : {self.config.gpaw_xc}")
    
        self._log(f"Temperature        : {self.config.temperature} K")
        self._log(f"Adaptive operators : {self.config.adaptive_operators}")
        self._log(f"Operators          : {[op[0] for op in self.config.operators]}")
        if self.simulation_box:
            self._log(f"Simulation box     : {self.simulation_box}")

        # Write worker log header
        if self.worker_log_file:
            with open(self.worker_log_file, 'w') as f:
                f.write("BHMC Worker Log\n")
                f.write(f"{n_processes} workers × {n_steps_per_worker} steps\n")
                # Write method/config info for each worker
                if self.config.backend == "psi4":
                    f.write(f"QM method         : {self.config.qm_method} / {self.config.qm_basis}\n")
                if self.config.backend == "xtb":
                    f.write(f"xTB method        : {self.config.xtb_method}\n")
                if self.config.backend == "gpaw":
                    f.write(f"GPAW mode         : {self.config.gpaw_mode}\n")
                    f.write(f"GPAW basis        : {self.config.gpaw_basis}\n")
                    f.write(f"GPAW xc           : {self.config.gpaw_xc}\n")
                f.write(f"Temperature        : {self.config.temperature} K\n")
                f.write(f"Adaptive operators : {self.config.adaptive_operators}\n")
                f.write(f"Operators          : {[op[0] for op in self.config.operators]}\n")
                if self.simulation_box:
                    f.write(f"Simulation box     : {self.simulation_box}\n")
                f.write("\n")

                        

        # SP timing estimate
        sp_time = self.estimate_single_point_time(initial_molecules[0])
        if sp_time is not None:
            est_min = sp_time * n_steps_per_worker * n_processes / 60
            self._log(f"Estimated run time : {est_min:.1f} min")
            print(f"Estimated run time : {est_min:.1f} min")

        config_dict = {
            'backend':            self.config.backend,
            'qm_method':          self.config.qm_method,
            'qm_basis':           self.config.qm_basis,
            'xtb_method':         self.config.xtb_method,
            'gpaw_mode':          self.config.gpaw_mode,
            'gpaw_basis':         self.config.gpaw_basis,
            'gpaw_xc':            self.config.gpaw_xc,
            'temperature':        self.config.temperature,
            'verbose':            self.config.verbose,
            'adaptive_operators': self.config.adaptive_operators,
            'adaptive_box':       self.config.adaptive_box,
            'box_update_interval':    self.config.box_update_interval,
            'box_target_acceptance':  self.config.box_target_acceptance,
            'box_acceptance_window':  self.config.box_acceptance_window,
            'box_growth_kp':          self.config.box_growth_kp,
            'box_growth_max':         self.config.box_growth_max,
            'box_max_scale':          self.config.box_max_scale,
            'box_stable_windows':     self.config.box_stable_windows,
        }

        sim_box_dict = self.simulation_box.to_dict() if self.simulation_box else None

        args_list = [
            (
                initial_molecules[wid],
                submolecule_indices,
                n_steps_per_worker,
                self.config.operators,
                config_dict,
                wid,
                sim_box_dict,
                self.worker_log_file,
            )
            for wid in range(n_processes)
        ]

        # ---- Profiling around the parallel section ----
        pr = cProfile.Profile()
        pr.enable()
        tracemalloc.start()

        try:
            self._log("Launching workers...")
            ctx = multiprocessing.get_context("spawn")
            results = [None] * n_processes
            with ProcessPoolExecutor(max_workers=n_processes, mp_context=ctx) as executor:
                future_to_wid = {
                    executor.submit(_bhmc_worker, args): args[5]
                    for args in args_list
                }
                for future in tqdm(as_completed(future_to_wid), total=n_processes, desc="BHMC Workers"):
                    wid = future_to_wid[future]
                    results[wid] = future.result()

            # ---- Collect results ----
            all_accepted: List[Tuple[Molecule, float]] = []
            self.worker_trajectories = {}
            self.worker_box_updates = {}
            self.worker_operator_acceptance = {}

            for r in results:
                wid = r['worker_id']
                accepted = r['accepted_structures']
                all_accepted.extend(accepted)
                self.worker_trajectories[wid] = r['energy_trajectory']
                self.worker_box_updates[wid] = r.get('box_updates', [])
                self.worker_operator_acceptance[wid] = r.get('operator_acceptance_rates', {})
                self._log(
                    f"  Worker {wid}: {len(accepted)} accepted, "
                    f"op rates: {r.get('operator_acceptance_rates', {})}"
                )

            self._log(f"Total accepted structures : {len(all_accepted)}")
            self._log(f"Memory (accepted list)    : {sys.getsizeof(all_accepted)/1e6:.2f} MB")

            if all_accepted:
                energies = [e for _, e in all_accepted]
                self._log(f"Energy (Ha) — min={min(energies):.6f}  max={max(energies):.6f}  "
                          f"mean={np.mean(energies):.6f}  std={np.std(energies):.6f}")

            self.accepted_structures = all_accepted

            # ---- Plots ----
            self.plot_worker_trajectories()
            self.plot_adaptive_box_updates()
            self.plot_worker_operator_acceptance()

        finally:
            pr.disable()
            current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            s = io.StringIO()
            pstats.Stats(pr, stream=s).sort_stats('cumulative').print_stats(10)
            self._log(f"Peak memory: {peak/1e6:.2f} MB")
            Path("profiling").mkdir(parents=True, exist_ok=True)
            with open("profiling/bhmc_profile.txt", 'w') as f:
                f.write(s.getvalue())
                f.write(f"\nPeak memory: {peak/1e6:.2f} MB\n")

        return self.accepted_structures

    # -------------------------------------------------------- structure analysis
    def analyze_train(
            self,
            reference_clustering: Clustering,
            structures: Optional[List[Tuple[Molecule, float]]] = None,
    ):
        """
        Analyzes the initial training phase of the BHMC sampling.

        Args:
            reference_clustering: The fitted Clustering instance returned by
                ClusterInitializer.initialize_from_xyz(). It carries the trained UMAP
                embedding/mapper (`reference_clustering.mapper`), the cluster labels and
                feature matrix it was trained on (`reference_clustering.labels`,
                `reference_clustering._fm()`), and the trained classifier
                (`reference_clustering.classifier_model`) — i.e. the whole reference model
                as one object, ready to be used for embedding/flagging/refitting the newly
                sampled structures (see Clustering.embed_new_structures,
                flag_novel_structures, refit_with_augmented_data).
            structures: Structures to analyze. Defaults to self.accepted_structures.
        """
        if structures is None:
            structures = self.accepted_structures
        # 1. Remove high energy outliers (Z-score > 3)
        mols = [m for m, _ in structures]
        energies = [e for _, e in structures]
        energy_arr = np.array(energies)
        z_e = (energy_arr - energy_arr.mean()) / (energy_arr.std() + 1e-12)
        mols     = [m for m, z in zip(mols, z_e)     if z <= 3.0]
        energies = [e for e, z in zip(energies, z_e) if z <= 3.0]
        self._log(f"After energy Z-score filter: {len(mols)} structures remain")

        # Log the energy distribution
        if energies:
            self._log(f"Energy (Ha) — min={min(energies):.6f}  max={max(energies):.6f}  "
                      f"mean={np.mean(energies):.6f}  std={np.std(energies):.6f}")
        else:
            self._log("No structures remain after energy Z-score filter.", level="warning")
            return
        # 2. Featurize with SOAP
        featurizer = Featurizer(FeaturizerConfig(descriptor_type="SOAP"))
        feature_mat = featurizer.build_feature_matrix(mols, energies=None, include_hbonds=False)
        self._log(f"SOAP features: {feature_mat.shape[0]} × {feature_mat.shape[1]}")
        
        # 3. Embedd the new structures using the reference Mapper and find novel structures
        novel_idx = reference_clustering.flag_novel_structures(feature_mat, threshold_percentile=25)
        self._log(f"Novel structures flagged: {len(novel_idx)} / {len(mols)}")
        
        

        



    def analyze_results(
        self,
        umap_model,
        classifier,
        structures: Optional[List[Tuple[Molecule, float]]] = None,
    ) -> List[Tuple[Molecule, float, int]]:
        """
        Run accepted structures through outlier removal, project into a pre-fitted
        UMAP space, and assign cluster labels via a pre-trained classifier.

        Returns list of (molecule, energy, cluster_label) for surviving structures.
        """
        if structures is None:
            structures = self.accepted_structures

        if not structures:
            self._log("No structures to analyse.", level="warning")
            return []

        mols = [m for m, _ in structures]
        energies = [e for _, e in structures]

        # Energy Z-score filter
        energy_arr = np.array(energies)
        z_e = (energy_arr - energy_arr.mean()) / (energy_arr.std() + 1e-12)
        mols     = [m for m, z in zip(mols, z_e)     if z <= 3.0]
        energies = [e for e, z in zip(energies, z_e) if z <= 3.0]
        self._log(f"After energy Z-score filter: {len(mols)} structures remain")

        if not mols:
            self._log("All structures filtered as high-energy outliers.", level="warning")
            return []

        featurizer = Featurizer(FeaturizerConfig(descriptor_type="SOAP"))
        feature_mat = featurizer.build_feature_matrix(mols, energies=None, include_hbonds=False)
        self._log(f"SOAP features: {feature_mat.shape[0]} × {feature_mat.shape[1]}")

        clustering = Clustering(
            feature_matrix=feature_mat, energies=energies,
            metric="cityblock", normalize=False, logger=self.logger,
        )

        # Feature Z-score filter
        idx = clustering.z_score_filtering()
        mols     = [mols[i]     for i in idx]
        energies = [energies[i] for i in idx]
        self._log(f"After feature Z-score filter: {len(mols)} structures remain")

        # Isolation forest
        idx = clustering.isolation_forest_outlier(contamination=0.4)
        mols     = [mols[i]     for i in idx]
        energies = [energies[i] for i in idx]
        self._log(f"After Isolation Forest: {len(mols)} structures remain")

        if not mols:
            self._log("No structures remain after outlier filtering.", level="warning")
            return []

        embedding = umap_model.transform(clustering._feature_matrix_normalized)
        labels = classifier.predict(embedding)
        labels_str = [f"Cluster {l}" for l in labels]
        self._log(f"Cluster distribution: {np.bincount(labels)}")

        if isinstance(classifier, KNeighborsClassifier):
            proba = clustering.KN_predict_probabilities(x_test=embedding, knn=classifier)
            clustering.plot_KN_probabilities(
                embedding=embedding, probabilities=proba, labels=labels_str,
                save_path="figures/knn_probabilities.png",
            )
        elif isinstance(classifier, SVC):
            proba = clustering.SVM_predict_probabilities(svm=classifier, x_test=embedding)
            clustering.plot_SVM_probabilities(
                embedding=embedding, probabilities=proba, labels=labels_str,
                save_path="figures/svm_probabilities.png",
            )

        clustering.plot_embedding(
            embedding=embedding, labels=labels_str,
            save_path="figures/bhmc_embedding.png",
        )

        return [(m, e, int(l)) for m, e, l in zip(mols, energies, labels)]

    # ----------------------------------------------------------------- plotting

    def plot_worker_trajectories(
        self,
        save_path: str = "figures/worker_energy_trajectories.png",
        trajectories: Optional[Dict] = None,
        title: str = "BHMC Worker Energy Trajectories",
    ):
        """Plot energy vs step for each worker."""
        import matplotlib.pyplot as plt

        if trajectories is None:
            trajectories = self.worker_trajectories
        if not trajectories:
            self._log("No trajectories to plot.", level="warning")
            return

        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        cmap = plt.get_cmap('tab20', max(len(trajectories), 1))

        plt.figure(figsize=(12, 6))
        for wid, traj in trajectories.items():
            steps, energies = zip(*traj)
            plt.plot(steps, energies, label=f"Worker {wid}", color=cmap(wid))
        plt.xlabel("Step")
        plt.ylabel("Energy (Hartree)")
        plt.title(title)
        plt.legend(ncols=4, loc="upper right", columnspacing=0.5)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()

    def plot_adaptive_box_updates(
        self,
        save_path: str = "figures/adaptive_box_trajectories.png",
        box_updates_by_worker: Optional[Dict] = None,
        title: str = "Adaptive Box Control",
    ):
        """Plot box size and window acceptance rate per worker."""
        import matplotlib.pyplot as plt

        if box_updates_by_worker is None:
            box_updates_by_worker = self.worker_box_updates

        valid = {wid: u for wid, u in box_updates_by_worker.items() if u}
        if not valid:
            self._log("No adaptive box updates to plot.", level="warning")
            return

        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig, (ax_size, ax_acc) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
        cmap = plt.get_cmap('tab20', max(len(valid), 1))

        for i, (wid, updates) in enumerate(sorted(valid.items())):
            color = cmap(i)
            steps = [u["step"] for u in updates]
            sizes = [u["box_size"] for u in updates]
            ax_size.plot(steps, sizes, marker='o', linewidth=1.8, markersize=4,
                         label=f"Worker {wid}", color=color)

            acc_steps = [u["step"] for u in updates if u["window_acceptance"] is not None]
            acc_vals  = [u["window_acceptance"] for u in updates if u["window_acceptance"] is not None]
            if acc_vals:
                ax_acc.plot(acc_steps, acc_vals, marker='o', linewidth=1.5, markersize=4,
                            label=f"Worker {wid}", color=color)

        target = self.config.box_target_acceptance
        window = self.config.box_acceptance_window
        ax_acc.axhline(target,          linestyle='--', linewidth=1.2, color='black', alpha=0.8, label='Target')
        ax_acc.axhline(target - window, linestyle=':',  linewidth=1.0, color='gray',  alpha=0.8, label='Window')
        ax_acc.axhline(target + window, linestyle=':',  linewidth=1.0, color='gray',  alpha=0.8)

        ax_size.set_ylabel("Box size (Å)")
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

    def plot_worker_operator_acceptance(
        self,
        save_path: str = "figures/worker_operator_acceptance.png",
    ):
        """Bar chart of operator acceptance rates per worker."""
        import matplotlib.pyplot as plt

        if not self.worker_operator_acceptance:
            return

        n = len(self.worker_operator_acceptance)
        ncols = max(1, min(4, n))
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharey=True)
        axes = axes.flatten() if n > 1 else [axes]

        for i, (wid, op_rates) in enumerate(sorted(self.worker_operator_acceptance.items())):
            ax = axes[i]
            ops   = list(op_rates.keys())
            rates = [op_rates[op] for op in ops]
            ax.bar(ops, rates, color='skyblue')
            ax.set_xticklabels(ops, rotation=45, ha='right')
            ax.set_title(f"Worker {wid}")
            ax.set_ylabel("Acceptance rate (%)")
            ax.set_ylim(0, 100)
            ax.grid(True, alpha=0.3)

        for j in range(i + 1, len(axes)):
            axes[j].axis('off')

        plt.tight_layout()
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=300)
        plt.close()

    # ---------------------------------------------------------- interpolation

    def generate_interpolated_candidates(
        self,
        representatives: List,
        submolecule_indices: Optional[List[List[int]]] = None,
        n_interpolations: int = 2,
        lambdas: Optional[List[float]] = None,
    ) -> List:
        """
        Generate interpolated structures between all pairs of representatives
        using the Kabsch algorithm.
        """
        from transformations import Transformation
        from cluster import StructureData

        transformer = Transformation()
        candidates = []
        n = len(representatives)
        self._log(f"Interpolating between {n} representatives ({n*(n-1)//2} pairs, {n_interpolations} points each)")

        for i in range(n):
            for j in range(i + 1, n):
                mols = transformer.kabsch_interpolate(
                    representatives[i].molecule,
                    representatives[j].molecule,
                    submolecule_indices=submolecule_indices,
                    lambdas=lambdas,
                    n_points=n_interpolations,
                )
                for mol in mols:
                    candidates.append(StructureData(
                        molecule=mol, energy=0.0, phase="interpolated",
                        metadata={"parent_i": i, "parent_j": j},
                    ))

        self._log(f"Generated {len(candidates)} interpolated candidates")
        return candidates


# ------------------------------------------------------------------ quick test

if __name__ == "__main__":
    from modules.bhmc_config import BHMCConfig

    config = BHMCConfig(
        backend="xtb",
        xtb_method="GFN2-xTB",
        temperature=300.0,
        verbose=True,
        adaptive_box=False,
    )

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

    bhmc = BHMC(config=config)
    accepted = bhmc.run(
        initial_molecules=[mol, mol],
        submolecule_indices=submolecule_indices,
        n_steps_per_worker=10,
        n_processes=2,
    )
    print(f"Accepted structures: {len(accepted)}")
