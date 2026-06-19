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
from modules.init import InitializationConfig
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
    
    def analyze_training_results(
            self,
            reference_clustering: Clustering,
            feature_mat_init: np.ndarray,
            mols_init: Optional[List[Molecule]] = None,
            structures: Optional[List[Tuple[Molecule, float]]] = None,
        ) -> Tuple[Optional[Clustering], Optional[List[Molecule]]]:
        """
        Analyze the first `training_frac` fraction of the BHMC run and use it to improve
        the reference classifier model.

        Workflow:
          1. Filter high-energy outliers
          2. Featurize with SOAP
          3. Embed into the reference UMAP by INIT
          4. Predict labels/probabilities with the reference classifier
          5. Flag "novel" structures (the ones the classifier is least confident about)
          6. Stack novel feature rows onto the initial raw feature matrix
          7. Run back through the same cleaning + clustering + classifier pipeline used
                during initialization (low-variance pruning, Isolation Forest, UMAP, LOF,
                clustering, classifier fit/evaluation), producing an updated Clustering instance
                with a more representative classifier.

        Args:
            reference_clustering: The Clustering instance fitted during initialization.
            feature_mat_init:     The raw (pre-pruning) feature matrix from initialization
                                   (e.g. ClusterInitializer.initialize_from_xyz's
                                   raw_feature_matrix return value).
            mols_init:            Molecules aligned 1:1 with feature_mat_init's rows (e.g.
                                   the corresponding clustering.raw_molecules). Required to
                                   resolve cluster representatives back to actual Molecule
                                   objects; if omitted, representatives cannot be resolved
                                   and None is returned in their place.
            structures:            (molecule, energy) pairs to analyze; defaults to
                                   self.accepted_structures.

        Returns:
            (updated_clustering, representative_mols) — representative_mols is one
            Molecule per retained cluster, or None if there was nothing to analyze,
            no novel structures were found, or mols_init was not provided.
        """
        self.logger.section("Analysis of BHMC Training Run")
        if structures is None:
            structures = self.accepted_structures

        if not structures:
            self._log("No structures to analyse.", level="warning")
            return None, None

        mols = [m for m, _ in structures]
        energies = [e for _, e in structures]
        energy_arr = np.array(energies)
        z_e = (energy_arr - energy_arr.mean()) / (energy_arr.std() + 1e-12)
        mols     = [m for m, z in zip(mols, z_e)     if z <= 3.0]
        energies = [e for e, z in zip(energies, z_e) if z <= 3.0]
        self._log(f"After energy Z-score filter: {len(mols)} structures remain")
        if not mols:
            self._log("All structures filtered as high-energy outliers.", level="warning")
            return None, None

        # Featurize with SOAP
        featurizer = Featurizer(FeaturizerConfig(descriptor_type="SOAP"))
        feature_mat = featurizer.build_feature_matrix(mols, energies=None, include_hbonds=False)
        self._log(f"SOAP features: {feature_mat.shape[0]} × {feature_mat.shape[1]}")

        # Embed the new structures using the reference_clustering object
        new_structures_embedding = reference_clustering.embed_new_structures(x_new = feature_mat)
        # Predict cluster labels using the classifier model in the reference clustering
        predicted_labels = reference_clustering.predict_labels(x_test = new_structures_embedding)
        predicted_probabilities = reference_clustering.predict_probabilities(x_test = new_structures_embedding)

        # Plot the probabilities
        reference_clustering.plot_probabilities(
            embedding=new_structures_embedding,
            probabilities=predicted_probabilities,
            labels=predicted_labels,
            title="BHMC Training Structures Cluster Probabilities",
            save_path="figures/bhmc_training_cluster_probabilities.png",
        )

        # Flag Novel structures based on a 25% probability threshold
        novel_indices = reference_clustering.flag_novel_structures(
            embedding_new = new_structures_embedding,
            threshold_percentile = 25.0
        )
        if len(novel_indices) == 0:
            self._log("No novel structures detected — keeping the existing reference model.")
            return reference_clustering, None

        # Combine Raw Features and Feature Matrix for the new structures
        novel_feature_mat = feature_mat[novel_indices]
        self._log(f"Novel structures feature mat shape: {novel_feature_mat.shape}")
        self._log(f"Combining with initial feature matrix shape: {feature_mat_init.shape}")
        combined_feature_mat = np.vstack([feature_mat_init, novel_feature_mat])

        # Keep the molecule list aligned 1:1 with combined_feature_mat's rows, so
        # cluster representative row indices can be resolved back to actual
        # structures after retraining (instead of being left as bare row indices).
        combined_mols = None
        if mols_init is not None:
            if len(mols_init) != feature_mat_init.shape[0]:
                raise ValueError(
                    f"mols_init has {len(mols_init)} entries but feature_mat_init has "
                    f"{feature_mat_init.shape[0]} rows — they must be aligned 1:1."
                )
            novel_mols = [mols[i] for i in novel_indices]
            combined_mols = list(mols_init) + novel_mols
        else:
            self._log(
                "mols_init not provided — cluster representatives will not be "
                "resolvable to Molecule objects.",
                level="warning",
            )

        # Retrain on the combined pool using the same recipe as initialize_from_xyz.
        n_clusters = (
            len(np.unique(reference_clustering.labels))
            if reference_clustering.labels is not None
            else 1
        )
        clustering_method = self.config.clustering_method
        classifier_backend = self.config.classifier_backend
        classifier_kwargs = self.config.classifier_kwargs

        self._log(
            f"\nRetraining clustering + {classifier_backend.upper()} classifier on "
            f"{combined_feature_mat.shape[0]} combined structures..."
        )
        updated_clustering, _, _, _ = Clustering.run_clustering_and_classifier_pipeline(
            feature_matrix=combined_feature_mat,
            n_clusters=n_clusters,
            molecules=combined_mols,
            clustering_method=clustering_method,
            classifier_backend=classifier_backend,
            classifier_kwargs=classifier_kwargs,
            metric=reference_clustering.metric,
            compute_representatives=False,
            logger=self.logger,
            log_fn=self._log,
            embedding_plot_path="figures/bhmc_retrained_umap_clusters.png",
            embedding_plot_title="Retrained UMAP Embedding with Cluster Labels",
            classifier_eval_dir="figures/classifier_evaluation_retrained",
        )
        # Obtain new representatives for the updated clustering
        # TODO: The main loop does not consider the energies right now so the "lowest_energy" gets an error this needs fixing
        representatives = updated_clustering.get_cluster_representatives(method="centroid")

        if updated_clustering.molecules is None:
            self._log(
                "Updated clustering has no tracked molecules (mols_init was not "
                "provided) — returning representative row indices cannot be resolved.",
                level="warning",
            )
            return updated_clustering, None

        representative_mols = [updated_clustering.molecules[idx] for idx in representatives.values()]
        return updated_clustering, representative_mols


    def analyze_sampling_results(
            self,
            reference_clustering: Clustering,
            structures: Optional[List[Tuple[Molecule, float]]] = None,
        ) -> Tuple[Optional[List[Molecule]], Optional[np.ndarray], Optional[Clustering], Optional[List[Molecule]]]:
        """
        Analyze the structures collected during the BHMC sampling phase against the
        current reference model (either the retrained model from
        analyze_training_results, or the original init model if no improvement was
        found there).

        Workflow:
          1. Filter high-energy outliers
          2. Featurize with SOAP
          3. Embed into the reference UMAP space
          4. Predict labels/probabilities with the reference classifier
          5. Flag "novel" structures (the ones the classifier is least confident about)
          6. Structures the classifier confidently recognized are kept as-is under
             their predicted label — they're already understood by the reference
             model, so there is nothing to recluster for them.
          7. Novel structures are re-clustered from scratch with
             Clustering.run_clustering_pipeline (cleaning + UMAP + clustering,
             but no classifier fit — there's no reason to trust a freshly fit
             classifier on a handful of outlier structures), then the
             lowest-energy structure per new cluster is taken as its representative.

        Args:
            reference_clustering: The Clustering instance to embed/predict against.
            structures:           (molecule, energy) pairs to analyze; defaults to
                                   self.accepted_structures.

        Returns:
            (confirmed_mols, confirmed_labels, novel_clustering, novel_representative_mols)
            confirmed_mols/confirmed_labels are the non-novel structures and their
            classifier-predicted labels (None if there was nothing to analyze).
            novel_clustering/novel_representative_mols are the reclustered novel
            pool and its lowest-energy representatives (None if no novel
            structures were found).
        """
        self.logger.section("Analysis of BHMC Sampling Run")
        if structures is None:
            structures = self.accepted_structures

        if not structures:
            self._log("No structures to analyse.", level="warning")
            return None, None, None, None

        mols = [m for m, _ in structures]
        energies = [e for _, e in structures]
        energy_arr = np.array(energies)
        z_e = (energy_arr - energy_arr.mean()) / (energy_arr.std() + 1e-12)
        mols     = [m for m, z in zip(mols, z_e)     if z <= 3.0]
        energies = [e for e, z in zip(energies, z_e) if z <= 3.0]
        self._log(f"After energy Z-score filter: {len(mols)} structures remain")
        if not mols:
            self._log("All structures filtered as high-energy outliers.", level="warning")
            return None, None, None, None

        # Featurize with SOAP
        featurizer = Featurizer(FeaturizerConfig(descriptor_type="SOAP"))
        feature_mat = featurizer.build_feature_matrix(mols, energies=None, include_hbonds=False)
        self._log(f"SOAP features: {feature_mat.shape[0]} × {feature_mat.shape[1]}")

        # Embed the new structures using the reference_clustering object
        new_structures_embedding = reference_clustering.embed_new_structures(x_new=feature_mat)
        # Predict cluster labels using the classifier model in the reference clustering
        predicted_labels = reference_clustering.predict_labels(x_test=new_structures_embedding)
        predicted_probabilities = reference_clustering.predict_probabilities(x_test=new_structures_embedding)

        # Plot the probabilities
        reference_clustering.plot_probabilities(
            embedding=new_structures_embedding,
            probabilities=predicted_probabilities,
            labels=predicted_labels,
            title="BHMC Sampling Structures Cluster Probabilities",
            save_path="figures/bhmc_sampling_cluster_probabilities.png",
        )

        # Flag Novel structures based on a 25% probability threshold
        novel_indices = reference_clustering.flag_novel_structures(
            embedding_new=new_structures_embedding,
            threshold_percentile=25.0
        )

        # Everything the classifier confidently recognized stays under its
        # predicted label — no reclustering needed for what the model already understands.
        novel_set = set(novel_indices.tolist())
        confirmed_indices = [i for i in range(len(mols)) if i not in novel_set]
        confirmed_mols = [mols[i] for i in confirmed_indices]
        confirmed_labels = predicted_labels[confirmed_indices]
        self._log(f"{len(confirmed_mols)} structures confidently recognized by the reference classifier — kept as-is.")

        if len(novel_indices) == 0:
            self._log("No novel structures detected in sampling run.")
            return confirmed_mols, confirmed_labels, None, None

        self._log(f"{len(novel_indices)} novel structures detected — reclustering them independently.")
        novel_mols = [mols[i] for i in novel_indices]
        novel_energies = [energies[i] for i in novel_indices]
        novel_feature_mat = feature_mat[novel_indices]

        n_novel_clusters = (
            len(np.unique(reference_clustering.labels)) if reference_clustering.labels is not None else 1
        )

        # A small novel pool doesn't warrant one cluster per reference cluster (or
        # worse, one cluster per point) — fall back to a fixed ratio of the novel
        # pool instead, so points aren't spread too thin across clusters.
        if len(novel_mols) < n_novel_clusters:
            ratio_clusters = max(1, round(len(novel_mols) * self.config.novel_cluster_ratio))
            self._log(
                f"Not enough novel structures ({len(novel_mols)}) to form "
                f"{n_novel_clusters} clusters — using a {self.config.novel_cluster_ratio:.0%} ratio "
                f"fallback instead: {ratio_clusters} cluster(s)."
            )
            n_novel_clusters = ratio_clusters



        novel_clustering, _, _, novel_rep_indices = Clustering.run_clustering_pipeline(
            feature_matrix=novel_feature_mat,
            n_clusters=n_novel_clusters,
            molecules=novel_mols,
            energies=novel_energies,
            clustering_method=self.config.clustering_method,
            metric=reference_clustering.metric,
            compute_representatives=True,
            representative_method="lowest_energy",
            logger=self.logger,
            log_fn=self._log,
            embedding_plot_path="figures/bhmc_sampling_novel_umap_clusters.png",
            embedding_plot_title="Novel Sampling Structures — UMAP Embedding with Cluster Labels",
        )
        novel_representative_mols = [novel_clustering.molecules[idx] for idx in novel_rep_indices.values()]
        self._log(f"Identified {len(novel_representative_mols)} lowest-energy representative(s) among the novel structures.")

        return confirmed_mols, confirmed_labels, novel_clustering, novel_representative_mols



    # ----------------------------------------------------------------- plotting

    def plot_worker_trajectories(
        self,
        save_path: str = "figures/worker_energy_trajectories.png",
        trajectories: Optional[Dict] = None,
        title: str = "BHMC Worker Energy Trajectories",
    ):
        """Plot energy vs step for each worker."""
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        from pathlib import Path

        if trajectories is None:
            trajectories = self.worker_trajectories
        if not trajectories:
            self._log("No trajectories to plot.", level="warning")
            return

        Path(save_path).parent.mkdir(parents=True, exist_ok=True)

        # Use a high-contrast qualitative colormap
        num_workers = len(trajectories)
        colors = cm.get_cmap('tab20', num_workers)


        plt.figure(figsize=(12, 6))

        # Sort keys for consistent color assignment
        for i, (wid, traj) in enumerate(sorted(trajectories.items())):
            steps, energies = zip(*traj) if traj else ([], [])

            # Add transparency and slightly thinned lines
            plt.plot(
                steps, 
                energies,
                label=f"Worker {wid}",
                color=colors(i),
                alpha = 0.7,
                linewidth=1.5
            )

        plt.xlabel("Step", fontsize=12, fontweight='bold')
        plt.ylabel("Energy (Hartree)", fontsize=12, fontweight='bold')
        plt.title(title, fontsize=14, fontweight='bold', pad=15)

        # Legend placement
        plt.legend(
            bbox_to_anchor = (1.02, 1),
            loc = "upper left",
            borderaxespad = 0,
            fontsize = "small",
            ncol = 1 + (num_workers // 20)
        )

        
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

    def plot_adaptive_box_updates(
        self,
        save_path: str = "figures/adaptive_box_trajectories.png",
        box_updates_by_worker: Optional[Dict] = None,
        title: str = "Adaptive Box Control",
    ):
        """Plot box size and window acceptance rate per worker."""
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        from pathlib import Path


        if box_updates_by_worker is None:
            box_updates_by_worker = self.worker_box_updates

        valid = {wid: u for wid, u in box_updates_by_worker.items() if u}
        if not valid:
            self._log("No adaptive box updates to plot.", level="warning")
            return

        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        
        fig, (ax_size, ax_acc) = plt.subplots(2, 1, figsize=(12, 8), sharex=True, gridspec_kw={'height_ratios': [1, 1]})
        
        colors = cm.get_cmap('tab20', max(len(valid), 1))


        for i, (wid, updates) in enumerate(sorted(valid.items())):
            color = colors(i)

            steps = [u["step"] for u in updates]
            
            
            # Plot box size
            ax_size.plot(steps, [u["box_size"] for u in updates],
                         marker=".", markersize=3, alpha=0.6, label=f"Worker {wid}", color=color)
            
            # Plot Acceptance Rate
            acc_data = [(u["step"], u["window_acceptance"]) for u in updates if "window_acceptance" in u]
            if acc_data:
                s, a = zip(*acc_data)
                ax_acc.plot(s, a, marker=".", markersize=3, alpha=0.6, color=color)

        target = self.config.box_target_acceptance
        window = self.config.box_acceptance_window

        ax_acc.axhline(target, color="black", linestyle="--", linewidth=1.5, label="Target Acceptance")
        ax_acc.axhspan(target - window, target + window, color="gray", alpha=0.2, label="Acceptance Window")

        # Formatting the Axes 
        ax_size.set_ylabel("Box size (Å)", fontsize=12, fontweight='bold')
        ax_acc.set_ylabel("Window Acceptance Rate (%)", fontsize=12, fontweight='bold')
        ax_acc.set_xlabel("Step", fontsize=12, fontweight='bold')
        ax_size.set_title(f"{title} — Box Size", fontsize=14, fontweight='bold', pad=15)
        ax_acc.set_title(f"{title} — Acceptance Rate", fontsize=14, fontweight='bold', pad=15)

        # Legend Outside to prevent data overlap
        fig.legend(
            loc="center right",
            bbox_to_anchor=(1.15, 0.5),
            fontsize = "small",
            frameon=False,
            ncol = 1 + (len(valid) // 20)
        )



        plt.tight_layout(rect=[0,0,0.9,1]) # Leave space on the right for the legend
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

    def plot_worker_operator_acceptance(
        self,
        save_path: str = "figures/worker_operator_acceptance.png",
    ):
        """Bar chart of operator acceptance rates per worker."""
        import matplotlib.pyplot as plt
        from pathlib import Path

        if not self.worker_operator_acceptance:
            self._log("No operator acceptance data to plot.", level="warning")
            return
        
        # Prepare the data for plotting
        workers = sorted(self.worker_operator_acceptance.keys())
        operators = sorted(next(iter(self.worker_operator_acceptance.values())).keys())

        # Create a 2D array: rows = workers, columns = operators
        data = np.array([[self.worker_operator_acceptance[wid].get(op, 0) for op in operators] for wid in workers]) 

        # Plot as a heatmap for clarity
        plt.figure(figsize = (max(8, len(operators) * 1.5), max(6, len(workers) * 0.5)))

        plt.imshow(data, cmap='viridis', aspect='auto', vmin=0, vmax=100)

        # Add text annotations inside cells for precision
        for i in range(len(workers)):
            for j in range(len(operators)):
                plt.text(j, i, f"{data[i, j]:.1f}%", ha='center', va='center', color='white' if data[i, j] < 50 else 'black')

        plt.colorbar(label='Acceptance Rate (%)')

        # Axis formatting 
        plt.xticks(range(len(operators)), operators, rotation=45, ha='right', fontsize=10)
        plt.yticks(range(len(workers)), [f"Worker {wid}" for wid in workers], fontsize=10)
        plt.title("Operator Acceptance Rates per Worker", fontsize=14, fontweight='bold', pad=15)
        plt.tight_layout()
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
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
