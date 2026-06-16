"""Initialization of molecular clusters for BHMC sampling.

The module provides a complete initialization workflow:
1. Load structures from XYZ files
2. Fragment into connected submolecules
3. Optionally optimize each submolecule
4. Build a simulation box
5. Generate and rank candidate initial configurations
"""

import datetime
import time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import qmc

from modules.calculator import EnergyEvaluator
from modules.cluster import Clustering
from modules.cluster_algs import kmeans_clustering
from modules.featurizer import Featurizer, FeaturizerConfig
from modules.box import SimulationBox
from modules.geometry import GeometryOps
from modules.logger import Logger
from modules.molecule_class import Molecule


@dataclass
class InitializationConfig:
    """Configuration for cluster initialization and prescreening.

    Energy backend settings control both submolecule optimization and the
    prescreening step in initialize_from_xyz.  Box settings control how the
    simulation sphere/cube is sized.  Clustering settings control how the
    candidate pool is reduced to n_workers representatives.
    """

    # --- Energy backend ---
    backend: str = "psi4"          # "psi4" | "xtb" | "ase_emt" | "gpaw"
    qm_method: str = "hf"          # QM method (psi4)
    qm_basis: str = "sto-3g"       # Basis set (psi4)
    xtb_method: str = "GFN2-xTB"  # xTB method
    gpaw_mode: str = "lcao"
    gpaw_basis: str = "dzp"
    gpaw_xc: str = "B3LYP"

    # --- Clustering / classifier ---
    clustering_method: str = "kmeans"       # "kmeans" | "dbscan" | "hdbscan"
    clustering_kwargs: Optional[dict] = None  # Extra kwargs forwarded to the clustering function

    classifier_backend: str = "knn"          # "knn" | "svm"
    classifier_kwargs: Optional[dict] = None  # Extra kwargs forwarded to the classifier

    # --- Simulation box ---
    box_type: str = "sphere"              # "sphere" | "cube"
    box_scale_factor: float = 2.5        # Legacy scale factor (covalent-radius mode)
    eta_factor: float = 4.0             # Packing factor for vdW-radius box sizing
    min_distance: float = 2.0           # Minimum intermolecular distance (Å)
    max_placement_attempts: int = 1000  # Max attempts per submolecule placement

    # --- Workflow switches ---
    rmsd_threshold: float = 0.5                     # RMSD threshold for duplicate detection (Å)
    optimize_submolecules: bool = True              # Optimize each fragment before placement
    optimize_parallel: bool = True                  # Run fragment optimizations in parallel
    optimize_cluster_representatives: bool = False  # Local opt of selected representatives
    verbose: bool = True

    def rank_candidates(
        self,
        molecules: List[Molecule],
        n_keep: int,
    ) -> List[Tuple[float, Molecule]]:
        """Evaluate all candidates and return the lowest-energy n_keep as (energy, molecule) pairs."""
        evaluator = EnergyEvaluator(
            backend=self.backend,
            qm_method=self.qm_method,
            qm_basis=self.qm_basis,
            xtb_method=self.xtb_method,
            gpaw_mode=self.gpaw_mode,
            gpaw_basis=self.gpaw_basis,
            gpaw_xc=self.gpaw_xc,
        )
        scored: List[Tuple[Optional[float], Molecule]] = []
        for mol in molecules:
            try:
                e = evaluator.evaluate(mol)
                scored.append((e, mol))
            except Exception:
                continue
        scored.sort(key=lambda x: x[0])
        return scored[:n_keep]


def run_clustering_and_classifier_pipeline(
    feature_matrix: np.ndarray,
    n_clusters: int,
    molecules: Optional[List[Molecule]] = None,
    energies: Optional[List[float]] = None,
    clustering_method: str = "kmeans",
    classifier_backend: str = "knn",
    classifier_kwargs: Optional[dict] = None,
    metric: str = "cityblock",
    compute_representatives: bool = True,
    logger: Optional[Logger] = None,
    log_fn: Optional[Callable[[str], None]] = None,
    embedding_plot_path: Optional[str] = None,
    embedding_plot_title: str = "UMAP Embedding with Cluster Labels",
    classifier_eval_dir: Optional[str] = "figures/classifier_evaluation",
) -> Tuple[Clustering, np.ndarray, np.ndarray, Optional[Dict[int, int]]]:
    """Clean a feature matrix and fit a clustering + classifier pipeline on it.

    Steps: low-variance feature pruning, Isolation Forest outlier removal, 10D
    UMAP embedding, Local Outlier Factor cleanup on the embedding, the
    configured clustering algorithm, metric evaluation, perturbation
    stability, an optional labeled embedding plot, and training/evaluating
    the online structure-assignment classifier.

    This is shared by `ClusterInitializer.initialize_from_xyz` (first build,
    where `molecules`/`energies` are known and representatives are picked)
    and `BHMC.analyze_training_results` (retraining on initial + novel
    features after a BHMC run, where representatives aren't needed) so both
    paths build their `Clustering` instance identically.

    Args:
        feature_matrix:          Descriptor matrix to clean and cluster.
        n_clusters:               Target cluster count (clamped to the post-cleanup pool size).
        molecules:                 Molecules aligned to feature_matrix rows, if available.
        energies:                  Energies aligned to feature_matrix rows, if available.
        clustering_method:         "kmeans" | "agglomerative" | "dbscan" | "hdbscan".
        classifier_backend:        "knn" | "svm" | "random_forest".
        classifier_kwargs:         Extra kwargs forwarded to the classifier constructor.
        metric:                    Distance metric for clustering and classification.
        compute_representatives:   Whether to pick a lowest-energy representative per
                                    cluster (requires `energies` to be set).
        logger:                    Optional Logger forwarded to the Clustering instance.
        log_fn:                    Optional logging callback (defaults to `print`).
        embedding_plot_path:       Where to save the labeled embedding plot, or None to skip.
        embedding_plot_title:      Title for the embedding plot.
        classifier_eval_dir:       Where to save classifier cross-validation diagnostics.

    Returns:
        (clustering, embedding_10d, labels, rep_indices) — rep_indices is None
        unless compute_representatives is True.
    """
    log = log_fn or print
    if classifier_backend.lower() not in {"knn", "svm", "random_forest"}:
        raise ValueError(
            f"Unknown classifier_backend: {classifier_backend!r}. "
            "Choose 'knn', 'svm', or 'random_forest'."
        )

    clustering = Clustering(
        feature_matrix=feature_matrix,
        energies=energies,
        molecules=molecules,
        metric=metric,
        normalize=False,
        logger=logger,
    )

    log("\n ----- Feature Filtering and Outlier Removal ----- \n")
    clustering.filter_low_variance_features(threshold=0.0005)
    clustering.detect_outliers_isolation_forest(contamination=0.4, n_estimators=100)

    embedding_10d = clustering.umap(n_components=10, n_neighbors=15, min_dist=0.1)
    clustering._feature_matrix_normalized = embedding_10d  # Cache the 10D embedding for downstream steps
    clustering.detect_outliers_local_outlier_factor(n_neighbors=20, contamination=0.3)

    # Outlier removal may have pruned rows — pull the cleaned pool back from clustering
    embedding_10d = clustering._feature_matrix_normalized
    log(f"After filtering and outlier removal, {embedding_10d.shape[0]} configurations remain for clustering.")

    n_clusters = min(n_clusters, embedding_10d.shape[0])
    labels, _ = clustering.run_clustering(method=clustering_method, n_clusters=n_clusters)
    clustering.evaluate_all_metrics(labels=labels, ignore_noise=True)

    def stability_clustering_fn(X, n_clusters=n_clusters):
        return kmeans_clustering(X, n_clusters=n_clusters, random_state=42)

    clustering.compute_pertubation_stability(
        clustering_fn=stability_clustering_fn,
        n_repeats=5,
        subsample_fraction=0.85,
        random_seed=42,
    )

    if embedding_plot_path:
        clustering.plot_embedding(
            embedding_10d,
            title=embedding_plot_title,
            labels=labels,
            save_path=embedding_plot_path,
        )

    rep_indices = None
    if compute_representatives:
        rep_indices = clustering.get_cluster_representatives(method="lowest_energy")

    kwargs = classifier_kwargs or {}
    log(f"\nTraining {classifier_backend.upper()} classifier on the 10D UMAP embedding...")
    clustering.train_classifier(
        model_type=classifier_backend.lower(),
        x_train=embedding_10d,
        y_train=labels,
        **kwargs,
    )
    clustering.evaluate_classifier(
        x=embedding_10d, y=labels, n_splits=5, save_dir=classifier_eval_dir
    )

    return clustering, embedding_10d, labels, rep_indices


class ClusterInitializer:
    """End-to-end initialization of molecular clusters for BHMC sampling.

    Workflow (see initialize_from_xyz):
      1. Load an XYZ structure and fragment it into connected submolecules.
      2. Optionally optimize each fragment with the configured backend.
      3. Build a simulation box sized by vdW radii.
      4. Sample a large pool of random/Sobol/grid candidate placements.
      5. Evaluate energies, filter outliers, cluster with UMAP+KMeans, and
         return one representative per cluster as the BHMC starting points.
    """

    def __init__(self, config: InitializationConfig, logger: Optional[Logger] = None):
        """Initialise with a configuration and an optional logger.

        An EnergyEvaluator is built immediately so that submolecule
        optimization and prescreening share the same backend settings.
        """
        self.config = config
        self.logger = logger
        self.energy_evaluator = EnergyEvaluator(
            backend=self.config.backend,
            qm_method=self.config.qm_method,
            qm_basis=self.config.qm_basis,
            xtb_method=self.config.xtb_method,
            gpaw_mode=self.config.gpaw_mode,
            gpaw_basis=self.config.gpaw_basis,
            gpaw_xc=self.config.gpaw_xc,
        )

    # ------------------------------------------------------------------
    # Logging
    # ------------------------------------------------------------------

    def _log(self, msg: str, level: str = "info") -> None:
        """Print msg to stdout (if verbose) and forward to the logger."""
        if self.config.verbose:
            print(msg)
        if self.logger:
            getattr(self.logger, level)(msg)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def initialize_from_xyz(
        self,
        xyz_file: str,
        n_workers: int = 1,
        n_configurations: int = 10000,
        n_sampling_workers: int = 1,
        placing_method: str = "random",
        grid_spacing: Optional[str] = "sobol",
        n_theta: Optional[int] = None,
        n_phi: Optional[int] = None,
        n_r: Optional[int] = None,
    ) -> Tuple[List[Molecule], List[List[int]], SimulationBox, "Clustering", np.ndarray]:
        """Run the full initialization workflow from an XYZ structure.

        Args:
            xyz_file:           Path to the XYZ file.
            n_workers:          Number of representatives to return (one per BHMC worker).
            n_configurations:   Candidates generated per sampling worker for prescreening.
            n_sampling_workers: Parallel workers for configuration generation and energy eval.
            placing_method:     Placement strategy — "random", "sobol", or "grid".
            grid_spacing:       Grid point distribution — "linear", "equal_volume_grid", or "sobol".
            n_theta:            Angular θ grid points (grid mode only).
            n_phi:              Angular φ grid points (grid mode only).
            n_r:                Radial grid points (grid mode only).

        Returns:
            Tuple of (selected_molecules, submolecule_indices, simulation_box,
            clustering, raw_feature_matrix).

            `clustering` is the fitted Clustering instance from the initialization
            pipeline — it carries the trained UMAP embedding/mapper, cluster labels,
            and the trained classifier (`clustering.classifier_model`), so the whole
            reference model can be passed around and reused (e.g. for incremental
            refits during BHMC) as a single object.

            `raw_feature_matrix` is the SOAP descriptor matrix built right after
            Z-score energy filtering but before any clustering-side feature pruning
            or outlier removal — keep it around to extend or retrain the model later
            without recomputing features from scratch.
        """
        time_start = time.time()
        self._log_run_configuration(
            xyz_file, n_workers, n_configurations, n_sampling_workers,
            placing_method, grid_spacing, n_theta, n_phi, n_r,
        )

        # Step 1: Load molecule
        molecule = self._load_molecule(xyz_file)

        # Step 2: Fragment into submolecules
        submolecules = self._fragment_molecule(molecule)

        # Step 2b: Record atom indices before any optimization
        submol_indices = [submol.get_index_in_parent() for submol in submolecules]

        # Step 3: Optimize submolecules (optional)
        if self.config.optimize_submolecules:
            submolecules = self._optimize_submolecules(submolecules)

        # Step 4: Create simulation box
        simulation_box = self._create_simulation_box(submolecules)

        # Step 5: Sample candidate configurations
        initial_molecules = self._generate_random_configurations(
            submolecules=submolecules,
            simulation_box=simulation_box,
            n_configurations=n_configurations,
            placing_method=placing_method,
            grid_spacing=grid_spacing,
            n_theta=n_theta,
            n_phi=n_phi,
            n_r=n_r,
            n_sampling_workers=n_sampling_workers,
        )

        # Step 6: Energy prescreening + Z-score filtering
        scored, failed = self._prescreen_candidates(initial_molecules, n_sampling_workers)
        mols_raw, energies_raw, e_min, e_max, n_workers = self._finalize_candidate_pool(
            scored, failed, n_workers
        )
        mols, energies = self._apply_outlier_filtering(mols_raw, energies_raw)

        # Step 7: Build the feature matrix and run the clustering pipeline
        if self.logger:
            self.logger.section("Initialization of Feature Matrix and Clustering")
        feature_mat_raw = self._build_feature_matrix(mols)
        clustering, embedding_10d, labels, rep_indices = run_clustering_and_classifier_pipeline(
            feature_matrix=feature_mat_raw,
            n_clusters=n_workers,
            molecules=mols,
            energies=energies,
            clustering_method=self.config.clustering_method,
            classifier_backend=self.config.classifier_backend,
            classifier_kwargs=self.config.classifier_kwargs,
            compute_representatives=True,
            logger=self.logger,
            log_fn=self._log,
            embedding_plot_path="figures/initialization_umap_clusters.png",
            embedding_plot_title="Initialization UMAP Embedding with Cluster Labels",
        )
        mols, energies = clustering.molecules, clustering.energies

        # Step 9: Pick one representative per cluster and resolve RMSD duplicates
        selected_molecules = self._select_unique_representatives(
            rep_indices, mols, energies, labels, embedding_10d
        )

        # Step 10: Optionally locally optimize the selected representatives
        if self.config.optimize_cluster_representatives:
            selected_molecules = self._optimize_representatives(selected_molecules)
        else:
            self._log("Skipping optimization of selected representatives (disabled in config).")

        self._log_initialization_summary(
            initial_molecules, submolecules, simulation_box, n_sampling_workers,
            n_configurations, scored, failed, selected_molecules, e_min, e_max, time_start,
        )
        self._check_clustering_artifacts(clustering)

        if self.logger:
            self.logger.separator()

        return selected_molecules, submol_indices, simulation_box, clustering, feature_mat_raw

    # ------------------------------------------------------------------
    # initialize_from_xyz: step helpers
    # ------------------------------------------------------------------

    def _log_run_configuration(
        self,
        xyz_file: str,
        n_workers: int,
        n_configurations: int,
        n_sampling_workers: int,
        placing_method: str,
        grid_spacing: Optional[str],
        n_theta: Optional[int],
        n_phi: Optional[int],
        n_r: Optional[int],
    ) -> None:
        """Log the initialization header and resolved run configuration (no-op without a logger)."""
        if not self.logger:
            return
        self.logger.header("Cluster Initialization")
        self.logger.info(
            f"Initialization started at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        )
        backend_str = (
            f"{self.config.backend} ({self.config.xtb_method})"
            if self.config.backend == "xtb"
            else self.config.backend
        )
        msg = (
            f"Configuration:\n"
            f"  XYZ file:              {xyz_file}\n"
            f"  n_workers:             {n_workers}\n"
            f"  n_configurations:      {n_configurations}\n"
            f"  Total candidates:      {n_sampling_workers * n_configurations}\n"
            f"  placing_method:        {placing_method}\n"
            f"  energy_backend:        {backend_str}\n"
            f"  classifier_backend:    {self.config.classifier_backend}\n"
        )
        if placing_method == "grid":
            msg += (
                f"  grid_spacing:          {grid_spacing}\n"
                f"  n_theta={n_theta}  n_phi={n_phi}  n_r={n_r}\n"
            )
        if self.config.classifier_kwargs:
            msg += f"  classifier_kwargs:     {self.config.classifier_kwargs}\n"
        if self.config.optimize_submolecules:
            msg += "  optimize_submolecules: True\n"
        if self.config.optimize_cluster_representatives:
            msg += "  optimize_cluster_representatives: True\n"
        self.logger.info(msg)

    def _prescreen_candidates(
        self,
        initial_molecules: List[Molecule],
        n_sampling_workers: int,
    ) -> Tuple[List[Tuple[float, Molecule]], int]:
        """Evaluate the energy of every candidate, in parallel if n_sampling_workers > 1.

        Returns:
            (scored, failed) — scored is a list of (energy, molecule) pairs for
            successful evaluations, failed is the count of failed evaluations.
        """
        scored: List[Tuple[float, Molecule]] = []
        failed = 0

        if n_sampling_workers > 1:
            self._log(
                f"\nEnergy prescreening: distributing {len(initial_molecules)} candidates "
                f"across {n_sampling_workers} workers..."
            )
            all_labels = [m.atom_labels.tolist() for m in initial_molecules]
            all_coords = [m.coordinates.tolist() for m in initial_molecules]

            indices    = list(range(len(initial_molecules)))
            chunk_size = max(1, len(indices) // n_sampling_workers)
            chunks     = [indices[i:i + chunk_size] for i in range(0, len(indices), chunk_size)]

            ctx = multiprocessing.get_context("spawn")
            energy_by_idx: dict = {}
            with ProcessPoolExecutor(max_workers=n_sampling_workers, mp_context=ctx) as executor:
                futures = {
                    executor.submit(
                        ClusterInitializer._evaluate_batch_static,
                        chunk,
                        [all_labels[i] for i in chunk],
                        [all_coords[i] for i in chunk],
                        self.config.backend,
                        self.config.qm_method,
                        self.config.qm_basis,
                        self.config.xtb_method,
                        self.config.gpaw_mode,
                        self.config.gpaw_basis,
                        self.config.gpaw_xc,
                        worker_id,
                    ): worker_id
                    for worker_id, chunk in enumerate(chunks)
                }
                for future in as_completed(futures):
                    worker_id = futures[future]
                    try:
                        batch_results = future.result()
                        ok = sum(1 for _, e in batch_results if e is not None)
                        self._log(f"  Worker {worker_id}: {ok}/{len(batch_results)} evaluations successful")
                        energy_by_idx.update({idx: e for idx, e in batch_results})
                    except Exception as exc:
                        self._log(f"  Worker {worker_id} batch failed: {exc}", level="error")

            for i, mol in enumerate(initial_molecules):
                e = energy_by_idx.get(i)
                if e is not None:
                    scored.append((e, mol))
                else:
                    failed += 1
        else:
            for mol in initial_molecules:
                try:
                    e = self.energy_evaluator.evaluate(mol)
                    scored.append((e, mol))
                except Exception as exc:
                    self._log(f"Energy evaluation failed for one candidate: {exc}", level="error")
                    failed += 1

        self._log(f"\nEnergy evaluation complete: {len(scored)} successful, {failed} failed")
        return scored, failed

    def _finalize_candidate_pool(
        self,
        scored: List[Tuple[float, Molecule]],
        failed: int,
        n_workers: int,
    ) -> Tuple[List[Molecule], List[float], float, float, int]:
        """Validate the scored pool, sort by energy, and write the initialization trajectory.

        Returns:
            (mols_raw, energies_raw, e_min, e_max, n_workers) — n_workers is
            clamped down to len(scored) if fewer valid configurations were found.
        """
        if len(scored) == 0:
            raise RuntimeError("All energy evaluations failed — cannot select initial configurations.")

        if len(scored) < n_workers:
            self._log(
                f"Warning: only {len(scored)} valid configurations found (n_workers={n_workers}). "
                "Returning all valid configurations."
            )
            n_workers = len(scored)

        scored.sort(key=lambda x: x[0])

        mols_raw     = [mol for _, mol in scored]
        energies_raw = [e for e, _ in scored]
        e_min, e_max = min(energies_raw), max(energies_raw)

        if self.logger:
            self.logger.write_xyz_trajectory(
                molecules=mols_raw,
                filepath="trajectories/initialization_trajectory.xyz",
                energies=energies_raw,
            )
            self.logger.separator()
            self.logger.section("Filtering and Data Preparation for Clustering")

        return mols_raw, energies_raw, e_min, e_max, n_workers

    def _apply_outlier_filtering(
        self,
        mols_raw: List[Molecule],
        energies_raw: List[float],
    ) -> Tuple[List[Molecule], List[float]]:
        """Remove high-energy outliers via Z-score filtering and plot the before/after distributions."""
        self._log("Performing Z-score filtering to remove high-energy outliers...")
        self._log(
            f"Energy statistics before filtering: "
            f"mean={np.mean(energies_raw):.6f}, std={np.std(energies_raw):.6f}, "
            f"min={min(energies_raw):.6f}, max={max(energies_raw):.6f} Hartree"
        )
        mols, energies = self.filter_high_energy_outliers(mols_raw, energies_raw, Z_threshold=3.0)
        self._log(
            f"Energy statistics after filtering: "
            f"mean={np.mean(energies):.6f}, std={np.std(energies):.6f}, "
            f"min={min(energies):.6f}, max={max(energies):.6f} Hartree"
        )
        self._log(f"Filtered out {len(mols_raw) - len(mols)} high-energy outliers.")

        self.plot_energy_distribution_filtering(
            energies_before=energies_raw,
            energies_after=energies,
            save_path="figures/energy_distribution_filtering.png",
        )
        return mols, energies

    def _build_feature_matrix(self, mols: List[Molecule]) -> np.ndarray:
        """Build the SOAP descriptor feature matrix for the given molecules."""
        featurizer = Featurizer(FeaturizerConfig(descriptor_type="SOAP"))
        return featurizer.build_feature_matrix(mols, energies=None, include_hbonds=False)

    def _select_unique_representatives(
        self,
        rep_indices: Dict[int, int],
        mols: List[Molecule],
        energies: List[float],
        labels: np.ndarray,
        embedding_10d: np.ndarray,
    ) -> List[Molecule]:
        """Pick the lowest-energy representative per cluster and resolve RMSD near-duplicates.

        If two representatives are within `rmsd_threshold` Å of each other, the
        higher-energy one is discarded and replaced with the structure in its
        cluster closest to the UMAP centroid (excluding the discarded index).
        """
        self._log("RMSD Based Checking for Duplicate Representatives...")
        cluster_order = list(rep_indices.keys())
        cluster_labels_arr = np.asarray(labels)
        discarded_clusters = set()

        for i, cluster_i in enumerate(cluster_order):
            if cluster_i in discarded_clusters:
                continue
            idx_i = rep_indices[cluster_i]
            mol_i = mols[idx_i]

            for cluster_j in cluster_order[i + 1:]:
                if cluster_j in discarded_clusters:
                    continue
                idx_j = rep_indices[cluster_j]
                mol_j = mols[idx_j]

                rmsd = self._calculate_rmsd(mol_i.coordinates, mol_j.coordinates)
                if rmsd < self.config.rmsd_threshold:
                    # Eliminate the one with higher energy
                    if energies[idx_i] < energies[idx_j]:
                        discarded_clusters.add(cluster_j)
                        self._log(f"Discarding cluster {cluster_j} (RMSD={rmsd:.3f} Å) in favor of cluster {cluster_i}")
                    else:
                        discarded_clusters.add(cluster_i)
                        self._log(f"Discarding cluster {cluster_i} (RMSD={rmsd:.3f} Å) in favor of cluster {cluster_j}")
                        break  # No need to compare cluster_i with others if it's discarded

        selected_molecules = []
        for cluster_id in cluster_order:
            if cluster_id not in discarded_clusters:
                # Keep the original globally lowest-energy cluster configuration
                selected_molecules.append(mols[rep_indices[cluster_id]])
                continue

            # Fall back to picking closest to centroid if the representative was discarded
            self._log(f"   Applying centroid fallback for cluster {cluster_id}...")
            member_indices = np.where(cluster_labels_arr == cluster_id)[0]

            # Exclude the duplicate configuration we just threw out
            fallback_pool = [idx for idx in member_indices if idx != rep_indices[cluster_id]]

            if not fallback_pool:
                self._log(f"   No valid fallback candidates for cluster {cluster_id}. Skipping this cluster.", level="warning")
                selected_molecules.append(mols[rep_indices[cluster_id]])  # Keep the original representative
                continue

            # Compute centroid in feature space
            cluster_embeddings = embedding_10d[cluster_labels_arr == cluster_id]
            umap_centroid = np.mean(cluster_embeddings, axis=0)

            best_fallback_idx = fallback_pool[0]
            min_distance = float("inf")
            for idx in fallback_pool:
                dist = np.linalg.norm(embedding_10d[idx] - umap_centroid)
                if dist < min_distance:
                    min_distance = dist
                    best_fallback_idx = idx
            selected_molecules.append(mols[best_fallback_idx])
            self._log(f"   Selected fallback representative for cluster {cluster_id} with UMAP distance {min_distance:.3f}")

        return selected_molecules

    def _optimize_representatives(self, selected_molecules: List[Molecule]) -> List[Molecule]:
        """Locally re-optimize each selected representative geometry."""
        self._log(f"\nOptimizing {len(selected_molecules)} selected representatives...")
        optimized_molecules = list(selected_molecules)
        for i, mol in enumerate(optimized_molecules):
            self._log(f"  Representative {i + 1}/{len(optimized_molecules)}...")
            try:
                optimized_molecules[i] = self.energy_evaluator.optimize_geometry(
                    mol, optimizer="LBFGS", trajectory_fp=None
                )
                self._log("    Optimization successful")
            except Exception as exc:
                self._log(f"    Optimization failed: {exc}", level="error")
        return optimized_molecules

    def _log_initialization_summary(
        self,
        initial_molecules: List[Molecule],
        submolecules: List[Molecule],
        simulation_box: SimulationBox,
        n_sampling_workers: int,
        n_configurations: int,
        scored: List[Tuple[float, Molecule]],
        failed: int,
        selected_molecules: List[Molecule],
        e_min: float,
        e_max: float,
        time_start: float,
    ) -> None:
        """Log a summary of the completed initialization run."""
        summary_lines = [
            f"{'=' * 60}",
            "Initialization Complete!",
            f"  Total atoms:      {len(initial_molecules[0].coordinates)}",
            f"  Submolecules:     {len(submolecules)}",
            f"  Box type:         {simulation_box.box_type}",
            f"  Candidates:       {n_sampling_workers * n_configurations} "
            f"({n_sampling_workers}w × {n_configurations}), valid: {len(scored)}, failed: {failed}",
            f"  Selected:         {len(selected_molecules)}",
            f"  Energy range:     {e_min:.6f} to {e_max:.6f} Hartree",
        ]
        if simulation_box.box_type == "sphere":
            summary_lines.append(f"  Box radius:       {simulation_box.radius:.2f} Å")
        else:
            summary_lines.append(f"  Box dimensions:   {simulation_box.box_dimensions}")
        summary_lines.append(f"  Box volume:       {simulation_box.get_volume():.2f} Å³")

        self._log("\n".join(summary_lines))
        self._log(f"Total initialization time: {time.time() - time_start:.2f} s")

    def _check_clustering_artifacts(self, clustering: Clustering) -> None:
        """Warn if the returned Clustering instance is missing its trained classifier or UMAP model."""
        if getattr(clustering, "classifier_model", None) is None:
            self._log("Warning: Classifier model is not set in the clustering object.", level="warning")
        if getattr(clustering, "umap_model", None) is None:
            self._log("Warning: UMAP model is not set in the clustering object.", level="warning")

    # ------------------------------------------------------------------
    # Filtering helpers
    # ------------------------------------------------------------------

    def filter_high_energy_outliers(
        self,
        molecules: List[Molecule],
        energies: List[float],
        Z_threshold: float = 3.0,
    ) -> Tuple[List[Molecule], List[float]]:
        """Remove structures whose Z-score exceeds Z_threshold on the high-energy side.

        Args:
            molecules:   Candidate molecules.
            energies:    Corresponding energies.
            Z_threshold: Structures with Z > Z_threshold are discarded.

        Returns:
            (filtered_molecules, filtered_energies)
        """
        energy_array = np.array(energies)
        mean_e = np.mean(energy_array)
        std_e  = np.std(energy_array)
        z_scores = (energy_array - mean_e) / std_e
        filtered_molecules, filtered_energies = [], []
        for mol, energy, z in zip(molecules, energies, z_scores):
            if z <= Z_threshold:
                filtered_molecules.append(mol)
                filtered_energies.append(energy)
        return filtered_molecules, filtered_energies

    def filter_duplicates_rmsd(
        self,
        molecules: List[Molecule],
        rmsd_threshold: float = 0.5,
    ) -> Tuple[List[Molecule], List[int]]:
        """Remove structurally duplicate molecules using pairwise RMSD.

        The first occurrence of a structure is kept; subsequent structures within
        rmsd_threshold Å of any already-kept structure are discarded.

        Args:
            molecules:      Candidate molecules.
            rmsd_threshold: RMSD cutoff (Å) for duplicate detection.

        Returns:
            (unique_molecules, duplicate_indices) where duplicate_indices lists
            the positions of the removed structures in the input list.
        """
        unique_molecules = []
        duplicate_indices = []
        for i, mol in enumerate(molecules):
            is_duplicate = False
            for umol in unique_molecules:
                if ClusterInitializer._calculate_rmsd(mol.coordinates, umol.coordinates) < rmsd_threshold:
                    is_duplicate = True
                    duplicate_indices.append(i)
                    break
            if not is_duplicate:
                unique_molecules.append(mol)
        if duplicate_indices:
            self._log(
                f"Filtered out {len(duplicate_indices)} duplicate structures (RMSD < {rmsd_threshold} Å)"
            )
        return unique_molecules, duplicate_indices

    def plot_energy_distribution_filtering(
        self,
        energies_before: List[float],
        energies_after: List[float],
        save_path: Optional[str] = None,
    ) -> None:
        """Plot energy distributions before and after Z-score filtering side-by-side.

        Args:
            energies_before: Energies of all candidates before filtering.
            energies_after:  Energies after outlier removal.
            save_path:       If given, save the figure here instead of displaying.
        """
        import seaborn as sns

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for ax, energies, color, title in (
            (axes[0], energies_before, "blue",  "Before Z-Score Filtering"),
            (axes[1], energies_after,  "green", "After Z-Score Filtering"),
        ):
            sns.histplot(energies, bins=30, kde=False, ax=ax, color=color)
            sns.kdeplot(energies, ax=ax, color=color, label="KDE")
            mean_e = np.mean(energies)
            std_e  = np.std(energies)
            ax.axvline(mean_e,         color="red",    linestyle="--", label=f"Mean: {mean_e:.2f}")
            ax.axvline(mean_e + std_e, color="orange", linestyle="--", label=f"Mean+σ: {mean_e + std_e:.2f}")
            ax.axvline(mean_e - std_e, color="orange", linestyle="--", label=f"Mean−σ: {mean_e - std_e:.2f}")
            ax.set_title(f"Energy Distribution — Initialization ({title})")
            ax.set_xlabel("Energy (Hartree)")
            ax.set_ylabel("Count")
            ax.legend()

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        plt.close()

    # ------------------------------------------------------------------
    # Private workflow steps
    # ------------------------------------------------------------------

    def _load_molecule(self, xyz_file: str) -> Molecule:
        """Load and return a Molecule from an XYZ file."""
        self._log(f"\n[1/5] Loading molecule from {xyz_file}")
        with open(xyz_file, "r") as f:
            xyz_content = f.read()
        molecule = Molecule.from_xyz(xyz_content)
        self._log(f"      {len(molecule.coordinates)} atoms — types: {set(molecule.atom_labels)}")
        return molecule

    def _fragment_molecule(self, molecule: Molecule) -> List[Molecule]:
        """Fragment the molecule into covalently connected submolecules."""
        self._log("\n[2/5] Fragmenting molecule into submolecules")
        molecule.compute_bonds()
        self._log(
            f"      {len(molecule._covalent_bonds)} covalent bonds, "
            f"{len(molecule._hydrogen_bonds)} hydrogen bonds"
        )
        submolecules = molecule.fragment_by_connectivity()
        self._log(f"      {len(submolecules)} submolecule(s):")
        for i, submol in enumerate(submolecules):
            self._log(f"        [{i + 1}] {len(submol.coordinates)} atoms")
        return submolecules

    def _optimize_submolecules(self, submolecules: List[Molecule]) -> List[Molecule]:
        """Optimize all submolecules, in parallel if configured."""
        self._log(f"\n[3/5] Optimizing {len(submolecules)} submolecule(s)...")
        if self.config.optimize_parallel and len(submolecules) > 1:
            return self._optimize_submolecules_parallel(submolecules)

        optimized = []
        for i, submol in enumerate(submolecules):
            self._log(f"      Submolecule {i + 1}: {len(submol.coordinates)} atoms...")
            try:
                optimized_submol = self.energy_evaluator.optimize_geometry(
                    submol, optimizer="BFGS", write_trajectory=False
                )
                optimized.append(optimized_submol)
                self._log("        Done")
            except Exception as exc:
                self._log(f"        Failed: {exc} — using original geometry", level="error")
                optimized.append(submol)
        return optimized

    def _optimize_submolecules_parallel(self, submolecules: List[Molecule]) -> List[Molecule]:
        """Optimize each submolecule in an isolated worker process.

        Uses a 'spawn' context to avoid fork-related issues with some backends
        (e.g. psi4, GPAW).  Results arrive out-of-order and are re-sorted by
        original index so the returned list order is preserved.
        """
        n_workers = min(len(submolecules), multiprocessing.cpu_count())
        self._log(f"      Parallel optimization on {n_workers} worker(s)...")
        optimized: List[Optional[Molecule]] = [None] * len(submolecules)

        submol_data = [
            {"atom_labels": submol.atom_labels, "coordinates": submol.coordinates, "index": i}
            for i, submol in enumerate(submolecules)
        ]
        ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as executor:
            future_to_idx = {
                executor.submit(
                    self._optimize_single_static,
                    data["atom_labels"].tolist(),
                    np.asarray(data["coordinates"], dtype=np.float64),
                    self.config.backend,
                    self.config.qm_method,
                    self.config.qm_basis,
                    self.config.xtb_method,
                    self.config.gpaw_mode,
                    self.config.gpaw_basis,
                    self.config.gpaw_xc,
                    data["index"],
                ): data["index"]
                for data in submol_data
            }
            for future in as_completed(future_to_idx):
                idx = future_to_idx[future]
                try:
                    optimized[idx] = future.result()
                    self._log(f"        Submolecule {idx + 1} optimized successfully")
                except Exception as exc:
                    self._log(f"        Submolecule {idx + 1} failed: {exc}", level="error")
                    optimized[idx] = submolecules[idx]
        return optimized

    def _create_simulation_box(self, submolecules: List[Molecule]) -> SimulationBox:
        """Build a simulation box sized from the collective vdW radii of all submolecules."""
        self._log("\n[4/5] Creating simulation box")
        all_vdw_radii = []
        total_atoms   = 0
        for submol in submolecules:
            all_vdw_radii.extend(submol.vdw_radii)
            total_atoms += len(submol.coordinates)

        simulation_box = SimulationBox.from_vdw_radii(
            vdw_radii=all_vdw_radii,
            n_atoms=total_atoms,
            box_type=self.config.box_type,
            eta_factor=self.config.eta_factor,
        )
        self._log(f"      Box type: {simulation_box.box_type}")
        if simulation_box.box_type == "sphere":
            self._log(f"      Radius:   {simulation_box.radius:.2f} Å")
            self._log(f"      Volume:   {simulation_box.get_volume():.2f} Å³")
        else:
            self._log(f"      Dims:     {simulation_box.box_dimensions}")
            self._log(f"      Volume:   {simulation_box.get_volume():.2f} Å³")
        return simulation_box

    def _generate_random_configurations(
        self,
        submolecules: List[Molecule],
        simulation_box: SimulationBox,
        n_configurations: int = 1,
        placing_method: str = "random",
        grid_spacing: Optional[str] = "sobol",
        n_theta: Optional[int] = None,
        n_phi: Optional[int] = None,
        n_r: Optional[int] = None,
        n_sampling_workers: int = 1,
    ) -> List[Molecule]:
        """Generate a pool of candidate configurations inside the simulation box.

        When n_sampling_workers > 1 each worker independently generates
        n_configurations placements, so the total pool size is
        n_sampling_workers × n_configurations.

        Args:
            submolecules:       Fragments to place.
            simulation_box:     Box used for placement.
            n_configurations:   Candidates per worker.
            placing_method:     "random", "sobol", or "grid".
            grid_spacing:       Grid distribution strategy (grid mode only).
            n_theta:            θ grid points (grid mode only).
            n_phi:              φ grid points (grid mode only).
            n_r:                Radial grid points (grid mode only).
            n_sampling_workers: Number of parallel sampling workers.

        Returns:
            Flat list of all generated Molecule objects.
        """
        if n_configurations <= 0:
            raise ValueError("n_configurations must be > 0")

        method = (placing_method or "random").lower().strip()
        if method not in {"random", "sobol", "grid"}:
            raise ValueError(
                f"Invalid placing method: {placing_method!r}. Choose from 'random', 'sobol', or 'grid'."
            )

        spacing = (grid_spacing or "sobol").lower().strip()
        if method == "grid":
            if spacing not in {"linear", "equal_volume_grid", "sobol"}:
                raise ValueError(
                    f"Invalid grid spacing: {grid_spacing!r}. "
                    "Choose from 'linear', 'equal_volume_grid', or 'sobol'."
                )
            if n_theta is None or n_phi is None or n_r is None:
                n_per_submol = max(1, n_configurations // len(submolecules))
                n_theta = n_phi = n_r = int(np.ceil(n_per_submol ** (1 / 3)))
                self._log(
                    f"      Grid auto-sizing: n_theta = n_phi = n_r = {n_theta} "
                    f"(target ≥ {n_configurations} configurations)"
                )
            self._log(f"      Grid spacing: {spacing}, n_theta={n_theta}, n_phi={n_phi}, n_r={n_r}")
        else:
            spacing = n_theta = n_phi = n_r = None

        if n_sampling_workers > 1:
            self._log(
                f"\n[5/5] Sampling {n_configurations} configurations on each of "
                f"{n_sampling_workers} workers (total: {n_sampling_workers * n_configurations})..."
            )
            submol_labels = [s.atom_labels.tolist() for s in submolecules]
            submol_coords = [s.coordinates.tolist() for s in submolecules]

            ctx = multiprocessing.get_context("spawn")
            all_configs: List[Molecule] = []
            with ProcessPoolExecutor(max_workers=n_sampling_workers, mp_context=ctx) as executor:
                futures = {
                    executor.submit(
                        ClusterInitializer._generate_batch_static,
                        submol_labels,
                        submol_coords,
                        simulation_box,
                        n_configurations,
                        method,
                        spacing,
                        n_theta,
                        n_phi,
                        n_r,
                        self.config.min_distance,
                        self.config.max_placement_attempts,
                        worker_id,
                    ): worker_id
                    for worker_id in range(n_sampling_workers)
                }
                for future in as_completed(futures):
                    worker_id = futures[future]
                    try:
                        batch = future.result()
                        all_configs.extend(batch)
                        self._log(f"  Worker {worker_id}: {len(batch)} configurations generated")
                    except Exception as exc:
                        self._log(f"  Worker {worker_id} failed: {exc}", level="error")
            self._log(f"  Total configurations collected: {len(all_configs)}")
            return all_configs

        self._log(f"\n[5/5] Sampling {n_configurations} configurations sequentially...")
        configurations: List[Molecule] = []
        for _ in range(n_configurations):
            mol = self._generate_configuration(
                submolecules=submolecules,
                simulation_box=simulation_box,
                placing_method=method,
                grid_spacing=spacing,
                n_theta=n_theta,
                n_phi=n_phi,
                n_r=n_r,
            )
            configurations.append(mol)
        return configurations

    def _generate_configuration(
        self,
        submolecules: List[Molecule],
        simulation_box: SimulationBox,
        placing_method: str = "random",
        grid_spacing: Optional[str] = None,
        n_theta: Optional[int] = None,
        n_phi: Optional[int] = None,
        n_r: Optional[int] = None,
    ) -> Molecule:
        """Place all submolecules inside the simulation box and return the combined Molecule.

        Placement strategies:
          - "random" — uniform random sampling with inverse-CDF.
          - "sobol"  — low-discrepancy quasi-random Sobol sequence.
          - "grid"   — systematic partition of the sphere into angular sectors.

        Raises:
            RuntimeError: if any submolecule cannot be placed within max_placement_attempts.
        """
        if not submolecules:
            raise ValueError("No submolecules to place")

        method = (placing_method or "random").lower().strip()
        if method not in {"random", "grid", "sobol"}:
            raise ValueError(
                f"Invalid placing method: {placing_method!r}. Choose from 'random', 'grid', or 'sobol'."
            )

        sobol_engine = sobol_engine_rotation = None
        if method == "sobol":
            sobol_engine          = qmc.Sobol(d=3, scramble=True)
            sobol_engine_rotation = qmc.Sobol(d=3, scramble=True)

        partition_points = partition_perm = partition_ptr = None
        if method == "grid":
            if simulation_box.box_type != "sphere":
                raise ValueError("Grid placement is currently only implemented for spherical boxes.")
            spacing = (grid_spacing or "sobol").lower().strip()
            if spacing not in {"linear", "equal_volume_grid", "sobol"}:
                raise ValueError(
                    f"Invalid grid spacing: {grid_spacing!r}. "
                    "Choose from 'linear', 'equal_volume_grid', or 'sobol'."
                )
            n_theta = 10 if n_theta is None else int(n_theta)
            n_phi   = 10 if n_phi   is None else int(n_phi)
            n_r     = 5  if n_r     is None else int(n_r)
            if n_theta <= 0 or n_phi <= 0 or n_r <= 0:
                raise ValueError("n_theta, n_phi, and n_r must be positive integers.")
            partition_points = self.generate_partition_points(
                center=simulation_box.center,
                radius=simulation_box.radius,
                n_partitions=len(submolecules),
                n_theta=n_theta,
                n_phi=n_phi,
                n_r=n_r,
                spacing=spacing,
                sobol_scramble=True,
                sobol_seed=None,
            )
            partition_perm = [np.random.permutation(len(pts)) for pts in partition_points]
            partition_ptr  = [0] * len(submolecules)

        placed_coords = np.empty((0, 3))
        placed_atoms: List[str] = []

        for i, submol in enumerate(submolecules):
            submol_com      = np.average(submol.coordinates, axis=0, weights=submol.masses)
            centered_coords = submol.coordinates - submol_com

            placed = False
            for _ in range(self.config.max_placement_attempts):
                if method == "random":
                    position = simulation_box.get_random_position(n_points=1)

                elif method == "sobol":
                    u = sobol_engine.random(n=1)[0]
                    if simulation_box.box_type == "cube":
                        half_dims = simulation_box.box_dimensions / 2.0
                        position  = simulation_box.center + (2.0 * u - 1.0) * half_dims
                    else:  # sphere
                        r         = simulation_box.radius * (u[0] ** (1.0 / 3.0))
                        cos_theta = 2.0 * u[1] - 1.0
                        sin_theta = np.sqrt(max(0.0, 1.0 - cos_theta ** 2))
                        phi       = 2.0 * np.pi * u[2]
                        position  = simulation_box.center + r * np.array(
                            [sin_theta * np.cos(phi), sin_theta * np.sin(phi), cos_theta]
                        )

                else:  # grid
                    perm = partition_perm[i]
                    ptr  = partition_ptr[i]
                    if ptr >= len(perm):
                        break
                    position         = partition_points[i][perm[ptr]]
                    partition_ptr[i] = ptr + 1

                rotation_matrix = (
                    self._sobol_rotation_matrix(sobol_engine_rotation)
                    if method == "sobol"
                    else self._random_rotation_matrix()
                )
                new_coords = (rotation_matrix @ centered_coords.T).T + position

                if i == 0 or self._check_min_distance(new_coords, placed_coords):
                    placed_coords = np.vstack([placed_coords, new_coords])
                    placed_atoms.extend(submol.atom_labels.tolist())
                    placed = True
                    break

            if not placed:
                raise RuntimeError(
                    f"Could not place submolecule {i + 1} after "
                    f"{self.config.max_placement_attempts} attempts. "
                    "Try increasing box size or reducing min_distance."
                )

        return Molecule.from_labels_and_coords(
            atom_labels=placed_atoms,
            coordinates=placed_coords.copy(),
        )

    def _check_min_distance(self, new_coords: np.ndarray, existing_coords: np.ndarray) -> bool:
        """Return True if every atom in new_coords is ≥ min_distance away from all existing atoms."""
        if len(existing_coords) == 0:
            return True
        diff = new_coords[:, np.newaxis, :] - existing_coords[np.newaxis, :, :]
        return bool(np.all(np.linalg.norm(diff, axis=2) >= self.config.min_distance))

    # ------------------------------------------------------------------
    # Static helpers (must be picklable for ProcessPoolExecutor)
    # ------------------------------------------------------------------

    @staticmethod
    def _calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
        """Return RMSD after optimal atom-correspondence alignment (Hungarian Algorithm)."""
        return GeometryOps.compute_optimal_correspondence_rmsd(coords1, coords2)

    @staticmethod
    def _random_rotation_matrix() -> np.ndarray:
        """Return a uniformly random 3D rotation matrix via Shoemake's quaternion method."""
        u = np.random.rand(3)
        q = np.array([
            np.sqrt(1 - u[0]) * np.sin(2 * np.pi * u[1]),
            np.sqrt(1 - u[0]) * np.cos(2 * np.pi * u[1]),
            np.sqrt(u[0])     * np.sin(2 * np.pi * u[2]),
            np.sqrt(u[0])     * np.cos(2 * np.pi * u[2]),
        ])
        q0, q1, q2, q3 = q
        return np.array([
            [1 - 2 * (q2 ** 2 + q3 ** 2), 2 * (q1 * q2 - q0 * q3),     2 * (q1 * q3 + q0 * q2)],
            [2 * (q1 * q2 + q0 * q3),     1 - 2 * (q1 ** 2 + q3 ** 2), 2 * (q2 * q3 - q0 * q1)],
            [2 * (q1 * q3 - q0 * q2),     2 * (q2 * q3 + q0 * q1),     1 - 2 * (q1 ** 2 + q2 ** 2)],
        ])

    def _sobol_rotation_matrix(self, sobol_engine) -> np.ndarray:
        """Return a quasi-random rotation matrix drawn from a Sobol engine."""
        u = sobol_engine.random(n=1)[0]
        q = np.array([
            np.sqrt(1 - u[0]) * np.sin(2 * np.pi * u[1]),
            np.sqrt(1 - u[0]) * np.cos(2 * np.pi * u[1]),
            np.sqrt(u[0])     * np.sin(2 * np.pi * u[2]),
            np.sqrt(u[0])     * np.cos(2 * np.pi * u[2]),
        ])
        q0, q1, q2, q3 = q
        return np.array([
            [1 - 2 * (q2 ** 2 + q3 ** 2), 2 * (q1 * q2 - q0 * q3),     2 * (q1 * q3 + q0 * q2)],
            [2 * (q1 * q2 + q0 * q3),     1 - 2 * (q1 ** 2 + q3 ** 2), 2 * (q2 * q3 - q0 * q1)],
            [2 * (q1 * q3 - q0 * q2),     2 * (q2 * q3 + q0 * q1),     1 - 2 * (q1 ** 2 + q2 ** 2)],
        ])

    @staticmethod
    def _generate_batch_static(
        submol_labels: List[List[str]],
        submol_coords: List[List],
        simulation_box,
        n_configurations: int,
        placing_method: str,
        grid_spacing: Optional[str],
        n_theta: Optional[int],
        n_phi: Optional[int],
        n_r: Optional[int],
        min_distance: float,
        max_placement_attempts: int,
        seed: int,
    ) -> List[Molecule]:
        """Generate a batch of configurations in an isolated worker process.

        Fully self-contained (no self) so it can be pickled by ProcessPoolExecutor.
        Failed placements are silently skipped; the returned list may therefore
        be shorter than n_configurations in degenerate cases.
        """
        import numpy as np
        from scipy.stats import qmc
        from modules.molecule_class import Molecule
        from modules.init import ClusterInitializer

        np.random.seed(seed)

        submolecules = [
            Molecule.from_labels_and_coords(labels, np.asarray(coords, dtype=np.float64))
            for labels, coords in zip(submol_labels, submol_coords)
        ]

        method  = placing_method
        spacing = grid_spacing

        batch_partition_points = None
        if method == "grid":
            batch_partition_points = ClusterInitializer.generate_partition_points(
                center=simulation_box.center,
                radius=simulation_box.radius,
                n_partitions=len(submolecules),
                n_theta=n_theta,
                n_phi=n_phi,
                n_r=n_r,
                spacing=spacing,
                sobol_scramble=True,
                sobol_seed=None,
            )

        def _rotation_matrix_from_u(u):
            q = np.array([
                np.sqrt(1 - u[0]) * np.sin(2 * np.pi * u[1]),
                np.sqrt(1 - u[0]) * np.cos(2 * np.pi * u[1]),
                np.sqrt(u[0])     * np.sin(2 * np.pi * u[2]),
                np.sqrt(u[0])     * np.cos(2 * np.pi * u[2]),
            ])
            q0, q1, q2, q3 = q
            return np.array([
                [1 - 2 * (q2 ** 2 + q3 ** 2), 2 * (q1 * q2 - q0 * q3),     2 * (q1 * q3 + q0 * q2)],
                [2 * (q1 * q2 + q0 * q3),     1 - 2 * (q1 ** 2 + q3 ** 2), 2 * (q2 * q3 - q0 * q1)],
                [2 * (q1 * q3 - q0 * q2),     2 * (q2 * q3 + q0 * q1),     1 - 2 * (q1 ** 2 + q2 ** 2)],
            ])

        def _check_dist(new_coords, placed_coords):
            diff = new_coords[:, np.newaxis, :] - placed_coords[np.newaxis, :, :]
            return np.all(np.linalg.norm(diff, axis=2) >= min_distance)

        configurations: List[Molecule] = []

        for _ in range(n_configurations):
            sobol_pos = qmc.Sobol(d=3, scramble=True) if method == "sobol" else None
            sobol_rot = qmc.Sobol(d=3, scramble=True) if method == "sobol" else None

            partition_perm = partition_ptr = None
            if method == "grid":
                partition_perm = [np.random.permutation(len(pts)) for pts in batch_partition_points]
                partition_ptr  = [0] * len(submolecules)

            placed_coords = np.empty((0, 3))
            placed_atoms: List[str] = []
            success = True

            for i, submol in enumerate(submolecules):
                com      = np.average(submol.coordinates, axis=0, weights=submol.masses)
                centered = submol.coordinates - com

                placed = False
                for _ in range(max_placement_attempts):
                    if method == "random":
                        position = simulation_box.get_random_position(n_points=1)
                    elif method == "sobol":
                        u = sobol_pos.random(n=1)[0]
                        if simulation_box.box_type == "cube":
                            half     = simulation_box.box_dimensions / 2.0
                            position = simulation_box.center + (2.0 * u - 1.0) * half
                        else:
                            r        = simulation_box.radius * (u[0] ** (1.0 / 3.0))
                            cos_t    = 2.0 * u[1] - 1.0
                            sin_t    = np.sqrt(max(0.0, 1.0 - cos_t ** 2))
                            phi      = 2.0 * np.pi * u[2]
                            position = simulation_box.center + r * np.array(
                                [sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t]
                            )
                    else:  # grid
                        perm = partition_perm[i]
                        ptr  = partition_ptr[i]
                        if ptr >= len(perm):
                            break
                        position         = batch_partition_points[i][perm[ptr]]
                        partition_ptr[i] = ptr + 1

                    rot = _rotation_matrix_from_u(
                        sobol_rot.random(n=1)[0] if method == "sobol" else np.random.rand(3)
                    )
                    new_coords = (rot @ centered.T).T + position

                    if i == 0 or _check_dist(new_coords, placed_coords):
                        placed_coords = np.vstack([placed_coords, new_coords])
                        placed_atoms.extend(submol.atom_labels.tolist())
                        placed = True
                        break

                if not placed:
                    success = False
                    break

            if success:
                configurations.append(
                    Molecule.from_labels_and_coords(
                        atom_labels=placed_atoms,
                        coordinates=placed_coords.copy(),
                    )
                )

        return configurations

    @staticmethod
    def _evaluate_batch_static(
        mol_indices: List[int],
        mol_labels: List[List[str]],
        mol_coords: List[List],
        backend: str,
        qm_method: str,
        qm_basis: str,
        xtb_method: str,
        gpaw_mode: str,
        gpaw_basis: str,
        gpaw_xc: str,
        worker_id: int,
    ) -> List[Tuple[int, Optional[float]]]:
        """Evaluate energies for a batch of molecules in an isolated worker process.

        Returns a list of (original_index, energy) pairs where energy is None on failure.
        """
        import os
        import tempfile
        from pathlib import Path
        import numpy as np
        from modules.molecule_class import Molecule
        from modules.calculator import EnergyEvaluator

        os.environ["OMP_NUM_THREADS"]      = "1"
        os.environ["MKL_NUM_THREADS"]      = "1"
        os.environ["OPENBLAS_NUM_THREADS"] = "1"

        base = Path.cwd() / ".psi4_scratch"
        base.mkdir(parents=True, exist_ok=True)
        task_scratch = Path(tempfile.mkdtemp(prefix=f"eval_{worker_id}_", dir=str(base)))
        os.environ["PSI_SCRATCH"] = str(task_scratch)

        evaluator = EnergyEvaluator(
            backend=backend,
            qm_method=qm_method,
            qm_basis=qm_basis,
            xtb_method=xtb_method,
            gpaw_mode=gpaw_mode,
            gpaw_basis=gpaw_basis,
            gpaw_xc=gpaw_xc,
        )

        results: List[Tuple[int, Optional[float]]] = []
        for idx, labels, coords in zip(mol_indices, mol_labels, mol_coords):
            try:
                mol = Molecule.from_labels_and_coords(labels, np.asarray(coords, dtype=np.float64))
                results.append((idx, evaluator.evaluate(mol)))
            except Exception:
                results.append((idx, None))

        try:
            import shutil
            if task_scratch.exists():
                shutil.rmtree(task_scratch, ignore_errors=True)
        except Exception:
            pass

        return results

    @staticmethod
    def _optimize_single_static(
        atom_labels: List[str],
        coordinates: np.ndarray,
        backend: str,
        method: str,
        basis: str,
        xtb_method: str,
        gpaw_mode: str,
        gpaw_basis: str,
        gpaw_xc: str,
        task_id: int,
    ) -> Molecule:
        """Optimize one submolecule in an isolated worker process.

        Args:
            atom_labels: Atom symbols.
            coordinates: Atomic coordinates (Å).
            backend:     Calculator backend identifier.
            method:      QM method string.
            basis:       Basis set string.
            xtb_method:  xTB method string.
            gpaw_mode:   GPAW mode.
            gpaw_basis:  GPAW basis set.
            gpaw_xc:     GPAW exchange-correlation functional.
            task_id:     Worker index used for scratch-directory naming.

        Returns:
            Optimized Molecule.

        Raises:
            RuntimeError: wrapping any exception from the optimizer.
        """
        import os
        import traceback
        import tempfile
        from pathlib import Path
        from modules.molecule_class import Molecule
        from modules.calculator import EnergyEvaluator

        os.environ["OMP_NUM_THREADS"]      = "1"
        os.environ["MKL_NUM_THREADS"]      = "1"
        os.environ["OPENBLAS_NUM_THREADS"] = "1"

        base = Path.cwd() / ".psi4_scratch"
        base.mkdir(parents=True, exist_ok=True)
        task_scratch = Path(tempfile.mkdtemp(prefix=f"task_{task_id}_", dir=str(base)))
        os.environ["PSI_SCRATCH"] = str(task_scratch)

        try:
            submol = Molecule.from_labels_and_coords(
                list(atom_labels), np.asarray(coordinates, dtype=np.float64).copy()
            )
            calc = EnergyEvaluator(
                backend=backend,
                qm_method=method,
                qm_basis=basis,
                xtb_method=xtb_method,
                gpaw_mode=gpaw_mode,
                gpaw_basis=gpaw_basis,
                gpaw_xc=gpaw_xc,
            )
            return calc.optimize_geometry(submol, optimizer="BFGS")
        except Exception as exc:
            traceback.print_exc()
            raise RuntimeError(f"Optimization failed for task {task_id}: {exc}")
        finally:
            try:
                import shutil
                if task_scratch.exists():
                    shutil.rmtree(task_scratch, ignore_errors=True)
            except Exception:
                pass

    # ------------------------------------------------------------------
    # Box partitioning utilities (used by grid placement)
    # ------------------------------------------------------------------

    @staticmethod
    def partition_sphere(center, radius, n: int):
        """Split the azimuth [0, 2π) into n equal sectors.

        Args:
            center: Sphere centre (unused; kept for API compatibility).
            radius: Sphere radius (unused; kept for API compatibility).
            n:      Number of sectors.

        Returns:
            List of (theta_start, theta_end) tuples in radians.
        """
        if n <= 0:
            raise ValueError("n must be a positive integer")
        angle_step = (2 * np.pi) / n
        return [(i * angle_step, (i + 1) * angle_step) for i in range(n)]

    @staticmethod
    def generate_partition_points(
        center,
        radius: float,
        n_partitions: int,
        n_theta: int = 10,
        n_phi: int = 10,
        n_r: int = 5,
        spacing: Optional[str] = "linear",
        sobol_scramble: bool = True,
        sobol_seed: Optional[int] = None,
    ):
        """Generate candidate placement points for each azimuthal sector of a sphere.

        Each sector covers an equal slice of [0, 2π).  Points within a sector
        are distributed according to spacing:

          - "linear"            — regular grid in (r, φ, θ).
          - "equal_volume_grid" — deterministic grid with equal-volume shells and
                                  equal-area latitude bands.
          - "sobol"             — quasi-random sampling uniform in volume per sector.

        Args:
            center:         Sphere centre (3-vector).
            radius:         Sphere radius (Å).
            n_partitions:   Number of azimuthal sectors (= number of submolecules).
            n_theta:        θ points per sector.
            n_phi:          φ points per sector.
            n_r:            Radial shells per sector.
            spacing:        Distribution strategy (see above).
            sobol_scramble: Whether to scramble the Sobol engine.
            sobol_seed:     Base seed for the Sobol engine (None = random).

        Returns:
            List of length n_partitions; each element is an (N, 3) float64 array.
        """
        center  = np.asarray(center, dtype=np.float64)
        spacing = (spacing or "linear").lower().strip()

        if radius <= 0:
            raise ValueError("radius must be positive")
        if n_partitions <= 0:
            raise ValueError("n_partitions must be > 0")
        if spacing not in {"linear", "equal_volume_grid", "sobol"}:
            raise ValueError(
                f"Unknown spacing method: {spacing!r}. "
                "Choose from 'linear', 'equal_volume_grid', or 'sobol'."
            )

        slices = ClusterInitializer.partition_sphere(center, radius, n_partitions)
        all_partition_points = []

        for p_idx, (theta_start, theta_end) in enumerate(slices):
            if spacing == "sobol":
                n_points = max(1, n_theta * n_phi * n_r)
                seed     = None if sobol_seed is None else sobol_seed + p_idx
                engine   = qmc.Sobol(d=3, scramble=sobol_scramble, seed=seed)
                u        = engine.random(n=n_points)

                r       = radius * np.cbrt(u[:, 0])
                cos_phi = 2.0 * u[:, 1] - 1.0
                sin_phi = np.sqrt(np.clip(1.0 - cos_phi ** 2, 0.0, 1.0))
                theta   = theta_start + (theta_end - theta_start) * u[:, 2]

                x = center[0] + r * sin_phi * np.cos(theta)
                y = center[1] + r * sin_phi * np.sin(theta)
                z = center[2] + r * cos_phi
                all_partition_points.append(np.column_stack((x, y, z)))
                continue

            theta_values = np.linspace(theta_start, theta_end, n_theta, endpoint=False)

            if spacing == "linear":
                phi_values = np.linspace(0, np.pi, n_phi)
                r_values   = np.linspace(0, radius, n_r)
            else:  # equal_volume_grid
                j = np.arange(n_phi)
                k = np.arange(n_r)
                cos_phi_values = -1.0 + 2.0 * (j + 0.5) / n_phi
                phi_values     = np.arccos(np.clip(cos_phi_values, -1.0, 1.0))
                r_values       = radius * np.cbrt((k + 0.5) / n_r)

            partition_points = [
                [
                    center[0] + r * np.sin(phi) * np.cos(theta),
                    center[1] + r * np.sin(phi) * np.sin(theta),
                    center[2] + r * np.cos(phi),
                ]
                for r in r_values
                for phi in phi_values
                for theta in theta_values
            ]
            all_partition_points.append(np.array(partition_points))

        return all_partition_points


# =============================== Debugging and Testing ===============================

def test_initializer(
    xyz_file: str,
    method: str = "hf",
    basis: str = "cc-pvdz",
    optimize: bool = True,
    box_type: str = "sphere",
    box_scale: float = 2.5,
    min_distance: float = 2.0,
    placing_method: str = "random",
    n_configurations: int = 1,
    save_output: bool = True,
) -> None:
    """Run an end-to-end initialization sanity check and print diagnostics.

    Args:
        xyz_file:         Path to the input XYZ file.
        method:           QM method for submolecule optimization.
        basis:            Basis set for submolecule optimization.
        optimize:         Whether to optimize submolecules.
        box_type:         Simulation box shape ("sphere" or "cube").
        box_scale:        Box scale factor.
        min_distance:     Minimum intermolecular distance (Å).
        placing_method:   Placement strategy ("random", "sobol", or "grid").
        n_configurations: Number of candidate configurations to generate.
        save_output:      Write generated configurations to initial_configurations.xyz.
    """
    print(f"\n{'=' * 80}")
    print(f"Testing ClusterInitializer with {xyz_file}")
    print(f"{'=' * 80}\n")
    print(f"Method: {method}, Basis: {basis}, Optimize: {optimize}")
    print(f"Box type: {box_type}, Box scale: {box_scale}, Min distance: {min_distance} Å")
    print(f"Placement method: {placing_method}")

    config = InitializationConfig(
        qm_method=method,
        qm_basis=basis,
        optimize_submolecules=optimize,
        box_type=box_type,
        box_scale_factor=box_scale,
        min_distance=min_distance,
        verbose=True,
    )
    initializer = ClusterInitializer(config=config)
    try:
        initial_molecules, submol_indices, simulation_box, _, _ = initializer.initialize_from_xyz(
            xyz_file,
            n_configurations=n_configurations,
            placing_method=placing_method,
        )
        initial_molecule = initial_molecules[0]

        print(f"\nInitialization successful! {len(initial_molecules)} configuration(s) generated.")
        print(f"\nMolecule Information:")
        print(f"  Total atoms:    {len(initial_molecule.coordinates)}")
        print(f"  Total mass:     {np.sum(initial_molecule.masses):.2f} amu")
        print(f"  Centre of mass: "
              f"{np.average(initial_molecule.coordinates, axis=0, weights=initial_molecule.masses)}")

        print(f"\nSubmolecule Information:  count={len(submol_indices)}")
        for i, indices in enumerate(submol_indices):
            print(f"    [{i + 1}] {len(indices)} atoms")

        print(f"\nSimulation Box:")
        print(f"  Type: {simulation_box.box_type}")
        if simulation_box.box_type == "sphere":
            print(f"  Radius: {simulation_box.radius:.2f} Å")
            print(f"  Volume: {simulation_box.get_volume():.2f} Å³")
        else:
            print(f"  Dims:   {simulation_box.box_dimensions}")
            print(f"  Volume: {simulation_box.get_volume():.2f} Å³")
        all_inside = simulation_box.is_inside(initial_molecule.coordinates)
        print(f"  All atoms inside box: {'Yes' if all_inside else 'No'}")

        coords = initial_molecule.coordinates
        dists  = [
            np.linalg.norm(coords[i] - coords[j])
            for i in range(len(coords))
            for j in range(i + 1, len(coords))
        ]
        print(f"\nDistance Check:")
        print(f"  Min interatomic distance: {min(dists):.2f} Å")
        print(f"  Max interatomic distance: {max(dists):.2f} Å")

        if len(submol_indices) > 1:
            print("\nSubmolecule COM Distances:")
            coms = []
            for indices in submol_indices:
                sub_coords = initial_molecule.coordinates[indices]
                sub_masses = initial_molecule.masses[indices]
                coms.append(np.average(sub_coords, axis=0, weights=sub_masses))
            for i in range(len(coms)):
                for j in range(i + 1, len(coms)):
                    print(f"  [{i + 1}] ↔ [{j + 1}]: {np.linalg.norm(coms[i] - coms[j]):.2f} Å")

        if save_output:
            with open("initial_configurations.xyz", "w") as f:
                for idx, mol in enumerate(initial_molecules):
                    f.write(f"{len(mol.coordinates)}\nConfiguration {idx + 1}\n")
                    for label, coord in zip(mol.atom_labels, mol.coordinates):
                        f.write(f"{label} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")
            print("\nInitial configurations saved to initial_configurations.xyz")

    except Exception as e:
        print(f"\nInitialization failed: {e}")


def test_submolecule_partition_mapping(
    initial_molecules: List[Molecule],
    submol_indices: List[List[int]],
    simulation_box: SimulationBox,
    atol: float = 1e-10,
    verbose: bool = True,
) -> bool:
    """Verify that each submolecule COM lies inside its assigned azimuthal sphere sector.

    Args:
        initial_molecules: Generated cluster configurations.
        submol_indices:    Atom index groups for each submolecule.
        simulation_box:    The spherical simulation box.
        atol:              Angular tolerance (radians) at sector boundaries.
        verbose:           Print pass/fail details.

    Returns:
        True if all configurations pass, False otherwise.
    """
    if simulation_box.box_type != "sphere":
        raise ValueError("Partition mapping test is only applicable for spherical boxes.")
    if not submol_indices:
        raise ValueError("No submolecules to test")

    center     = np.asarray(simulation_box.center, dtype=np.float64)
    n_submols  = len(submol_indices)
    partitions = ClusterInitializer.partition_sphere(center, simulation_box.radius, n_submols)
    failures   = []

    for cfg_idx, mol in enumerate(initial_molecules):
        for i, atom_idx in enumerate(submol_indices):
            idx        = np.asarray(atom_idx, dtype=int)
            sub_coords = mol.coordinates[idx]
            sub_masses = mol.masses[idx]
            com        = np.average(sub_coords, axis=0, weights=sub_masses)
            rel        = com - center

            theta = np.arctan2(rel[1], rel[0])
            if theta < 0:
                theta += 2 * np.pi

            start_angle, end_angle = partitions[i]
            if i == n_submols - 1:
                ok = (theta >= start_angle - atol) and (theta <= end_angle + atol)
            else:
                ok = (theta >= start_angle - atol) and (theta < end_angle + atol)

            if not ok:
                failures.append((cfg_idx, i, theta, start_angle, end_angle))

    if verbose:
        if not failures:
            print(f"Partition mapping test passed for all {len(initial_molecules)} configurations.")
        else:
            print(f"Partition mapping test FAILED for {len(failures)} case(s):")
            for cfg_idx, sub_idx, theta, start, end in failures:
                print(
                    f"  Config {cfg_idx + 1}, Submolecule {sub_idx + 1}: "
                    f"theta={theta:.4f} not in [{start:.4f}, {end:.4f})"
                )
    return len(failures) == 0


def plot_placement_comparison(
    R: float = 5.0,
    N: int = 500,
    save_path: str = "figures/placement_comparison.png",
) -> None:
    """Produce a 2×3 figure comparing the three placement strategies on a spherical box.

    Rows    — XY projection (top view) and XZ projection (side view).
    Columns — pseudo-random, Sobol, equal-volume grid.

    Args:
        R:         Sphere radius (Å).
        N:         Number of sample points per method.
        save_path: Output figure path.
    """
    from scipy.stats import qmc

    rng = np.random.default_rng(42)

    # Pseudo-random (inverse-CDF)
    u     = rng.uniform(0, 1, N)
    phi   = rng.uniform(0, 2 * np.pi, N)
    cos_t = rng.uniform(-1, 1, N)
    sin_t = np.sqrt(1.0 - cos_t ** 2)
    r_rnd = R * np.cbrt(u)
    rand_pts = np.column_stack([r_rnd * sin_t * np.cos(phi),
                                 r_rnd * sin_t * np.sin(phi),
                                 r_rnd * cos_t])

    # Sobol quasi-random
    engine  = qmc.Sobol(d=3, scramble=True, seed=42)
    sv      = engine.random(N)
    r_sob   = R * np.cbrt(sv[:, 0])
    cos_t_s = 2.0 * sv[:, 1] - 1.0
    sin_t_s = np.sqrt(np.clip(1.0 - cos_t_s ** 2, 0.0, 1.0))
    phi_s   = 2.0 * np.pi * sv[:, 2]
    sobol_pts = np.column_stack([r_sob * sin_t_s * np.cos(phi_s),
                                  r_sob * sin_t_s * np.sin(phi_s),
                                  r_sob * cos_t_s])

    # Equal-volume grid
    n_side  = max(1, int(np.round(N ** (1.0 / 3.0))))
    grid_pts = []
    for k in range(n_side):
        r_g = R * ((k + 0.5) / n_side) ** (1.0 / 3.0)
        for i in range(n_side):
            cos_phi_g = -1.0 + 2.0 * (i + 0.5) / n_side
            sin_phi_g = np.sqrt(max(0.0, 1.0 - cos_phi_g ** 2))
            for j in range(n_side):
                theta_g = 2.0 * np.pi * j / n_side
                grid_pts.append([r_g * sin_phi_g * np.cos(theta_g),
                                  r_g * sin_phi_g * np.sin(theta_g),
                                  r_g * cos_phi_g])
    grid_pts = np.array(grid_pts)

    methods = [
        (rand_pts,  "Random Uniform\n(Inverse-CDF)",     "#2878B5"),
        (sobol_pts, "Sobol Sequence\n(Low-Discrepancy)", "#E07B39"),
        (grid_pts,  "Grid-Based\n(Deterministic)",        "#3BAA4A"),
    ]
    proj_pairs = [
        (0, 1, "x (Å)", "y (Å)", "XY projection"),
        (0, 2, "x (Å)", "z (Å)", "XZ projection"),
    ]

    fig, axes  = plt.subplots(2, 3, figsize=(13, 8.5))
    circle_t   = np.linspace(0, 2 * np.pi, 300)

    for row_idx, (xi, yi, xlabel, ylabel, proj_label) in enumerate(proj_pairs):
        for col_idx, (pts, title, color) in enumerate(methods):
            ax = axes[row_idx, col_idx]
            ax.scatter(pts[:, xi], pts[:, yi], s=7, alpha=0.55, color=color, linewidths=0)
            ax.plot(R * np.cos(circle_t), R * np.sin(circle_t),
                    color="black", linewidth=1.2, zorder=5)
            ax.set_aspect("equal")
            ax.set_xlim(-R * 1.15, R * 1.15)
            ax.set_ylim(-R * 1.15, R * 1.15)
            ax.set_xlabel(xlabel, fontsize=10)
            ax.set_ylabel(ylabel, fontsize=10)
            ax.tick_params(labelsize=9)
            if row_idx == 0:
                ax.set_title(title, fontsize=11, fontweight="bold", pad=8)
            ax.text(0.03, 0.97, proj_label,
                    transform=ax.transAxes, va="top", ha="left", fontsize=8, color="dimgray")
            ax.text(0.97, 0.97, f"N = {len(pts)}",
                    transform=ax.transAxes, va="top", ha="right", fontsize=8, color="dimgray")

    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Placement comparison figure saved to {save_path}")


if __name__ == "__main__":
    xyz_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/co2_h2o.xyz"
    plot_placement_comparison(R=5.0, N=500, save_path="figures/placement_comparison.png")
