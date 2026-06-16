"""  
Module for various clustering methods and analysis of the sampled molecular configurations
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
from typing import List, Tuple, Optional, Dict, Any, Union, Callable
from dataclasses import dataclass, field
from collections import defaultdict



# Imports for profiling and memory tracking
import cProfile 
import pstats
import io 
import tracemalloc
import sys 
import time
from sklearn.metrics import adjusted_rand_score

# Import the logger
from modules.logger import Logger

# Import the Submodules 
import modules.feature_analysis as fa
import modules.outlier_detection as od
import modules.cluster_vis as vis
import modules.cluster_algs as algs
import modules.dim_reduction as dr
import modules.classifiers as clfs
import modules.cluster_metrics as met





from scipy.stats import zscore
from multiprocessing import Pool, cpu_count
from sklearn.preprocessing import StandardScaler

import time


# UMAP Import 
try:
    from umap import UMAP
except ImportError:
    UMAP = None
try:
    from hdbscan import HDBSCAN
except ImportError:
    HDBSCAN = None





class Clustering:
    """ 
    Generalized Data Container for Clustering
    Manages the raw feature configuration, scaling transformation, tracking labels and centralized login
    """ 


    def __init__(self,
                 feature_matrix: np.ndarray,
                 energies: Optional[List[float]] = None,
                 molecules: Optional[List[Any]] = None,
                 metric: str = "cityblock",
                 logger: Optional[Logger] = None,
                 normalize: bool = False,
    ):
        if not isinstance(feature_matrix, np.ndarray):
            raise ValueError("Feature matrix must be a numpy array")

        # Core Data Attributes
        self.feature_matrix: np.ndarray = feature_matrix
        self.molecules: Optional[List[Any]] = molecules
        self.energies: Optional[np.ndarray] = energies
        self.logger: Optional[Logger] = logger


        # Shared pipline context attributes (dynamically populated by various methods)
        self.metric: str = metric  # Distance metric for clustering and dimensionality reduction
        self.labels: Optional[np.ndarray] = None  # Cluster labels assigned to each structure
        self.corr_matrix: Optional[np.ndarray] = None  # Correlation matrix for feature analysis
        self.feature_stats: Dict[str, Any] = {}  # Dictionary to store feature statistics like mean and std
        self.outlier_models: Dict[str, Any] = {}  # Store fitted outlier detection models for reference
        self.representatives: Dict[int,int] = None  # Indices of representative structures for each cluster
        self.embedding: Optional[np.ndarray] = None  # Dimensionality reduction embedding for visualization and analysis
        self.umap_model: Optional[Any] = None  # Fitted UMAP model for embedding new structures

        # Unsupervised State Variables 
        self.labels: Optional[np.ndarray] = None
        self.scaler: Optional[StandardScaler] = None 
        self._feature_matrix_normalized: Optional[np.ndarray] = None

        # Supervised State Variables
        self.classifier_model: Optional[Any] = None
        self.classifier_probabilities: Optional[np.ndarray] = None
        

        if normalize:
            self.normalize()

    @property
    def feature_matrix(self) -> np.ndarray:
        """ 
        Public property interface replacing the _fm() helper function.
        Always exposes the operation feature matrix (normalized if available, otherwise raw) and ensures that the normalization step is applied if requested during initialization.
        """
        if self._feature_matrix_normalized is not None:
            return self._feature_matrix_normalized
        return self.feature_matrix_raw

    @feature_matrix.setter
    def feature_matrix(self, value: np.ndarray) -> None:
        """  
        Allows direct override of the core feature matrix space.
        """
        if not isinstance(value, np.ndarray):
            raise ValueError("Feature matrix must be a numpy array")
        self.feature_matrix_raw = value

    @property
    def umap_metric(self) -> str:
        """  
        Translates the metric names to UMAP-compatible naming conventions
        """
        return "manhattan" if self.metric == "cityblock" else self.metric
    
    def log(self, msg: str, level: str = "info") -> None:
        """ 
        Centralized logging utility routing messages to the configured logger.
        Falls back to regular printing if no logger instance is provided.
        """
        if self.logger and hasattr(self.logger, level):
            log_func = getattr(self.logger, level)
            log_func(msg)
        else:
            prefix = f"[{level.upper()}]" if level else ""
            print(f"{prefix} {msg}")

    def normalize(self) -> None:
        """ 
        Fits a standard scaler to the raw feature matrix, transforms it and
        safely handles any infinite values or NaNs resulting from zero-variance features
        """
        self.scaler = StandardScaler()
        normed = self.scaler.fit_transform(self.feature_matrix_raw)
        # Replace NaNs with 0.0 and inf with large finite numbers to ensure downstream clustering algorithms can operate without errors
        self._feature_matrix_normalized = np.nan_to_num(normed)

        self.log(
            f"Feature matrix normalized: mean ~ {np.mean(self._feature_matrix_normalized):.4f}, std ~ {np.std(self._feature_matrix_normalized):.4f}"
        )

    # =========================================================
    # FEATURE ANALYSIS INTERFACE --> modules/feature_analysis.py
    # =========================================================

    def compute_feature_stats(self) -> Dict[str, Dict[str, float]]:
        """ 
        Computes summary statistics for each active descriptor column
        """
        return fa.feature_statistics(self)
    
    def compute_spearman(self, plot: bool = False, save_dir: str = "figures") -> np.ndarray:
        """  
        Computes and caches the Spearman rank correlation matrix for the active feature matrix
        """
        return fa.spearman_correlation(self, plot=plot, save_dir=save_dir)

    def prune_features(self, threshold: float = 0.8) -> List[int]:
        """  
        Eliminates highly collinear descriptor columns using upper-triangular screenshot of the Spearman correlation
        """
        return fa.filter_correlation_spearman(self, threshold=threshold)

    def filter_low_variance_features(self, threshold: float = 0.0005) -> List[int]:
        """  
        Identifies and removes features with variance below the specified threshold to reduce noise and dimensionality
        """
        return fa.filter_low_variance_features(self, threshold=threshold)

    # =========================================================
    # OUTLIER DETECTION INTERFACE --> modules/outlier_detection.py
    # =========================================================

    def detect_outliers_zscore(self, threshold: float = 3.0, prune: bool = True) -> List[int]:
        """  
        Wrapper for z-score outlier detection method in outlier_detection module
        """
        return od.z_score_outlier_detection(self, threshold=threshold, prune=prune)
    
    def detect_outliers_isolation_forest(
            self,
            contamination: float = 0.05,
            n_estimators: int = 100,
            bootstrap: bool = False,
            prune: bool = True,
            random_state: Optional[int] = None
        ) -> np.ndarray:
        """
        Wrapper for Isolation Forest outlier detection method in outlier_detection module
        """
        return od.run_isolation_forest(
            self,
            contamination=contamination,
            n_estimators=n_estimators,
            bootstrap=bootstrap,
            prune=prune,
            random_state=random_state
        )
    
    def detect_outliers_local_outlier_factor(
            self,
            contamination: float = 0.05,
            n_neighbors: int = 20,
            prune: bool = True,
        ) -> np.ndarray:
        """  
        Wrapper for Local Outlier Factor detection method in outlier_detection module
        """
        return od.run_local_outlier_factor(
            self,
            contamination=contamination,
            n_neighbors=n_neighbors,
            prune=prune
        )
    
    # =====================================================
    # DIMENSIONALITY REDUCTION INTERFACE --> modules/dim_reduction.py
    # =====================================================

    def pca(self, n_components: int = 2, **kwargs) -> np.ndarray:
        """  
        Wrapper for PCA dimensionality reduction method in dim_reduction module
        """
        return dr.pca(self, n_components=n_components, **kwargs)
    
    def kernel_based_pca(self, n_components: int = 2, kernel: str = "rbf", gamma: Optional[float] = None, n_jobs: int = -1) -> np.ndarray:
        """  
        Wrapper for Kernel PCA dimensionality reduction method in dim_reduction module
        """
        return dr.kernel_based_pca(self, n_components=n_components, kernel=kernel, gamma=gamma, n_jobs=n_jobs)
    
    def tsne(self, n_components: int = 2, perplexity: float = 30.0, n_iter: int = 1000, random_state: Optional[int] = None) -> np.ndarray:
        """  
        Wrapper for t-SNE dimensionality reduction method in dim_reduction module
        """
        return dr.tsne(self, n_components=n_components, perplexity=perplexity, n_iter=n_iter, random_state=random_state)
    
    def umap(self, n_components: int = 2, n_neighbors: int = 15, min_dist: float = 0.1, random_state: Optional[int] = None) -> np.ndarray:
        """  
        Wrapper for UMAP dimensionality reduction method in dim_reduction module
        """
        if UMAP is None:
            raise ImportError("UMAP is not installed. Run `pip install umap-learn` to use this method.")
        return dr.umap(self, n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist, random_state=random_state)
    
    def umap_densmap(self, n_components: int = 2, n_neighbors: int = 15, min_dist: float = 0.1, random_state: Optional[int] = None) -> np.ndarray:
        """  
        Wrapper for UMAP DensMAP dimensionality reduction method in dim_reduction module
        """
        if UMAP is None:
            raise ImportError("UMAP is not installed. Run `pip install umap-learn` to use this method.")
        return dr.umap_densmap(self, n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist, random_state=random_state)
    
    def umap_metric_learning(self, n_components: int, x_train: np.ndarray, y_train: np.ndarray, n_neighbors: int) -> np.ndarray:
        """  
        Wrapper for UMAP metric learning dimensionality reduction method in dim_reduction module
        """
        if UMAP is None:
            raise ImportError("UMAP is not installed. Run `pip install umap-learn` to use this method.")
        return dr.umap_metric_learning(self, n_components=n_components, x_train=x_train, y_train=y_train, n_neighbors=n_neighbors)
    
    def embed_new_structures(self, x_new: np.ndarray) -> np.ndarray:
        """
        Wrapper for embedding new structures into an existing umap model
        """
        if UMAP is None:
            raise ImportError("UMAP is not installed. Run `pip install umap-learn` to use this method.")
        return dr.embed_new_structures(self, x_new)
    
    def flag_novel_structures(
        self,
        embedding_new: np.ndarray,
        classifier: Optional[Any] = None,
        threshold_percentile: float = 25.0
        ) -> np.ndarray:
        """ 
        Wrapper for the flag_novel_structures function in the dim_reduction module, which identifies novel structures based on their distance in the embedding space to the nearest cluster centers or classifier decision boundaries.
        """
        return dr.flag_novel_structures(self, embedding_new, classifier=classifier, threshold_percentile=threshold_percentile)
    

    # ==================================================
    # VISUALIZATION INTERFACE --> modules/cluster_vis.py
    # ==================================================
    def plot_embedding(self, embedding: np.ndarray, title: str = "Embedding", labels: Optional[np.ndarray] = None, save_path: Optional[str] = None) -> None:
        """ 
        Wrapper for the plot_embedding function in the cluster_vis module
        """
        vis.plot_embedding(self, embedding, title=title, labels=labels, save_path=save_path)

    def plot_probabilities(self, embedding: np.ndarray, probabilities: np.ndarray, labels: Optional[np.ndarray] = None, title: str = "Classifier Assignment Probabilities", save_path: Optional[str] = None) -> None:
        """  
        Wrapper for the plot_probabilities function in the cluster_vis module
        """
        vis.plot_probabilities(self, embedding, probabilities, labels=labels, title=title, save_path=save_path)
        

    # =========================================================
    # CLUSTERING ALGORITHMS INTERFACE --> modules/cluster_algs.py
    # =========================================================
       
    def run_clustering(
            self,
            method: str = "kmeans",
            **kwargs
        ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """ 
        Wrapper for the run_clustering function in the cluster_algs module

        Kwargs are passed directly to the underlying clustering algorithm function. Supported methods include 'kmeans', 'agglomerative', 'dbscan', and 'hdbscan'.
        """
        # Save labels to cache
        self.labels, metadata = algs.run_clustering(self, method=method, **kwargs)
        return self.labels, metadata

    def get_cluster_representatives(self, labels: Optional[np.ndarray] = None, method: str = "lowest_energy") -> Dict[int, int]:
        """  
        Wrapper method that extracts one representative molecular configuration per partition group
        """
        return algs.get_cluster_representatives(self, labels=labels, method=method)


    # ========================================================
    # CLASSIFICATION INTERFACE --> modules/classifiers.py
    # ========================================================

    def train_classifier(
            self,
            model_type: str,
            x_train: np.ndarray,
            y_train: np.ndarray,
            save_path: Optional[str] = None,
            **kwargs
        ) -> Any:
        """ 
        Wrapper method that trains a classification model, saves it to the active Clustering context
        """ 
        return clfs.train_classifier(self, model_type=model_type, x_train=x_train, y_train=y_train, save_path=save_path, **kwargs)
    
    def predict_labels(self, x_test: np.ndarray) -> np.ndarray:
        """ 
        Wrapper method that uses the trained classifier model to predict labels for new data points
        """
        return clfs.predict_labels(self, x_test)
    
    def predict_probabilities(self, x_test: np.ndarray) -> np.ndarray:
        """  
        Wrapper method that uses the trained classifier model to predict class probabilities for new data points
        """
        return clfs.predict_probabilities(self, x_test)
    
    def evaluate_classifier(self, x: np.ndarray, y: np.ndarray, n_splits: int = 5, save_dir: Optional[str] = None) -> Dict[str, Any]:
        """ 
        Wrapper method that evaluates the classifier model using stratified k-fold cross-validation
        """
        return clfs.evaluate_classifier(self, x=x, y=y, n_splits=n_splits, save_dir=save_dir)
    
    # =========================================================
    # CLUSTER EVALUATION INTERFACE --> modules/cluster_metrics.py
    # =========================================================
    def evaluate_all_metrics(self, labels: Optional[np.ndarray] = None, ignore_noise: bool = True) -> Dict[str, float]:
        """  
        Wrapper method that evaluates all implemented clustering metrics in the cluster_metrics module and returns a dictionary of metric_name -> score
        """
        return met.evaluate_all_metrics(self, labels=labels, ignore_noise=ignore_noise)


    def calculate_wccs(self, labels: Optional[np.ndarray] = None) -> float:
        """ 
        Wrapper method that calculates the WCSS for each cluster
        """
        return met.calculate_wcss_per_cluster(self, labels=labels)
    
    def compute_pertubation_stability(self, clustering_fn: Callable, n_repeats: int = 10, subsample_fraction: float = 0.85, random_seed: int = 42, **clustering_kwargs) -> Dict[str, Any]:
        """  
        Wrapper method that computes the perturbation stability of the clustering using the function in cluster_metrics module
        """
        return met.compute_pertubation_stability(self, clustering_fn=clustering_fn, n_repeats=n_repeats, subsample_fraction=subsample_fraction, random_seed=random_seed, **clustering_kwargs)

    
    # ----- Evaluation and Assessment -----
    def _resolve_labels(self, labels: Optional[np.ndarray]) -> np.ndarray:
        """  
        Helper function to resolve labels for evaluation. If labels are provided, use them. Otherwise, use the stored self.labels. If neither is available, raise an error.
        """
        lbl = labels if labels is not None else self.labels
        if lbl is None:
            raise ValueError("No labels provided for evaluation. Please run clustering first or provide labels.")
        return lbl
    
    def _mask_noise(self, labels: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """  
        Helper function to return (feature_matrix, labels) with noise points (label=-1) removed for evaluation metrics that require it
        """
        mask = labels != -1
        return self._fm()[mask], labels[mask]
    

    

    
    def elbow_analysis(
        self,
        k_range: range = range(2, 11),
        n_init: int = 10,
        random_state: int = 42,
        save_path: str = "figures/elbow_analysis.png"
        ) -> Dict[int, float]:
        """
        Run k-means for each k in k_range and recort the within cluster sum of squares (inertia) to perform an elbow analysis for determining the optimal number of clusters. 
        The 'elbow' point in the curve is a heuristic for the optimal k 
        Returns a dict mapping of k->inertia for each k in k_range
        """
        inertias: Dict[int, float] = {}
        for k in k_range:
            km = KMeans(n_clusters=k, n_init=n_init, random_state=random_state)
            km.fit(self._fm())
            inertias[k] = float(km.inertia_)
            self._log(f"Elbow analysis: k={k}, inertia={km.inertia_:.4f}")
        # Plot the elbow curve
        plt.figure(figsize=(8, 5))
        plt.plot(list(inertias.keys()), list(inertias.values()), marker='o')
        plt.title('Elbow Analysis for K-Means Clustering')
        plt.xlabel('Number of Clusters (k)')
        plt.ylabel('Inertia (Within-Cluster Sum of Squares)')
        plt.xticks(list(inertias.keys()))
        plt.grid(True)
        plt.savefig(save_path)
        plt.close()
        # Return the inertia values for each k
        return inertias


    def kmeans_noise_robustness(self, n_clusters: int = 5, noise_levels = (0.0, 0.01, 0.02, 0.05), n_runs = 5, random_state: int = 42) -> Dict[float, List[float]]:
        """  
        Evaluate the robustness of the K-Means clustering with the
        effect of noise by adding different levels of Gaussian noise to the feature matrix
        """
        fm = self._fm()
        rng = np.random.default_rng(random_state)
        base_model = KMeans(
            n_clusters=n_clusters, n_init = 10, random_state=random_state
        )
        base_labels = base_model.fit_predict(fm)
        results = {}
        feature_scale = fm.std(axis=0, keepdims=True) # scale the noise to the feature scale
        feature_scale[feature_scale == 0] = 1.0 # avoid division by zero for constant features

        for sigma in noise_levels:
            aris = []
            for run in range(n_runs):
                noise = rng.normal(loc=0.0, scale=sigma, size=fm.shape) * feature_scale
                fm_noisy = fm + noise
                noisy_model = KMeans(
                    n_clusters=n_clusters, n_init = 10, random_state=random_state + run + 1
                )
                noisy_labels = noisy_model.fit_predict(fm_noisy)
                aris.append(adjusted_rand_score(base_labels, noisy_labels))
            results[sigma] = {
                "mean_ari": float(np.mean(aris)),
                "std_ari": float(np.std(aris))
            }
        
        return base_labels, results

    
    def silhouette_analysis(
        self,
        k_range: range = range(2, 11),
        n_init: int = 10,
        random_state: int = 42,
        save_path: str = "figures/silhouette_analysis.png"
        )  -> int:
        """
        Run k-means for each k in k_range and record the silhouette score
        Plot score vs k, then marks the best k in red and returns that best k
        Use this alongside the elbow_analysis() to choose k confidently.
        """
        from sklearn.metrics import silhouette_score as _sk_silhouette_score
        scores: Dict[int, float] = {}
        for k in k_range:
            km = KMeans(n_clusters=k, n_init=n_init, random_state=random_state)
            lbl = km.fit_predict(self._fm())
            scores[k] = _sk_silhouette_score(self._fm(), lbl)
            self._log(f"Silhouette analysis: k={k}, silhouette_score={scores[k]:.4f}") 
        
        best_k = max(scores, key=scores.__getitem__)
        plt.figure(figsize=(8, 5))
        plt.plot(list(scores.keys()), list(scores.values()), marker='o', linewidth=2, color='blue')
        plt.scatter(best_k, scores[best_k], color='red', s=100, label=f'Best k={best_k}')
        plt.title('Silhouette Analysis for K-Means Clustering')
        plt.xlabel('Number of Clusters (k)')
        plt.ylabel('Silhouette Score')
        plt.xticks(list(scores.keys()))
        plt.grid(True)
        plt.legend()
        plt.savefig(save_path)
        plt.close()
        return best_k, scores 
    
    # ------ Hyperparameter Estimation for Density-Based Clustering Methods ----
    def estimate_dbscan_eps(self, min_samples, metric: Optional[str] = None):
        """ 
        Calculates the k = min_samples nearest neighbor distances and plots this distances from each point to its k-th nearest neighbor sorted in ascending order. The 'elbow' point in this curve is a heuristic for the optimal eps parameter for DBSCAN.
        """
        from sklearn.neighbors import NearestNeighbors
        from kneed import KneeLocator
        neigh = NearestNeighbors(n_neighbors=min_samples, metric=metric or self.metric)
        # Fit the model to the feature matrix and compute the distances to the nearest neighbors
        neighbors_fit = neigh.fit(self._fm())

        # Get the distances to the nearest neighbors
        distances, indices = neighbors_fit.kneighbors(self._fm())
        # Sort the distances to find the nee
        # Take the last column which corresponds to the distance to the k-th nearest neighbor
        k_distances = np.sort(distances[:, -1])
        # Plot the k-distance graph

        # Determine the knee
        kl = KneeLocator(range(len(k_distances)), k_distances, curve='convex', direction='increasing')
        optimal_eps = k_distances[kl.knee] if kl.knee is not None else None
        self._log(f"Estimated optimal eps for DBSCAN: {optimal_eps:.4f} (knee at index {kl.knee})")


        plt.figure(figsize=(8, 5))
        plt.plot(k_distances, marker='o', linestyle='-', color='blue')
        plt.axvline(x=kl.knee, color='red', linestyle='--', label=f'Knee at index {kl.knee}, eps={optimal_eps:.4f}')
        plt.title(f'k-Distance Graph for DBSCAN (k={min_samples})')
        plt.xlabel('Points sorted by distance to k-th nearest neighbor')
        plt.ylabel(f'Distance to {min_samples}-th nearest neighbor')
        plt.grid(True)
        plt.savefig("figures/dbscan_k_distance_graph.png")
        plt.close()
        self.optimal_dbscan_eps = optimal_eps
        return optimal_eps

    def _algorithm_complexity_registry(self) -> Dict[str, Any]:
        """ 
        Returns the complexity functions O(n,m,k) for supported clustering algorithms.

        Values are rough operation-count models (not wall-clock time)

        Args used here:
            n = number of structures
            m = number of features
            k = number of clusters (for methods that require it)
        """
        def agglomerative(n: int, m:int, k: int, **kwargs) -> float:
            """ 
            Typically O(n^3) or O(n^2 * m) for m features

            With simple linkage one can achieve O(n^2)
            """
            return float(n * n * m)

        def kmeans(n: int, m: int, k: int, n_iter: int = 300, **kwargs) -> float:
            """ 
            O(n * m * k * n_iter) where n_iter is the number of iterations until convergence
            """
            return float(n * m * k * n_iter)
        
        def dbscan_avg(n: int, m:int, k:int, **kwargs) -> float:
            """ 
            Average case with spatial indexing asa rough model
            
            O(n * log(n) * m) for feature distance calculations with efficient indexing
            """
            return float(n * np.log2(max(n,2)) * m)
        
        def dbscan_worst(n: int, m: int, k: int, **kwargs) -> float:
            """ 
            Worst case (no spatial indexing or all points in one cluster)
            O(n^2 * m) for pairwise distance calculations
            """
            return float(n * n * m)
        
        def hdbscan_avg(n: int, m:int, k:int, **kwargs) -> float:
            """ 
            Average case pretty similar to DBSCAN with efficient spatial indexing
            """
            return float(n * np.log2(max(n,2)) * m)

        return {
            "agglomerative": agglomerative,
            "kmeans": kmeans,
            "dbscan_avg": dbscan_avg,
            "dbscan_worst": dbscan_worst,
            "hdbscan_avg": hdbscan_avg
            }

    def benchmark_operations_per_sample(
        self,
        k_values: Tuple[int, ...] = (2, 3, 5, 10, 20, 50),
        algorithms: Optional[List[str]] = None,
        n_iter_kmeans: int = 300,) -> Dict[int, List[Dict[str,float]]]:
        """
        Compute and rank operations per sample for each algorithm and each k:
        Similar Implementation see:

        https://doi.org/10.48550/arXiv.2106.12792


        Routine:
        1) Determine dataset size (n,m)
        2) For each k and algorutm compute O(n,m,k)/n
        3) Sort ascending by operations per sample
        """
        self._log("--- Benchmarking Clustering Algorithms Complexity ---")

        feature_mat = self._fm()
        n, m = feature_mat.shape
        if n == 0:
            self._log("No structures in feature matrix.")
            return {}

        registry = self._algorithm_complexity_registry()
        selected_algorithms = algorithms if algorithms is not None else list(registry.keys())

        # Validate algorithm names
        unknown = [a for a in selected_algorithms if a not in registry]
        if unknown:
            raise ValueError(f"Unknown algorithms specified: {unknown}. Valid options are: {list(registry.keys())}")

        results_by_k: Dict[int, List[Dict[str, float]]] = {}

        self._log(f"Operation Benchmark start (n={n}, m={m}, k_values={k_values})")

        for k in k_values:
            steps: List[Dict[str, float]] = []

            for algo_name in selected_algorithms:
                # Obtain the complexity function
                complexity_fn = registry[algo_name]
                ops_total = complexity_fn(n=n, m=m, k=k, n_iter=n_iter_kmeans)
                ops_per_sample = ops_total / float(n)

                steps.append({
                    "algorithm": algo_name,
                    "k": float(k),
                    "ops_total": float(ops_total),
                    "ops_per_sample": float(ops_per_sample)
                })

            steps.sort(key=lambda x: x["ops_per_sample"])
            results_by_k[int(k)] = steps

            # Write to the logger
            self._log(f"--- k={k} ---")
            for step in steps:
                self._log(f"{step['algorithm']:<15} : Total Ops = {step['ops_total']:.2e}, Ops/Sample = {step['ops_per_sample']:.2e}")

            # Search for algorithm with lowest ops/sample per k and log it
            for step in steps:
                if step["ops_per_sample"] == steps[0]["ops_per_sample"]:
                    self._log(f"Lowest Ops/Sample for k={k}: {step['algorithm']} with {step['ops_per_sample']:.2e} ops/sample")
                    break

            self._log(f"--- End of Clustering Operation Benchmark ---\n")

        # Use the highest k value to determine the most suited algorithm
        best_algorithms = []
        for k in sorted(results_by_k.keys()):
            best_algorithms.append(results_by_k[k][0]["algorithm"])
        best_algorithm_overall = max(set(best_algorithms), key=best_algorithms.count)
        self._log(f"Most suited algorithm across k values: {best_algorithm_overall}")

        return results_by_k, best_algorithm_overall

# ------ Testing and Validation of Algorithm Design for Phase A

def algorithm_test_phase_a_clustering():
    """ 
    Tests the proposed algorithm for the clustering of the phase a structures
    For this we perform the following statistical analysis

    1) Normalize the Feature Matrix --> done in init
    2) Use Spearman correlation to identify highly correlated features
    """
    from sklearn.datasets import load_digits
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, homogeneity_completeness_v_measure
    X, true_labels = load_digits(return_X_y=True)
    # This "feature matrix" is used as a test for the algorithm
    print(f"Testing algorithm on Digits datatest feature matrix with shape: {X.shape}")
    # Plot the histograms of the 5 first features to visuaize their distributions
    plt.figure(figsize=(12, 8))
    for i in range(5):
        plt.subplot(2, 3, i+1)
        sns.histplot(X[:, i], bins=20, kde=True)
        plt.title(f"Feature {i} Distribution")
    plt.tight_layout()
    plt.savefig("figures/feature_distributions.png")
    plt.close()
    Cluster = Clustering(feature_matrix=X)
    Cluster.spearman_correlation()
    # Determine which cluster algorithm to use
    results = Cluster.benchmark_operations_per_sample(k_values=(2, 3, 5, 10), algorithms=["agglomerative", "kmeans", "dbscan_avg", "hdbscan_avg"])
    print("Clustering Algorithm Complexity Benchmark:")
    for k, algos in results[0].items():
        print(f"  k={k}:")
        for algo in algos:
            print(f"    {algo['algorithm']:<15} : Total Ops = {algo['ops_total']:.2e}, Ops/Sample = {algo['ops_per_sample']:.2e}")
    print(f"Most suited algorithm across k values: {results[1]}")
 
    # Filter According to Spearman correlation with a threshold of 0.9
    selected_indices = Cluster.filter_correlation_spearman(threshold=0.9)
    print(f"Length of selected features after filtering: {len(selected_indices)}")

    # Reduce to 2D with UMAP first — DBSCAN works in this space where clusters are visible
    embedding = Cluster.umap(n_components=2, n_neighbors=15, min_dist=0.1, random_state=42)

    # Swap the feature matrix to the UMAP embedding so all subsequent methods operate on it
    Cluster._feature_matrix_normalized = embedding

    # Run benchmark again
    results = Cluster.benchmark_operations_per_sample(k_values=(2, 3, 5, 10), algorithms=["agglomerative", "kmeans", "dbscan_avg", "hdbscan_avg"])
    print("Clustering Algorithm Complexity Benchmark after UMAP embedding:") 
    for k, algos in results[0].items():
        print(f"  k={k}:")
        for algo in algos:
            print(f"    {algo['algorithm']:<15} : Total Ops = {algo['ops_total']:.2e}, Ops/Sample = {algo['ops_per_sample']:.2e}")
    print(f"Most suited algorithm across k values after UMAP embedding: {results[1]}")


    # Hyperparameter estimation in 2D embedding space (curse of dimensionality avoided)
    min_poins = 5  # standard heuristic for 2D: ~2 * n_dims + 1
    optimal_eps = Cluster.estimate_dbscan_eps(min_samples=min_poins, metric='euclidean')
    print(f"Estimated optimal eps for DBSCAN (on UMAP embedding): {optimal_eps:.4f} with min_samples={min_poins}")

    # Run DBSCAN on the embedding
    pred_labels = Cluster.dbscan_clustering(eps=optimal_eps, min_samples=min_poins)

    # Evaluate clustering against ground truth (noise points labeled -1 are excluded)
    mask = pred_labels != -1
    ari = adjusted_rand_score(true_labels[mask], pred_labels[mask])
    nmi = normalized_mutual_info_score(true_labels[mask], pred_labels[mask])
    hom, com, vmes = homogeneity_completeness_v_measure(true_labels[mask], pred_labels[mask])
    n_noise = np.sum(~mask)
    n_clusters = len(set(pred_labels[mask]))
    print(f"\n--- Clustering Evaluation (noise points excluded: {n_noise}/{len(pred_labels)}) ---")
    print(f"  Clusters found:      {n_clusters}  (true: 10)")
    print(f"  Adjusted Rand Index: {ari:.4f}  (1.0 = perfect)")
    print(f"  Norm. Mutual Info:   {nmi:.4f}  (1.0 = perfect)")
    print(f"  Homogeneity:         {hom:.4f}")
    print(f"  Completeness:        {com:.4f}")
    print(f"  V-measure:           {vmes:.4f}")

    # Plot UMAP embedding colored by DBSCAN cluster labels
    Cluster.plot_embedding(embedding, title="UMAP Embedding with DBSCAN Clusters", save_path="figures/umap_embedding_dbscan.png")

if __name__ == "__main__":
    # Initialize a full test class and try all the methods
    mock_feature_matrix = np.random.rand(100, 20)  # 100 samples, 20 features
    mock_cluster = Clustering(feature_matrix=mock_feature_matrix, energies=None, normalize=True)
    # Run the feature analysis methods
    mock_cluster.compute_feature_stats()
    mock_cluster.compute_spearman(plot=True)
    mock_cluster.prune_features(threshold=0.8)
    # Run outlier detection methods
    mock_cluster.detect_outliers_zscore(threshold=3.0)
    mock_cluster.detect_outliers_isolation_forest(contamination=0.1, prune=True)
    mock_cluster.detect_outliers_local_outlier_factor(contamination=0.1, prune=True)
    # Run clustering algorithms
    mock_cluster.run_clustering(method="kmeans", n_clusters=5, n_init=10, random_state=42)
    mock_cluster.run_clustering(method="agglomerative", n_clusters=5, linkage='ward')
    mock_cluster.run_clustering(method="dbscan", eps=0.5, min_samples=5)
    if HDBSCAN is not None:
        mock_cluster.run_clustering(method="hdbscan", min_cluster_size=5, min_samples=5)
    # Run dimensionality reduction methods
    mock_cluster.pca(n_components=2)
    mock_cluster.kernel_based_pca(n_components=2, kernel='rbf', gamma=None)
    mock_cluster.tsne(n_components=2, perplexity=30.0, n_iter=1000, random_state=42)
    if UMAP is not None:
        mock_cluster.umap(n_components=2, n_neighbors=15, min_dist=0.1, random_state=42)
        mock_cluster.umap_densmap(n_components=2, n_neighbors=15, min_dist=0.1, random_state=42)
    # Run visualization methods
    embedding = mock_cluster.umap(n_components=2, n_neighbors=15, min_dist=0.1, random_state=42)
    mock_cluster.plot_embedding(embedding, title="UMAP Embedding", save_path="figures/umap_embedding.png")
    mock_cluster.plot_probabilities(embedding, probabilities=np.random.rand(100, 5), title="Mock Classifier Probabilities", save_path="figures/mock_probabilities.png")
    # Run classification methods
    x_train = np.random.rand(80, 20)
    y_train = np.random.randint(0, 5, size=80)
    mock_cluster.train_classifier(model_type="rf", x_train=x_train, y_train=y_train)
    x_test = np.random.rand(20, 20)
    mock_cluster.predict_labels(x_test)
    mock_cluster.predict_probabilities(x_test)
    # Run cluster evaluation methods
    mock_cluster.evaluate_all_metrics(labels=np.random.randint(0, 5, size=100))
    
    # Check all vars in the class
    print("\n--- Clustering Class Variables ---")
    print(vars(mock_cluster))