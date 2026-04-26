"""  
Module for various clustering methods and analysis of the sampled molecular configurations
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
from typing import List, Tuple, Optional, Dict, Any, Union
from dataclasses import dataclass, field
from collections import defaultdict



# Imports for profiling and memory tracking
import cProfile 
import pstats
import io 
import tracemalloc
import sys 
import time

from sklearn.cluster import AgglomerativeClustering as SklearnAgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
from scipy.spatial.distance import cdist, pdist
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from multiprocessing import Pool, cpu_count

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

from modules.molecule_class import Molecule
from modules.geometry import GeometryOps
from modules.box import SimulationBox
from modules.logger import Logger
from modules.descriptors import FeatureExtractor

def _compute_hbond_single(args: Tuple) -> Tuple[int,float]:
    """  
    Worker function to compute the H-bond configs, top-level for pickling

    Args:
        args: Tuple of (atom_labels, coordinates, angle_threshold)
            - atom_labels: List[str]
            - coordinates: np.ndarray of shape (n_atoms, 3)
            - angle_threshold: float in degrees
            - max_distance: float in Angstroms the donor-acceptor distance threshold
    """
    atom_labels, coordinates, angle_threshold, max_distance  = args 

    # Reconstruct molecule
    mol = Molecule(name="worker_mol")
    mol.add_atoms_batch(atom_labels, coordinates)
    mol.compute_bonds()
    valid_configs = mol.get_valid_hbond_configurations(angle_threshold=angle_threshold, max_distance=max_distance)
    
    # Distance based filtering
    valid_configs = [c for c in valid_configs if c.donor_acceptor_distance is not None and c.donor_acceptor_distance <= max_distance]
    
    n_hbonds = len(valid_configs)
    avg_angle = float(np.mean([c.angle for c in valid_configs])) if valid_configs else 0.0
    avg_da_distance = float(np.mean([c.donor_acceptor_distance for c in valid_configs])) if valid_configs else 0.0
    return n_hbonds, avg_angle, avg_da_distance 


class Clustering:
    """ 
    General Clustering Class for Clustering of Molecular Structures

    Takes as an input a feature matrix of the type
    (n_structures, n_features) and performs clustering using various algorithms.
    """
    def __init__(self,
                 feature_matrix: np.ndarray,
                 logger: Optional[Logger] = None,
                 normalize: bool = True,
    ):
        if not isinstance(feature_matrix, np.ndarray):
            raise ValueError("Feature matrix must be a numpy array")
        self.feature_matrix_raw = feature_matrix
        self.logger = logger
        self.labels: Optional[np.ndarray] = None
        self._feature_matrix_normalized = None
        self.scaler: Optional[StandardScaler] = None
        if normalize:
            self._normalize()

    def _log(self, msg: str)-> None:
        if self.logger:
            self.logger.info(msg) 

    def _normalize(self) -> None:
        self._scaler = StandardScaler()
        normed = self._scaler.fit_transform(self.feature_matrix_raw)
        self._feature_matrix_normalized = np.nan_to_num(normed)  # Replace NaNs with 0.0 and inf with large finite numbers
        self._log(f"Feature matrix normalized: mean ~ {np.mean(self._feature_matrix_normalized):.4f}, std ~ {np.std(self._feature_matrix_normalized):.4f}")

    def _fm(self) -> np.ndarray:
        """ 
        Helper function to get the feature matrix (normalized if available, otherwise raw)
        """
        return self._feature_matrix_normalized if self._feature_matrix_normalized is not None else self.feature_matrix_raw
    
    # ---- Clustering ------------

    def cluster(self,method: str="kmeans", n_clusters: Optional[int]= None, **kwargs) -> None:
        """  
        Unified cluster dispatcher to call the different clustering methods
        """
        dispatch = {
            "kmeans": self.kmeans_clustering,
            "agglomerative": self.agglomerative_clustering,
            "dbscan": self.dbscan_clustering,
            "hdbscan": self.hdbscan_clustering,
        }
        if method not in dispatch:
            raise ValueError(f"Invalid clustering method: {method}. Choose from {list(dispatch.keys())}")
        fn = dispatch[method]
        if n_clusters is not None and method in ["kmeans", "agglomerative"]:
            self._log(f"Clustering with method: {method}, n_clusters: {n_clusters}")
            self.labels = fn(n_clusters=n_clusters, **kwargs)
        else:
            self._log(f"Clustering with method: {method} (no n_clusters needed)")
            self.labels = fn(**kwargs)

    def kmeans_clustering(self, n_clusters: int = 5, n_init: int = 10,
                          max_iter: int=300, random_state: int = 42) -> np.ndarray:
        """ 
        Perform K-Means clustering on the feature matrix
        
        Args:
            n_clusters: Number of clusters to form
            n_init: Number of time the k-means algorithm will be run with different centroid seeds. The final results will be the best output of n_init consecutive runs in terms of inertia.
            max_iter: Maximum number of iterations of the k-means algorithm for a single run.
            random_state: Determines random number generation for centroid initialization. Use an int to make the randomness deterministic.
        """
        km = KMeans(n_clusters=n_clusters, n_init=n_init, max_iter=max_iter, random_state=random_state)
        self.labels = km.fit_predict(self._fm())
        self._log(f"K-Means clustering completed: n_clusters={n_clusters}, inertia={km.inertia_:.4f}")
        return self.labels  
    
    def agglomerative_clustering(self, n_clusters: int = 5, linkage: str = "ward") -> np.ndarray:
        """ 
        Perform Agglomerative Clustering on the feature matrix
        
        Args:
            n_clusters: Number of clusters to form
            linkage: Which linkage criterion to use. The linkage criterion determines which distance to use between sets of observation. The algorithm will merge the pairs of cluster that minimize this criterion. 
                     - "ward": minimizes the variance of the clusters being merged (only works with Euclidean distance)
                     - "complete": uses the maximum distances between all observations of the two sets
                     - "average": uses the average of the distances of each observation of the two sets
                     - "single": uses the minimum of the distances between all observations of the two sets
        """
        ac = SklearnAgglomerativeClustering(n_clusters=n_clusters, linkage=linkage)
        self.labels = ac.fit_predict(self._fm())
        self._log(f"Agglomerative clustering completed: n_clusters={n_clusters}, linkage={linkage}")
        return self.labels

    def dbscan_clustering(self, eps: float = 0.5, min_samples: int = 5) -> np.ndarray:
        """ 
        Perform DBSCAN clustering on the feature matrix
        
        Args:
            eps: The maximum distance between two samples for one to be considered as in the neighborhood of the other. This is not a maximum bound on the distances of points within a cluster. This is the most important DBSCAN parameter to choose appropriately for your data set and distance function.
            min_samples: The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. This includes the point itself.
        """
        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        self.labels = dbscan.fit_predict(self._fm())
        n_clusters = len(set(self.labels)) - (1 if -1 in self.labels else 0)
        n_noise = list(self.labels).count(-1)
        self._log(f"DBSCAN clustering completed: eps={eps}, min_samples={min_samples}, n_clusters={n_clusters}, n_noise={n_noise}")
        return self.labels
        
    def hdbscan_clustering(self, min_cluster_size: int = 5, min_samples: Optional[int] = None) -> np.ndarray:
        """ 
        Perform HDBSCAN clustering on the feature matrix
        
        Args:
            min_cluster_size: The minimum size of clusters; single linkage splits that contain fewer points than this will be considered points "falling out" of a cluster instead of a cluster splitting into two new clusters.
            min_samples: The number of samples in a neighborhood for a point to be considered as a core point. This includes the point itself. If None, the default value is the same as min_cluster_size.
        """
        if HDBSCAN is None:
            raise ImportError("HDBSCAN is not installed. Please install it with `pip install hdbscan` to use this method.")
        hdbscan = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples)
        self.labels = hdbscan.fit_predict(self._fm())
        n_clusters = len(set(self.labels)) - (1 if -1 in self.labels else 0)
        n_noise = list(self.labels).count(-1)
        self._log(f"HDBSCAN clustering completed: min_cluster_size={min_cluster_size}, min_samples={min_samples}, n_clusters={n_clusters}, n_noise={n_noise}")
        return self.labels

    # ---- Dimensionality reduction ------------

    def pca(self, n_components: int=2) -> Tuple[np.ndarray, np.ndarray]:
        """   
        Applies a PCA dimensionality reduction

        returns (embedding, explained_variance_ratio)
        """
        p = PCA(n_components=n_components)
        embedding = p.fit_transform(self._fm())
        self._log(f"PCA completed: n_components={n_components}, explained_variance_ratio={p.explained_variance_ratio_}")
        return embedding, p.explained_variance_ratio_

    def tsne(self, n_components: int=2, perpelxity: float = 30.0, n_iter: int=1000, random_state: int = 42) -> np.ndarray:
        """  
        Applies a T-SNE Dimensionality Reduction
        """
        ts = TSNE(n_components=n_components, perplexity=perpelxity, n_iter=n_iter, random_state=random_state)
        embedding = ts.fit_transform(self._fm())
        self._log(f"T-SNE completed: n_components={n_components}, perplexity={perpelxity}, n_iter={n_iter}")
        return embedding
    
    def umap(self, n_components: int=2, n_neighbors: int=15, min_dist: float=0.1, random_state: int = 42) -> np.ndarray:
        """  
        Applies a UMAP Dimensionality Reduction
        """
        if UMAP is None:
            raise ImportError("UMAP is not installed. Please install it with `pip install umap-learn` to use this method.")
        umap = UMAP(n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist, random_state=random_state)
        embedding = umap.fit_transform(self._fm())
        self._log(f"UMAP completed: n_components={n_components}, n_neighbors={n_neighbors}, min_dist={min_dist}")
        return embedding

    # ---- Plotting of Dimensionality Reduction Methods -----
    
    def plot_embedding(self, embedding: np.ndarray, title: str = "Embedding", save_path: Optional[str] = None) -> None:
        """  
        Helper function to plot a 2D embedding with cluster labels if available
        """
        plt.figure(figsize=(8, 6))
        if self.labels is not None:
            sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=self.labels, palette='tab10', s=100, alpha=0.7)
            plt.legend(title="Cluster")
        else:
            sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], color='blue', s=100, alpha=0.7)
        plt.title(title)
        plt.xlabel("Component 1")
        plt.ylabel("Component 2")
        plt.grid(True)
        if save_path:
            plt.savefig(save_path)
        else:
            plt.savefig(f"figures/{title.replace(' ', '_').lower()}.png")
        plt.close()
    

    
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
    
    def silhouette_score(self, labels: Optional[np.ndarray] = None) -> float:
        """ 
        Determine the Sillhouette coefficient.
        This is a measrue how similar each point is to its own cluster compared vs neighboring clusters
        Objects with a high silhouette coefficient are well matched to their own cluster and poorly matched to neighboring clusters.
        """
        from sklearn.metrics import silhouette_score as _sk_silhouette_score
        lbl = self._resolve_labels(labels)
        fm, lbl = self._mask_noise(lbl)
        if len(set(lbl)) < 2:
            self._log("Silhouette score cannot be computed with less than 2 clusters (after masking noise). Returning 0.0", level="warning")
            return 0.0
        score = _sk_silhouette_score(fm, lbl)
        self._log(f"Silhouette Score: {score:.4f}")
        return score
    
    def davies_bouldin_score(self, labels: Optional[np.ndarray] = None) -> float:
        """ 
        DB = 1/n sum(max((sigma_i + sigma_j) / d(c_i, c_j))) for i != j
        where n is the number of clusters, sigma_i is the average distance of all points in cluster i to the centroid of cluster i, and d(c_i, c_j) is the distance between the centroids of clusters i and j.

        Since algorithms that produce clusters with low intra-cluster distances (high intra-cluster similarity) and high inter-cluster distances will habve a low
        Davies-Bouldin index, the clustering algorithm that produces a collection of clusters with the smallest Davies-Bouldin index is considered the best.
        """
        from sklearn.metrics import davies_bouldin_score as _sk_davies_bouldin_score
        lbl = self._resolve_labels(labels)
        fm, lbl = self._mask_noise(lbl)
        if len(set(lbl)) < 2:
            self._log("Davies-Bouldin score cannot be computed with less than 2 clusters (after masking noise). Returning inf", level="warning")
            return float('inf')
        score = _sk_davies_bouldin_score(fm, lbl)
        self._log(f"Davies-Bouldin Score: {score:.4f}")
        return score
    
    def calinski_harabasz_score(self, labels: Optional[np.ndarray] = None) -> float:
        """  
        Calinski-Harabasz Index = (Tr / (k - 1)) / (Tw / (n - k))
        where Tr is the between-clusters dispersion matrix, Tw is the within-cluster dispersion matrix, n is the number of samples, and k is the number of clusters.
        """
        from sklearn.metrics import calinski_harabasz_score as _sk_calinski_harabasz_score
        lbl = self._resolve_labels(labels)
        fm, lbl = self._mask_noise(lbl)
        if len(set(lbl)) < 2:
            self._log("Calinski-Harabasz score cannot be computed with less than 2 clusters (after masking noise). Returning 0.0", level="warning")
            return 0.0
        score = _sk_calinski_harabasz_score(fm, lbl)
        self._log(f"Calinski-Harabasz Score: {score:.4f}")
        return score
    
    def evaluate(self, labels: Optional[np.ndarray] = None) -> Dict[str, float]:
        """  
        Evaluate the clustering using multiple metrics and return a dictionary of scores
        """
        lbl = self._resolve_labels(labels)
        scores = {
            "silhouette_score": self.silhouette_score(lbl),
            "davies_bouldin_score": self.davies_bouldin_score(lbl),
            "calinski_harabasz_score": self.calinski_harabasz_score(lbl)
        }
        return scores
    
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

    # ----- Obtain cluster representatives -----
    def get_cluster_representatives(self, labels: Optional[np.ndarray] = None, method: str = "centroid") -> Dict[int, int]:
        """
        Return one representative index per cluster.

        Args:
            labels: Cluster label array. If None uses self.labels (run clustering first).
            method: "centroid" picks the point closest to the cluster mean in feature space;
                    "medoid" picks the point with the smallest average pairwise distance to all
                    cluster members (slower but more robust to non-convex clusters).

        Returns:
            Dict mapping cluster_label -> index into the original data (rows of feature_matrix).
            Noise points (label == -1 from DBSCAN/HDBSCAN) are skipped.
        """
        lbl = labels if labels is not None else self.labels
        if lbl is None:
            raise ValueError("No cluster labels available. Run clustering first or pass labels.")

        fm = self._fm()
        representatives: Dict[int, int] = {}

        for label in np.unique(lbl):
            if label == -1:
                continue
            mask = lbl == label
            cluster_features = fm[mask]
            cluster_indices = np.where(mask)[0]

            if method == "centroid":
                centroid = cluster_features.mean(axis=0)
                local_idx = int(np.argmin(np.linalg.norm(cluster_features - centroid, axis=1)))
            elif method == "medoid":
                dists = cdist(cluster_features, cluster_features)
                local_idx = int(np.argmin(np.mean(dists, axis=1)))
            else:
                raise ValueError(f"Unknown method '{method}'. Choose 'centroid' or 'medoid'.")

            representatives[int(label)] = int(cluster_indices[local_idx])

        self._log(f"Cluster representatives ({method}): {representatives}")
        return representatives

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
        k_values: Tuple[int, ...] = (2, 5, 10, 20, 50, 100),
        phase: Optional[str] = None,
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

        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        n = len(structures)
        if n == 0:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return {}
        # m = feature dimension used by clustering
        feature_mat = self.feature_matrix(normalize=True)
        m = feature_mat.shape[1]

        registry = self._algorithm_complexity_registry()
        selected_algorithms = algorithms if algorithms is not None else list(registry.keys())

        # Validate algorithm names
        unknown = [a for a in selected_algorithms if a not in registry]
        if unknown:
            raise ValueError(f"Unknown algorithms specified: {unknown}. Valid options are: {list(registry.keys())}")
        
        results_by_k: Dict[int, List[Dict[str, float]]] = {}

        self._log(f"Operation Benchmark start (phase = {phase if phase else 'all'}, n={n}, m={m}, k_values={k_values})")
        
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
        
        return results_by_k


@dataclass
class StructureData:
    """  
    Container for Molecular structure and associated properties
    """
    molecule: Molecule
    energy: float
    phase: str = "unknown"
    step: int = -1
    worker_id: int = -1
    accepted: bool = True
    dipole_vector: Optional[List[float]] = None 
    dipole_magnitude: float = 0.0
    mulliken_charges: Optional[List[float]] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

class BHMCAnalyzer:
    """   
    Analyzer for the Basin Hopping Monte Carlo results
    """
    def __init__(self, name: str = "BHMC_Analysis", submolecule_indices: Optional[List[int]] = None, logger: Optional[Logger] = None):
        self.name = name 
        self.structures: List[StructureData] = []
        self.phases: Dict[str, List[StructureData]] = defaultdict(list)
        self.submolecule_indices = submolecule_indices
        self.logger = logger
        self.labels = None # Cluster labels

        # Storage of feature matrix
        self._feature_matrix_raw: Optional[np.ndarray] = None 
        self._feature_matrix_normalized: Optional[np.ndarray] = None 
        self._non_zero_var_indices: Optional[List[int]] = None  # Track which features are kept after filtering
        self.used_features: Optional[List[str]] = None  # To keep track of which features are used


    def _log(self, msg: str, level: str = "info"):
        """ 
        Log a message if the logger is available
        """
        if self.logger:
            getattr(self.logger, level)(msg)




    def add_structure(self,
                      molecule: Molecule,
                      energy: float, 
                      phase: str = "unknown",
                      dipole_vector: Optional[List[float]] = None,
                      dipole_magnitude: float = 0.0,
                      mulliken_charges: Optional[List[float]] = None,
                      **kwargs):
        """   
        Add a structure to the analysis

        Args:
            molecule: Molecule object representing the structure
            energy: Energy of the structure
            phase: Phase of the BHMC process (e.g., "phase_a", "phase_b", "phase_c")
            dipole_vector: Optional list of 3 floats representing the dipole vector
            dipole_magnitude: Optional float representing the magnitude of the dipole moment
            **kwargs: Additional properties (e.g., step, worker_id, accepted)
        """
        structure_data = StructureData(
            molecule=molecule,
            energy=energy,
            phase=phase,
            dipole_vector = dipole_vector,
            dipole_magnitude=dipole_magnitude,
            mulliken_charges=mulliken_charges,
            metadata=kwargs
        )
        self.structures.append(structure_data)
        self.phases[phase].append(structure_data)

    def add_structures_batch(self,
                             structures: List[Tuple[Molecule, float]],
                             phase: str = "unknown"):
        """   
        Add multiple structures at once

        Args:
            structures: List of tuples containing (Molecule, energy) or (Molecule, energy, dipole_vec, dipole_magnitude)
            phase: Phase of the BHMC process for all structures
        """
        for item in structures:
            if len(item) == 4:
                mol, energy, dipole_vec, dipole_mag = item 
                self.add_structure(molecule=mol,
                                   energy=energy,
                                   phase=phase,
                                   dipole_vector=dipole_vec,
                                   dipole_magnitude=dipole_mag) 
                
            elif len(item) == 5:
                mol, energy, dipole_vec, dipole_mag, mulliken_charges = item 
                self.add_structure(molecule=mol,
                                   energy=energy,
                                   phase=phase,
                                   dipole_vector=dipole_vec,
                                   dipole_magnitude=dipole_mag,
                                   mulliken_charges=mulliken_charges)
            
            elif len(item) == 2:
                mol, energy = item
                self.add_structure(molecule=mol,
                                      energy=energy,
                                      phase=phase)
                

            else:
                raise ValueError(f"Invalid structure tuple length: {len(item)}. Expected 2 or 4.")
        self._log(f"Added batch of {len(structures)} structures to phase: {phase}. Total structures in this phase: {len(self.phases[phase])}")
        # Log Data usage of the new structures
        self._log(f"Size of new structures added: {sys.getsizeof(self.structures[-1]) * len(structures) / 1e6:.2f} MB")




    # ==================== Energy Statistics =====================
    def get_energy_statistics(self, phase: Optional[str] = None) -> Dict[str, float]:
        """   
        Compute energy statistics.
        
        Args:
            phase: Specified phase to compute statistics for. If None, computes for all structures.

        Returns:
            Dictionary with Statistics
        """
        if phase:
            energies = [s.energy for s in self.phases[phase]]
        else:
            energies = [s.energy for s in self.structures]

        if not energies:
            self._log(f"No energy statistics!", level="warning")
            return {}
        
        energies = np.array(energies)
        stats = {
                "n_structures": len(energies),
                "min": np.min(energies),
                "max": np.max(energies),
                "mean": np.mean(energies),
                "median": np.median(energies),
                "std": np.std(energies),
                "q25": np.percentile(energies, 25),
                "q75": np.percentile(energies, 75),
                "iqr": np.percentile(energies, 75) - np.percentile(energies, 25)
        }
        self._log(f"Energy Statistics (phase: {phase if phase else 'all'}):")
        self._log(f"  N Structures: {stats['n_structures']}")
        self._log(f"  Min: {stats['min']:.6f}")
        self._log(f"  Max: {stats['max']:.6f}")
        self._log(f"  Mean: {stats['mean']:.6f}")
        self._log(f"  Median: {stats['median']:.6f}")
        self._log(f"  Std: {stats['std']:.6f}")
        self._log(f"  Q25: {stats['q25']:.6f}")
        self._log(f"  Q75: {stats['q75']:.6f}")
        self._log(f"  IQR: {stats['iqr']:.6f}")
        return stats


    def compute_rmsd_matrix(self, phase: Optional[str] = None) -> np.ndarray:
        """   
        Compute pairwise RMSD Matrix for the structures in the specified phase
        
        Returns:
            RMSD Matrix (n_structures x n_structures)
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        n = len(structures)
        self._log(f"Computing {n}x{n} RMSD matrix for phase: {phase if phase else 'all'}")
        coords_list = [s.molecule.coordinates for s in structures]
        rmsd_matrix = np.zeros((n, n))
        for i in range(n):
            ci = coords_list[i]
            for j in range(i+1, n):
                rmsd = self._calculate_rmsd(ci, coords_list[j])
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd
        self._log(f"RMSD matrix computed: Range: [{rmsd_matrix[rmsd_matrix > 0].min():.4f}, {rmsd_matrix.max():.4f} ]")
        return rmsd_matrix


    def plot_rmsd_heatmap(self, phase: Optional[str] = None):
        """  
        Helper function to plot the pairwise rmsd matrix as a heatmap
        """
        rmsd_matrix = self.compute_rmsd_matrix(phase=phase)
        plt.figure(figsize=(10, 8))
        sns.heatmap(rmsd_matrix, cmap='viridis')
        plt.title(f'RMSD Matrix Heatmap - Phase: {phase if phase else "All"}')
        plt.xlabel('Structure Index')
        plt.ylabel('Structure Index')
        plt.savefig(f"figures/rmsd_heatmap_{phase if phase else 'all'}.png")
        plt.close()


        
    def plot_int_d_vs_e(self, phase: Optional[str] = None):
        """   
        Plots the intermolecular distance vs energy for the specified phase
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures

        energies = [s.energy for s in structures]
        distances = self.intermolecular_distance(phase=phase)

        plt.figure(figsize=(10, 6))
        plt.scatter(distances, energies, color='purple', alpha=0.7)
        plt.title(f'Intermolecular Distance vs Energy - Phase: {phase if phase else "All"}')
        plt.xlabel('Intermolecular Distance')
        plt.ylabel('Energy')
        plt.grid(True)
        plt.savefig(f"figures/int_d_vs_energy_{phase if phase else 'all'}.png")
        plt.close()

    @staticmethod 
    def _calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
        """   
        Calculate optimal RMSD between two sets of coordinate points.
        For this we use the respective function of Geometry Ops
        No centering is required because we use the optimal correspondence function
        """
        coords_2_opt = GeometryOps.find_optimal_correspondence(coords1, coords2)
        diff = coords1 - coords_2_opt
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        return rmsd

    def plot_rmsd_matrix(self, phase: Optional[str] = None):
        """   
        Plot the RMSD matrix as a heatmap
        """
        rmsd_matrix = self.compute_rmsd_matrix(phase=phase)
        plt.figure(figsize=(10, 8))
        sns.heatmap(rmsd_matrix, cmap='viridis')
        plt.title(f'RMSD Matrix Heatmap - Phase: {phase if phase else "All"}')
        plt.xlabel('Structure Index')
        plt.ylabel('Structure Index')
        plt.savefig(f"figures/rmsd_matrix_{phase if phase else 'all'}.png")
        plt.close()

    def plot_energy_vs_rmsd(self, reference_structure: Optional[StructureData] = None, phase: Optional[str] = None):
        """   
        Plots Energy vs RMSD to a reference structure for the specified phase

        Args:
            reference_structure: Optional reference structure to compute RMSD against. If None, it uses to lowest energy structure as reference.
            phase: Specified phase to plot for. If None, plots for all structures.
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        
        if not structures:
            print(f"No structures found for phase: {phase}")
            return
        
        if reference_structure is None:
            reference_structure = self.get_lowest_energy_structure(phase=phase)

        energies = [s.energy for s in structures]
        rmsd_values = [self._calculate_rmsd(s.molecule.coordinates, reference_structure.molecule.coordinates) for s in structures]

        plt.figure(figsize=(10, 6))
        plt.scatter(rmsd_values, energies, color='blue', alpha=0.7)
        plt.title(f'Energy vs RMSD to Reference - Phase: {phase if phase else "All"}')
        plt.xlabel('RMSD to Reference Structure')
        plt.ylabel('Energy')
        plt.grid(True)
        plt.savefig(f"figures/energy_vs_rmsd_{phase if phase else 'all'}.png")
        plt.close()

    def plot_energy_distribution(self, phase: Optional[str] = None, bins: int = 50, save_path: Optional[str] = None):
        """
        Plots the energy distribution for the specified phase
        Args:
            phase: Specified phase to plot for. If None, plots for all structures.
            bins: Number of bins for the histogram
        """
        if phase:
            energies = [s.energy for s in self.phases[phase]]
        else:
            energies = [s.energy for s in self.structures]

        plt.figure(figsize=(10, 6))
        plt.hist(energies, bins=bins, color='green', alpha=0.7)
        plt.title(f'Energy Distribution - Phase: {phase if phase else "All"}')
        plt.xlabel('Energy')
        plt.ylabel('Frequency')
        plt.grid(True)
        if save_path:
            plt.savefig(save_path)
        else:
            plt.savefig(f"figures/energy_distribution_{phase if phase else 'all'}.png")
        plt.close()

    def plot_rg_vs_energy(self, phase: Optional[str] = None):
        """   
        Plots the Radius of Gyration vs Energy for the specified phase 

        Args:
            phase: Specified phase to plot for. If None, plots for all structures.
        """ 
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            print(f"No structures found for phase: {phase}")
            return
        energies = [s.energy for s in structures]
        rg_values = [self.radius_of_gyration(s.molecule.coordinates, s.molecule.masses) for s in structures]
        plt.figure(figsize=(10, 6))
        plt.scatter(rg_values, energies, color='red', alpha=0.7)
        plt.title(f'Radius of Gyration vs Energy - Phase: {phase if phase else "All"}')
        plt.xlabel('Radius of Gyration')
        plt.ylabel('Energy')
        plt.grid(True)
        plt.savefig(f"figures/rg_vs_energy_{phase if phase else 'all'}.png")
        plt.close()

    

    def get_lowest_energy_structure(self,phase: Optional[str] = None) -> Optional[StructureData]:
        """   
        Returns the structure with the lowest energy in the specified phase
        
        Args:
            phase: Specified phase to search in. If None, searches in all structures.
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        
        if not structures:
            return None
        
        lowest_energy_structure = min(structures, key=lambda s: s.energy)
        return lowest_energy_structure

    
    
    

    @staticmethod
    def compute_interatomic_distance(coords: np.ndarray) -> np.ndarray:
        """   
        Computes the interatomic distance matrix for a given set of coordinates

        Returns:
            Interatomic distance matrix (n_atoms x n_atoms)
        """
        distance_matrix = cdist(coords, coords)
        # Flatten the upper triangle of the distance matrix to get a descriptor vector
        return distance_matrix[np.triu_indices_from(distance_matrix, k=1)]
    
    def compute_interatomic_distance_matrix(self, phase: Optional[str] = None) -> List[np.ndarray]:
        """    
        Compute the interatomic distance matrix for all the structures
        
        Each row is one structure's descriptor vecotr
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures

        descriptors = []
        for s in structures:
            distance_matrix = self.compute_interatomic_distance(s.molecule.coordinates)
            descriptors.append(distance_matrix)

        return descriptors



    @staticmethod 
    def determine_rotational_constants(coords: np.ndarray, masses: np.ndarray) -> Tuple[float, float, float]:
        """   
        Determines the rotational constants (A,B,C) for a given set of coordinates and masses

        Algorithm:
        1. Center coordinates at the center of mass
        2. Compute the inertia tensor
        3. Diagonalize the inertia tensor to get the principal moments of inertia (I_A, I_B, I_C)
        4. Convert to Rotational Constants using the formula:
           A = h / (8 * π^2 * I_A)
           B = h / (8 * π^2 * I_B)
           C = h / (8 * π^2 * I_C)
        """
        h = 6.62607015e-34  # Planck's constant in J*s
        com = GeometryOps.center_of_mass(coords, masses)
        coords_centered = coords - com
        I = GeometryOps.inertia_tensor(coords_centered, masses)
        I_A, I_B, I_C = np.linalg.eigvalsh(I)  # Get principal moments of inertia
        A = h / (8 * np.pi**2 * I_A) if I_A > 0 else 0.0
        B = h / (8 * np.pi**2 * I_B) if I_B > 0 else 0.0
        C = h / (8 * np.pi**2 * I_C) if I_C > 0 else 0.0
        return A, B, C
    

    # ====================== Data Cleanup ========================

    """ 
    We should first of perform a rmsd filtering to remove duplicate structures and then perform a energy-based filtering to remove
    both high energy outliers and low energy duplicates that are very similar to the lowest energy structure.
    This will ensure we have a diverse set of representative structures
    """
    def rmsd_filtering(self, threshold: float = 0.2, phase: Optional[str] = None):
        """   
        Filters out structures that are within a certain RMSD threshold of each other

        Args:
            threshold: RMSD threshold for filtering
            phase: Specified phase to filter. If None, filters all structures.
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return
        n_before = len(structures)
        self._log(f"RMSD filtering (threshold={threshold:.3f} Å) - Phase: {phase if phase else 'all'} - Starting with {n_before} structures")

        # Extract all coordinates
        all_coords = np.array([s.molecule.coordinates for s in structures]) 
        n_atoms = all_coords[0].shape[0]
        selected_indices = []
        selected_flat = np.empty((len(structures), n_atoms * 3), dtype=all_coords[0].dtype)
        n_selected = 0
        for i, coords in enumerate(all_coords):
            candidate_flat = coords.ravel()
            if n_selected == 0:
                selected_indices.append(i)
                selected_flat[n_selected] = candidate_flat
                n_selected += 1
                continue 
            diff = selected_flat[:n_selected] - candidate_flat
            rmsd_values = np.sqrt(np.mean(diff.reshape(n_selected, n_atoms, 3)**2, axis=(1,2)))

            if np.all(rmsd_values >= threshold):
                selected_indices.append(i)
                selected_flat[n_selected] = candidate_flat
                n_selected += 1

        filtered_structures = [structures[i] for i in selected_indices]
        if phase:
            self.phases[phase] = filtered_structures
        else:
            self.structures = filtered_structures
        n_after = len(filtered_structures)
        self._log(f"RMSD filtering completed - Phase: {phase if phase else 'all'} - {n_after} structures remaining (removed {n_before - n_after})")  

    def filter_high_energy_outliers(self, method: str = "iqr", factor: float = 1.5, phase: Optional[str] = None) -> None:
        """
        Remove high energy outliers based on the chosen method:

        IQR: Removes strctures with energy > Q3 + factor * IQR
        STD: Removes structures with energy > mean + factor * std
        ABSOLUTE: Removes structures with energy > absolute threshold (factor is the threshold value in this case)
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures 
        energies = np.array([s.energy for s in structures])
        if method == "iqr":
            q75, q25 = np.percentile(energies, [75, 25])
            iqr = q75 - q25
            threshold = q75 + factor * iqr
        elif method == "std":
            mean = np.mean(energies)
            std = np.std(energies)
            threshold = mean + factor * std
        elif method == "absolute":
            threshold = factor
        else:
            raise ValueError(f"Invalid method: {method}. Choose from 'iqr', 'std', 'absolute'.")
    
        filtered_structures = [s for s in structures if s.energy <= threshold]
        if phase:
            self.phases[phase] = filtered_structures
        else:
            self.structures = filtered_structures
        self._log(f"High energy outlier filtering (method={method}, threshold={threshold:.6f}) - Phase: {phase if phase else 'all'} - {len(filtered_structures)} structures remaining (removed {len(structures) - len(filtered_structures)})")


                   

    def _compute_com_trajectories(self, structures: List[StructureData]) -> List[np.ndarray]:
        """
        Compute geometric center (centroid) trajectories for each submolecule.
        
        Note: Using geometric mean (centroid) instead of mass-weighted COM.
        
        Args:
            structures: List of StructureData objects
            
        Returns:
            List of numpy arrays, each of shape (n_frames, 3) containing centroid coordinates
        """
        n_submols = len(self.submolecule_indices)
        n_frames = len(structures)
        
        centroid_trajectories = []
        
        for submol_idx in self.submolecule_indices:
            centroid_traj = np.zeros((n_frames, 3))
            
            for frame_idx, structure_data in enumerate(structures):
                mol = structure_data.molecule
                
                # Ensure coordinates are numpy array
                if isinstance(mol.coordinates, list):
                    mol_coords = np.array(mol.coordinates)
                else:
                    mol_coords = mol.coordinates
                
                # Extract submolecule coordinates
                submol_coords = mol_coords[submol_idx]
                
                # Calculate geometric center (simple mean of coordinates)
                centroid = np.mean(submol_coords, axis=0)
                centroid_traj[frame_idx] = centroid
            
            centroid_trajectories.append(centroid_traj)
        
        return centroid_trajectories


    @staticmethod
    def remove_string_digits(s: str) -> str:
        """   
        Utility function to remove digits from a string (e.g., "H2O" -> "HO")
        """
        return ''.join(filter(lambda x: not x.isdigit(), s))


    # =========================== Hydrogen Bond Analysis ===========================


    def compute_hbond_features(self,
                               angle_threshold: float = 150.0,
                               max_distance: float = 3.5,
                               phase: Optional[str] = None,
                               n_processes: Optional[int] = None) -> np.ndarray:
        """ 
        Computes hydrogen bond features for all structures

        Computes:
            + Number of valid Hydrogen bond configuration 
            + Average Hydrogen bond angle (0.0 if no H-Bonds found)
            + Average Hydrogen bond distance (0.0 if no H-Bonds found)

        Uses multiprocessing to parallelize the computation across structures
            
        Args:
            angle_treshold: Minimum angle (in degrees) for a valid hydrogen bond
            max_distance: Maximum distance (in Å) for a valid hydrogen bond
            phase: Specified phase to compute for. If None, computes for all structures.
            n_processes: Number of parallel processes to use for computation. If None, uses all available cores.

        Returns:
            Numpy array of shape (n_structures, 3) with columns [n_hbonds, avg_angle, avg_da_distance]
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        n_structures = len(structures)
        self._log(f"Computing H-bond features for {n_structures} (angle_threshold={angle_threshold}°) - Phase: {phase if phase else 'all'}")

        if n_structures == 0:
            return np.zeros((0, 3))

        # Number of processes
        if n_processes is None:
            n_processes = min(cpu_count() - 2, n_structures)
        n_processes = max(1, min(n_processes, n_structures))
        self._log(f"Using {n_processes} parallel processes for H-bond computation")

        worker_args = []
        for structure_data in structures:
            mol = structure_data.molecule
            worker_args.append((
                list(mol.atom_labels),
                np.array(mol.coordinates),
                angle_threshold,
                max_distance
            ))
        
        # Determine the chunksie
        chunksize = max(1, min(20, math.ceil(n_structures / (n_processes * 4))))

        if n_processes > 1:
            self._log("Starting multiprocessing pool for H-bond feature computation...")
            with Pool(processes=n_processes) as pool:
                results = pool.map(_compute_hbond_single, worker_args, chunksize=chunksize)
        else:
            self._log("Computing H-bond features sequentially (n_processes=1)...")
            results = list(map(_compute_hbond_single, worker_args))
        # Collect results
        features = np.zeros((n_structures, 3))
        total_hbonds = 0
        structures_with_hbonds = 0
        for i ,(n_hbonds, avg_angle, avg_da_distance) in enumerate(results):
            features[i, 0] = n_hbonds
            features[i, 1] = avg_angle
            features[i, 2] = avg_da_distance
            total_hbonds += n_hbonds
            if n_hbonds > 0:
                structures_with_hbonds += 1
        self._log(f"Structures with H-bonds: {structures_with_hbonds}/{n_structures}")
        self._log(f"Total H-bond configs: {total_hbonds}")
        if structures_with_hbonds > 0:
            self._log(f"Avg H-bonds per struct (with H-bonds): {total_hbonds / structures_with_hbonds:.2f}")
            hbond_angles = features[features[:, 0] > 0, 1]
            hbond_distances = features[features[:, 0] > 0, 2]
            self._log(f"Avg H-bond angle: {np.mean(hbond_angles):.2f}° (std: {np.std(hbond_angles):.2f}°)")
            self._log(f"Min/Max H-bond angle: {np.min(hbond_angles):.2f}° / {np.max(hbond_angles):.2f}°")
            self._log(f"Avg D-A distance: {np.mean(hbond_distances):.2f} Å (std: {np.std(hbond_distances):.2f} Å)")

        return features


    def compute_mulliken_charge_features(self, phase: Optional[str] = None) -> np.ndarray:
        """  
        Small helper function to compute Mulliken charge features for all structures.
        The Mulliken charges are defined for each atom so we flatten it into a single vector for each structure.
        Also compute important statistics
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        n_structures = len(structures)
        self._log(f"Computing Mulliken charge features for {n_structures} - Phase: {phase if phase else 'all'}")

        if n_structures == 0:
            return np.zeros((0, 1))

        available_charge_lengths = [len(s.mulliken_charges) for s in structures if s.mulliken_charges is not None]
        if not available_charge_lengths:
            self._log("No Mulliken charges found. Returning zero placeholder feature.", level="warning")
            return np.zeros((n_structures, 1))

        max_atoms = max(available_charge_lengths)
        features = np.zeros((n_structures, max_atoms))
        for i, s in enumerate(structures):
            if s.mulliken_charges is not None:
                charges = np.array(s.mulliken_charges)
                features[i, :len(charges)] = charges
        # Log statistics
        all_charges = features[features != 0]
        self._log(f"Mulliken Charge Statistics (non-zero values):")
        self._log(f"  Count: {len(all_charges)}")
        self._log(f"  Mean: {np.mean(all_charges):.4f}")
        self._log(f"  Std: {np.std(all_charges):.4f}")
        self._log(f"  Min: {np.min(all_charges):.4f}")
        self._log(f"  Max: {np.max(all_charges):.4f}")
        return features


    def compute_coloumb_matrix_features(self, phase: Optional[str] = None) -> np.ndarray:
        """ 
        Computes the Coloumb-Matrix and and sorts the upper triangle matrix to get a descriptor.
        This descriptor is of length n_atoms * (n_atoms - 1) / 2 for each structure.

        The coloumb matrix is symmetric and the elements C_ij are defined by

        C_ij = 0.5 * Z_i^2.4 if i == j
        C_ij = Z_i * Z_j / |R_i - R_j| if i != j

        where Z_i is the atomic number of atom i and R_i is the position vector of atom i.
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures

        if not structures:
            return np.zeros((0, 0))
        
        # Compute the feature length
        max_atoms = max(len(s.molecule.atom_labels) for s in structures)
        feature_length = max_atoms * (max_atoms - 1) // 2

        features = np.zeros((len(structures), feature_length))
        for i, s in enumerate(structures):
            Z = np.asarray(s.molecule.atomic_numbers, dtype=np.float64)
            coords = s.molecule.coordinates
            N = len(Z)
            if len(Z) == 0:
                raise ValueError(f"Structure {i} has no atomic numbers defined.")
            if len(coords) == 0:
                raise ValueError(f"Structure {i} has no coordinates defined.")

            # Build Coulomb matrix: C_ii = 0.5*Z_i^2.4, C_ij = Z_i*Z_j/|R_i-R_j|
            dist = cdist(coords, coords)
            zz = np.outer(Z, Z)
            with np.errstate(divide='ignore', invalid='ignore'):
                mat = np.divide(zz, dist, out=np.zeros_like(zz), where=dist > 1e-12)
            np.fill_diagonal(mat, 0.5 * (Z ** 2.4))

            # Keep permutation-robust descriptor style analogous to interatomic distances.
            upper_tri = mat[np.triu_indices(N, k=1)]
            features[i, :len(upper_tri)] = np.sort(upper_tri)
        return features

    def eigenspectrum_coloumb_matrix(self, phase: Optional[str] = None) -> np.ndarray:
        """ 
        Computes the eigenspectrum of the coloumb matrix for each structure and uses the sorted eigenvalues as features.
        """ 
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures

        if not structures:
            return np.zeros((0, 0))
        
        max_atoms = max(len(s.molecule.atom_labels) for s in structures)
        features = np.zeros((len(structures), max_atoms))
        for i, s in enumerate(structures):
            Z = np.asarray(s.molecule.atomic_numbers, dtype=np.float64)
            coords = s.molecule.coordinates
            N = len(Z)
            if len(Z) == 0:
                raise ValueError(f"Structure {i} has no atomic numbers defined.")
            if len(coords) == 0:
                raise ValueError(f"Structure {i} has no coordinates defined.")

            dist = cdist(coords, coords)
            zz = np.outer(Z, Z)
            with np.errstate(divide='ignore', invalid='ignore'):
                mat = np.divide(zz, dist, out=np.zeros_like(zz), where=dist > 1e-12)
            np.fill_diagonal(mat, 0.5 * (Z ** 2.4))

            eigenvalues = np.linalg.eigvalsh(mat)
            features[i, :len(eigenvalues)] = np.sort(eigenvalues)[::-1]  # Sort in descending order
        return features

    

    # ====================== Clustering Methods =======================

    def _ensure_feature_matrix(self, remove_zero_variance: bool = True):
        """   
        Compute and cache the feature matrix if not already done.

        Further removes zero variance features to avoid issues in clustering algorithms and to reduce dimensionality.
        """
        # Start Profiling and Memory Tracking
        start_time = time.time()
        pr = cProfile.Profile()
        pr.enable()
        tracemalloc.start()
        try:
            if self._feature_matrix_raw is not None:
                return 
            self._log("Computing and caching feature matrix...")

            # Base Energy and Dipoles
            lowest_energy = self.get_lowest_energy_structure()
            delta_e = np.array([s.energy - lowest_energy.energy for s in self.structures]).reshape(-1, 1)
            dipole_mags = np.array([s.dipole_magnitude for s in self.structures]).reshape(-1, 1)

            # Get Mulliken charge features
            mulliken_features = self.compute_mulliken_charge_features()
            


            # H-Bond Features
            hbond_features = self.compute_hbond_features()

            # Extractor for Geometric Features
            extractor = FeatureExtractor()
            geom_features = extractor.extract_fast_features(self.structures, self.submolecule_indices)

            rmsd_values = np.array([self._calculate_rmsd(s.molecule.coordinates, self.get_lowest_energy_structure().molecule.coordinates) for s in self.structures])
            rotational_constants = np.array([self.determine_rotational_constants(s.molecule.coordinates, s.molecule.masses) for s in self.structures])
            intermolecular_distances = np.array(self.compute_interatomic_distance_matrix())
            coloumb_matrix_features = self.compute_coloumb_matrix_features()

            

            self._feature_matrix_raw = np.hstack((
                delta_e, 
                rmsd_values.reshape(-1, 1), 
                geom_features[:, :1],  # Rg
                rotational_constants, 
                hbond_features, 
                intermolecular_distances,
                coloumb_matrix_features,
                geom_features[:, 1:4],  # Asphericity, Acylindricity, Kappa^2
                dipole_mags, 
                mulliken_features
            )


            )

            n_iad = intermolecular_distances.shape[1] if intermolecular_distances.ndim == 2 else 0
            n_coloumb = coloumb_matrix_features.shape[1] if coloumb_matrix_features.ndim == 2 else 0
            n_mulliken = mulliken_features.shape[1] if mulliken_features.ndim == 2 else 0
            all_feature_names = [
                'Delta E', 'RMSD', 'Rg', 'Rot A', 'Rot B', 'Rot C',
                'Num H-bonds', 'Avg H-bond Angle', 'Avg D-A Distance'
            ]
            all_feature_names.extend([f'Dist_{i}' for i in range(n_iad)])
            all_feature_names.extend([f'Coulomb_{i}' for i in range(n_coloumb)])
            all_feature_names.extend(['Asphericity', 'Acylindricity', 'Kappa^2', 'Dipole Magnitude'])
            all_feature_names.extend([f'Mulliken_{i}' for i in range(n_mulliken)])
            self.used_features = all_feature_names

            self._log(f"Feature matrix shape: {self._feature_matrix_raw.shape}")
            self._log(f"Features: ['Delta E', 'RMSD', 'Rg', 'Rot A', 'Rot B', 'Rot C', 'Num H-bonds', 'Avg H-bond Angle', 'Avg D-A Distance', 'Intermolecular Distances...', 'Coulomb Matrix...', 'Asphericity', 'Acylindricity', 'Kappa^2', 'Dipole Magnitude', 'Mulliken Charges...']")
            

            feature_variances = []
            for i in range(self._feature_matrix_raw.shape[1]):
                var = np.var(self._feature_matrix_raw[:, i])
                feature_variances.append(var)
                self._log(f"Feature {i} variance: {var:.6f}")
            if remove_zero_variance:
                non_zero_var_indices = [i for i, var in enumerate(feature_variances) if var > 1e-6]
                self._non_zero_var_indices = non_zero_var_indices  # Store which features are kept
                self._feature_matrix_raw = self._feature_matrix_raw[:, non_zero_var_indices]
                if self.used_features is not None:
                    self.used_features = [self.used_features[i] for i in non_zero_var_indices]
                self._log(f"Removed {len(feature_variances) - len(non_zero_var_indices)} zero-variance features. New shape: {self._feature_matrix_raw.shape}")
             
            # Log Statistics of remaining features
            for i in range(self._feature_matrix_raw.shape[1]):
                col = self._feature_matrix_raw[:, i]
                self._log(f"Feature {i} stats: mean={np.mean(col):.4f}, std={np.std(col):.4f}, min={np.min(col):.4f}, max={np.max(col):.4f}")
            
            
            scaler = StandardScaler()
            
            self._feature_matrix_normalized = scaler.fit_transform(self._feature_matrix_raw)
            self._log("Feature matrix computed and cached.")

            # Log memory usage of the feature matrix
            self._log(f"Size of raw feature matrix: {self._feature_matrix_raw.nbytes / 1e6:.2f} MB")
            self._log(f"Size of normalized feature matrix: {self._feature_matrix_normalized.nbytes / 1e6:.2f} MB")
        finally:
            pr.disable()
            s = io.StringIO()
            ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
            ps.print_stats(10)  # Print top 10 functions
            current, peak = tracemalloc.get_traced_memory()
            # Log Memory Usage
            self._log(f"Current memory usage: {current / 1e6:.2f} MB; Peak memory usage: {peak / 1e6:.2f} MB")
            tracemalloc.stop()
            total_time = time.time() - start_time
            self._log(f"Time taken to compute feature matrix: {total_time:.2f} seconds")
            # Write CPU usage to profiling folder
            with open("profiling/feature_matrix_cpu.txt", "w") as f:
                f.write(f"--- Profiling Feature Matrix Computation ---\n")
                f.write(s.getvalue())



    def feature_matrix(self, normalize = True) -> np.ndarray:
        """   
        Current Feature Matrix wifth Columns

        1. Delta E to lowest energy structure
        2. RMSD to lowest energy structure
        3. Radius of Gyration
        4. Rotational Constants (A,B,C)
        5. Intermolecular Distance Matrix (flattened upper triangle)
        6. Number of valid H-bond configurations
        7. Average H-bond angle
        8. Average D-A distance
        9. Gyration Tensor Shape Descriptors (Asphericity, Acylindricity, Kappa^2)
        10. Coulomb Matrix Descriptor (sorted upper triangle)
        11. Dipole Magnitude (Debye)
        12. Mulliken Charges (flattened)
        """
        self._ensure_feature_matrix()
        if normalize:
            return self._feature_matrix_normalized
        else:
            return self._feature_matrix_raw
        
    
    # ================================= Clustering Algorithms =================================

    def cluster(self,
                method: str = "agglomerative",
                n_clusters: Optional[int] = None,
                phase: Optional[str] = None,
                **kwargs) -> np.ndarray:
        """  
        Clustering Interface 

        Args:
            method: Clustering method to use ("agglomerative", "kmeans")
            n_clusters: Number of clusters to form
            phase: Phase to cluster. If None, clusters all structures.
            **kwargs: Additional keyword arguments:
                 - For Agglomerative Clustering: linkage (default "ward")
                 - For kmeans: n_init (int), max_iter (int), random_state (int)
                 - For DBSCAN: eps (float), min_samples (int)
                 - For HDBSCAN: min_cluster_size (int), min_samples (int)
        Returns:
            Cluster labels array
        """
        if method == "agglomerative":
            if n_clusters is None:
                raise ValueError("n_clusters must be specified for agglomerative clustering")
            return self.agglomerative_clustering(n_clusters=n_clusters, phase=phase, **kwargs)
        elif method == "kmeans":
            if n_clusters is None:
                raise ValueError("n_clusters must be specified for kmeans clustering")
            return self.kmeans_clustering(n_clusters=n_clusters, phase=phase, **kwargs)
        elif method == "dbscan":
            return self.dbscan_clustering(**kwargs, phase=phase)
        elif method == "hdbscan":
            return self.hdbscan_clustering(**kwargs, phase=phase)
        
        else:
            raise ValueError(f"Unsupported clustering method: {method}. Supported methods: 'agglomerative', 'kmeans'")



    def dbscan_clustering(self,
                          eps: float = 0.5,
                          min_samples: int = 5,
                          phase: Optional[str] = None) -> np.ndarray:
        """  
        Perform DBSCAN clustering on the structures

        DBSCAN is density-based and does not require n_clusters. It discovers clusters from the
        data based on eps (neighborhood radius) and min_samples
        Noise Points are labeled as -1

        Args:
            eps: The maximum distance between two samples for them to be considered as in the same neighborhood.
            min_samples: The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. This includes the point itself.
            phase: Specified phase to cluster. If None, clusters all structures.
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return np.array([])
        self._log(f"Performing DBSCAN Clustering (eps={eps}, min_samples={min_samples}) for phase: {phase if phase else 'all'} - n_structures: {len(structures)}")
        feature_matrix = self.feature_matrix()
        model = DBSCAN(eps=eps, min_samples=min_samples)
        labels = model.fit_predict(feature_matrix)
        self.labels = labels
        n_clusters = len(set(labels) - {-1}) # Exclude noise label (-1)
        n_noise = np.sum(labels == -1)
        self._log(f"DBSCAN found {n_clusters} clusters and {n_noise} noise points")
        return labels

    def hdbscan_clustering(self,
                           min_cluster_size: int = 5,
                           min_samples: Optional[int] = None,
                           phase: Optional[str] = None) -> np.ndarray:
        """ 
        Perform HDBSCAN clustering on the structures

        HDBSCAN is hierarchical density-based clustering that automatically determines the number of clusters.
        It handles varying density better than DBSCAN
        Noise points are labeled as -1
        """
        if HDBSCAN is None:
            raise ImportError("HBDSCAN is not installed")

        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return np.array([])
        if min_samples is None:
            min_samples = min_cluster_size

        self._log(f"Performing HDBSCAN Clustering (min_cluster_size={min_cluster_size}, min_samples={min_samples}) for phase: {phase if phase else 'all'} - n_structures: {len(structures)}")
        feature_matrix = self.feature_matrix()
        model = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples)
        labels = model.fit_predict(feature_matrix)
        self.labels = labels
        n_clusters = len(set(labels) - {-1}) # Exclude noise label (-1)
        n_noise = np.sum(labels == -1)
        self._log(f"HDBSCAN found {n_clusters} clusters and {n_noise} noise points")
        return labels
    

    def agglomerative_clustering(self, n_clusters: int = 5, phase: Optional[str] = None, linkage: str = "ward") -> np.ndarray:
        """   
        Perform Agglomerative Clustering on the structures in the specified phase

        Args:
            n_clusters: Number of clusters to form
            phase: Specified phase to cluster. If None, clusters all structures.
            linkage: Linkage criterion for agglomerative clustering ("ward", "complete", "average", "single")
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            print(f"No structures found for phase: {phase}")
            return np.array([])

        self._log(f"Performing Agglomerative Clustering (n_clusters={n_clusters}, linkage='{linkage}') for phase: {phase if phase else 'all'} - n_structures: {len(structures)}")
        feature_matrix = self.feature_matrix()
        clustering_model = SklearnAgglomerativeClustering(n_clusters=n_clusters, linkage=linkage)
        labels = clustering_model.fit_predict(feature_matrix)
        self.labels = labels  # Cache labels for later use
        return labels
    
    def kmeans_clustering(self,
                          n_clusters: int = 5,
                          phase: Optional[str] = None,
                          n_init: int = 10,
                          max_iter: int = 300,
                          random_state: Optional[int] = 42) -> np.ndarray:
        """ 
        Perform KMeans clustering on the structures

        Args:
            n_clusters: Number of clusters to form
            phase: Specified phase to cluster. If None, clusters all structures.
            n_init: Number of time the k-means algorithm will be run with different centroid seeds. The final results will be the best output of n_init consecutive runs in terms of inertia.
            max_iter: Maximum number of iterations of the k-means algorithm for a single run.
            random_state: Determines random number generation for centroid initialization. Use an int to make the randomness deterministic.
        
        Returns:
            Cluster labels array
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return np.array([])
        
        self._log(f"Performing KMeans Clustering (n_clusters={n_clusters}, n_init={n_init}, max_iter={max_iter}) for phase: {phase if phase else 'all'} - n_structures: {len(structures)}")
        feature_matrix = self.feature_matrix()
        model = KMeans(n_clusters=n_clusters,
                       n_init=n_init,
                       max_iter=max_iter,
                       random_state=random_state)
        labels = model.fit_predict(feature_matrix)
        self.labels = labels  # Cache labels for later use

        return labels


    # ================================ Obtain the Cluster Representatives ================================
    def get_cluster_representatives(self, method: str = "centroid", labels: Optional[np.ndarray] = None) -> List[StructureData]:
        """  
        Get the representative structures for each cluster

        
            +   Method == "centroid": For each cluster computes the centroid in the normalized feature space and selects the structure closest to the centroid (Euclidean distance)
            +   Method == "medoid": Selects the point with the smallest average distance to all points in the cluster
            +   Method == "density": Alias for "medoid" (robust for DBSCAN/HDBSCAN clusters)
            +   Method == "lowest_energy": Selects the lowest energy structure in each cluster as the representative
        Args:
            labels: Clusterer labels array. If None, uses cached labels from the last clustering run. If no cached labels, raises  ValueError
        
        Returns:
            List of StructureData objects representing the cluster centers
        """
        if labels is None:
            labels = self.labels
        if labels is None:
            raise ValueError("No cluster labels available. Run a clustering method first")

        feature_mat = self.feature_matrix(normalize=True)
        unique_labels = np.unique(labels)
        representatives = []

        def _medoid_index(cluster_features: np.ndarray) -> int:
            # Compute pairwise distances and choose the point with minimal mean distance.
            dists = cdist(cluster_features, cluster_features)
            return int(np.argmin(np.mean(dists, axis=1)))

        if method == 'centroid':
            for label in unique_labels:
                if label == -1:
                    # Skip noise points form DBSCAN/HDBSCAN
                    continue
                mask = labels == label
                cluster_features  = feature_mat[mask]
                cluster_indices = np.where(mask)[0]

                centroid = cluster_features.mean(axis=0)
                distances = np.linalg.norm(cluster_features - centroid, axis=1)
                closest_idx = cluster_indices[np.argmin(distances)]
                rep = self.structures[closest_idx]
                self._log(f"Cluster {label}: representative index {closest_idx}, distance to centroid: {distances.min():.3f}, E={rep.energy:.6f} Ha")
                representatives.append(rep)
        
            self._log(f"Total clusters (excluding noise): {len(unique_labels) - (1 if -1 in unique_labels else 0)}. Representatives obtained: {len(representatives)}")
        
        elif method == 'medoid' or method == 'density':
            for label in unique_labels:
                if label == -1:
                    # Skip noise points form DBSCAN/HDBSCAN
                    continue
                mask = labels == label
                cluster_features = feature_mat[mask]
                cluster_indices = np.where(mask)[0]
                medoid_local = _medoid_index(cluster_features)
                closest_idx = cluster_indices[medoid_local]
                rep = self.structures[closest_idx]
                self._log(f"Cluster {label}: representative index {closest_idx}, medoid avg dist: {np.mean(cdist(cluster_features, cluster_features)[medoid_local]):.3f}, E={rep.energy:.6f} Ha")
                representatives.append(rep)
            self._log(f"Total clusters (excluding noise): {len(unique_labels) - (1 if -1 in unique_labels else 0)}. Representatives obtained: {len(representatives)}")

        elif method == 'lowest_energy':
            for label in unique_labels:
                if label == -1:
                    # Skip noise points form DBSCAN/HDBSCAN
                    continue
                mask = labels == label
                cluster_structures = [self.structures[i] for i in np.where(mask)[0]]
                lowest_energy_struct = min(cluster_structures, key=lambda s: s.energy)
                representatives.append(lowest_energy_struct)
                self._log(f"Cluster {label}: representative index {self.structures.index(lowest_energy_struct)}, E={lowest_energy_struct.energy:.6f} Ha")
            self._log(f"Total clusters (excluding noise): {len(unique_labels) - (1 if -1 in unique_labels else 0)}. Representatives obtained: {len(representatives)}")
        else:
            raise ValueError(
                "Unsupported method for selecting cluster representatives: "
                f"{method}. Supported methods: 'centroid', 'medoid', 'density', 'lowest_energy'"
            )
        
        
        return representatives

    

    def pca(self, n_components: int = 2, phase: Optional[str] = None) -> np.ndarray:
        """   
        Performs PCA Analysis on the structures in the specified phase

        Args:
            n_components: Number of principal components to compute
            phase: Specified phase to analyze. If None, analyzes all structures.
        """ 
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            print(f"No structures found for phase: {phase}")
            return np.array([])
        feature_matrix = self.feature_matrix()
        pca_model = PCA(n_components=n_components)
        pca_result = pca_model.fit_transform(feature_matrix)
        explained_variance = pca_model.explained_variance_ratio_
        return pca_result, explained_variance
        


    def tsne(self, n_components: int = 2, perplexity: float = 30.0, n_iter: int = 1000, random_state: int = 42, phase: Optional[str] = None) -> np.ndarray:
        """   
        Performs t-SNE dimensionality reduction on the feature matrix.

        Args:
            n_components: Number of t-SNE dimensions (typically 2)
            perplexity: Perplexity parameter for t-SNE
            n_iter: Number of iterations for optimization
            random_state: Random seed for reproducibility
            phase: Specified phase to analyze. If None, analyzes all structures.

        Returns:
            t-SNE embedding array of shape (n_structures, n_components)
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return np.array([])
        
        feature_mat = self.feature_matrix()
        # Adjust perplexity if needed: must be less than n_samples
        effective_perplexity = min(perplexity, len(structures) - 1)
        tsne_model = TSNE(n_components=n_components, perplexity=effective_perplexity,
                          n_iter=n_iter, random_state=random_state)
        tsne_result = tsne_model.fit_transform(feature_mat)
        return tsne_result

    def umap(self, n_components: int = 2, n_neighbors: int = 15, min_dist: float = 0.1, random_state: int = 42, phase: Optional[str] = None) -> np.ndarray:
        """   
        Performs UMAP dimensionality reduction on the feature matrix.

        Args:
            n_components: Number of UMAP dimensions (typically 2)
            n_neighbors: Number of neighbors for UMAP
            min_dist: Minimum distance parameter for UMAP
            random_state: Random seed for reproducibility
            phase: Specified phase to analyze. If None, analyzes all structures.

        Returns:
            UMAP embedding array of shape (n_structures, n_components)
        """
        if UMAP is None:
            raise ImportError("UMAP is not installed. Install with `pip install umap-learn`.")
        
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return np.array([])
        
        feature_mat = self.feature_matrix()
        # Adjust n_neighbors if needed
        effective_neighbors = min(n_neighbors, len(structures) - 1)
        umap_model = UMAP(n_components=n_components, n_neighbors=effective_neighbors,
                          min_dist=min_dist, random_state=random_state)
        umap_result = umap_model.fit_transform(feature_mat)
        return umap_result

    def get_feature_names(self) -> List[str]:
        """ 
        Helper Method that returns a list of feature names corresponding to the columns of the feature matrix.
        Handles both filtered and unfiltered feature matrices.
        """
        self._ensure_feature_matrix()
        if self.used_features is None:
            return []
        return self.used_features
    
    def plot_pca_loadings(self, phase: Optional[str] = None, pc_index: int= 0, top_n_features: int = -1):
        """  
        Plots the PCA loadings to show which original features contribute most to a specific principal component

        Higher absolute loading values indicate that the feature has a stronger influence on that principal component.

        Args:
            phase: Specified phase to analyze. If None, analyzes all structures.
            pc_index: Index of the principal component to analyze (0-based)
            top_n_features: Number of top features to display based on absolute loading values, -1 is all features
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return

        feature_matrix = self.feature_matrix()
        feature_names = self.get_feature_names()

        if top_n_features == -1 or top_n_features > len(feature_names):
            top_n_features = len(feature_names)

        n_comps = min(feature_matrix.shape[1], 20)  # Limit to 20 components for loading analysis
        pca_model = PCA(n_components=n_comps)
        pca_model.fit(feature_matrix)

        if pc_index >= n_comps:
            self._log(f"Requested PC index {pc_index} exceeds number of computed components {n_comps}. Adjusting to {n_comps - 1}", level="warning")
            return 
        
        # Get the loadings for the specified principal component
        loadings = pca_model.components_[pc_index]

        # Sort features by their absolute contribution
        sorted_indices = np.argsort(np.abs(loadings))[::-1][:top_n_features]
        sorted_loadings = loadings[sorted_indices]
        sorted_names = [feature_names[i] for i in sorted_indices]

        plt.figure(figsize=(12, 6))

        # Positive Correlations in blue, negative in red
        colors = ["red" if l < 0 else "blue" for l in sorted_loadings]
        bars = plt.bar(range(top_n_features), sorted_loadings, color=colors, alpha=0.7, edgecolor='black')
        plt.xticks(range(top_n_features), sorted_names, rotation=45, ha='right',fontsize=9)
        explained_var = pca_model.explained_variance_ratio_[pc_index] * 100
        plt.title(f'PCA Loadings for PC{pc_index + 1} ({explained_var:.1f}% variance) - Phase: {phase if phase else "All"}', fontsize=12)
        plt.ylabel('Loading Value')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"figures/pca_loadings_pc{pc_index + 1}_{phase if phase else 'all'}.png", dpi=200)
        plt.close()
        self._log(f"Top {top_n_features} features for PC{pc_index + 1}:")
        for name, loading in zip(sorted_names, sorted_loadings):
            self._log(f" {name:<20} : {loading:.4f}")

    def plot_pca_3D_clustered(self,
                              n_clusters: Optional[int] = None,
                              phase: Optional[str] = None,
                              cluster_method: str = "agglomerative",
                              **cluster_kwargs):
        """ 
        Plots the PCA as a 3D dimensional scatter plot.
        The x-y dimensions are the first two principal components and the z-dimension is the energies
        
        Args:
            n_clusters: Number of clusters to form (only for agglomerative and kmeans)
            phase: Specified phase to analyze. If None, analyzes all structures.
            cluster_method: Clustering method ("agglomerative", "kmeans", "dbscan", "hdbscan")
            **cluster_kwargs: Passed to the clustering method
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            print(f"No structures found for phase: {phase}")
            return
        # Obtain the PCA
        pca_result, explained_variance = self.pca(phase=phase)
        if cluster_method in ["dbscan", "hdbscan"]:
            labels = self.cluster(method=cluster_method, phase=phase, **cluster_kwargs)
            title_suffix = f"clusters={len(set(labels) - {-1})}, noise={np.sum(labels == -1)}"
        else:
            if n_clusters is None:
                raise ValueError("n_clusters must be specified for agglomerative or kmeans clustering")
            labels = self.cluster(method=cluster_method, n_clusters=n_clusters, phase=phase, **cluster_kwargs)
            title_suffix = f"k={n_clusters}"
        energies = np.array([s.energy for s in (self.phases[phase] if phase else self.structures)])
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(pca_result[:, 0], pca_result[:, 1], energies, c=labels, cmap='tab10', alpha=0.7)
        ax.set_xlabel(f'PC1 ({explained_variance[0]*100:.1f}%)')
        ax.set_ylabel(f'PC2 ({explained_variance[1]*100:.1f}%)')
        ax.set_zlabel('Energy (Ha)')
        ax.set_title(f'3D PCA + {cluster_method.capitalize()} ({title_suffix}) - Phase: {phase if phase else "All"}', fontsize=12)
        plt.colorbar(scatter, label='Cluster Label')
        plt.grid(True, alpha=0.3)
        plt.savefig(f"figures/pca_3d_{cluster_method}_{phase if phase else 'all'}.png", dpi=200, bbox_inches='tight')
        plt.close()

    def plot_pca_explained_variance(self, n_components: int = 20, phase: Optional[str] = None):
        """  
        Plots the explained variance ratio for PCA components
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            print(f"No structures found for phase: {phase}")
            return
        pca_result, explained_variance = self.pca(n_components=n_components, phase=phase)
        plt.figure(figsize=(10, 6))
        plt.bar(range(1, n_components + 1), explained_variance * 100, color='skyblue')
        plt.xlabel('Principal Component')
        plt.ylabel('Explained Variance (%)')
        plt.title(f'PCA Explained Variance - Phase: {phase if phase else "All"}')
        plt.xticks(range(1, n_components + 1))
        plt.grid(True, alpha=0.3)
        plt.savefig(f"figures/pca_explained_variance_{phase if phase else 'all'}.png")
        plt.close()


    def plot_pca_clustered(self, 
                           n_clusters: Optional[int] = None,
                           phase: Optional[str] = None, 
                           cluster_method: str = "agglomerative",
                           **cluster_kwargs):
        """   
        Plots the PCA projection colored by cluster labels.

        Args:
            n_clusters: Number of clusters to form
            phase: Specified phase to analyze. If None, analyzes all structures.
            cluster_method: Clustering method ("agglomerative", "kmeans")
            **cluster_kwargs: Passed to the clustering method
        """
        pca_result, explained_variance = self.pca(phase=phase)


        if cluster_method in ["dbscan", "hdbscan"]:
            labels = self.cluster(method=cluster_method, phase=phase, **cluster_kwargs)
            title_suffix = f"clusters={len(set(labels) - {-1})}, noise={np.sum(labels == -1)}"
        else:
            if n_clusters is None:
                raise ValueError("n_clusters must be specified for agglomerative or kmeans clustering")
            labels = self.cluster(method=cluster_method, n_clusters=n_clusters, phase=phase, **cluster_kwargs)
            title_suffix = f"k={n_clusters}"




        plt.figure(figsize=(10, 6))
        scatter = plt.scatter(pca_result[:, 0], pca_result[:, 1], c=labels, cmap='tab10', alpha=0.7)
        plt.xlabel(f'PC1 ({explained_variance[0]*100:.1f}%)')
        plt.ylabel(f'PC2 ({explained_variance[1]*100:.1f}%)')
        plt.title(f'PCA + {cluster_method.capitalize()} ({title_suffix}) - Phase: {phase if phase else "All"}')
        plt.colorbar(scatter, label='Cluster Label')
        plt.grid(True, alpha=0.3)
        plt.savefig(f"figures/pca_{cluster_method}_{phase if phase else 'all'}.png", dpi=200, bbox_inches='tight')
        plt.close()

    def plot_tsne_clustered(self,
                            n_clusters: Optional[int] = None,
                            phase: Optional[str] = None,
                            cluster_method: str = "agglomerative",
                            n_iter: int = 1000,
                            random_state: int = 42,
                            colored_by: str = "cluster",
                            **cluster_kwargs) -> None:
        """
        Plots t-SNE projection colored by cluster labels or energy.
        4 subplots with perplexities [20, 40, 60, 80].
        """
        if cluster_method in ["dbscan", "hdbscan"]:
            labels = self.cluster(method=cluster_method, phase=phase, **cluster_kwargs)
            title_suffix = f"clusters={len(set(labels) - {-1})}, noise={np.sum(labels == -1)}"
        else:
            if n_clusters is None:
                raise ValueError("n_clusters must be specified for agglomerative or kmeans clustering")
            labels = self.cluster(method=cluster_method, n_clusters=n_clusters, phase=phase, **cluster_kwargs)
            title_suffix = f"k={n_clusters}"

        perplexities = [20, 40, 60, 80]
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        for ax, perplexity in zip(axes.flatten(), perplexities):
            tsne_result = self.tsne(n_components=2, perplexity=perplexity, n_iter=n_iter, random_state=random_state, phase=phase)
            if colored_by == "cluster":
                scatter = ax.scatter(tsne_result[:, 0], tsne_result[:, 1], c=labels, cmap='tab10', alpha=0.7)
                plt.colorbar(scatter, ax=ax, label='Cluster Label')
            elif colored_by == "energy":
                energies = np.array([s.energy for s in (self.phases[phase] if phase else self.structures)])
                scatter = ax.scatter(tsne_result[:, 0], tsne_result[:, 1], c=energies, cmap='viridis', alpha=0.7)
                plt.colorbar(scatter, ax=ax, label='Energy')
            ax.set_title(f't-SNE (perp={perplexity}) - {colored_by}', fontsize=12)
            ax.set_xlabel('t-SNE Dim 1')
            ax.set_ylabel('t-SNE Dim 2')
            ax.grid(True, alpha=0.3)
        plt.suptitle(f"t-SNE + {cluster_method.capitalize()} ({title_suffix}) - Phase: {phase if phase else 'All'}", fontsize=14, fontweight='bold')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(f"figures/tsne_{cluster_method}_{phase if phase else 'all'}_{colored_by}.png", dpi=200, bbox_inches='tight')
        plt.close()

    def plot_umap_clustered(self, 
                            n_clusters: Optional[int] = None,
                            phase: Optional[str] = None,
                            cluster_method: str = "agglomerative",
                            random_state: int = 42,
                            colored_by: str = "cluster",
                            n_neighbors: int = 15,
                            **cluster_kwargs) -> None:
        """
        Plots UMAP projection colored by cluster labels or energy.
        Uses a single plot with a specified n_neighbors.
        Noise Points (if present with label -1) are plottet in light gray with low opacitc
        """
        if UMAP is None:
            raise ImportError("UMAP is not installed.")
        
        if cluster_method in ["dbscan", "hdbscan"]:
            labels = self.cluster(method=cluster_method, phase=phase, **cluster_kwargs)
            title_suffix = f"clusters={len(set(labels) - {-1})}, noise={np.sum(labels == -1)}"
        else:
            if n_clusters is None:
                raise ValueError("n_clusters must be specified for agglomerative or kmeans clustering")
            labels = self.cluster(method=cluster_method, n_clusters=n_clusters, phase=phase, **cluster_kwargs)
            title_suffix = f"k={n_clusters}"
        
        plt.figure(figsize=(10, 6))
        umap_result = self.umap(n_components=2, n_neighbors=n_neighbors, random_state=random_state, phase=phase)
        
        if colored_by == "cluster": 
            # Mask noise points and core points
            noise_mask = labels == -1
            core_mask = ~noise_mask

            # Plot noise points first 
            if np.any(noise_mask):
                plt.scatter(umap_result[noise_mask, 0], umap_result[noise_mask, 1], c='lightgray',
                            alpha=0.8, label='Noise', edgecolor='none')
            if np.any(core_mask):
                scatter = plt.scatter(umap_result[core_mask, 0], umap_result[core_mask, 1], c=labels[core_mask],
                                      cmap='tab10', alpha=0.7, edgecolor='k', linewidth=0.5)
                plt.colorbar(scatter, label='Cluster Label')
        
        elif colored_by == "energy":
            energies = np.array([s.energy for s in (self.phases[phase] if phase else self.structures)])
            scatter = plt.scatter(umap_result[:, 0], umap_result[:, 1], c=energies, cmap='viridis', alpha=0.7)
            plt.colorbar(scatter, label='Energy')
        
        plt.title(f'UMAP + {cluster_method.capitalize()} ({title_suffix}) - Phase: {phase if phase else "All"}', fontsize=12)
        plt.xlabel('UMAP Dim 1')
        plt.ylabel('UMAP Dim 2')
        plt.grid(True, alpha=0.3)
        plt.savefig(f"figures/umap_{cluster_method}_{phase if phase else 'all'}_{colored_by}.png", dpi=200, bbox_inches='tight')
        plt.close()

    def plot_umap_3D_clustered(self,
                              n_clusters: Optional[int] = None,
                              phase: Optional[str] = None,
                              cluster_method: str = "agglomerative",
                              random_state: int = 42,
                              n_neighbors: int = 15,
                              **cluster_kwargs) -> None:
        """  
        Plots a 3D UMAP projection:
        x-y are the first two UMAP dimensions and z is the energy.
        The points are colored by cluster labels.

        Args:
            n_clusters: Number of clusters to form (only for agglomerative and kmeans)
            phase: Specified phase to analyze. If None, analyzes all structures.
            cluster_method: Clustering method ("agglomerative", "kmeans", "dbscan", "hdbscan")
            n_neighbors: Number of neighbors for UMAP
            **cluster_kwargs: Passed to the clustering method
        """
        if UMAP is None:
            raise ImportError("UMAP is not installed.")
        if cluster_method in ["dbscan", "hdbscan"]:
            labels = self.cluster(method=cluster_method, phase=phase, **cluster_kwargs)
            title_suffix = f"clusters={len(set(labels) - {-1})}, noise={np.sum(labels == -1)}"
        else:
            if n_clusters is None:
                raise ValueError("n_clusters must be specified for agglomerative or kmeans clustering")
            labels = self.cluster(method=cluster_method, n_clusters=n_clusters, phase=phase, **cluster_kwargs)
            title_suffix = f"k={n_clusters}"
        umap_result = self.umap(n_components=2, n_neighbors=n_neighbors, random_state=random_state, phase=phase)
        energies = np.array([s.energy for s in (self.phases[phase] if phase else self.structures)])
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(umap_result[:, 0], umap_result[:, 1], energies, c=labels, cmap='tab10', alpha=0.7)
        ax.set_xlabel('UMAP Dim 1')
        ax.set_ylabel('UMAP Dim 2')
        ax.set_zlabel('Energy (Ha)')
        ax.set_title(f'3D UMAP + {cluster_method.capitalize()} ({title_suffix}) - Phase: {phase if phase else "All"}', fontsize=12)
        plt.colorbar(scatter, label='Cluster Label')
        plt.grid(True, alpha=0.3)
        plt.savefig(f"figures/umap_3d_{cluster_method}_{phase if phase else 'all'}.png", dpi=200, bbox_inches='tight')
        plt.close()



    def log_representative_features(
            self,
            representatives: List[StructureData],
            labels: Optional[np.ndarray] = None,
            output_file: str = "cluster_features.log"
        ) -> None:
        """ 
        Log the feature matrix rows for the given representative structures

        Args:
            representatives: List of StructureData objects representing the cluster centers
            labels: Cluster labels array to identify which representative belongs to which cluster. If None, uses cached labels.
        """
        self._ensure_feature_matrix()
        raw = self._feature_matrix_raw
        normalized = self._feature_matrix_normalized

        rep_indices = []
        for rep in representatives:
            try:
                idx = self.structures.index(rep)
                rep_indices.append(idx)
            except ValueError:
                if self.logger:
                    self._log(f"Representative structure not found in structures list: {rep}", level="warning")
                continue
        
        if not rep_indices:
            if self.logger:
                self._log("No valid representative indices found. Skipping logging.", level="warning")
            return
        
        if labels is not None:
            row_labels = [f"Cluster {labels[idx]}" for idx in rep_indices]
        else:
            row_labels = [f"Rep {i}" for i in range(len(rep_indices))]

        col_labels = self.get_feature_names()
        
        rep_matrix_raw = raw[rep_indices]
        rep_matrix_norm = normalized[rep_indices]

        # Compute maximum variane
        variances = np.var(rep_matrix_raw, axis=0)
        max_var_idx = np.argmax(variances)
        most_changed_feature = col_labels[max_var_idx] if max_var_idx < len(col_labels) else f"Feature_{max_var_idx}"
        max_variance_val = variances[max_var_idx]

        if self.logger:
            self._log(f"Writing representative features to {output_file}")
            self._log(f"Most varied feature among representatives: {most_changed_feature} (variance: {max_variance_val:.4f})")

        with open(output_file, "w") as f:
            f.write("=== Cluster Representatives Feature Matrix ===\n")
            f.write(f"Most prominent discriminator: {most_changed_feature} (variance: {max_variance_val:.4f})\n")
            
            header = f"{'Cluster':<15} | " + " | ".join([f"{col:<15}" for col in col_labels])

            f.write("--- Raw Feature Matrix ---\n")
            f.write(header + "\n")
            f.write("-" * len(header) + "\n")
            for row_label, row_data in zip(row_labels, rep_matrix_raw):
                row_str = f"{row_label:<15} | " + " | ".join([f"{val:<15.4f}" for val in row_data])
                f.write(row_str + "\n")
            
            f.write("\n--- Normalized Feature Matrix ---\n")
            f.write(header + "\n")
            f.write("-" * len(header) + "\n")
            for row_label, row_data in zip(row_labels, rep_matrix_norm):
                row_str = f"{row_label:<15} | " + " | ".join([f"{val:<15.4f}" for val in row_data])
                f.write(row_str + "\n")

    # Benchmark on Speed on Clustering Methods:
    

    
    # =============================== Noise Assesment in Dataset ===============================

    def add_noise(
            self,
            data_set: np.ndarray,
            random_state: Optional[int] = 42
        ) -> np.ndarray:
        """ 
        Adds artificial uniform noise per dimension and appends it to the dataset.


        Mapping:
        - a = min(Dataset) per dimension
        - b = max(Dataset) per dimension
        - Noise = (b-a) * rand(n,m) + a
        - DataSet = Dataset.join(Noise)
        """
        if data_set.ndim != 2:
            raise ValueError("Input data_set must be a 2D array (n_samples, n_features)")
        n,m = data_set.shape 
        rng = np.random.default_rng(seed=random_state)
        a = np.min(data_set, axis=0)
        b = np.max(data_set, axis=0)
        noise = (b-a) * rng.random((n, m)) + a
        return np.vstack([data_set, noise])
    
    def get_distance_to_mean(self, data_set: np.ndarray) -> np.ndarray:
        """ 
        Computes sorted Euclidean distance of each sample to the dataset mean
        """
        if data_set.ndim != 2:
            raise ValueError("Input data_set must be a 2D array (n_samples, n_features)")
        if data_set.shape[0] == 0:
            return np.array([])
        mean_vec = np.mean(data_set, axis=0)
        dist_mean = np.linalg.norm(data_set - mean_vec, axis=1)
        return np.sort(dist_mean)
    

    def noise_assessment(
            self,
            data_set: Optional[np.ndarray] = None,
            bins: int = 50,
            random_state: Optional[int] = 42,
            save_path: str = "figures/noise_assessment_dist2mean.png"
        ) -> Dict[str, Any]:
        """ 
        Assessment of noise by comparing distance-to-mean distributions before and after adding
        synthetic uniform noise.
        
        Steps:
        1) z-score standardization
        2) add uniform noise and append to dataset
        3) compute sorted eucliden dist-to-mean for both
        4) compare histograms
        """
        self._log("---- Noise Assessment Start ----")

        # if no external dataset is provided we use the raw feature matrix from analyzer
        if data_set is None:
            x = self.feature_matrix(normalize=False)
            self._log("Using internal raw feature matrix for noise assessment.")
        else:
            x = np.asarray(data_set, dtype=np.float64)
            self._log("Using provided dataset for noise assessment.")
        if x.ndim != 2:
            raise ValueError("Input data_set must be a 2D array (n_samples, n_features)")
        if x.shape[0] < 2:
            self._log("Not enough samples for noise assessment. Need at least 2 samples.", level="warning")
            return {}
        
        # z score normalizaton
        mu = np.mean(x, axis=0)
        sigma = np.std(x, axis=0)
        sigma_safe = np.where(sigma < 1e-12, 1.0, sigma)
        x_norm = (x - mu) / sigma_safe

        # Add synthetic noise and compute distance
        x_noise = self.add_noise(x_norm, random_state=random_state)
        dist_clean = self.get_distance_to_mean(x_norm)
        dist_noise = self.get_distance_to_mean(x_noise)

        # Build comparable histrogram densities
        all_vals = np.concatenate([dist_clean, dist_noise])
        hist_range = float(np.min(all_vals)), float(np.max(all_vals))
        hist_clean, edges = np.histogram(dist_clean, bins=bins, range=hist_range, density=True)
        hist_noise, _ = np.histogram(dist_noise, bins=bins, range=hist_range, density=True)
        bin_width = edges[1] - edges[0] if len(edges) > 1 else 1.0
        overlap = float(np.sum(np.minimum(hist_clean, hist_noise)) * bin_width)

        # Simple impact metrics
        mean_shift = float(np.mean(dist_noise) - np.mean(dist_clean))
        std_ratio = float(np.std(dist_noise) / np.std(dist_clean)) if np.std(dist_clean) > 1e-12 else float('inf')
        q95_shift = float(np.percentile(dist_noise, 95) - np.percentile(dist_clean, 95))

        plt.figure(figsize=(10, 6))
        plt.hist(dist_clean, bins=bins, range=hist_range, density=True, alpha=0.6, label='Clean Data', color='blue')
        plt.hist(dist_noise, bins=bins, range=hist_range, density=True, alpha=0.6, label='With Noise', color='red')
        plt.xlabel('Distance to Mean')
        plt.ylabel('Density')
        plt.title('Noise Assessment: Distance to Mean Distribution')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(save_path, dpi=200)
        plt.close()

        self._log(f"Noise Assessment Results:")
        self._log(f"Mean Distance Shift: {mean_shift:.4f}")
        self._log(f"Std Ratio (Noise/Clean): {std_ratio:.4f}")
        self._log(f"95th Percentile Shift: {q95_shift:.4f}")
        self._log(f"Distribution Overlap: {overlap:.4f}")
        self._log("---- Noise Assessment End ----")

    # ================================ Testing functions ================================

    def test_chunksize_effect(self, angle_threshold=150.0, max_distance=3.5, n_processes=None):
        """  
        Testing function how different chunksizes affect the H-bond feature calculation
        """
        n_structures = len(self.structures)
        if n_structures == 0:
            self._log("No structures available for testing chunksize effect.", level="warning")
            return
        
        if n_processes is None:
            n_processes = max(1, cpu_count() - 2)
        n_processes = min(n_processes, n_structures)

        chunksizes = [1, 2, 4, 8, 16, 32, 64, max(1, math.ceil(n_structures / n_processes * 4))]
        results = []
        print(f"Testing chunksize effect on H-bond feature calculation with {n_structures} structures and {n_processes} processes...")
        worker_args = [
            (
                list(s.molecule.atom_labels),
                np.array(s.molecule.coordinates),
                angle_threshold,
                max_distance
            )
            for s in self.structures
        ]
        times = []
        for chunksize in chunksizes:
            start = time.time()
            with Pool(processes=n_processes) as pool:
                pool.map(_compute_hbond_single, worker_args, chunksize=chunksize)
            end = time.time()
            elapsed = end - start
            times.append(elapsed)

        # plot results
        plt.figure(figsize=(10, 6))
        plt.plot(chunksizes, times, marker='o')
        plt.xscale('log')
        plt.xlabel('Chunksize')
        plt.ylabel('Time (seconds)')
        plt.title(f'Effect of Chunksize on H-bond Feature Calculation (n_structures={n_structures}, processes={n_processes})')
        plt.grid(True, alpha=0.3)
        plt.savefig("figures/chunksize_effect.png", dpi=200)
        plt.close()

        