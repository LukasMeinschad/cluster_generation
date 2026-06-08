""" 
Module for unsupervised partitioning and density-based clustering algorithms
optimized for high-dimensional feature spaces and large datasets
"""
import numpy as np
from typing import Any, Optional, Dict, List, Tuple, Callable
from sklearn.cluster import KMeans, DBSCAN
from sklearn.cluster import AgglomerativeClustering as SklearnAgglomerativeClustering
import warnings
from scipy.spatial.distance import cdist


try:
    from hdbscan import HDBSCAN
except ImportError:
    HDBSCAN = None
    print("HDBSCAN is not installed. Please install it with 'pip install hdbscan' to use HDBSCAN clustering.")

def kmeans_clustering(
        feature_matrix: np.ndarray,
        n_clusters: int = 5,
        n_init: int = 10,
        max_iter: int = 300,
        random_state: Optional[int] = None,
        log_fn: Optional[Callable[[str], None]] = None
    ) -> KMeans:
    """ 
    Performs K-Means clustering on the provided feature matrix

    Args:
        feature_matrix: The input data for clustering
        n_clusters: The number of clusters to form
        n_init: Number of time the k-means algorithm will be run with different centroid seeds
        max_iter: Maximum number of iterations of the k-means algorithm for a single run
        random_state: Determines random number generation for centroid initialization
        log_fn: Optional logging function for progress updates
    
    Returns:
        Tuple containing:
              - labels (np.ndarray): Cluster labels for each point in the dataset
              - metadata (dict): Metrics
    """ 
    km = KMeans(n_clusters=n_clusters, n_init=n_init, max_iter=max_iter, random_state=random_state)
    labels = km.fit_predict(feature_matrix)
    metadata = {"inertia": km.inertia_, "n_clusters": n_clusters}
    if log_fn:
        log_fn(f"K-Means clustering completed with inertia: {metadata['inertia']}")
    return labels, metadata

def agglomerative_clustering(
    feature_matrix: np.ndarray,
    n_clusters: int = 5,
    linkage: str = "ward",
    log_fn: Optional[Callable[[str], None]] = None
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Perform Agglomerative Hierarchical Clustering on the feature matrix.
    
    Args:
        feature_matrix: Array of shape (n_structures, n_features).
        n_clusters: Number of clusters to form.
        linkage: Linkage criterion to use ('ward', 'complete', 'average', 'single').
        log_fn: Optional callback function for logging.
        
    Returns:
        Tuple of (labels, metadata_dict).
    """
    ac = SklearnAgglomerativeClustering(n_clusters=n_clusters, linkage=linkage)
    labels = ac.fit_predict(feature_matrix)
    
    metadata = {"n_clusters": n_clusters, "linkage": linkage}
    
    if log_fn:
        log_fn(f"Agglomerative clustering completed: n_clusters={n_clusters}, linkage={linkage}")
        
    return labels, metadata


def dbscan_clustering(
    feature_matrix: np.ndarray,
    eps: float = 0.5,
    min_samples: int = 5,
    metric: str = "cityblock",
    log_fn: Optional[Callable[[str], None]] = None
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Perform DBSCAN clustering on the feature matrix.
    
    Args:
        feature_matrix: Array of shape (n_structures, n_features).
        eps: Maximum distance between two samples for neighborhood consideration.
        min_samples: Number of samples in a neighborhood for a core point.
        metric: Distance metric to use.
        log_fn: Optional callback function for logging.
        
    Returns:
        Tuple of (labels, metadata_dict).
    """
    dbscan = DBSCAN(eps=eps, min_samples=min_samples, metric=metric)
    labels = dbscan.fit_predict(feature_matrix)
    
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)
    
    metadata = {"n_clusters": n_clusters, "n_noise": n_noise, "eps": eps, "metric": metric}
    
    if log_fn:
        log_fn(f"DBSCAN clustering completed: eps={eps}, min_samples={min_samples}, "
               f"metric={metric}, n_clusters={n_clusters}, n_noise={n_noise}")
               
    return labels, metadata


def hdbscan_clustering(
    feature_matrix: np.ndarray,
    min_cluster_size: int = 5,
    min_samples: Optional[int] = None,
    metric: str = "cityblock",
    log_fn: Optional[Callable[[str], None]] = None
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Perform HDBSCAN clustering on the feature matrix.
    
    Args:
        feature_matrix: Array of shape (n_structures, n_features).
        min_cluster_size: Minimum size of clusters.
        min_samples: Number of samples in a neighborhood for a core point.
        metric: Distance metric to use.
        log_fn: Optional callback function for logging.
        
    Returns:
        Tuple of (labels, metadata_dict).
    """
    if HDBSCAN is None:
        raise ImportError("HDBSCAN is not installed. Run `pip install hdbscan` to use this method.")
        
    hdbscan = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples, metric=metric)
    labels = hdbscan.fit_predict(feature_matrix)
    
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)
    
    metadata = {"n_clusters": n_clusters, "n_noise": n_noise, "min_cluster_size": min_cluster_size, "metric": metric}
    
    if log_fn:
        log_fn(f"HDBSCAN clustering completed: min_cluster_size={min_cluster_size}, "
               f"min_samples={min_samples}, metric={metric}, n_clusters={n_clusters}, n_noise={n_noise}")
               
    return labels, metadata


def run_clustering(
    method: str,
    feature_matrix: np.ndarray,
    n_clusters: Optional[int] = None,
    metric: str = "cityblock",
    log_fn: Optional[Callable[[str], None]] = None,
    **kwargs
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Unified cluster dispatcher to easily run any supported algorithm.
    
    Args:
        method: Choice of 'kmeans', 'agglomerative', 'dbscan', or 'hdbscan'.
        feature_matrix: Array of shape (n_structures, n_features).
        n_clusters: Targeted cluster count (used for kmeans/agglomerative).
        metric: Distance metric string (used for density-based methods).
        log_fn: Optional logging handler.
        **kwargs: Forwarded directly to the underlying algorithm function.
    """
    dispatch = {
        "kmeans": lambda fm, **kw: kmeans_clustering(fm, n_clusters=n_clusters or 5, log_fn=log_fn, **kw),
        "agglomerative": lambda fm, **kw: agglomerative_clustering(fm, n_clusters=n_clusters or 5, log_fn=log_fn, **kw),
        "dbscan": lambda fm, **kw: dbscan_clustering(fm, metric=metric, log_fn=log_fn, **kw),
        "hdbscan": lambda fm, **kw: hdbscan_clustering(fm, metric=metric, log_fn=log_fn, **kw),
    }
    
    method_lower = method.lower()
    if method_lower not in dispatch:
        raise ValueError(f"Invalid clustering method: {method}. Choose from {list(dispatch.keys())}")
        
    return dispatch[method_lower](feature_matrix, **kwargs)

if __name__ == "__main__":
    # Example usage with synthetic data
    from sklearn.datasets import make_blobs
    
    X, _ = make_blobs(n_samples=300, centers=4, n_features=10, random_state=42)
    
    labels, metadata = run_clustering("kmeans", X, n_clusters=4, log_fn=print)
    print(f"Cluster labels: {labels}")
    print(f"Metadata: {metadata}")

def get_cluster_representatives(
    context: Any, 
    labels: Optional[np.ndarray] = None, 
    method: str = "centroid"
) -> Dict[int, int]:
    """
    Extracts one representative molecular configuration per partition group.

    Args:
        context: The central Clustering context container holding the feature matrix and settings.
        labels: Optional explicit label array override. If None, reads context.labels.
        method: The optimization criterion:
                "centroid" -> picks the point closest to the cluster mean in feature space.
                "medoid" -> picks the structure with the smallest average pairwise distance to 
                            all other group members (robust to non-convex descriptors).
                "lowest_energy" -> picks the configuration with the absolute lowest potential energy.

    Returns:
        Dict[int, int]: Mapping of cluster_label -> absolute row index in the original feature matrix.
                        Noise structures (label == -1) are skipped automatically.
    """
    # 1. Resolve label configuration
    lbl = labels if labels is not None else getattr(context, "labels", None)
    if lbl is None:
        raise ValueError("No cluster labels available. Run a clustering algorithm first or pass labels explicitly.")

    # 2. Guard for energetic selections
    energies = getattr(context, "energies", None)
    if method == "lowest_energy" and energies is None:
        raise ValueError("method='lowest_energy' requires an energy array to be provided at context construction.")

    fm = context.feature_matrix
    metric = getattr(context, "metric", "cityblock")
    representatives: Dict[int, int] = {}

    # 3. Step through unique partition labels
    unique_labels = np.unique(lbl)
    for label in unique_labels:
        if label == -1:  # Filter density-based noise flags (DBSCAN/HDBSCAN)
            continue
            
        mask = (lbl == label)
        cluster_features = fm[mask]
        cluster_indices = np.where(mask)[0]

        # 4. Apply geometric or thermodynamic optimization slice
        if method == "centroid":
            centroid = cluster_features.mean(axis=0)
            local_idx = int(np.argmin(np.linalg.norm(cluster_features - centroid, axis=1)))
        elif method == "medoid":
            dists = cdist(cluster_features, cluster_features, metric=metric)
            local_idx = int(np.argmin(np.mean(dists, axis=1)))
        elif method == "lowest_energy":
            local_idx = int(np.argmin(energies[cluster_indices]))
        else:
            raise ValueError(f"Unknown method '{method}'. Choose 'centroid', 'medoid', or 'lowest_energy'.")

        # Map back to the absolute global matrix row coordinates
        representatives[int(label)] = int(cluster_indices[local_idx])

    # Cache target representative registers on the tracking context
    context.representatives = representatives
    context.log(f"Cluster representatives identified via '{method}' optimization: {representatives}")
    
    return representatives 