""" 
Clustering performance metrics and evaluation
Handles noise-aware computations for Silhouette Score, Davies-Bouldin Index, and Calinski-Harabasz Index
"""
import numpy as np
from typing import Optional, Any, Tuple, Dict
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score

def evaluate_all_metrics(
        cluster_class: Any,
        feature_matrix: Optional[np.ndarray] = None,
        labels: Optional[np.ndarray] = None,
        ignore_noise: bool = True
    ) -> Dict[str, float]:
    """ 
    Computes a comprehensive suite of unsupervised evaluation metrics
    Updates the active Clustering container and handles noise filtering for accurate performance assessment.
    
    Args:
        cluster_class: The Clustering instance containing the feature matrix and labels to be evaluated
        feature_matrix: Optional feature matrix to use for evaluation (defaults to cluster_class.feature_matrix)
        labels: Optional cluster labels to evaluate (defaults to cluster_class.labels)
        ignore_noise: Whether to exclude noise points (label -1) from metric calculations
    """
    X = feature_matrix if feature_matrix is not None else getattr(cluster_class, 'feature_matrix', None)
    y = labels if labels is not None else getattr(cluster_class, 'labels', None)
    if X is None or y is None:
        raise ValueError("Feature matrix and labels must be provided either as arguments or as attributes of the cluster_class.")

    # Get the distance metric
    metric = getattr(cluster_class, 'metric', 'euclidean')

    # Optionally filter density based noise points (label -1) from evaluation
    X_clean, y_clean, n_clusters = _prepare_and_filter_labels(X, y, ignore_noise)

    # Guard against insufficient clusters
    if n_clusters < 2:
        cluster_class.log("Not enough clusters to compute metrics. Returning NaN for all metrics.")
        return {"silhouette:": None, "davies_bouldin": None, "calinski_harabasz": None}
    
    cluster_class.log(f"Evaluating metrics on {X_clean.shape[0]} samples across {n_clusters} clusters using metric: {metric}")
    try:
        sil = float(silhouette_score(X_clean, y_clean, metric=metric))
    except Exception as e:
        cluster_class.log(f"Error computing Silhouette Score: {e}")
        sil = None
    try:
        db = float(davies_bouldin_score(X_clean, y_clean))
    except Exception as e:
        cluster_class.log(f"Error computing Davies-Bouldin Index: {e}")
        db = None
    try:
        ch = float(calinski_harabasz_score(X_clean, y_clean))
    except Exception as e:
        cluster_class.log(f"Error computing Calinski-Harabasz Index: {e}")
        ch = None

    # Log internal execution profile metrics
    cluster_class.log("Internal Cluster Validation Statistics:")
    if sil is not None:
        cluster_class.log(f"  -> Silhouette Score: {sil:.4f} (Higher is better, range [-1, 1])")
    else:
        cluster_class.log("  -> Silhouette Score: Not computable")
    if db is not None:
        cluster_class.log(f"  -> Davies-Bouldin Index: {db:.4f} (Lower is better, range [0, inf])")
    else:
        cluster_class.log("  -> Davies-Bouldin Index: Not computable")
    if ch is not None:
        cluster_class.log(f"  -> Calinski-Harabasz Index: {ch:.4f} (Higher is better, range [0, inf])")
    else:
        cluster_class.log("  -> Calinski-Harabasz Index: Not computable")

    scores = {"silhouette": sil, "davies_bouldin": db, "calinski_harabasz": ch}
    cluster_class.clustering_metrics = scores
    return scores


def calculate_wcss_per_cluster(cluster_class, labels: Optional[np.ndarray] = None) -> Dict[int, float]:
    """ 
    Calculates the Within-cluster-sum-of-squares which is defined via

    WCSS = sum(||x_i - c_k||^2) for all x_i in cluster k, where c_k is the centroid of cluster k
    """ 
    lbl = labels if labels is not None else getattr(cluster_class, 'labels', None)
    X = cluster_class.feature_matrix
    if lbl is None or X is None:
        raise ValueError("Labels and feature matrix must be provided either as arguments or as attributes of the cluster_class.")
    if X is None or X.size == 0:
        raise ValueError("Feature matrix is empty or not set in the cluster_class.")
    wcss_per_cluster = {}
    for cluster_id in np.unique(lbl):
        if cluster_id == -1:
            continue  # Skip noise points
        cluster_points = X[lbl == cluster_id]
        centroid = np.mean(cluster_points, axis=0)
        wcss = np.sum(np.linalg.norm(cluster_points - centroid, axis=1) ** 2)
        wcss_per_cluster[cluster_id] = wcss
        cluster_class.log(f"Cluster {cluster_id}: WCSS = {wcss:.4f} with {cluster_points.shape[0]} points.")
    cluster_class.wcss_per_cluster = wcss_per_cluster
    return wcss_per_cluster

def _prepare_and_filter_labels(
    X: np.ndarray, 
    y: np.ndarray, 
    ignore_noise: bool
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Internal data sanitation utility. Removes unclustered noise points (-1) 
    so standard mathematical validation metrics evaluate correctly.
    """
    unique_labels = np.unique(y)
    has_noise = -1 in unique_labels

    if ignore_noise and has_noise:
        non_noise_mask = (y != -1)
        X_filtered = X[non_noise_mask]
        y_filtered = y[non_noise_mask]
        n_clusters = len(unique_labels) - 1
        return X_filtered, y_filtered, n_clusters

    n_clusters = len(unique_labels)
    return X, y, n_clusters

if __name__ == "__main__":
    # Example usage with synthetic data
    from sklearn.datasets import make_blobs
    
    X, y = make_blobs(n_samples=300, centers=4, n_features=10, random_state=42)
    
    class DummyClusterClass:
        def log(self, message):
            print(message)
    
    cluster_class = DummyClusterClass()
    metrics = evaluate_all_metrics(cluster_class, feature_matrix=X, labels=y)
    print(f"Computed Metrics: {metrics}")