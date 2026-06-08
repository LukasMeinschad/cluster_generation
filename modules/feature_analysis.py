""" 
Module for statistical preprocessing, correlation analysis
"""
import os
from typing import List, Tuple, Optional, Dict, Any
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import rankdata
from scipy.spatial.distance import pdist, squareform



def feature_statistics(cluster_class: Any) -> Dict[str, Dict[str, float]]:
    """  
    Computes row/column summary statistics for the active feature matrix
    """
    X = cluster_class.feature_matrix
    n_samples, n_features = X.shape
    cluster_class.log(f"Computing feature statistics for {n_samples} samples and {n_features} features.")
    
    feature_stats = {}
    for i in range(n_features):
        col_data = X[:, i]
        mean_val = float(np.mean(col_data))
        std_val = float(np.std(col_data))
        feature_stats[f"feature_{i}"] = {"mean": mean_val, "std": std_val}
        cluster_class.log(f"  -> Feature {i}: mean={mean_val:.4f}, std={std_val:.4f}")

    cluster_class.feature_stats = feature_stats
    return feature_stats

def spearman_correlation(
        cluster_class: Any,
        plot: bool = False,
        save_dir: str = "figures"
    ) -> np.ndarray:
    """ 
    Computes the Spearman Rank Correlation matrix using a
    optimized rank transform approach
    """
    X = cluster_class.feature_matrix
    if X is None or X.size == 0:
        raise ValueError("Feature matrix is empty or not set in the cluster_class.")
    if X.ndim == 1:
        cluster_class.log("Feature matrix is 1-D; returning trivial correlation matrix of 1.0.")
        corr = np.array([[1.0]])
        cluster_class.corr_matrix = corr
        return corr
    
    cluster_class.log(f"Computing Spearman Rank Correlation for {X.shape[1]} features.")
    
    # Vectorized rank transform 
    ranked = np.apply_along_axis(rankdata, 0, X)
    corr = np.corrcoef(ranked, rowvar=False)
    corr = np.nan_to_num(corr, nan=0.0, posinf=1.0, neginf=-1.0)

    # Cache inside tracked context space
    cluster_class.corr_matrix = corr

    if plot:
        os.makedirs(save_dir, exist_ok=True)
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr, annot=True, fmt=".2f", cmap="coolwarm", vmin=-1, vmax=1)
        plt.title("Spearman Rank Correlation    Matrix")
        plt.tight_layout()
        save_path = os.path.join(save_dir, "spearman_correlation_matrix.png")
        plt.savefig(save_path)
        plt.close()
        cluster_class.log(f"Spearman correlation matrix plot saved to: {save_path}")

def filter_correlation_spearman(
        cluster_class: Any,
        threshold: float = 0.8
    ) -> List[int]:
    """ 
    Identifies and removes redundant, highly correlated descriptor columns to combat collinearity and reduce dimensionality
    """
    corr_matrix = getattr(cluster_class, 'corr_matrix', None)
    if corr_matrix is None:
        corr_matrix = spearman_correlation(cluster_class, plot=False)
    n_features = corr_matrix.shape[0]
    upper_tri_indices = np.triu_indices(n_features, k=1)

    high_corr_mask = np.abs(corr_matrix[upper_tri_indices]) > threshold
    cols_to_drop = np.unique(upper_tri_indices[1][high_corr_mask])
    selected_features = [i for i in range(n_features) if i not in cols_to_drop]

    # Mutate the context feature matric in place
    cluster_class.feature_matrix = cluster_class.feature_matrix[:, selected_features]
    cluster_class.corr_matrix = corr_matrix[np.ix_(selected_features, selected_features)]
    cluster_class.log(f"Filtered out {len(cols_to_drop)} highly correlated features (threshold={threshold}). Remaining features: {len(selected_features)}.")
    return selected_features

def determine_distance_metric(context: Any) -> str:
    """
    Heuristic decision tree that samples feature space overlaps to select 
    the optimal distance metric topology for downstream clustering.
    """
    corr_matrix = getattr(context, "corr_matrix", None)
    if corr_matrix is None:
        corr_matrix = spearman_correlation(context, plot=False)

    if corr_matrix.shape[0] <= 1:
        context.metric = "euclidean"
        return "euclidean"

    # Extract upper triangular values to measure overall feature interdependence
    triu_vals = corr_matrix[np.triu_indices_from(corr_matrix, k=1)]
    avg_cor = float(np.mean(np.abs(triu_vals))) if triu_vals.size > 0 else 0.0

    # Branch 1: Dense global collinearity requires covariance normalization
    if avg_cor > 0.5:
        context.log(f"High feature correlation detected (avg_cor={avg_cor:.4f}). Routing to Mahalanobis metric.")
        context.metric = "mahalanobis"
        return "mahalanobis"

    # Branch 2: Low redundancy -> find metric maximizing spatial separation contrast
    X = context.feature_matrix
    results = {}
    for metric_candidate in ["cityblock", "euclidean", "cosine"]:
        try:
            dists = pdist(X, metric=metric_candidate)
            if dists.size > 0 and np.mean(dists) > 0:
                results[metric_candidate] = (np.max(dists) - np.min(dists)) / np.mean(dists)
            else:
                results[metric_candidate] = -1.0
        except Exception:
            results[metric_candidate] = -1.0

    best_metric = max(results, key=results.get)
    context.log(
        f"Low feature correlation (avg_cor={avg_cor:.4f}). "
        f"Selected '{best_metric}' based on contrast ratio maximization."
    )
    
    context.metric = best_metric
    return best_metric

if __name__ == "__main__":
    import modules.cluster as cluster_module
    # Example usage and testing of feature analysis functions
    mock_data = np.random.rand(100, 20)  # 100 samples, 20 features
    mock_cluster = Clustering(feature_matrix=mock_data, energies=None, normalize=True)
    stats = feature_statistics(mock_cluster)
    corr_matrix = spearman_correlation(mock_cluster, plot=True)
    selected_features = filter_correlation_spearman(mock_cluster, threshold=0.8)
    best_metric = determine_distance_metric(mock_cluster)
    print(f"Selected distance metric: {best_metric}")
