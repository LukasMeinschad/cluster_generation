""" 
Module for identifying and filtering anomalous or unphysical molecular configurations
using different outlier detection methods
"""

from typing import Any, Optional, Tuple, Dict, List
import numpy as np
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor




def z_score_outlier_detection(
        cluster_class: Any,
        threshold: float = 3.0,
        prune: bool = True
    ) -> List[int]:
    """  
    Filter out samples that have any feature with a z-score above the given threshold.
    This is used for the first data cleanup before clustering
    """
    X = cluster_class.feature_matrix
    cluster_class.log(f"Running Z-score outlier detection with threshold={threshold} on {X.shape[0]} samples and {X.shape[1]} features.")
    mean = np.mean(X, axis=0)
    std = np.std(X, axis=0)
    z_scores = np.abs((X - mean) / std)
    outlier_mask = np.any(z_scores > threshold, axis=1)
    n_outliers = np.sum(outlier_mask)
    cluster_class.log(f"Z-score method detected {n_outliers} outliers out of {X.shape[0]} samples.")
    cluster_class.outlier_mask = outlier_mask
    if prune:
        _prune_context(cluster_class)
    
    return outlier_mask


def run_isolation_forest(
        cluster_class: Any,
        contamination: float = 0.05,
        n_estimators: int = 100,
        bootstrap: bool = False,
        prune: bool = True,
        random_state: Optional[int] = None
    ) -> np.ndarray:
    """ 
    Detects anomalies using an Isolation Forest by isolating features through random partitioning.
    Good for catching global outliers in high-dimensional spaces.

    Args:
        cluster_class: The clustering context object containing the feature matrix and logging method
        contamination: The proportion of outliers in the data set
        n_estimators: The number of base estimators in the ensemble
        bootstrap: Whether to use bootstrap samples when building trees
        prune: Whether to remove detected outliers from the feature matrix and update context
        random_state: Controls the randomness of the estimator
    """
    X = cluster_class.feature_matrix
    cluster_class.parameters("Isolation Forest", [
        ("contamination", contamination),
        ("n_estimators", n_estimators),
        ("bootstrap", bootstrap),
        ("random_state", random_state),
    ])
    iso_forest = IsolationForest(
        n_estimators=n_estimators,
        contamination=contamination,
        bootstrap=bootstrap,
        random_state=random_state,
        n_jobs=-1
    )
    preds = iso_forest.fit_predict(X)
    scores = iso_forest.decision_function(X)

    # Initialize or update model registry in cluster context
    if not hasattr(cluster_class, "outlier_models"):
        cluster_class.outlier_models = {}
    cluster_class.outlier_models['isolation_forest'] = iso_forest

    # Map to boolean mask where True = Outlier
    outlier_mask = (preds == -1)
    n_outliers = np.sum(outlier_mask)
    cluster_class.substep(f"Detected {n_outliers} outliers out of {X.shape[0]} samples")

    cluster_class.outlier_mask = outlier_mask
    cluster_class.anomaly_scores = scores

    if prune:
        _prune_context(cluster_class)

    return outlier_mask

def run_local_outlier_factor(
        cluster_class: Any,
        contamination: float = 0.05,
        n_neighbors: int = 20,
        prune: bool = True
    ) -> np.ndarray:
    """ 
    Detects anomalies using Local Outlier Factor (LOF) by measuring local density deviations.
    Compares the local denisty of a point to that of its neighbors, making it effective for identifying local outliers in clustered data.
    
    Args:
        cluster_class: The clustering context object containing the feature matrix and logging method
        contamination: The proportion of outliers in the data set
        n_neighbors: Number of neighbors to use for local density estimation
        prune: Whether to remove detected outliers from the feature matrix and update context
    """

    X = cluster_class.feature_matrix
    metric = getattr(cluster_class, "metric", "cityblock")
    cluster_class.parameters("Local Outlier Factor", [
        ("contamination", contamination),
        ("n_neighbors", n_neighbors),
        ("metric", metric),
    ])
    lof = LocalOutlierFactor(
        n_neighbors=n_neighbors,
        contamination=contamination,
        metric=metric,
        n_jobs=-1
    )
    preds = lof.fit_predict(X)
    scores = lof.negative_outlier_factor_

    if not hasattr(cluster_class, "outlier_models"):
        cluster_class.outlier_models = {}
    cluster_class.outlier_models['local_outlier_factor'] = lof

    # Map to boolean mask where True = Outlier
    outlier_mask = (preds == -1)
    n_outliers = np.sum(outlier_mask)
    cluster_class.substep(f"Detected {n_outliers} outliers out of {X.shape[0]} samples")
    cluster_class.outlier_mask = outlier_mask
    cluster_class.anomaly_scores = scores
    if prune:
        _prune_context(cluster_class)
    return outlier_mask
    
def _prune_context(context: Any) -> None:
    """
    Internal utility to slice context variables in place, erasing flagged anomalies.
    """
    outlier_mask = context.outlier_mask
    inlier_mask = ~outlier_mask
    
    old_shape = context.feature_matrix.shape

    context.feature_matrix = context.feature_matrix[inlier_mask]

    if hasattr(context, "_feature_matrix_normalized") and context._feature_matrix_normalized is not None:
        context._feature_matrix_normalized = context._feature_matrix_normalized[inlier_mask]

    if hasattr(context, "molecules") and context.molecules is not None:
        # Cleanly synchronize also the molecule list if it exists
        # Inlier mask is boolean, so we can use it to filter the list of molecules
        context.molecules = [mol for i, mol in enumerate(context.molecules) if inlier_mask[i]]



    if hasattr(context, "energies") and context.energies is not None:
        if isinstance(context.energies, np.ndarray):
            context.energies = context.energies[inlier_mask]
        else:
            context.energies = [e for i, e in enumerate(context.energies) if inlier_mask[i]]
    
    # Cleanly synchronize secondary parallel tracking attributes if they exist
    if hasattr(context, "labels") and context.labels is not None:
        context.labels = context.labels[inlier_mask]
        
    if hasattr(context, "structures") and context.structures is not None:
        context.structures = [s for i, s in enumerate(context.structures) if inlier_mask[i]]

    context.substep(
        f"Context pruning complete: {old_shape[0]} → {context.feature_matrix.shape[0]} structures",
        level=2,
    )

if __name__ == "__main__":
    from modules.cluster import Clustering
    # Example usage and testing of outlier detection functions
    mock_data = np.random.rand(100, 10)  # 100 samples, 10 features
    mock_cluster = Clustering(feature_matrix=mock_data, energies=None, normalize=True)
    outliers_iso = run_isolation_forest(cluster_class=mock_cluster, contamination=0.1, prune=True)
    outliers_lof = run_local_outlier_factor(cluster_class=mock_cluster, contamination=0.1, prune=True)

