"""
Unified dimensionality reduction modul
Handles PCA, Kernel PCA, UMAP etc...
"""
from typing import Optional, Any, Tuple, Dict, List
import numpy as np
from sklearn.decomposition import PCA, KernelPCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, AgglomerativeClustering




try:
    from umap import UMAP
except ImportError:
    UMAP = None
    print("UMAP is not installed. Please install it with 'pip install umap-learn' to use UMAP dimensionality reduction.")

def pca(cluster_class: Any, n_components: int = 2, **kwargs) -> PCA:
    """ 
    Applies linear principal component analysis to the active feature matrix

    Args:
        cluster_class: The Clustering instance containing the feature matrix to be reduced
        n_components: The number of principal components to compute
        **kwargs: Additional keyword arguments for PCA constructor

    Returns:
        Tuple containing:
           - embedding (np.ndarray): The reduced dimensionality representation of the data
           - explained_variance_ratio (np.ndarray): The variance explained by each of the selected components
    """
    cluster_class.log(f"Applying PCA with n_components={n_components} and parameters: {kwargs}")

    transfomer = PCA(n_components=n_components, **kwargs)
    embedding = transfomer.fit_transform(cluster_class.feature_matrix)
    explained_variance_ratio = transfomer.explained_variance_ratio_
    cluster_class.log(f"PCA completed. Explained variance ratio: {explained_variance_ratio}")
    return embedding, explained_variance_ratio

def kernel_based_pca(
        cluster_class: Any,
        n_components: int = 2,
        kernel: str = 'rbf',
        gamma: Optional[float] = None,
        n_jobs: int = -1
    ) -> Tuple[np.ndarray, np.ndarray]:
    """  
    Performs non-linear Kernel PCA to project configurations into a higher-dimensional space before applying linear PCA
    
    Returns:
        Tuple containing:
           - embedding (np.ndarray): The reduced dimensionality representation of the data
           - explained_variance_ratio (np.ndarray): The variance explained by each of the selected components
    """
    cluster_class.log(f"Applying Kernel PCA with n_components={n_components}, kernel={kernel}, gamma={gamma}, n_jobs={n_jobs}")
    kpca = KernelPCA(n_components=n_components, kernel=kernel, gamma=gamma, n_jobs=n_jobs)
    embedding = kpca.fit_transform(cluster_class.feature_matrix)
    
    # Approximated explained variance ratio from eigenvalues
    if kpca.eigenvalues_ is not None:
        explained_variance_ratio = kpca.eigenvalues_ / np.sum(kpca.eigenvalues_)
        cluster_class.log(f"Kernel PCA completed. Approximated explained variance ratio: {explained_variance_ratio}")
    else:
        explained_variance_ratio = np.array([])
        cluster_class.log("Kernel PCA completed. Eigenvalues not available, cannot compute explained variance ratio.")
    return embedding, explained_variance_ratio

def tsne(cluster_class: Any, n_components: int = 2, perplexity: float = 30.0, n_iter: int = 1000, random_state: Optional[int] = None)-> np.ndarray:
    """  
    Applies t-Distributed Stochastic Neighbor Embedding for non-linear dimensionality reduction
    """
    cluster_class.log(f"Applying t-SNE with n_components={n_components}, perplexity={perplexity}, n_iter={n_iter}, random_state={random_state}")
    tsne_model = TSNE(n_components=n_components, perplexity=perplexity, n_iter=n_iter, random_state=random_state)
    embedding = tsne_model.fit_transform(cluster_class.feature_matrix)
    cluster_class.log("t-SNE completed.")
    return embedding

def umap(cluster_class: Any, n_components: int = 2, n_neighbors: int = 15, min_dist: float = 0.1, random_state: Optional[int] = None) -> np.ndarray:
    """  
    Applies unsupervised Uniform Manifold Approximation and Projection (UMAP)
    Stores the fitted model inside the context for out-of-sample mapping and inverse transformation
    """
    _check_umap_availability()
    cluster_class.log(f"Applying UMAP with n_components={n_components}, n_neighbors={n_neighbors}, min_dist={min_dist}, random_state={random_state}")
    umap_model = UMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=cluster_class.umap_metric,
        random_state=random_state
    )
    embedding = umap_model.fit_transform(cluster_class.feature_matrix)
    # Cache the fitted umap model
    cluster_class.umap_model = umap_model
    return embedding

def umap_densmap(cluster_class: Any, n_neighbors: int = 15, min_dist: float= 0.1, dens_lambda: float = 0.8, random_state: Optional[int] = None) -> np.ndarray:
    """ 
    Runs DensMAP a density-preserving variant of UMAP that preserves local density
    variances from a high dimensional featuere matrix inside the low-dimensional embedding
    """
    _check_umap_availability()
    cluster_class.log(f"Applying DensMAP with n_neighbors={n_neighbors}, min_dist={min_dist}, dens_lambda={dens_lambda}, random_state={random_state}")
    umap_model = UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        dens_lambda=dens_lambda,
        metric=cluster_class.umap_metric,
        random_state=random_state
    )
    embedding = umap_model.fit_transform(cluster_class.feature_matrix)
    # Cache the fitted umap model
    cluster_class.umap_model = umap_model
    return embedding

def umap_metric_learning(cluster_class: Any, x_train: np.ndarray, y_train: np.ndarray, n_neighbors: int = 15) -> np.ndarray:
    """ 
    Fits a supervised UMAP space using categorical structural labels to force metric learning and enhance class separation in the embedding
    """
    _check_umap_availability()
    cluster_class.log(f"Applying supervised UMAP Metric Learning with (n_neighbors={n_neighbors}) on training data with {x_train.shape[0]} samples and {x_train.shape[1]} features.")
    umap_model = UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        metric=cluster_class.umap_metric,
    ).fit(x_train, y=y_train)
    
    # Update the internal tracked properties inside shared context space
    cluster_class.supervised_mapper = umap_model
    cluster_class.trained_embedding = umap_model.embedding_
    cluster_class.log("Supervised UMAP Metric Learning completed and model cached in context.")
    return umap_model.embedding_

def umap_transform(cluster_class: Any, x_test:np.ndarray) -> np.ndarray:
    """  
    Projects raw, out-of-sample data points directly into an existing UMAP manifold
    using whichever pre-fitted mapper is available in the context
    """
    umap_model = getattr(cluster_class, 'supervised_mapper', None) or getattr(cluster_class, 'umap_model', None)
    if umap_model is None:
        raise ValueError("No fitted UMAP model found in context. Please fit a UMAP model before attempting to transform new data.")
    
    embedding = umap_model.transform(x_test)
    cluster_class.log(f"Transformed {x_test.shape[0]} out-of-sample data points into UMAP space.")
    return embedding

def embed_new_structures(cluster_class: Any, x_new: np.ndarray) -> np.ndarray:
    """ 
    Helper wrapping umap_transform for checking/embedding new structures into an existing UMAP space
    """
    embedding = umap_transform(cluster_class, x_new)
    cluster_class.log(f"Embedded {x_new.shape[0]} new structures into existing UMAP space.")
    return embedding

def flag_novel_structures(
    context: Any, 
    x_new: np.ndarray, 
    classifier: Optional[Any] = None, 
    threshold_percentile: float = 25.0
) -> np.ndarray:
    """
    Embeds new structures and flags indices whose classification margins fall below a confidence 
    percentile threshold, signaling highly unique/novel microstates.
    """
    active_classifier = classifier if classifier is not None else getattr(context, "classifier_model", None)
    if active_classifier is None:
        raise ValueError("No valid classification model found in context. Train or provide one explicitly.")

    # Project coordinates into the persistent low-dimensional space
    embedding = embed_new_structures(context, x_new)
    
    # Analyze output boundary profiles
    probabilities = active_classifier.predict_proba(embedding)
    max_probabilities = np.max(probabilities, axis=1)
    
    threshold = np.percentile(max_probabilities, threshold_percentile)
    novel_indices = np.where(max_probabilities <= threshold)[0]
    
    context.log(
        f"Novelty analysis: Flagged {len(novel_indices)}/{len(x_new)} arrays as novel configurations "
        f"(Max Prediction confidence <= {threshold:.4f} at the {threshold_percentile}th percentile)"
    )
    return novel_indices

def _check_umap_availability() -> None:
    """Helper function guarding runtime imports"""
    if UMAP is None:
        raise ImportError("UMAP is not available. Please install it with 'pip install umap-learn' to use UMAP dimensionality reduction.")
    
if __name__ == "__main__":
    from modules.cluster import Clustering
    # Example usage and testing of dimensionality reduction functions
    mock_data = np.random.rand(100, 50)  # 100 samples, 50 features
    mock_cluster = Clustering(feature_matrix=mock_data, energies=None, normalize=True)
    pca_embedding, pca_variance = pca(mock_cluster, n_components=5)
    kernel_embedding, kernel_variance = kernel_based_pca(mock_cluster, n_components=5, kernel='rbf', gamma=0.1)
    tsne_embedding = tsne(mock_cluster, n_components=2, perplexity=30, n_iter=1000, random_state=42)
    if UMAP is not None:
        umap_embedding = umap(mock_cluster, n_components=2, n_neighbors=15, min_dist=0.1, random_state=42)
    print("Dimensionality reduction tests completed successfully.")