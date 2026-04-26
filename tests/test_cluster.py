import numpy as np
import pytest

# import iris dataset for testing
from sklearn.datasets import load_iris
from modules.cluster import Clustering

@pytest.fixture
def iris_data():
    """ 
    Provides the Iris feature matrix and the labels
    """
    data = load_iris()
    return data.data, data.target

@pytest.fixture
def clustering_instance(iris_data):
    """ 
    Provides a instance of the Clustering class initialized
    """    
    X, _ = iris_data
    return Clustering(feature_matrix=X, normalize=True)

def test_clustering_initialization(clustering_instance):
    """  
    Test that the feature matrix is correctly assigned and normalized in the Clustering class
    """
    clustering = clustering_instance
    # Check that the feature matrix is assigned
    assert hasattr(clustering, 'feature_matrix_raw'), "Clustering instance should have a feature_matrix attribute"
    assert hasattr(clustering, '_feature_matrix_normalized'), "Clustering instance should have a feature_matrix_normalized attribute"
    # Check that the raw feature matrix is the same as the input
    assert np.array_equal(clustering.feature_matrix_raw, clustering._feature_matrix_normalized) == False, "Raw and normalized feature matrices should not be the same"
    # Check that the normalized feature matrix has mean ~0 and std ~1
    mean = np.mean(clustering._feature_matrix_normalized, axis=0)
    std = np.std(clustering._feature_matrix_normalized, axis=0)
    assert np.allclose(mean, 0, atol=1e-6), "Normalized feature matrix should have mean close to 0"
    assert np.allclose(std, 1, atol=1e-6), "Normalized feature matrix should have std close to 1"

def test_clustering_kmeans(clustering_instance):
    """  
    Test the KMeans clustering method on the Iris dataset
    """
    clustering = clustering_instance
    n_clusters = 3
    labels = clustering.kmeans_clustering(n_clusters=n_clusters, random_state=42)
    # Check that the number of unique labels is equal to n_clusters
    unique_labels = np.unique(labels)
    assert len(unique_labels) == n_clusters, f"KMeans should produce {n_clusters} clusters, but got {len(unique_labels)}"
    # Check that the labels are integers
    assert issubclass(labels.dtype.type, np.integer), "KMeans labels should be integers"

def test_clustering_agglomerative(clustering_instance):
    """  
    Test the Agglomerative clustering method on the Iris dataset
    """
    clustering = clustering_instance
    n_clusters = 3 
    labels = clustering.agglomerative_clustering(n_clusters=n_clusters)
    # Check that the number of unique labels is equal to n_clusters
    unique_labels = np.unique(labels)
    assert len(unique_labels) == n_clusters, f"Agglomerative should produce {n_clusters} clusters, but got {len(unique_labels)}"
    # Check that the labels are integers
    assert issubclass(labels.dtype.type, np.integer), "Agglomerative labels should be integers"

def test_clustering_dbscan(clustering_instance):
    """ 
    Test the DBSCAN clustering method on the Iris dataset
    """
    clustering = clustering_instance 
    eps = 0.5
    min_samples = 5
    labels = clustering.dbscan_clustering(eps=eps, min_samples=min_samples)
    # Check that the labels are integers
    assert issubclass(labels.dtype.type, np.integer), "DBSCAN labels should be integers"
    # Check that there are some noise points (label = -1)
    assert np.any(labels == -1), "DBSCAN should identify some noise points with label = -1"

def test_clustering_hdbscan(clustering_instance):
    """  
    Test the HDBSCAN clustering method on the Iris dataset
    """
    clustering = clustering_instance 
    min_cluster_size = 5
    labels = clustering.hdbscan_clustering(min_cluster_size=min_cluster_size)
    # Check that the labels are integers
    assert issubclass(labels.dtype.type, np.integer), "HDBSCAN labels should be integers"
    # Check that there are some noise points (label = -1)
    assert np.any(labels == -1), "HDBSCAN should identify some noise points with label = -1"

def test_pca(clustering_instance):
    """  
    Test the PCA dimensionality reduction method on the Iris dataset
    """
    clustering = clustering_instance 
    n_components = 2
    reduced_data, variance_explained_ratio = clustering.pca(n_components=n_components)
    # Check that the reduced data has the correct shape
    assert reduced_data.shape[1] == n_components, f"PCA should reduce to {n_components} components, but got {reduced_data.shape[1]}"
    # Check that the reduced data is not all zeros
    assert not np.allclose(reduced_data, 0), "PCA reduced data should not be all zeros"



def test_tsne(clustering_instance):
    """ 
    Test the t-SNE dimensionality reduction method on the Iris dataset
    """
    clustering = clustering_instance 
    n_components = 2
    reduced_data = clustering.tsne(n_components=n_components, random_state=42)
    # Check that the reduced data has the correct shape
    assert reduced_data.shape[1] == n_components, f"t-SNE should reduce to {n_components} components, but got {reduced_data.shape[1]}"
    # Check that the reduced data is not all zeros
    assert not np.allclose(reduced_data, 0), "t-SNE reduced data should not be all zeros"

def test_umap(clustering_instance):
    """  
    Test the UMAP dimensionality reduction method on the Iris dataset
    """
    clustering = clustering_instance 
    n_components = 2
    reduced_data = clustering.umap(n_components=n_components, random_state=42)
    # Check that the reduced data has the correct shape
    assert reduced_data.shape[1] == n_components, f"UMAP should reduce to {n_components} components, but got {reduced_data.shape[1]}"
    # Check that the reduced data is not all zeros
    assert not np.allclose(reduced_data, 0), "UMAP reduced data should not be all zeros"

def test_embedding_plotting(clustering_instance):
    """ 
    Test that the plotting for PCA, t-SNE, and UMAP embeddings works without errors
    """
    clustering = clustering_instance
    n_components = 3 
    # Test by plotting for PCA
    reduced_data_pca, _ = clustering.pca(n_components=n_components)
    clustering.plot_embedding(reduced_data_pca, title="PCA Embedding")
    # Test by plotting for t-SNE
    reduced_data_tsne = clustering.tsne(n_components=n_components, random_state=42)
    clustering.plot_embedding(reduced_data_tsne, title="t-SNE Embedding")
    # Test by plotting for UMAP
    reduced_data_umap = clustering.umap(n_components=n_components, random_state=42)
    clustering.plot_embedding(reduced_data_umap, title="UMAP Embedding")

def test_evaluate(clustering_instance):
    """  
    Test the evaluation of multiple metrics
    """
    clustering = clustering_instance 
    n_clusters = 3
    labels_kmeans = clustering.kmeans_clustering(n_clusters=n_clusters, random_state=42)
    labels_agglomerative = clustering.agglomerative_clustering(n_clusters=n_clusters)
    labels_dbscan = clustering.dbscan_clustering(eps=0.5, min_samples=5)
    labels_hdbscan = clustering.hdbscan_clustering(min_cluster_size=5)

    metrics_kmeans = clustering.evaluate(labels_kmeans)
    metrics_agglomerative = clustering.evaluate(labels_agglomerative)
    metrics_dbscan = clustering.evaluate(labels_dbscan)
    metrics_hdbscan = clustering.evaluate(labels_hdbscan)

    # Check that the metrics are dictionaries with expected keys
    expected_keys = {"silhouette_score", "calinski_harabasz_score", "davies_bouldin_score"}
    assert isinstance(metrics_kmeans, dict) and expected_keys.issubset(metrics_kmeans.keys()), "KMeans evaluation should return a dictionary with the expected metric keys"
    assert isinstance(metrics_agglomerative, dict) and expected_keys.issubset(metrics_agglomerative.keys()), "Agglomerative evaluation should return a dictionary with the expected metric keys"
    assert isinstance(metrics_dbscan, dict) and expected_keys.issubset(metrics_dbscan.keys()), "DBSCAN evaluation should return a dictionary with the expected metric keys"
    assert isinstance(metrics_hdbscan, dict) and expected_keys.issubset(metrics_hdbscan.keys()), "HDBSCAN evaluation should return a dictionary with the expected metric keys"

def test_elbow_analysis(clustering_instance):
    """ 
    Test the implementation of the elbow analysis for the clustering method
    """
    clustering = clustering_instance
    k_range = range(1, 10)
    inertia_values = clustering.elbow_analysis(k_range=k_range, random_state=42)
    # Inertia values are dictionary with k as keys and inertia as values
    assert isinstance(inertia_values, dict), "Elbow analysis should return a dictionary of inertia values with k as keys and inertia as values"
    # Check that inertia values are returned for each k in the range
    assert len(inertia_values) == len(k_range), f"Elbow analysis should return inertia values for each k in the range, but got {len(inertia_values)} inertia values for {len(k_range)} k values"
    # Check that inertia values are non-negative
    assert all(inertia >= 0 for inertia in inertia_values.values()), "Inertia values should be non-negative"

def test_silhouette_analysis(clustering_instance):
    """  
    Test the implementation of silhouette analysis for the clustering method
    """
    clustering = clustering_instance
    k_range = range(2, 10)  # Silhouette score is not defined for k=1
    best_k, scores = clustering.silhouette_analysis(k_range=k_range, random_state=42)
    # Check that best_k is integer and scores a dictionary with k as keys and silhouette scores as values
    assert isinstance(best_k, int), "Best k from silhouette analysis should be an integer"
    assert isinstance(scores, dict), "Silhouette analysis should return a dictionary of silhouette scores with k as keys and silhouette scores as values"
    # Check that silhouette scores are returned for each k in the range
    assert len(scores) == len(k_range), f"Silhouette analysis should return silhouette scores for each k in the range, but got {len(scores)} silhouette scores for {len(k_range)} k values"
    # Check that silhouette scores are between -1 and 1
    assert all(-1 <= score <= 1 for score in scores.values()), "Silhouette scores should be between -1 and 1"