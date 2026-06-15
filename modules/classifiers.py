"""
Unified classification engine for molecular configurations
Handles, training, prediction, serialization and uncertainty profiling
"""
import os 
import joblib
import numpy as np
from typing import Optional, Any
from scipy.stats import entropy
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.base import clone
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import classification_report, confusion_matrix




def train_classifier(
        cluster_class: Any,
        model_type: str,
        x_train: np.ndarray,
        y_train: np.ndarray,
        save_path: Optional[str] = None,
        **kwargs
    ) -> Any:
    """ 
    Trains a classification model, saves it to the active Clustering context
    and optionally serializes the trained checkpoint to disk.

    Args:
        cluster_class: The Clustering instance to which the trained model will be attached
        model_type: The type of classifier to train (e.g., 'knn', 'svm', 'rf')
        x_train: The training feature matrix
        y_train: The training labels
        save_path: Optional path to save the trained model checkpoint
        **kwargs: Additional keyword arguments for the classifier constructor
    """
    model_type = model_type.lower()
    cluster_class.log(f"Initializing classifier training for model type: {model_type}")

    if model_type == 'knn':
        # Fall back to containers distance metric choice if not explicitly overridden
        if 'metric' not in kwargs:
            kwargs['metric'] = cluster_class.metric
        model = KNeighborsClassifier(**kwargs)
    elif model_type == 'svm':
        # Enable internal probability calibration 
        if "probability" not in kwargs:
            kwargs['probability'] = True
        model = SVC(**kwargs)
    elif model_type == 'rf':
        model = RandomForestClassifier(**kwargs)
    else:
        raise ValueError(f"Unsupported model type: {model_type}. Supported types are: 'knn', 'svm', 'rf'.")
    
    # Train model on the provided data space
    model.fit(x_train, y_train)
    cluster_class.log(f"Sucsessfully trained {model_type} classifier on {x_train.shape[0]} samples with {x_train.shape[1]} features.")

    # mutate model to store in cluster class context
    cluster_class.classifier_model = model

    if save_path:
        raise NotImplementedError("Model checkpoint saving is not yet implemented. Please implement this functionality to enable model persistence.")

    return model

def predict_labels(cluster_class: Any, x_test: np.ndarray) -> np.ndarray:
    """  
    Predicts precise, discrete cluster assignments using the stored classifier model in the Clustering context.
    """
    _ensure_model_is_loaded(cluster_class)
    predictions = cluster_class.classifier_model.predict(x_test)
    
    cluster_class.log(f"Completed cluster label assignment for {x_test.shape[0]} samples.")
    return predictions

def predict_probabilities(cluster_class: Any, x_test: np.ndarray) -> np.ndarray:
    """  
    Calculates assignment probabilities across clusters and computes Shannon Entropy to flag ambiguous configurations.
    """
    _ensure_model_is_loaded(cluster_class)

    # Generate probabilities
    probabilities = cluster_class.classifier_model.predict_proba(x_test)

    # Store matrix natively into the Clustering container context
    cluster_class.classifier_probabilities = probabilities

    # Informational entropy calculation
    molecule_entropies = entropy(probabilities, axis=1)
    avg_entropy = np.mean(molecule_entropies)

    # Statistical distribution metric of high-confidence vs. ambiguous predictions
    max_probabilities = np.max(probabilities, axis=1)
    mean_max_prob = np.mean(max_probabilities)
    std_max_prob = np.std(max_probabilities)
    cluster_class.log("Model Uncertainty Profiling & Confidence Analysis:")
    cluster_class.log(f"   -> Average Shannon Entropy across predictions: {avg_entropy:.4f}")
    cluster_class.log(f"   -> Mean of maximum predicted probabilities: {mean_max_prob:.4f}")
    cluster_class.log(f"   -> Standard deviation of maximum predicted probabilities: {std_max_prob:.4f}")

    return probabilities



    
def _ensure_model_is_loaded(cluster_class: Any) -> None:
    """ 
    Internal sanity check to ensure that a classifier model is loaded
    """
    if not hasattr(cluster_class, 'classifier_model') or cluster_class.classifier_model is None:
        raise ValueError("No classifier model is currently loaded in the Clustering context. Please train or load a model before making predictions.")

if __name__ == "__main__":
    # Import cluster class for testing
    from modules.cluster import Clustering

    # Testing of the classifier functons
    mock_features = np.random.rand(100, 10)  # 100 samples, 10 features
    mock_labels = np.random.randint(0, 5, size=100)  # 5 clusters
    # Create a mock Clustering instance for testing
    mock_cluster = Clustering(feature_matrix=mock_features, energies=None, normalize=True)
    # Train a KNN classifier
    train_classifier(cluster_class=mock_cluster, model_type='knn', x_train=mock_features, y_train=mock_labels, n_neighbors=3)
    # Predict labels and probabilities
    predicted_labels = predict_labels(cluster_class=mock_cluster, x_test=mock_features)
    predicted_probabilities = predict_probabilities(cluster_class=mock_cluster, x_test=mock_features)
    print("Predicted Labels:", predicted_labels)
    print("Predicted Probabilities:", predicted_probabilities)
