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

# Validation and evaluation utilities
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
    cluster_class.substep(f"Trained {model_type.upper()} classifier on {x_train.shape[0]} samples with {x_train.shape[1]} features")

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
    cluster_class.parameters("Model Uncertainty Profiling & Confidence Analysis", [
        ("Average Shannon entropy", f"{avg_entropy:.4f}"),
        ("Mean max predicted probability", f"{mean_max_prob:.4f}"),
        ("Std max predicted probability", f"{std_max_prob:.4f}"),
    ])

    return probabilities


def evaluate_classifier(
        cluster_class: Any,
        x: np.ndarray,
        y: np.ndarray,
        n_splits: int = 5,
        save_dir: Optional[str] = None
    ) -> None:
    """  
    Evaluates the performance of the stored classifier model using stratified k-fold cross-validation
    
    Args:        
        cluster_class: The Clustering instance containing the classifier model to evaluate
        x: The feature matrix for evaluation
        y: The true labels for evaluation
        n_splits: The number of folds for cross-validation
        save_dir: Optional directory to save evaluation reports and confusion matrices
    """
    _ensure_model_is_loaded(cluster_class)
    prod_model = cluster_class.classifier_model

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    # Arrays to store out-of-fold predictions and probabilities
    oof_preds = np.zeros_like(y)
    unique_classes = np.unique(y)
    oof_probs = np.zeros((len(y), len(unique_classes)))

    cluster_class.substep(f"Starting {n_splits}-fold cross-validation for classifier evaluation...")

    for fold, (train_idx, val_idx) in enumerate(skf.split(x,y)):
        x_tr, x_val = x[train_idx], x[val_idx]
        y_tr, y_val = y[train_idx], y[val_idx]

        # clone() copies the hyperparameters (KNN Metrics or SVM kernels)
        # without copying fitted weights
        fold_model = clone(prod_model)
        fold_model.fit(x_tr, y_tr)

        oof_preds[val_idx] = fold_model.predict(x_val)
        if hasattr(fold_model, "predict_proba"):
            oof_probs[val_idx] = fold_model.predict_proba(x_val)
        
    # Compile scikit-learn performane dictionaries and printout text
    report_dict = classification_report(y, oof_preds, output_dict=True)
    report_text = classification_report(y, oof_preds)

    cluster_class.log("Classifier evaluation report (full per-class table at DEBUG level):")
    cluster_class.log(report_text, level="debug")

    # Profile boundary confidence margins if probabilities are available
    if hasattr(prod_model, "predict_proba"):
        max_probs = np.max(oof_probs, axis=1)
        cluster_class.parameters("Classifier Evaluation Summary", [
            ("Avg out-of-fold confidence", f"{np.mean(max_probs):.4f}"),
            ("Below 60% confidence", f"{np.sum(max_probs < 0.6)} / {len(max_probs)} ({np.mean(max_probs < 0.6) * 100:.2f}%)"),
        ])

    # Generate and Save Confusion matrix plot
    cm = confusion_matrix(y, oof_preds)
    plt.figure(figsize=(8,6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=unique_classes, yticklabels=unique_classes)
    plt.title("Confusion Matrix")
    plt.xlabel("Predicted Labels")
    plt.ylabel("True Labels")

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        cm_path = os.path.join(save_dir, "confusion_matrix.png")
        plt.savefig(cm_path)
        cluster_class.log(f"Confusion matrix saved to: {cm_path}")
    plt.close()

    return report_dict



    
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
