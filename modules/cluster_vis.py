""" 
Helper Module for Plotting of the Clustering Class
"""
import os
from typing import Any, Optional, Dict, List, Tuple
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_embedding(
        cluster_class: Any,
        embedding: np.ndarray,
        title: str = "Embedding",
        labels: Optional[np.ndarray] = None,
        save_path: Optional[str] = None,
    ) -> None:
    """ 
    Plots a 2D dimensionality reduction space with cluster markeers
    If the incoming array has > 2 dimensions, the first two columns will be used
    """
    x, y = embedding[:, 0], embedding[:, 1]
    plt.figure(figsize=(8, 6))
    
    # Resolbe labels
    lbl = labels if labels is not None else getattr(cluster_class, 'labels', None)

    if lbl is not None:
        sns.scatterplot(x=x, y=y, hue=lbl, palette='tab10', legend='full', s=50)
        plt.legend(title="Cluster")
    else:
        sns.scatterplot(x=x, y=y, color='blue', s=50)

    plt.title(title)
    plt.xlabel("Component 1")
    plt.ylabel("Component 2")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_path:
        out_path = save_path
    else:
        os.makedirs("figures", exist_ok=True)
        out_path = f"figures/{title.replace(' ', '_').lower()}.png"
    plt.savefig(out_path, dpi=300)
    cluster_class.log(f"Saved embedding plot to: {out_path}")
    plt.close()

def plot_probabilities(
        cluster_class: Any,
        embedding: np.ndarray,
        probabilities: np.ndarray,
        labels: Optional[np.ndarray] = None,
        title: str = "Classifier Assignment Probabilities",
        save_path: Optional[str] = None
    ) -> None:
    """
    Plots a 2D embedding space where point marker scale and alpha transparency 
    reflect the assignment confidence (maximum probability) of each structure.

    Args:
        context: The central Clustering context container.
        embedding: Array of shape (n_samples, 2) matching the coordinate space.
        probabilities: Array of shape (n_samples, n_clusters) containing model probabilities.
        labels: Optional explicit label array override. If None, reads context.labels.
        title: Title of the generated plot.
        save_path: Optional exact export path. Defaults to 'figures/<title>.png'.
    """
    x, y = embedding[:, 0], embedding[:, 1]
    
    # Extract the highest probability score per structure as our certainty metric
    confidences = np.max(probabilities, axis=1)
    sizes = 10 + 90 * confidences  # Map confidence smoothly from scale 10 to 100

    fig, ax = plt.subplots(figsize=(10, 8))
    lbl = labels if labels is not None else getattr(cluster_class, "labels", None)

    if lbl is not None:
        unique_labels = sorted(set(lbl))
        palette = sns.color_palette('tab10', len(unique_labels))
        color_map = {l: palette[i] for i, l in enumerate(unique_labels)}

        # Build per-point RGBA sequences where assignment certainty shapes the alpha channel
        rgba = np.array([(*color_map[l], float(conf)) for l, conf in zip(lbl, confidences)])
        ax.scatter(x, y, c=rgba, s=sizes)

        handles = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[l], markersize=8, label=l)
            for l in unique_labels
        ]
        ax.legend(handles=handles, title="Cluster Assignment")
    else:
        # Fallback to standard color matrix modulated by sample classification weights
        rgba = np.array([(0.12, 0.47, 0.71, float(c)) for c in confidences])
        ax.scatter(x, y, c=rgba, s=sizes)

    ax.set_title(f"{title}\n(Size & Alpha reflect prediction confidence)")
    ax.set_xlabel("Component 1")
    ax.set_ylabel("Component 2")
    ax.grid(True, alpha=0.3)

    if save_path:
        out_path = save_path
    else:
        os.makedirs("figures", exist_ok=True)
        out_path = f"figures/{title.replace(' ', '_').lower()}.png"

    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    cluster_class.log(f"Classification confidence map saved to: {out_path}")
    plt.close(fig)


