"""
Helper Module for Plotting of the Clustering Class
"""
import os
from typing import Any, Optional, Dict, List, Tuple
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# =========================================================
# Shared scientific-plot helpers
# =========================================================

def _resolve_labels(cluster_class: Any, labels: Optional[np.ndarray], n_points: int) -> Optional[np.ndarray]:
    """
    Resolves the label array used for color-coding, validating its length against
    the number of plotted points. Falls back to the cluster_class's stored labels
    if none are explicitly provided, and to None if neither is usable.
    """
    if labels is not None:
        if len(labels) != n_points:
            raise ValueError(f"Length of labels ({len(labels)}) does not match number of points ({n_points}).")
        return labels

    internal_labels = getattr(cluster_class, "labels", None)
    if internal_labels is not None and len(internal_labels) == n_points:
        return internal_labels
    return None


def _resolve_palette(n_labels: int) -> List:
    """
    Resolves a categorical color palette sized to the number of cluster labels,
    switching to a higher-capacity HSV palette once the qualitative palette runs
    out of visually distinct colors (>10 clusters).
    """
    if n_labels > 10:
        return sns.color_palette("gist_stern", n_colors=n_labels)
    return sns.color_palette("tab10", n_labels)


def _resolve_save_path(save_path: Optional[str], title: str) -> str:
    """
    Resolves the output path for a figure, defaulting to a slugified title under figures/.
    """
    if save_path:
        return save_path
    os.makedirs("figures", exist_ok=True)
    return f"figures/{title.replace(' ', '_').lower()}.png"


def _style_axes(ax: plt.Axes) -> None:
    """
    Applies the shared scientific/paper styling: light dashed gridlines kept
    below the data, and despined top/right axes.
    """
    ax.grid(True, linestyle="--", alpha=0.5, color="#e0e0e0", zorder=0)
    ax.set_axisbelow(True)  # Keep grid lines below scatter points
    sns.despine(ax=ax, top=True, right=True)


def _style_title_and_labels(ax: plt.Axes, title: str) -> None:
    """Applies the shared scientific typography for the title and axis labels."""
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel("Component 1", fontsize=12, labelpad=8)
    ax.set_ylabel("Component 2", fontsize=12, labelpad=8)


def _place_legend(ax: plt.Axes, title: str, n_labels: int, handles: Optional[List] = None) -> None:
    """Places a legend outside the axes, wrapping into two columns for many clusters."""
    ncol = 2 if n_labels > 12 else 1
    kwargs = dict(
        title=title,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        ncol=ncol,
        frameon=True,
        facecolor="white",
        edgecolor="black",
    )
    if handles is not None:
        ax.legend(handles=handles, **kwargs)
    else:
        ax.legend(**kwargs)


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
    # Use Seaborn theme
    sns.set_theme(style="ticks", context="paper")

    x, y = embedding[:, 0], embedding[:, 1]
    fig, ax = plt.subplots(figsize=(9, 6))

    lbl = _resolve_labels(cluster_class, labels, len(x))

    # Render plot based on resolved labels
    if lbl is not None:
        unique_labels = len(np.unique(lbl))
        palette = _resolve_palette(unique_labels)

        if unique_labels > 10:
            base_markers = ["o", "s", "^", "D", "H"]
            markers = [base_markers[i % len(base_markers)] for i in range(unique_labels)]
        else:
            markers = ["o"] * unique_labels  # Use circle markers for fewer clusters

        sns.scatterplot(x=x, y=y, hue=lbl, style=lbl, markers=markers,
                        palette=palette, legend='full', alpha=0.8, s=25, edgecolor="white", linewidth=0.3, ax=ax)

        _place_legend(ax, title="Clusters", n_labels=unique_labels)

    else:
        sns.scatterplot(x=x, y=y, color='blue', alpha=0.6, s=20, edgecolor="none", ax=ax)

    # Nicer Typography and Styling
    _style_title_and_labels(ax, title)
    _style_axes(ax)
    plt.tight_layout()

    out_path = _resolve_save_path(save_path, title)
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
        cluster_class: The central Clustering context container.
        embedding: Array of shape (n_samples, 2) matching the coordinate space.
        probabilities: Array of shape (n_samples, n_clusters) containing model probabilities.
        labels: Optional explicit label array override. If None, reads cluster_class.labels.
        title: Title of the generated plot.
        save_path: Optional exact export path. Defaults to 'figures/<title>.png'.
    """
    # Use the same scientific theme as plot_embedding
    sns.set_theme(style="ticks", context="paper")

    x, y = embedding[:, 0], embedding[:, 1]
    lbl = _resolve_labels(cluster_class, labels, len(x))

    # Extract the highest probability score per structure as our certainty metric
    confidences = np.max(probabilities, axis=1)
    sizes = 10 + 90 * confidences  # Map confidence smoothly from scale 10 to 100

    fig, ax = plt.subplots(figsize=(9, 6))

    if lbl is not None:
        unique_labels = sorted(set(lbl))
        palette = _resolve_palette(len(unique_labels))
        color_map = {l: palette[i] for i, l in enumerate(unique_labels)}

        # Build per-point RGBA sequences where assignment certainty shapes the alpha channel
        rgba = np.array([(*color_map[l], float(conf)) for l, conf in zip(lbl, confidences)])
        ax.scatter(x, y, c=rgba, s=sizes, edgecolor="white", linewidth=0.3)

        handles = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[l], markersize=8, label=l)
            for l in unique_labels
        ]
        _place_legend(ax, title="Clusters", n_labels=len(unique_labels), handles=handles)
    else:
        # Fallback to standard color matrix modulated by sample classification weights
        rgba = np.array([(0.12, 0.47, 0.71, float(c)) for c in confidences])
        ax.scatter(x, y, c=rgba, s=sizes, edgecolor="white", linewidth=0.3)

    _style_title_and_labels(ax, f"{title}\n(marker size & alpha reflect prediction confidence)")
    _style_axes(ax)
    plt.tight_layout()

    out_path = _resolve_save_path(save_path, title)
    fig.savefig(out_path, dpi=300)
    cluster_class.log(f"Classification confidence map saved to: {out_path}")
    plt.close(fig)
