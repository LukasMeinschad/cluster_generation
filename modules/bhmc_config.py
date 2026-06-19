from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict, Any

@dataclass
class BHMCConfig:
    """Configuration for Basin Hopping Monte Carlo sampling."""

    # Energy backend
    backend: str = "psi4"
    qm_method: str = "hf"
    qm_basis: str = "sto-3g"
    xtb_method: str = "GFN2-xTB"
    gpaw_mode: str = "lcao"
    gpaw_basis: str = "dzp"
    gpaw_xc: str = "B3LYP"

    # Sampling
    temperature: float = 400.0      # Metropolis temperature in Kelvin
    verbose: bool = False
    adaptive_operators: bool = True  # Pass adaptive=True to operators that support it

    # Operator pool: list of (name, weight) for both local and non-local operators.
    # Weights are relative; they are normalised inside the worker.
    operators: List[Tuple[str, float]] = field(default_factory=lambda: [
        # --- non-local (basin-hop moves) ---
        ("twist",                    0.6),
        ("large_displacement",       0.5),
        ("mirror",                   0.4),
        ("roto_reflection",          0.3),
        ("exchange",                 0.3),
        ("random_so3",               0.2),
        ("principal_axis_rotation",  0.3),
        ("com_com_approach",         0.3),
        ("com_com_separation",       0.1),
        # --- local (basin-exploration moves) ---
        ("random_displacement",      3.0),
        ("random_rotation",          2.5),
        ("local_twist",              2.0),
        ("correlated_displacement",  2.0),
        ("small_principal_axis_rotation", 2.0),
    ])

    # Adaptive box control
    adaptive_box: bool = True
    box_update_interval: int = 10
    box_target_acceptance: float = 0.6
    box_acceptance_window: float = 0.05
    box_growth_kp: float = 0.6
    box_growth_max: float = 1.15
    box_max_scale: float = 2.0
    box_stable_windows: int = 3

    # Clustering and Classifier Retraining
    clustering_method: str = "kmeans"
    classifier_backend: str = "knn"
    classifier_kwargs: Optional[Dict[str, Any]] = None

    # Fraction of a novel-structure pool to use as its re-clustering target,
    # so a handful of novel points doesn't get split into one cluster each.
    novel_cluster_ratio: float = 0.1
