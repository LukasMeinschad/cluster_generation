import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional 

def _ensure_dir(filepath: str):
    """  
    Ensures that the directory for the plot file exists
    """
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)

def plot_energy_distribution(
        energies: np.ndarray,
        bins: int = 70,
        phase: str = "All",
        save_path: str = "figures/energy_distribution.png"):
    """  
    Plots the energy distribution for a given set of structures
    """
    _ensure_dir(save_path)
    plt.figure(figsize=(10, 6))
    sns.histplot(energies, bins=bins, kde=True, color="skyblue", alpha=0.7)
    plt.title(f"Energy Distribution - Phase {phase}")
    plt.xlabel("Energy (Hartree)")
    plt.ylabel("Count")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()


