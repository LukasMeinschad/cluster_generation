import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional
from scipy.spatial.distance import pdist, squareform

""" 
A various collection of visualization functions for analysis of molecules and BHMC results
"""

from molecule_class import Molecule

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

def plot_contour_vdw_radii(molecule: Molecule, save_path: Optional[str] = None):
    """ 
    Makes three contour plots showing the following:

    The VDW radius times 0.1, ... 1 as contours around each atom in the molecule.
    This is then done for the x-y, x-z and y-z planes in a three column subplot
    """ 
    coords = molecule.coordinates
    vdw_radii = molecule.vdw_radii

    planes = [
        ("X","Y", (0,1)),
        ("X","Z", (0,2)),
        ("Y","Z", (1,2))
    ]
    levels = np.linspace(0, 1.0, 10)  # Contour levels from 0 to 1 times the VDW radius
    padding = vdw_radii.max() + 1.5
    res = 150 

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    def get_plane_data(idx_a, idx_b):
        # Helper function to calculate grid and sitance field for a specified projecton
        a_min, a_max = coords[:, idx_a].min() - padding, coords[:, idx_a].max() + padding
        b_min, b_max = coords[:, idx_b].min() - padding, coords[:, idx_b].max() + padding
        
        grid_a = np.linspace(a_min, a_max, res)
        grid_b = np.linspace(b_min, b_max, res)
        AA, BB = np.meshgrid(grid_a, grid_b)

        # Initialize field with infinty
        Z = np.full(AA.shape, np.inf)

        # Calculate scaled distance field: dist / VDW radius
        for i in range(len(coords)):
            dist = np.sqrt((AA - coords[i, idx_a])**2 + (BB - coords[i, idx_b])**2)
            scaled_dist = dist / vdw_radii[i]
            Z = np.minimum(Z, scaled_dist)  # Take the minimum across all atoms
        return AA, BB, Z
    
    # Plotting loop
    for i, (label_a, label_b, indices) in enumerate(planes):
        ax = axes[i]
        AA, BB, Z = get_plane_data(*indices)

        cp = ax.contour(AA, BB, Z, levels=levels, cmap='magma', linewidths=0.8)

        # Plot nuclei positions
        ax.scatter(coords[:, indices[0]], coords[:, indices[1]], color='black', s=15, label='Nuclei')
        ax.set_title(f"{label_a}-{label_b} Plane")
        ax.set_xlabel(f"{label_a} (Å)")
        ax.set_ylabel(f"{label_b} (Å)")
        ax.legend()
        ax.set_aspect('equal', adjustable='box')
    
    # Add colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
    fig.colorbar(cp, cax=cbar_ax, label='Scaled Distance (dist / VDW radius)')
    plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout to make room for colorbar
    if save_path:
        _ensure_dir(save_path)
        plt.savefig(save_path, dpi=300)
    plt.close()

def plot_pairwise_distance_heatmap(molecule: Molecule, save_path: Optional[str] = None):
    """ 
    Plots a heatmap of the pairwise distances between atoms in the molecule
    """
    coords = molecule.coordinates
    dists = squareform(pdist(coords))
    vdw_radii = molecule.vdw_radii
    # Print vdw radii * 0.5, *0.8 * 1.0 for debugging
    print("VDW Radii:", vdw_radii)
    print("VDW Radii * 0.5:", vdw_radii * 0.5)
    print("VDW Radii * 0.8:", vdw_radii * 0.8)
    print("VDW Radii * 1.0:", vdw_radii * 1.0)


    plt.figure(figsize=(8, 6))
    sns.heatmap(dists, cmap='viridis', xticklabels=False, yticklabels=False)
    # Add text annotations for the distances
    for i in range(dists.shape[0]):
        for j in range(dists.shape[1]):
            if i < j:  # Only annotate upper triangle to avoid clutter
                plt.text(j + 0.5, i + 0.5, f"{dists[i, j]:.2f}", ha='center', va='center', color='white', fontsize=6)


    plt.title("Pairwise Distance Heatmap")
    plt.xlabel("Atom Index")
    plt.ylabel("Atom Index")
    plt.tight_layout()
    
    if save_path:
        _ensure_dir(save_path)
        plt.savefig(save_path, dpi=300)
    plt.close()

    

