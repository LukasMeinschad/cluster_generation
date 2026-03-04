"""  
Module for various clustering methods and analysis of the sampled molecular configurations
"""
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional, Dict, Any, Union
from dataclasses import dataclass, field
from collections import defaultdict
import seaborn as sns
from sklearn.cluster import AgglomerativeClustering as SklearnAgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import cdist, pdist
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from multiprocessing import Pool, cpu_count
import time


# UMAP Import 
try:
    from umap import UMAP
except ImportError:
    UMAP = None
    print("UMAP is not installed. Install with `pip install umap-learn` to use UMAP dimensionality reduction.")




from molecule_class import Molecule
from transformations import GeometryOps
from box import SimulationBox
from logger import Logger

def _compute_hbond_single(args: Tuple) -> Tuple[int,float]:
    """  
    Worker function to compute the H-bond configs, top-level for pickling

    Args:
        args: Tuple of (atom_labels, coordinates, angle_threshold)
            - atom_labels: List[str]
            - coordinates: np.ndarray of shape (n_atoms, 3)
            - angle_threshold: float in degrees
            - max_distance: float in Angstroms the donor-acceptor distance threshold
    """
    atom_labels, coordinates, angle_threshold, max_distance  = args 

    # Reconstruct molecule
    mol = Molecule(name="worker_mol")
    mol.add_atoms_batch(atom_labels, coordinates)
    mol.compute_bonds()
    valid_configs = mol.get_valid_hbond_configurations(angle_threshold=angle_threshold, max_distance=max_distance)
    
    # Distance based filtering
    valid_configs = [c for c in valid_configs if c.donor_acceptor_distance is not None and c.donor_acceptor_distance <= max_distance]
    
    n_hbonds = len(valid_configs)
    avg_angle = float(np.mean([c.angle for c in valid_configs])) if valid_configs else 0.0
    avg_da_distance = float(np.mean([c.donor_acceptor_distance for c in valid_configs])) if valid_configs else 0.0
    return n_hbonds, avg_angle, avg_da_distance 



@dataclass
class StructureData:
    """  
    Container for Molecular structure and associated properties
    """
    molecule: Molecule
    energy: float
    phase: str = "unknown"
    step: int = -1
    worker_id: int = -1
    accepted: bool = True
    dipole_vector: Optional[List[float]] = None 
    dipole_magnitude: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)

class BHMCAnalyzer:
    """   
    Analyzer for the Basin Hopping Monte Carlo results
    """
    def __init__(self, name: str = "BHMC_Analysis", submolecule_indices: Optional[List[int]] = None, logger: Optional[Logger] = None):
        self.name = name 
        self.structures: List[StructureData] = []
        self.phases: Dict[str, List[StructureData]] = defaultdict(list)
        self.submolecule_indices = submolecule_indices
        self.logger = logger
        self.labels = None # Cluster labels

        # Storage of feature matrix
        self._feature_matrix_raw: Optional[np.ndarray] = None 
        self._feature_matrix_normalized: Optional[np.ndarray] = None 




    def _log(self, msg: str, level: str = "info"):
        """ 
        Log a message if the logger is available
        """
        if self.logger:
            getattr(self.logger, level)(msg)




    def add_structure(self,
                      molecule: Molecule,
                      energy: float, 
                      phase: str = "unknown",
                      dipole_vector: Optional[List[float]] = None,
                      dipole_magnitude: float = 0.0,
                      **kwargs):
        """   
        Add a structure to the analysis

        Args:
            molecule: Molecule object representing the structure
            energy: Energy of the structure
            phase: Phase of the BHMC process (e.g., "phase_a", "phase_b", "phase_c")
            dipole_vector: Optional list of 3 floats representing the dipole vector
            dipole_magnitude: Optional float representing the magnitude of the dipole moment
            **kwargs: Additional properties (e.g., step, worker_id, accepted)
        """
        structure_data = StructureData(
            molecule=molecule,
            energy=energy,
            phase=phase,
            dipole_vector = dipole_vector,
            dipole_magnitude=dipole_magnitude,
            metadata=kwargs
        )
        self.structures.append(structure_data)
        self.phases[phase].append(structure_data)

    def add_structures_batch(self,
                             structures: List[Tuple[Molecule, float]],
                             phase: str = "unknown"):
        """   
        Add multiple structures at once

        Args:
            structures: List of tuples containing (Molecule, energy) or (Molecule, energy, dipole_vec, dipole_magnitude)
            phase: Phase of the BHMC process for all structures
        """
        for item in structures:
            if len(item) == 4:
                mol, energy, dipole_vec, dipole_mag = item 
                self.add_structure(molecule=mol,
                                   energy=energy,
                                   phase=phase,
                                   dipole_vector=dipole_vec,
                                   dipole_magnitude=dipole_mag) 
                
            elif len(item) == 2:
                mol, energy = item
                self.add_structure(molecule=mol,
                                      energy=energy,
                                      phase=phase)
            else:
                raise ValueError(f"Invalid structure tuple length: {len(item)}. Expected 2 or 4.")
        self._log(f"Added batch of {len(structures)} structures to phase: {phase}. Total structures in this phase: {len(self.phases[phase])}")




    # ==================== Energy Statistics =====================
    def get_energy_statistics(self, phase: Optional[str] = None) -> Dict[str, float]:
        """   
        Compute energy statistics.
        
        Args:
            phase: Specified phase to compute statistics for. If None, computes for all structures.

        Returns:
            Dictionary with Statistics
        """
        if phase:
            energies = [s.energy for s in self.phases[phase]]
        else:
            energies = [s.energy for s in self.structures]

        if not energies:
            self._log(f"No energy statistics!", level="warning")
            return {}
        
        energies = np.array(energies)
        stats = {
                "n_structures": len(energies),
                "min": np.min(energies),
                "max": np.max(energies),
                "mean": np.mean(energies),
                "median": np.median(energies),
                "std": np.std(energies),
                "q25": np.percentile(energies, 25),
                "q75": np.percentile(energies, 75),
                "iqr": np.percentile(energies, 75) - np.percentile(energies, 25)
        }
        self._log(f"Energy Statistics (phase: {phase if phase else 'all'}):")
        self._log(f"  N Structures: {stats['n_structures']}")
        self._log(f"  Min: {stats['min']:.6f}")
        self._log(f"  Max: {stats['max']:.6f}")
        self._log(f"  Mean: {stats['mean']:.6f}")
        self._log(f"  Median: {stats['median']:.6f}")
        self._log(f"  Std: {stats['std']:.6f}")
        self._log(f"  Q25: {stats['q25']:.6f}")
        self._log(f"  Q75: {stats['q75']:.6f}")
        self._log(f"  IQR: {stats['iqr']:.6f}")
        return stats


    def compute_rmsd_matrix(self, phase: Optional[str] = None) -> np.ndarray:
        """   
        Compute pairwise RMSD Matrix for the structures in the specified phase
        
        Returns:
            RMSD Matrix (n_structures x n_structures)
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        n = len(structures)
        self._log(f"Computing {n}x{n} RMSD matrix for phase: {phase if phase else 'all'}")
        coords_list = [s.molecule.coordinates for s in structures]
        rmsd_matrix = np.zeros((n, n))
        for i in range(n):
            ci = coords_list[i]
            for j in range(i+1, n):
                rmsd = self._calculate_rmsd(ci, coords_list[j])
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd
        self._log(f"RMSD matrix computed: Range: [{rmsd_matrix[rmsd_matrix > 0].min():.4f}, {rmsd_matrix.max():.4f} ]")
        return rmsd_matrix


    def plot_rmsd_heatmap(self, phase: Optional[str] = None):
        """  
        Helper function to plot the pairwise rmsd matrix as a heatmap
        """
        rmsd_matrix = self.compute_rmsd_matrix(phase=phase)
        plt.figure(figsize=(10, 8))
        sns.heatmap(rmsd_matrix, cmap='viridis')
        plt.title(f'RMSD Matrix Heatmap - Phase: {phase if phase else "All"}')
        plt.xlabel('Structure Index')
        plt.ylabel('Structure Index')
        plt.savefig(f"figures/rmsd_heatmap_{phase if phase else 'all'}.png")
        plt.close()




    # TODO: Restructure this function  
    def intermolecular_distance(self, phase: Optional[str] = None) -> np.ndarray:
        """   
        Computes the intermolecular distance for each structure in the specified phase
        For this we use the submolecule indices and compute the distance between the respective
        geometric means
        TODO: Maybe add com distance as well
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures

        distances = []

        for structure in structures:
            coords = structure.molecule.coordinates
            submolecule_coords = [coords[indices] for indices in self.submolecule_indices]
            submolecule_means = [np.mean(sub_coords, axis=0) for sub_coords in submolecule_coords]
            if len(submolecule_means) == 2:
                distance = np.linalg.norm(submolecule_means[0] - submolecule_means[1])
                distances.append(distance)
            else:
                distances.append(0.0)  # TODO this has to be a matrix with 3 submolecules?
        
        return np.array(distances)
        
    def plot_int_d_vs_e(self, phase: Optional[str] = None):
        """   
        Plots the intermolecular distance vs energy for the specified phase
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures

        energies = [s.energy for s in structures]
        distances = self.intermolecular_distance(phase=phase)

        plt.figure(figsize=(10, 6))
        plt.scatter(distances, energies, color='purple', alpha=0.7)
        plt.title(f'Intermolecular Distance vs Energy - Phase: {phase if phase else "All"}')
        plt.xlabel('Intermolecular Distance')
        plt.ylabel('Energy')
        plt.grid(True)
        plt.savefig(f"figures/int_d_vs_energy_{phase if phase else 'all'}.png")
        plt.close()

    @staticmethod 
    def _calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
        """   
        Calculate RMSD between two sets of coordinates

        Employs the Kabsch Algorithm to find the Optimal Rotation and then Computes the RMSD

        Steps:
        1. Center the Coordinates at the Centroid
        2. Compute the Matrix H = P^T * Q where P and Q are the centered coordinates
        3. Compute the SVD of H = U * S * V^T 
        4. See if Orthogonal Matrix have Reflections d = det(U * V^T) = det(U) * det(V^T)
        5. Calculate the Rotationam Matrix

        R = U (1 1 d) V^T
        """
        # Center the coordinates
        coords1_centered = coords1 - np.mean(coords1, axis=0)
        coords2_centered = coords2 - np.mean(coords2, axis=0)
        # Covariance Matrix
        H = coords1_centered.T @ coords2_centered
        U, S, Vt = np.linalg.svd(H)
        d = np.linalg.det(U) * np.linalg.det(Vt)
        R = U @ np.diag([1, 1, d]) @ Vt
        coords1_rotated = coords1_centered @ R
        rmsd = np.sqrt(np.mean(np.sum((coords1_rotated - coords2_centered)**2, axis=1)))
        return rmsd

    def plot_rmsd_matrix(self, phase: Optional[str] = None):
        """   
        Plot the RMSD matrix as a heatmap
        """
        rmsd_matrix = self.compute_rmsd_matrix(phase=phase)
        plt.figure(figsize=(10, 8))
        sns.heatmap(rmsd_matrix, cmap='viridis')
        plt.title(f'RMSD Matrix Heatmap - Phase: {phase if phase else "All"}')
        plt.xlabel('Structure Index')
        plt.ylabel('Structure Index')
        plt.savefig(f"figures/rmsd_matrix_{phase if phase else 'all'}.png")
        plt.close()

    def plot_energy_vs_rmsd(self, reference_structure: Optional[StructureData] = None, phase: Optional[str] = None):
        """   
        Plots Energy vs RMSD to a reference structure for the specified phase

        Args:
            reference_structure: Optional reference structure to compute RMSD against. If None, it uses to lowest energy structure as reference.
            phase: Specified phase to plot for. If None, plots for all structures.
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        
        if not structures:
            print(f"No structures found for phase: {phase}")
            return
        
        if reference_structure is None:
            reference_structure = self.get_lowest_energy_structure(phase=phase)

        energies = [s.energy for s in structures]
        rmsd_values = [self._calculate_rmsd(s.molecule.coordinates, reference_structure.molecule.coordinates) for s in structures]

        plt.figure(figsize=(10, 6))
        plt.scatter(rmsd_values, energies, color='blue', alpha=0.7)
        plt.title(f'Energy vs RMSD to Reference - Phase: {phase if phase else "All"}')
        plt.xlabel('RMSD to Reference Structure')
        plt.ylabel('Energy')
        plt.grid(True)
        plt.savefig(f"figures/energy_vs_rmsd_{phase if phase else 'all'}.png")
        plt.close()

    def plot_energy_distribution(self, phase: Optional[str] = None, bins: int = 50):
        """
        Plots the energy distribution for the specified phase
        Args:
            phase: Specified phase to plot for. If None, plots for all structures.
            bins: Number of bins for the histogram
        """
        if phase:
            energies = [s.energy for s in self.phases[phase]]
        else:
            energies = [s.energy for s in self.structures]

        plt.figure(figsize=(10, 6))
        plt.hist(energies, bins=bins, color='green', alpha=0.7)
        plt.title(f'Energy Distribution - Phase: {phase if phase else "All"}')
        plt.xlabel('Energy')
        plt.ylabel('Frequency')
        plt.grid(True)
        plt.savefig(f"figures/energy_distribution_{phase if phase else 'all'}.png")
        plt.close()

    def plot_rg_vs_energy(self, phase: Optional[str] = None):
        """   
        Plots the Radius of Gyration vs Energy for the specified phase 

        Args:
            phase: Specified phase to plot for. If None, plots for all structures.
        """ 
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            print(f"No structures found for phase: {phase}")
            return
        energies = [s.energy for s in structures]
        rg_values = [self.radius_of_gyration(s.molecule.coordinates, s.molecule.masses) for s in structures]
        plt.figure(figsize=(10, 6))
        plt.scatter(rg_values, energies, color='red', alpha=0.7)
        plt.title(f'Radius of Gyration vs Energy - Phase: {phase if phase else "All"}')
        plt.xlabel('Radius of Gyration')
        plt.ylabel('Energy')
        plt.grid(True)
        plt.savefig(f"figures/rg_vs_energy_{phase if phase else 'all'}.png")
        plt.close()

    

    def get_lowest_energy_structure(self,phase: Optional[str] = None) -> Optional[StructureData]:
        """   
        Returns the structure with the lowest energy in the specified phase
        
        Args:
            phase: Specified phase to search in. If None, searches in all structures.
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        
        if not structures:
            return None
        
        lowest_energy_structure = min(structures, key=lambda s: s.energy)
        return lowest_energy_structure

    
    
    

    @staticmethod
    def compute_interatomic_distance(coords: np.ndarray) -> np.ndarray:
        """   
        Computes the interatomic distance matrix for a given set of coordinates

        Returns:
            Interatomic distance matrix (n_atoms x n_atoms)
        """
        distance_matrix = cdist(coords, coords)
        # Flatten the upper triangle of the distance matrix to get a descriptor vector
        return distance_matrix[np.triu_indices_from(distance_matrix, k=1)]
    
    def compute_interatomic_distance_matrix(self, phase: Optional[str] = None) -> List[np.ndarray]:
        """    
        Compute the interatomic distance matrix for all the structures
        
        Each row is one structure's descriptor vecotr
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures

        descriptors = []
        for s in structures:
            distance_matrix = self.compute_interatomic_distance(s.molecule.coordinates)
            descriptors.append(distance_matrix)

        return descriptors



    @staticmethod 
    def determine_rotational_constants(coords: np.ndarray, masses: np.ndarray) -> Tuple[float, float, float]:
        """   
        Determines the rotational constants (A,B,C) for a given set of coordinates and masses

        Algorithm:
        1. Center coordinates at the center of mass
        2. Compute the inertia tensor
        3. Diagonalize the inertia tensor to get the principal moments of inertia (I_A, I_B, I_C)
        4. Convert to Rotational Constants using the formula:
           A = h / (8 * π^2 * I_A)
           B = h / (8 * π^2 * I_B)
           C = h / (8 * π^2 * I_C)
        """
        h = 6.62607015e-34  # Planck's constant in J*s
        com = GeometryOps.center_of_mass(coords, masses)
        coords_centered = coords - com
        I = GeometryOps.inertia_tensor(coords_centered, masses)
        I_A, I_B, I_C = np.linalg.eigvalsh(I)  # Get principal moments of inertia
        A = h / (8 * np.pi**2 * I_A) if I_A > 0 else 0.0
        B = h / (8 * np.pi**2 * I_B) if I_B > 0 else 0.0
        C = h / (8 * np.pi**2 * I_C) if I_C > 0 else 0.0
        return A, B, C
    


    @staticmethod
    def radius_of_gyration(coords: np.ndarray, masses: np.ndarray) -> float:
        """   
        Computes the radius of gyration (R_g) for a given set of coordinates and masses

        R_g = sqrt( (1/N) * sum((r_i - r_cm)^2) )

        where r_i is the position vector of atom i, r_cm is the center of mass position vector, and N is the number of atoms.
        """
        com = GeometryOps.center_of_mass(coords, masses)
        N = coords.shape[0]
        diff = coords - com
        rg = np.sqrt(np.sum(masses * np.sum(diff**2, axis=1)) / np.sum(masses))
        return rg

    # ====================== Data Cleanup ========================

    """ 
    We should first of perform a rmsd filtering to remove duplicate structures and then perform a energy-based filtering to remove
    both high energy outliers and low energy duplicates that are very similar to the lowest energy structure.
    This will ensure we have a diverse set of representative structures
    """
    def rmsd_filtering(self, threshold: float = 0.5, phase: Optional[str] = None):
        """   
        Filters out structures that are within a certain RMSD threshold of each other

        Args:
            threshold: RMSD threshold for filtering
            phase: Specified phase to filter. If None, filters all structures.
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return
        n_before = len(structures)
        self._log(f"RMSD filtering (threshold={threshold:.3f} Å) - Phase: {phase if phase else 'all'} - Starting with {n_before} structures")

        # Extract all coordinates
        all_coords = np.array([s.molecule.coordinates for s in structures]) 
        n_atoms = all_coords[0].shape[0]
        selected_indices = []
        selected_flat = np.empty((len(structures), n_atoms * 3), dtype=all_coords[0].dtype)
        n_selected = 0
        for i, coords in enumerate(all_coords):
            candidate_flat = coords.ravel()
            if n_selected == 0:
                selected_indices.append(i)
                selected_flat[n_selected] = candidate_flat
                n_selected += 1
                continue 
            diff = selected_flat[:n_selected] - candidate_flat
            rmsd_values = np.sqrt(np.mean(diff.reshape(n_selected, n_atoms, 3)**2, axis=(1,2)))

            if np.all(rmsd_values >= threshold):
                selected_indices.append(i)
                selected_flat[n_selected] = candidate_flat
                n_selected += 1

        filtered_structures = [structures[i] for i in selected_indices]
        if phase:
            self.phases[phase] = filtered_structures
        else:
            self.structures = filtered_structures
        n_after = len(filtered_structures)
        self._log(f"RMSD filtering completed - Phase: {phase if phase else 'all'} - {n_after} structures remaining (removed {n_before - n_after})")  




    def plot_com_trajectory_2d_projection(self,
                                          phase: Optional[str] = None,
                                          save_path: Optional[str] = None,
                                          simulation_box: Optional[SimulationBox] = None,
                                          show_box: bool = True,
                                          plot_type: str = "density",
                                          cmap: str = "viridis",
                                          gridsize: int = 10,
                                          show_trajectory: bool = False,
                                          alpha_trajectory: float = 0.1,
                                          separate_submolecules: bool = True) -> None:
        """   
        Plots 2D projections of submolecule center trajectories

        Args:
            phase: Which phase to plot ("A", "B", "all")
            save_path: Path to save the figure
            show_box: Whether to show the simulation box boundaries
            plot_type: "density" for hexbin density plot, "scatter" for scatter plot
                        "contour" for contour plot "kde" for kernel density estimate
            cmap: Colormap for density plot
            gridsize: Grid size for hexbin plot
            show_trajectory: Whether to show the trajectory lines between points
            alpha_trajectory: Alpha value for trajectory lines
            separate_submolecules: If True, plot each submolecule in separate row
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            print(f"No structures found for phase: {phase}")
            return
        
        n_submols = len(self.submolecule_indices)
        
        # Collect all COM positions per submolecules
        com_trajectories = [[] for _ in range(n_submols)]
        for mol in structures:
            for i, indices in enumerate(self.submolecule_indices):
                coords = mol.molecule.coordinates[indices]
                com = np.mean(coords, axis=0)
                com_trajectories[i].append(com)

        for i in range(n_submols):
            com_trajectories[i] = np.array(com_trajectories[i])
        
        projections = [
            ("X", "Y", 0, 1),
            ("X", "Z", 0, 2),
            ("Y", "Z", 1, 2)
        ]
        
        # Get box limits for consistent scaling
        if simulation_box is not None:
            if simulation_box.box_type == "sphere":
                lim = simulation_box.radius * 1.1
            else:
                lim = np.max(simulation_box.box_dimensions) / 2 * 1.1
        else:
            all_coords = np.vstack(com_trajectories)
            lim = np.max(np.abs(all_coords)) * 1.1

        # Different colormaps for each submolecule
        submol_cmaps = ['Blues', 'Oranges', 'Greens', 'Reds', 'Purples'][:n_submols]
        submol_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'][:n_submols]
        start_markers = ['o', 's', '^', 'D', 'v'][:n_submols]
        end_markers = ['*', 'P', 'X', 'h', '8'][:n_submols]

        if separate_submolecules and n_submols > 1:
            # Create grid: n_submols rows x 3 columns
            fig, axes = plt.subplots(n_submols, 3, figsize=(16, 5 * n_submols))
            if n_submols == 1:
                axes = axes.reshape(1, -1)
            
            hb_list = []  # Store hexbin objects for colorbars
            
            for submol_i in range(n_submols):
                traj = com_trajectories[submol_i]
                
                for col, (xlabel, ylabel, xi, yi) in enumerate(projections):
                    ax = axes[submol_i, col]
                    
                    ax.set_xlabel(f'{xlabel} (Å)', fontsize=11)
                    ax.set_ylabel(f'{ylabel} (Å)', fontsize=11)
                    ax.set_title(f"Submolecule {submol_i + 1} | {xlabel}-{ylabel}", 
                                fontsize=12, fontweight='bold')
                    ax.set_xlim(-lim, lim)
                    ax.set_ylim(-lim, lim)
                    ax.set_aspect('equal')
                    ax.grid(True, alpha=0.3, linestyle='--')

                    # Draw simulation box
                    if show_box and simulation_box is not None:
                        if simulation_box.box_type == "sphere":
                            circle = plt.Circle(
                                (simulation_box.center[xi], simulation_box.center[yi]),
                                simulation_box.radius,
                                fill=False, color="gray", linestyle='--', linewidth=1.5, alpha=0.7
                            )
                            ax.add_patch(circle)
                        elif simulation_box.box_type == "cube":
                            half = simulation_box.box_dimensions / 2
                            rec = plt.Rectangle(
                                (simulation_box.center[xi] - half[xi],
                                 simulation_box.center[yi] - half[yi]),
                                2 * half[xi], 2 * half[yi],
                                fill=False, color="gray", linestyle='--', linewidth=1.5, alpha=0.7
                            )
                            ax.add_patch(rec)
                    
                    x = traj[:, xi]
                    y = traj[:, yi]
                    
                    if plot_type == "density":
                        hb = ax.hexbin(x, y, gridsize=gridsize, cmap=submol_cmaps[submol_i], 
                                      mincnt=1, alpha=0.9, edgecolors='none')
                        if col == 2:  # Store last column hexbin for colorbar
                            hb_list.append(hb)
                    
                    elif plot_type == "scatter":
                        ax.scatter(x, y, color=submol_colors[submol_i], alpha=0.5, 
                                  s=20, edgecolors='none')
                    
                    elif plot_type == "contour":
                        H, xedges, yedges = np.histogram2d(x, y, bins=gridsize,
                                                          range=[[-lim, lim], [-lim, lim]])
                        X, Y = np.meshgrid((xedges[:-1] + xedges[1:]) / 2,
                                          (yedges[:-1] + yedges[1:]) / 2)
                        ax.contourf(X, Y, H.T, levels=10, cmap=submol_cmaps[submol_i], alpha=0.7)
                        ax.contour(X, Y, H.T, levels=5, colors='black', alpha=0.3, linewidths=0.5)
                    
                    if show_trajectory:
                        ax.plot(x, y, color=submol_colors[submol_i], alpha=alpha_trajectory, linewidth=0.5)
                    
                    # Mark start and end points
                    ax.scatter(x[0], y[0], color=submol_colors[submol_i], marker=start_markers[submol_i], 
                              s=150, edgecolors='black', linewidths=2, zorder=10, label='Start')
                    ax.scatter(x[-1], y[-1], color=submol_colors[submol_i], marker=end_markers[submol_i], 
                              s=200, edgecolors='black', linewidths=2, zorder=10, label='End')
                    
                    # Add legend only to last column
                    if col == 2:
                        ax.legend(loc='upper left', fontsize=9)
            
            # Add colorbars for each submolecule row
            plt.tight_layout()
            fig.subplots_adjust(right=0.92)
            
            if plot_type == "density" and hb_list:
                for i, hb in enumerate(hb_list):
                    # Position colorbar on the right side of each row
                    cbar_ax = fig.add_axes([0.94, 0.1 + (n_submols - 1 - i) * (0.8 / n_submols), 
                                           0.02, 0.7 / n_submols])
                    cbar = fig.colorbar(hb, cax=cbar_ax)
                    cbar.set_label(f'Submol {i+1} Visits', fontsize=10)
        
        else:
            # Original combined plot with overlay
            fig, axes = plt.subplots(1, 3, figsize=(16, 6))
            
            hb = None
            
            for ax, (xlabel, ylabel, xi, yi) in zip(axes, projections):
                ax.set_xlabel(f'{xlabel} (Å)', fontsize=11)
                ax.set_ylabel(f'{ylabel} (Å)', fontsize=11)
                ax.set_title(f"{xlabel}-{ylabel} Projection", fontsize=12, fontweight='bold')
                ax.set_xlim(-lim, lim)
                ax.set_ylim(-lim, lim)
                ax.set_aspect('equal')
                ax.grid(True, alpha=0.3, linestyle='--')

                # Draw simulation box
                if show_box and simulation_box is not None:
                    if simulation_box.box_type == "sphere":
                        circle = plt.Circle(
                            (simulation_box.center[xi], simulation_box.center[yi]),
                            simulation_box.radius,
                            fill=False, color="gray", linestyle='--', linewidth=1.5, alpha=0.7
                        )
                        ax.add_patch(circle)
                    elif simulation_box.box_type == "cube":
                        half = simulation_box.box_dimensions / 2
                        rec = plt.Rectangle(
                            (simulation_box.center[xi] - half[xi],
                             simulation_box.center[yi] - half[yi]),
                            2 * half[xi], 2 * half[yi],
                            fill=False, color="gray", linestyle='--', linewidth=1.5, alpha=0.7
                        )
                        ax.add_patch(rec)
                
                # Plot each submolecule with different appearance
                for i, traj in enumerate(com_trajectories):
                    x = traj[:, xi]
                    y = traj[:, yi]
                    
                    if plot_type == "density":
                        hb = ax.hexbin(x, y, gridsize=gridsize, cmap=submol_cmaps[i], 
                                      mincnt=1, alpha=0.6, edgecolors='none')
                    
                    elif plot_type == "scatter":
                        ax.scatter(x, y, color=submol_colors[i], alpha=0.5, 
                                  s=20, edgecolors='none', label=f'Submol {i+1}')
                    
                    if show_trajectory:
                        ax.plot(x, y, color=submol_colors[i], alpha=alpha_trajectory, linewidth=0.5)
                    
                    # Mark start and end points with different markers per submolecule
                    ax.scatter(x[0], y[0], color=submol_colors[i], marker=start_markers[i], 
                              s=150, edgecolors='black', linewidths=2, zorder=10)
                    ax.scatter(x[-1], y[-1], color=submol_colors[i], marker=end_markers[i], 
                              s=200, edgecolors='black', linewidths=2, zorder=10)

            # Create legend
            legend_elements = []
            for i in range(n_submols):
                legend_elements.append(plt.scatter([], [], c=submol_colors[i], s=100, 
                                                  marker=start_markers[i], edgecolors='black',
                                                  label=f'Start {i+1}'))
                legend_elements.append(plt.scatter([], [], c=submol_colors[i], s=120, 
                                                  marker=end_markers[i], edgecolors='black',
                                                  label=f'End {i+1}'))
            
            axes[2].legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), 
                          borderaxespad=0., title="Submolecules", fontsize=9)

            plt.tight_layout()
            
            # Add combined colorbar at bottom (note: only shows last submolecule's scale)
            if plot_type == "density" and hb is not None:
                fig.subplots_adjust(bottom=0.18)
                cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.025])
                cbar = fig.colorbar(hb, cax=cbar_ax, orientation='horizontal')
                cbar.set_label('Visit Count', fontsize=11)

        # Title
        fig.suptitle(f'Submolecule Center of Mass Distribution\n'
                     f'Phase: {phase.upper() if phase else "All"} | {len(structures)} structures | '
                     f'{n_submols} submolecules',
                     fontsize=13, fontweight='bold', y=1.02)
    
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
            print(f"Saved trajectory density plot to {save_path}")
    
        plt.close()
                   

    def _compute_com_trajectories(self, structures: List[StructureData]) -> List[np.ndarray]:
        """
        Compute geometric center (centroid) trajectories for each submolecule.
        
        Note: Using geometric mean (centroid) instead of mass-weighted COM.
        
        Args:
            structures: List of StructureData objects
            
        Returns:
            List of numpy arrays, each of shape (n_frames, 3) containing centroid coordinates
        """
        n_submols = len(self.submolecule_indices)
        n_frames = len(structures)
        
        centroid_trajectories = []
        
        for submol_idx in self.submolecule_indices:
            centroid_traj = np.zeros((n_frames, 3))
            
            for frame_idx, structure_data in enumerate(structures):
                mol = structure_data.molecule
                
                # Ensure coordinates are numpy array
                if isinstance(mol.coordinates, list):
                    mol_coords = np.array(mol.coordinates)
                else:
                    mol_coords = mol.coordinates
                
                # Extract submolecule coordinates
                submol_coords = mol_coords[submol_idx]
                
                # Calculate geometric center (simple mean of coordinates)
                centroid = np.mean(submol_coords, axis=0)
                centroid_traj[frame_idx] = centroid
            
            centroid_trajectories.append(centroid_traj)
        
        return centroid_trajectories


    @staticmethod
    def remove_string_digits(s: str) -> str:
        """   
        Utility function to remove digits from a string (e.g., "H2O" -> "HO")
        """
        return ''.join(filter(lambda x: not x.isdigit(), s))

    def compute_gyration_tensor_features(self, phase: Optional[str] = None) -> np.ndarray:
        """  
        Computes shape descriptors based on the eigenvalues of the gyration tensor for each structure

        + Asphericity: Measures the deviation of the molecule's shape from a perfect spherical symmetry. Values near zero indicate a perfectly spherical distribution
        + Acylindrity: Measures the deviation of the molecule's shape from cylindrical symmetry (rod-like). Values near zero indicate a perfect cylindrical distribution
        + Relative Shape Anisotropy: kappa^2 = 0 is spherically symmetric, kappa^2 = 1 is perfect rod shaped
        
        Args:
            phase: Specified phase to compute for. If None, computes for all structures.
        """ 
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        
        self._log(f"Computing gyration tensor based shape descriptors for {len(structures)} structures - Phase: {phase if phase else 'all'}")
        features = np.zeros((len(structures), 3))  # Columns: [asphericity, acylindricity, kappa^2]
        for i, s in enumerate(structures):
            coords = s.molecule.coordinates
            masses = s.molecule.masses
            com = GeometryOps.center_of_mass(coords, masses)
            centered = coords - com 

            # Gyration Tensor S_alpha,beta = 1/sum(m_i) * sum(m_i * r_i,alpha * r_i,beta)
            total_mass = np.sum(masses)
            S = np.zeros((3, 3))
            for j in range(len(masses)):
                S += masses[j] * np.outer(centered[j], centered[j])
            S /= total_mass
            eigenvalues = np.sort(np.linalg.eigvalsh(S))  # Sort eigenvalues in ascending order: lambda_1 <= lambda_2 <= lambda_3
            l1,l2,l3 = eigenvalues

            asphericity = l3 - 0.5 * (l1 + l2)
            acylindricity = l2 - l1
            trace = l1 + l2 + l3
            if trace > 1e-12:
                kappa_squared = (3/2) * (l1**2 + l2**2 + l3**2) / (trace**2) - 0.5
            else:
                kappa_squared = 0.0

            features[i] = [asphericity, acylindricity, kappa_squared]
        
        self._log(f"Gyration tensor features computed - Asphericity range: {features[:,0].min():.3f} to {features[:,0].max():.3f}, "
                    f"Acylindricity range: {features[:,1].min():.3f} to {features[:,1].max():.3f}, "
                    f"Kappa^2 range: {features[:,2].min():.3f} to {features[:,2].max():.3f}")
        return features


    # =========================== Hydrogen Bond Analysis ===========================



    def compute_hbond_features(self,
                               angle_threshold: float = 150.0,
                               max_distance: float = 3.5,
                               phase: Optional[str] = None,
                               n_processes: Optional[int] = None) -> np.ndarray:
        """ 
        Computes hydrogen bond features for all structures

        Computes:
            + Number of valid Hydrogen bond configuration 
            + Average Hydrogen bond angle (0.0 if no H-Bonds found)
            + Average Hydrogen bond distance (0.0 if no H-Bonds found)

        Uses multiprocessing to parallelize the computation across structures
            
        Args:
            angle_treshold: Minimum angle (in degrees) for a valid hydrogen bond
            max_distance: Maximum distance (in Å) for a valid hydrogen bond
            phase: Specified phase to compute for. If None, computes for all structures.
            n_processes: Number of parallel processes to use for computation. If None, uses all available cores.

        Returns:
            Numpy array of shape (n_structures, 3) with columns [n_hbonds, avg_angle, avg_da_distance]
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        n_structures = len(structures)
        self._log(f"Computing H-bond features for {n_structures} (angle_threshold={angle_threshold}°) - Phase: {phase if phase else 'all'}")

        if n_structures == 0:
            return np.zeros((0, 3))

        # Number of processes
        if n_processes is None:
            n_processes = min(cpu_count() - 2, n_structures)
        n_processes = max(1, min(n_processes, n_structures))
        self._log(f"Using {n_processes} parallel processes for H-bond computation")

        worker_args = []
        for structure_data in structures:
            mol = structure_data.molecule
            worker_args.append((
                list(mol.atom_labels),
                np.array(mol.coordinates),
                angle_threshold,
                max_distance
            ))
        if n_processes > 1:
            self._log("Starting multiprocessing pool for H-bond feature computation...")
            with Pool(processes=n_processes) as pool:
                results = pool.map(_compute_hbond_single, worker_args, chunksize=max(1, n_structures // (n_processes * 4)))
        else:
            self._log("Computing H-bond features sequentially (n_processes=1)...")
            results = list(map(_compute_hbond_single, worker_args))
        # Collect results
        features = np.zeros((n_structures, 3))
        total_hbonds = 0
        structures_with_hbonds = 0
        for i ,(n_hbonds, avg_angle, avg_da_distance) in enumerate(results):
            features[i, 0] = n_hbonds
            features[i, 1] = avg_angle
            features[i, 2] = avg_da_distance
            total_hbonds += n_hbonds
            if n_hbonds > 0:
                structures_with_hbonds += 1
        self._log(f"Structures with H-bonds: {structures_with_hbonds}/{n_structures}")
        self._log(f"Total H-bond configs: {total_hbonds}")
        if structures_with_hbonds > 0:
            self._log(f"Avg H-bonds per struct (with H-bonds): {total_hbonds / structures_with_hbonds:.2f}")
            hbond_angles = features[features[:, 0] > 0, 1]
            hbond_distances = features[features[:, 0] > 0, 2]
            self._log(f"Avg H-bond angle: {np.mean(hbond_angles):.2f}° (std: {np.std(hbond_angles):.2f}°)")
            self._log(f"Min/Max H-bond angle: {np.min(hbond_angles):.2f}° / {np.max(hbond_angles):.2f}°")
            self._log(f"Avg D-A distance: {np.mean(hbond_distances):.2f} Å (std: {np.std(hbond_distances):.2f} Å)")

        return features



    # ====================== Clustering Methods =======================

    def _ensure_feature_matrix(self):
        """   
        Compute and cache the feature matrix if not already done.
        """
        if self._feature_matrix_raw is not None:
            return 
        self._log("Computing and caching feature matrix...")
        delta_e = np.array([s.energy - self.get_lowest_energy_structure().energy for s in self.structures])
        rmsd_values = np.array([self._calculate_rmsd(s.molecule.coordinates, self.get_lowest_energy_structure().molecule.coordinates) for s in self.structures])
        rg_values = np.array([self.radius_of_gyration(s.molecule.coordinates, s.molecule.masses) for s in self.structures])
        rotational_constants = np.array([self.determine_rotational_constants(s.molecule.coordinates, s.molecule.masses) for s in self.structures])
        intermolecular_distances = np.array(self.compute_interatomic_distance_matrix())
        hbond_features = self.compute_hbond_features()
        gyration_tensor_features = self.compute_gyration_tensor_features()
        dipole_magnitudes = np.array([s.dipole_magnitude for s in self.structures])



        self._feature_matrix_raw = np.hstack((delta_e.reshape(-1, 1),
                                            rmsd_values.reshape(-1, 1), 
                                            rg_values.reshape(-1, 1), 
                                            rotational_constants, 
                                            intermolecular_distances,
                                            hbond_features,
                                            gyration_tensor_features,
                                            dipole_magnitudes.reshape(-1, 1)))
        self._log(f"Feature matrix shape: {self._feature_matrix_raw.shape}")
        self._log(f"Features: ['Delta E', 'RMSD', 'Rg', 'Rot A', 'Rot B', 'Rot C', 'Num H-bonds', 'Avg H-bond Angle', 'Avg D-A Distance', 'Intermolecular Distances...', 'Asphericity', 'Acylindricity', 'Kappa^2', 'Dipole Magnitude']")
        scaler = StandardScaler()
        self._feature_matrix_normalized = scaler.fit_transform(self._feature_matrix_raw)
        self._log("Feature matrix computed and cached.")




    def feature_matrix(self, normalize = True) -> np.ndarray:
        """   
        Current Feature Matrix wifth Columns

        1. Delta E to lowest energy structure
        2. RMSD to lowest energy structure
        3. Radius of Gyration
        4. Rotational Constants (A,B,C)
        5. Intermolecular Distance Matrix (flattened upper triangle)
        6. Number of valid H-bond configurations
        7. Average H-bond angle
        8. Average D-A distance
        9. Gyration Tensor Shape Descriptors (Asphericity, Acylindricity, Kappa^2)
        10. Dipole Magnitude (Debye)
        """
        self._ensure_feature_matrix()
        if normalize:
            return self._feature_matrix_normalized
        else:
            return self._feature_matrix_raw
        
        
    
    def agglomerative_clustering(self, n_clusters: int = 5, phase: Optional[str] = None, linkage: str = "ward") -> np.ndarray:
        """   
        Perform Agglomerative Clustering on the structures in the specified phase

        Args:
            n_clusters: Number of clusters to form
            phase: Specified phase to cluster. If None, clusters all structures.
            linkage: Linkage criterion for agglomerative clustering ("ward", "complete", "average", "single")
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            print(f"No structures found for phase: {phase}")
            return np.array([])

        self._log(f"Performing Agglomerative Clustering (n_clusters={n_clusters}, linkage='{linkage}') for phase: {phase if phase else 'all'} - n_structures: {len(structures)}")
        
        feature_matrix = self.feature_matrix()
        clustering_model = SklearnAgglomerativeClustering(n_clusters=n_clusters, linkage=linkage)
        labels = clustering_model.fit_predict(feature_matrix)
        self.labels = labels

        unique, counts = np.unique(labels, return_counts=True)
        self._log(f"Cluster distribution:")
        for cluster_id, count in zip(unique, counts):
            cluster_energies = [s.energy for s,l in zip(structures, labels) if l == cluster_id]
            self._log(f" Cluster {cluster_id}: {count} structures | Energy range: {min(cluster_energies):.2f} to {max(cluster_energies):.2f}")
    

        return labels
    
    def pca(self, n_components: int = 2, phase: Optional[str] = None) -> np.ndarray:
        """   
        Performs PCA Analysis on the structures in the specified phase

        Args:
            n_components: Number of principal components to compute
            phase: Specified phase to analyze. If None, analyzes all structures.
        """ 
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            print(f"No structures found for phase: {phase}")
            return np.array([])
        feature_matrix = self.feature_matrix()
        pca_model = PCA(n_components=n_components)
        pca_result = pca_model.fit_transform(feature_matrix)
        explained_variance = pca_model.explained_variance_ratio_
        return pca_result, explained_variance

    def plot_pca_explained_variance(self, n_components: int = 20, phase: Optional[str] = None):
        """  
        Plots the explained variance ratio for PCA components
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            print(f"No structures found for phase: {phase}")
            return
        pca_result, explained_variance = self.pca(n_components=n_components, phase=phase)
        plt.figure(figsize=(10, 6))
        plt.bar(range(1, n_components + 1), explained_variance * 100, color='skyblue')
        plt.xlabel('Principal Component')
        plt.ylabel('Explained Variance (%)')
        plt.title(f'PCA Explained Variance - Phase: {phase if phase else "All"}')
        plt.xticks(range(1, n_components + 1))
        plt.grid(True, alpha=0.3)
        plt.savefig(f"figures/pca_explained_variance_{phase if phase else 'all'}.png")
        plt.close()


    def plot_pca_agglomerative(self, n_clusters: int = 5, phase: Optional[str] = None, linkage: str = "ward"):
        """   
        Plots the PCA projection of the structures colored by Agglomerative Clustering labels

        Args:
            n_clusters: Number of clusters to form
            phase: Specified phase to analyze. If None, analyzes all structures.
            linkage: Linkage criterion for agglomerative clustering ("ward", "complete", "average", "single")
        """
        pca_result, explained_variance = self.pca(phase=phase)
        # Use cached labels if available with matching n_clusters, otherwise compute
        if self.labels is None or len(np.unique(self.labels)) != n_clusters:
            labels = self.agglomerative_clustering(n_clusters=n_clusters, phase=phase, linkage=linkage)
        else:
            labels = self.labels
        plt.figure(figsize=(10, 6))
        scatter = plt.scatter(pca_result[:, 0], pca_result[:, 1], c=labels, cmap='tab10', alpha=0.7)
        plt.xlabel(f'PC1 ({explained_variance[0]*100:.1f}%)')
        plt.ylabel(f'PC2 ({explained_variance[1]*100:.1f}%)')
        plt.title(f'PCA Projection with Agglomerative Clustering - Phase: {phase if phase else "All"}')
        plt.colorbar(scatter, label='Cluster Label')
        plt.grid(True, alpha=0.3)
        plt.savefig(f"figures/pca_agglomerative_{phase if phase else 'all'}.png")
        plt.close()

    def get_cluster_representatives(self, phase: Optional[str] = None, method = "lowest_energy") -> Dict[int, StructureData]:
        """   
        Gets the representative structure for each cluster in the specified phase

        Methods:
            "closest_to_centroid": Returns the structure closest to the cluster centroid
            "lowest_energy": Returns the structure with the lowest energy in each cluster
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return {}
        if self.labels is None:
            print("Clustering has not been performed yet. Please run Clustering method first.")
            return {}
        cluster_representatives = {}
        for cluster_id in np.unique(self.labels):
            cluster_structures = [s for s, label in zip(structures, self.labels) if label == cluster_id]
            if method == "lowest_energy":
                representative = min(cluster_structures, key=lambda s: s.energy)
                cluster_representatives[cluster_id] = representative
        
        # Return a list of molecules 
        cluster_representatives_mol = [rep.molecule for rep in cluster_representatives.values()]
        return cluster_representatives_mol

    def tsne(self, n_components: int = 2, perplexity: float = 30.0, n_iter: int = 1000, random_state: int = 42, phase: Optional[str] = None) -> np.ndarray:
        """   
        Performs t-SNE dimensionality reduction on the feature matrix

        t-SNE (t-distributed Stochastic Neighbor Embedding) is a non-linear dimensionality reduction technique for visualization of high-dimenional data

        Args:
            n_components: Number of dimensions to reduce to (typically 2 or 3 for visualization)
            perplexity: Perplexity parameter for t-SNE (related to number of nearest neighbors)
            n_iter: Number of iterations for optimization
            random_state: Random seed for reproducibility
            phase: Specified phase to analyze. If None, analyzes all structures.
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            return np.array([])
        feature_mat = self.feature_matrix()

        # Adjust perplexity < n_samples
        n_samples = feature_mat.shape[0]
        effective_perplexity = min(perplexity, max(5.0, n_samples - 1))
        if effective_perplexity != perplexity:
            perplexity = effective_perplexity
        tsne_model = TSNE(n_components=n_components, perplexity=perplexity, n_iter=n_iter, random_state=random_state)
        tsne_result = tsne_model.fit_transform(feature_mat)
        return tsne_result

    def umap(self,
             n_components: int = 2,
             n_neighbors: int = 15,
             min_dist: float = 0.1,
             metric: str = 'euclidean',
             random_state: int = 42,
             phase: Optional[str] = None) -> np.ndarray:
        """  
        Performs UMAP dimensionality reduction on the feature matrix
        
        Args:
            n_components: Number of dimensions to reduce to (typically 2 or 3 for visualization)
            n_neighbors: Number of neighbors to consider for local structure
            min_dist: Minimum distance between points in the embedding space
            metric: Distance metric to use (e.g., 'euclidean', 'manhattan')
            random_state: Random seed for reproducibility
            phase: Specified phase to analyze. If None, analyzes all structures.
        """
        if UMAP is None:
            raise ImportError("UMAP is not installed. Please install umap-learn to use this method.")
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            return np.array([])
        feature_mat = self.feature_matrix()
        n_samples = feature_mat.shape[0]
        effective_n_neighbors = min(n_neighbors, n_samples - 1)
        if effective_n_neighbors != n_neighbors:
            self._log(f"Adjusting n_neighbors from {n_neighbors} to {effective_n_neighbors} due to number of samples ({n_samples})", level="warning")
            n_neighbors = effective_n_neighbors
        
        umap_model = UMAP(
            n_components = n_components,
            n_neighbors = n_neighbors,
            min_dist = min_dist,
            metric = metric,
            random_state = random_state
        )
        umap_result = umap_model.fit_transform(feature_mat) 
        return umap_result

    def plot_umap_agglomerative(self, n_clusters: int = 5, 
                                phase: Optional[str] = None,
                                linkage: str = "ward", 
                                random_state: int = 42, 
                                colored_by: str = "cluster") -> None:
        """  
        Plots the UMAP projection of the structures colored by Agglomerative Clustering labels

        Creates a 2x2 grid with different n_neighbors values [10, 15, 30, 50]
        
        Args:
            n_clusters: Number of clusters to form
            phase: Specified phase to analyze. If None, analyzes all structures.
            linkage: Linkage criterion for agglomerative clustering ("ward", "complete", "average", "single")
            random_state: Random seed for reproducibility
            colored_by: Whether to color points by "cluster" labels or "energy" values
        """
        if UMAP is None:
            raise ImportError("UMAP is not installed. Please install umap-learn to use this method.")
    
        # Compute clustering once before the loop
        if self.labels is None or len(np.unique(self.labels)) != n_clusters:
            labels = self.agglomerative_clustering(n_clusters=n_clusters, phase=phase, linkage=linkage)
        else:
            labels = self.labels
        
        n_neighbors_list = [10, 15, 30, 50]
        fig,axes = plt.subplots(2, 2, figsize=(16, 12))
        for ax, n_neighbors in zip(axes.flatten(), n_neighbors_list):
            umap_result = self.umap(n_components=2, n_neighbors=n_neighbors, random_state=random_state, phase=phase)
            if colored_by == "cluster":
                scatter = ax.scatter(umap_result[:, 0], umap_result[:, 1], c=labels, cmap='tab10', alpha=0.7)
                ax.set_title(f'UMAP (n_neighbors={n_neighbors}) - Colored by Cluster', fontsize=12)
                plt.colorbar(scatter, ax=ax, label='Cluster Label')
            elif colored_by == "energy":
                energies = np.array([s.energy for s in (self.phases[phase] if phase else self.structures)])
                scatter = ax.scatter(umap_result[:, 0], umap_result[:, 1], c=energies, cmap='viridis', alpha=0.7)
                ax.set_title(f'UMAP (n_neighbors={n_neighbors}) - Colored by Energy', fontsize=12)
                plt.colorbar(scatter, ax=ax, label='Energy')
            ax.set_xlabel('UMAP Dimension 1')
            ax.set_ylabel('UMAP Dimension 2')
            ax.grid(True, alpha=0.3)
        plt.suptitle(f'UMAP Projection with Agglomerative Clustering - Phase: {phase if phase else "All"}', fontsize=14, fontweight='bold')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(f"figures/umap_agglomerative_{phase if phase else 'all'}_{colored_by}.png")
        plt.close()



    def plot_tsne_agglomerative(self,
                                n_clusters: int = 5,
                                phase: Optional[str] = None,
                                linkage: str = "ward",
                                n_iter: int = 1000,
                                random_state: int = 42,
                                colored_by: str = "cluster") -> None:
        """
        Plots the t-SNE projection of the structures colored by Agglomerative Clustering labels

        This makes 4x4 subplots with perplexities [20, 40, 60, 80]
        Args:
            n_clusters: Number of clusters to form
            phase: Specified phase to analyze. If None, analyzes all structures.
            linkage: Linkage criterion for agglomerative clustering ("ward", "complete", "average", "single")
            n_iter: Number of iterations for t-SNE optimization
            random_state: Random seed for reproducibility
            colored_by: Whether to color points by "cluster" labels or "energy" values
        """
        # Compute clustering once before the loop
        if self.labels is None or len(np.unique(self.labels)) != n_clusters:
            labels = self.agglomerative_clustering(n_clusters=n_clusters, phase=phase, linkage=linkage)
        else:
            labels = self.labels

        perplexities = [20, 40, 60, 80]
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        for ax, perplexity in zip(axes.flatten(), perplexities):
            tsne_result = self.tsne(n_components=2, perplexity=perplexity, n_iter=n_iter, random_state=random_state, phase=phase)
            if colored_by == "cluster":
                scatter = ax.scatter(tsne_result[:, 0], tsne_result[:, 1], c=labels, cmap='tab10', alpha=0.7)
                ax.set_title(f't-SNE (perplexity={perplexity}) - Colored by Cluster', fontsize=12)
                plt.colorbar(scatter, ax=ax, label='Cluster Label')
            elif colored_by == "energy":
                energies = np.array([s.energy for s in (self.phases[phase] if phase else self.structures)])
                scatter = ax.scatter(tsne_result[:, 0], tsne_result[:, 1], c=energies, cmap='viridis', alpha=0.7)
                ax.set_title(f't-SNE (perplexity={perplexity}) - Colored by Energy', fontsize=12)
                plt.colorbar(scatter, ax=ax, label='Energy')
            ax.set_xlabel('t-SNE Dimension 1')
            ax.set_ylabel('t-SNE Dimension 2')
            ax.grid(True, alpha=0.3)
        plt.suptitle(f't-SNE Projection with Agglomerative Clustering - Phase: {phase if phase else "All"}', fontsize=14, fontweight='bold')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(f"figures/tsne_agglomerative_{phase if phase else 'all'}_{colored_by}.png")
        plt.close()









if __name__ == "__main__":
    import os
    import sys 
    def parse_multi_xyz(filepath: str) -> List[Tuple[Molecule, float]]:
        """  
        Parse a multi-structure XYZ file into a list of (Molecule, energy)
        """
        structures = []
        with open(filepath, 'r') as f:
            lines = f.readlines()

        i = 0
        sample_idx = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line:
                i += 1 
                continue
            n_atoms = int(line)
            comment = lines[i+1].strip()
            mol = Molecule(name=f"Sample_{sample_idx}")
            for j in range(n_atoms):
                parts = lines[i+2+j].split()
                label = parts[0]
                coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                mol.add_atom(label, coords)
            # Use sample index as a dummy energy
            energy = float(sample_idx) * 0.01
            structures.append((mol,energy))
            i += 2 + n_atoms
            sample_idx += 1
        return structures
    
    xyz_path = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/test_structures.xyz"
    submolecule_indices = [[0,1,2], [3,4,5]]  # Example submolecule indices for two water molecule
    structures = parse_multi_xyz(xyz_path)
    print(f"Loading Structures from {xyz_path}...")
    print(f"Parsed {len(structures)} structures.")
    

    multiplies = 10 # More structures for more effective speed comparison
    structures_large = structures * multiplies
    n_large = len(structures_large)
    print(f"Testing H-bond feature computation on {n_large} structures with multiprocessing...")
    analyzer = BHMCAnalyzer(submolecule_indices=submolecule_indices)
    analyzer.add_structures_batch(structures_large)

    # Test donor acceptor dist
    features = analyzer.compute_hbond_features(angle_threshold=150.0, max_distance=3.5, n_processes=4)
    print(f"Computed H-bond features for {features.shape[0]} structures. Sample output (first 5 rows):")
    print("n_hbonds | avg_angle | avg_donor-acceptor_distance")
    print(features[:5]) 


#    # Serial
#    t0 = time.time()
#    features_serial = analyzer.compute_hbond_features(n_processes=1)
#    t_serial = time.time() - t0
#    print(f"Serial computation time: {t_serial:.2f} seconds")
#
#    procs = [2, 4, 8, 16]
#    for n_procs in procs:
#        t0 = time.time()
#        features_parallel = analyzer.compute_hbond_features(n_processes=n_procs)
#        t_parallel = time.time() - t0
#        print(f"Parallel computation time with {n_procs} processes: {t_parallel:.2f} seconds | Speedup: {t_serial / t_parallel:.2f}x")


