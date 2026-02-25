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

from molecule_class import Molecule
from transformations import GeometryOps
from box import SimulationBox

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
    metadata: Dict[str, Any] = field(default_factory=dict)

class BHMCAnalyzer:
    """   
    Analyzer for the Basin Hopping Monte Carlo results
    """
    def __init__(self, name: str = "BHMC_Analysis", submolecule_indices: Optional[List[int]] = None):
        self.name = name 
        self.structures: List[StructureData] = []
        self.phases: Dict[str, List[StructureData]] = defaultdict(list)
        self.submolecule_indices = submolecule_indices



    def add_structure(self,
                      molecule: Molecule,
                      energy: float, 
                      phase: str = "unknown",
                      **kwargs):
        """   
        Add a structure to the analysis

        Args:
            molecule: Molecule object representing the structure
            energy: Energy of the structure
            phase: Phase of the BHMC process (e.g., "phase_a", "phase_b", "phase_c")
            **kwargs: Additional properties (e.g., step, worker_id, accepted)
        """
        structure_data = StructureData(
            molecule=molecule,
            energy=energy,
            phase=phase,
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
            structures: List of tuples containing (Molecule, energy)
            phase: Phase of the BHMC process for all structures
        """
        for mol, energy in structures:
            self.add_structure(mol, energy, phase=phase)

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
            return {}
        
        energies = np.array(energies)
        return {
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
        coords_list = [s.molecule.coordinates for s in structures]
        rmsd_matrix = np.zeros((n, n))
        for i in range(n):
            ci = coords_list[i]
            for j in range(i+1, n):
                rmsd = self._calculate_rmsd(ci, coords_list[j])
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd
        return rmsd_matrix



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
            print(f"No structures found for phase: {phase}")
            return
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


    def plot_com_trajectory_2d_projections(
            self,
            simulation_box: Optional[SimulationBox] = None,
            save_path: str = "figures/centroid_trajectory_2d_projections.png",
            phase: Optional[str] = None
        ):
        """
        Plot 2D projections of geometric center (centroid) trajectories onto XY, XZ, and YZ planes.
        """
        # Get structures
        if phase:
            structures = self.phases.get(phase, [])
        else:
            structures = self.structures
        
        if not structures:
            print(f"No structures available for phase: {phase}")
            return
        
        if not self.submolecule_indices:
            print("No submolecule indices available.")
            return
        
        n_frames = len(structures)
        n_submols = len(self.submolecule_indices)
        
        print(f"\nPlotting 2D projections for {n_frames} structures, {n_submols} submolecules...")
        
        # Compute centroid trajectories (geometric center)
        centroid_trajectories = self._compute_com_trajectories(structures)
        
        # Debug output - check if coordinates are actually changing
        for i, traj in enumerate(centroid_trajectories):
            print(f"\n  Submol {i+1} trajectory:")
            print(f"    Shape: {traj.shape}")
            print(f"    X range: [{traj[:, 0].min():.4f}, {traj[:, 0].max():.4f}] (spread: {traj[:, 0].max() - traj[:, 0].min():.4f})")
            print(f"    Y range: [{traj[:, 1].min():.4f}, {traj[:, 1].max():.4f}] (spread: {traj[:, 1].max() - traj[:, 1].min():.4f})")
            print(f"    Z range: [{traj[:, 2].min():.4f}, {traj[:, 2].max():.4f}] (spread: {traj[:, 2].max() - traj[:, 2].min():.4f})")
            
            # Check if trajectory is actually moving
            movement = np.sqrt(np.sum((traj[-1] - traj[0])**2))
            print(f"    Total displacement (start to end): {movement:.4f} Å")
            
            # Check variance in each dimension
            print(f"    X variance: {np.var(traj[:, 0]):.6f}")
            print(f"    Y variance: {np.var(traj[:, 1]):.6f}")
            print(f"    Z variance: {np.var(traj[:, 2]):.6f}")
        
        # Create figure with 3 subplots
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        # Define projections
        projections = [
            ('X', 'Y', 0, 1),
            ('X', 'Z', 0, 2),
            ('Y', 'Z', 1, 2)
        ]
        
        colors = plt.cm.tab10.colors
        
        # Plot each projection
        for ax, (label1, label2, idx1, idx2) in zip(axes, projections):
            print(f"\nPlotting {label1}-{label2} projection...")
            
            # Plot trajectory for each submolecule
            for i, centroid_traj in enumerate(centroid_trajectories):
                color = colors[i % len(colors)]
                
                # Extract coordinates for this projection
                x_coords = centroid_traj[:, idx1]
                y_coords = centroid_traj[:, idx2]
                
                print(f"  Submol {i+1}: plotting {len(x_coords)} points")
                print(f"    {label1} coords: min={x_coords.min():.4f}, max={x_coords.max():.4f}")
                print(f"    {label2} coords: min={y_coords.min():.4f}, max={y_coords.max():.4f}")
                
                # Plot trajectory line with increased visibility
                line = ax.plot(x_coords, y_coords,
                              color=color, alpha=0.8, linewidth=3,
                              label=f'Submol {i+1}', 
                              marker='o', markersize=2, markevery=10)[0]
                
                print(f"    Line plotted: {line}")
                
                # Mark start point (large circle)
                start = ax.scatter(x_coords[0], y_coords[0],
                                  color=color, marker='o', s=150,
                                  edgecolors='black', linewidths=2,
                                  zorder=10, label=f'Start {i+1}')
                
                # Mark end point (star)
                end = ax.scatter(x_coords[-1], y_coords[-1],
                                color=color, marker='*', s=200,
                                edgecolors='black', linewidths=2,
                                zorder=10, label=f'End {i+1}')
                
                print(f"    Start marker at: ({x_coords[0]:.4f}, {y_coords[0]:.4f})")
                print(f"    End marker at: ({x_coords[-1]:.4f}, {y_coords[-1]:.4f})")
            
            # Draw simulation box projection
            if simulation_box is not None:
                self._draw_box_projection_2d(ax, simulation_box)
            
            # Formatting
            ax.set_xlabel(f'{label1} (Å)', fontsize=12)
            ax.set_ylabel(f'{label2} (Å)', fontsize=12)
            ax.set_title(f'{label1}-{label2} Projection', fontsize=13, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='--')
            ax.set_aspect('equal', adjustable='box')
            
            # Show legend
            ax.legend(fontsize=9, loc='best', framealpha=0.9)
            
            print(f"  Axis limits: x=[{ax.get_xlim()}], y=[{ax.get_ylim()}]")
        
        # Main title
        plt.suptitle(
            f'Geometric Center Trajectory 2D Projections\n'
            f'Phase: {phase if phase else "All"} | {n_frames} structures | {n_submols} submolecules',
            fontsize=14, fontweight='bold', y=1.02
        )
        
        plt.tight_layout()
        
        # Save figure
        from pathlib import Path
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"\n✓ Saved 2D projections to: {save_path}")
        
        # Also show the plot to verify
        # plt.show()
        
        plt.close()
    
    @staticmethod
    def _draw_box_projection_2d(ax, simulation_box: SimulationBox):
        """Draw 2D projection of simulation box."""
        if simulation_box.box_type == "sphere":
            # Draw circle
            circle = plt.Circle(
                (0, 0),
                simulation_box.radius,
                fill=False,
                edgecolor='gray',
                linewidth=2,
                linestyle='--',
                alpha=0.5,
                zorder=1
            )
            ax.add_patch(circle)
        
        elif simulation_box.box_type == "cube":
            # Draw square
            L = simulation_box.box_length
            rect = plt.Rectangle(
                (-L/2, -L/2),
                L,
                L,
                fill=False,
                edgecolor='gray',
                linewidth=2,
                linestyle='--',
                alpha=0.5,
                zorder=1
            )
            ax.add_patch(rect)
    
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

    # ====================== Clustering Methods =======================

    def feature_matrix(self, normalize = True) -> np.ndarray:
        """   
        Current Feature Matrix wifth Columns

        1. Delta E to lowest energy structure
        2. RMSD to lowest energy structure
        3. Radius of Gyration
        4. Rotational Constants (A,B,C)
        5. Intermolecular Distance Matrix (flattened upper triangle)
        """
        delta_e = np.array([s.energy - self.get_lowest_energy_structure().energy for s in self.structures])
        rmsd_values = np.array([self._calculate_rmsd(s.molecule.coordinates, self.get_lowest_energy_structure().molecule.coordinates) for s in self.structures])
        rg_values = np.array([self.radius_of_gyration(s.molecule.coordinates, s.molecule.masses) for s in self.structures])
        rotational_constants = np.array([self.determine_rotational_constants(s.molecule.coordinates, s.molecule.masses) for s in self.structures])
        intermolecular_distances = np.array(self.compute_interatomic_distance_matrix())

        feature_matrix = np.hstack((delta_e.reshape(-1, 1), rmsd_values.reshape(-1, 1), rg_values.reshape(-1, 1), rotational_constants, intermolecular_distances)) 

        # Normalize using StandardScaler
        scaler = StandardScaler()
        if normalize:
            feature_matrix = scaler.fit_transform(feature_matrix)


        return feature_matrix
    
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
        feature_matrix = self.feature_matrix()
        clustering_model = SklearnAgglomerativeClustering(n_clusters=n_clusters, linkage=linkage)
        labels = clustering_model.fit_predict(feature_matrix)
        self.labels = labels
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
        return pca_result
    
    def plot_pca_agglomerative(self, n_clusters: int = 5, phase: Optional[str] = None, linkage: str = "ward"):
        """   
        Plots the PCA projection of the structures colored by Agglomerative Clustering labels

        Args:
            n_clusters: Number of clusters to form
            phase: Specified phase to analyze. If None, analyzes all structures.
            linkage: Linkage criterion for agglomerative clustering ("ward", "complete", "average", "single")
        """
        pca_result = self.pca(phase=phase)
        labels = self.agglomerative_clustering(n_clusters=n_clusters, phase=phase, linkage=linkage)
        plt.figure(figsize=(10, 6))
        sns.scatterplot(x=pca_result[:, 0], y=pca_result[:, 1], hue=labels, palette='Set2', alpha=0.7)
        plt.title(f'PCA Projection with Agglomerative Clustering - Phase: {phase if phase else "All"}')
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.legend(title='Cluster')
        plt.grid(True)
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
            print(f"No structures found for phase: {phase}")
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
            
     
    






