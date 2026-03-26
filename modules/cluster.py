"""  
Module for various clustering methods and analysis of the sampled molecular configurations
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from typing import List, Tuple, Optional, Dict, Any, Union
from dataclasses import dataclass, field
from collections import defaultdict


from sklearn.cluster import AgglomerativeClustering as SklearnAgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN
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
try:
    from hdbscan import HDBSCAN
except ImportError:
    HDBSCAN = None

from molecule_class import Molecule
from geometry import GeometryOps
from box import SimulationBox
from logger import Logger
from descriptors import FeatureExtractor

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
        Calculate optimal RMSD between two sets of coordinate points.
        For this we use the respective function of Geometry Ops
        No centering is required because we use the optimal correspondence function
        """
        coords_2_opt = GeometryOps.find_optimal_correspondence(coords1, coords2)
        diff = coords1 - coords_2_opt
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
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
    

    # ====================== Data Cleanup ========================

    """ 
    We should first of perform a rmsd filtering to remove duplicate structures and then perform a energy-based filtering to remove
    both high energy outliers and low energy duplicates that are very similar to the lowest energy structure.
    This will ensure we have a diverse set of representative structures
    """
    def rmsd_filtering(self, threshold: float = 0.2, phase: Optional[str] = None):
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

        # Base Energy and Dipoles
        lowest_energy = self.get_lowest_energy_structure()
        delta_e = np.array([s.energy - lowest_energy.energy for s in self.structures]).reshape(-1, 1)
        dipole_mags = np.array([s.dipole_magnitude for s in self.structures]).reshape(-1, 1)
        
        # H-Bond Features
        hbond_features = self.compute_hbond_features()

        # Extractor for Geometric Features
        extractor = FeatureExtractor()
        geom_features = extractor.extract_fast_features(self.structures, self.submolecule_indices)

        rmsd_values = np.array([self._calculate_rmsd(s.molecule.coordinates, self.get_lowest_energy_structure().molecule.coordinates) for s in self.structures])
        rotational_constants = np.array([self.determine_rotational_constants(s.molecule.coordinates, s.molecule.masses) for s in self.structures])
        intermolecular_distances = np.array(self.compute_interatomic_distance_matrix())



        self._feature_matrix_raw = np.hstack((
            delta_e, 
            rmsd_values.reshape(-1, 1), 
            geom_features[:, :1],  # Rg
            rotational_constants, 
            hbond_features, 
            intermolecular_distances,
            geom_features[:, 1:4],  # Asphericity, Acylindricity, Kappa^2
            dipole_mags
        )

        )
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
        
    
    # ================================= Clustering Algorithms =================================

    def cluster(self,
                method: str = "agglomerative",
                n_clusters: Optional[int] = None,
                phase: Optional[str] = None,
                **kwargs) -> np.ndarray:
        """  
        Clustering Interface 

        Args:
            method: Clustering method to use ("agglomerative", "kmeans")
            n_clusters: Number of clusters to form
            phase: Phase to cluster. If None, clusters all structures.
            **kwargs: Additional keyword arguments:
                 - For Agglomerative Clustering: linkage (default "ward")
                 - For kmeans: n_init (int), max_iter (int), random_state (int)
                 - For DBSCAN: eps (float), min_samples (int)
                 - For HDBSCAN: min_cluster_size (int), min_samples (int)
        Returns:
            Cluster labels array
        """
        if method == "agglomerative":
            if n_clusters is None:
                raise ValueError("n_clusters must be specified for agglomerative clustering")
            return self.agglomerative_clustering(n_clusters=n_clusters, phase=phase, **kwargs)
        elif method == "kmeans":
            if n_clusters is None:
                raise ValueError("n_clusters must be specified for kmeans clustering")
            return self.kmeans_clustering(n_clusters=n_clusters, phase=phase, **kwargs)
        elif method == "dbscan":
            return self.dbscan_clustering(**kwargs, phase=phase)
        elif method == "hdbscan":
            return self.hdbscan_clustering(**kwargs, phase=phase)
        
        else:
            raise ValueError(f"Unsupported clustering method: {method}. Supported methods: 'agglomerative', 'kmeans'")



    def dbscan_clustering(self,
                          eps: float = 0.5,
                          min_samples: int = 5,
                          phase: Optional[str] = None) -> np.ndarray:
        """  
        Perform DBSCAN clustering on the structures

        DBSCAN is density-based and does not require n_clusters. It discovers clusters from the
        data based on eps (neighborhood radius) and min_samples
        Noise Points are labeled as -1

        Args:
            eps: The maximum distance between two samples for them to be considered as in the same neighborhood.
            min_samples: The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. This includes the point itself.
            phase: Specified phase to cluster. If None, clusters all structures.
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return np.array([])
        self._log(f"Performing DBSCAN Clustering (eps={eps}, min_samples={min_samples}) for phase: {phase if phase else 'all'} - n_structures: {len(structures)}")
        feature_matrix = self.feature_matrix()
        model = DBSCAN(eps=eps, min_samples=min_samples)
        labels = model.fit_predict(feature_matrix)
        self.labels = labels
        n_clusters = len(set(labels) - {-1}) # Exclude noise label (-1)
        n_noise = np.sum(labels == -1)
        self._log(f"DBSCAN found {n_clusters} clusters and {n_noise} noise points")
        return labels

    def hdbscan_clustering(self,
                           min_cluster_size: int = 5,
                           min_samples: Optional[int] = None,
                           phase: Optional[str] = None) -> np.ndarray:
        """ 
        Perform HDBSCAN clustering on the structures

        HDBSCAN is hierarchical density-based clustering that automatically determines the number of clusters.
        It handles varying density better than DBSCAN
        Noise points are labeled as -1
        """
        if HDBSCAN is None:
            raise ImportError("HBDSCAN is not installed")

        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return np.array([])
        if min_samples is None:
            min_samples = min_cluster_size

        self._log(f"Performing HDBSCAN Clustering (min_cluster_size={min_cluster_size}, min_samples={min_samples}) for phase: {phase if phase else 'all'} - n_structures: {len(structures)}")
        feature_matrix = self.feature_matrix()
        model = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples)
        labels = model.fit_predict(feature_matrix)
        self.labels = labels
        n_clusters = len(set(labels) - {-1}) # Exclude noise label (-1)
        n_noise = np.sum(labels == -1)
        self._log(f"HDBSCAN found {n_clusters} clusters and {n_noise} noise points")
        return labels
    

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
        self.labels = labels  # Cache labels for later use
        return labels
    
    def kmeans_clustering(self,
                          n_clusters: int = 5,
                          phase: Optional[str] = None,
                          n_init: int = 10,
                          max_iter: int = 300,
                          random_state: Optional[int] = 42) -> np.ndarray:
        """ 
        Perform KMeans clustering on the structures

        Args:
            n_clusters: Number of clusters to form
            phase: Specified phase to cluster. If None, clusters all structures.
            n_init: Number of time the k-means algorithm will be run with different centroid seeds. The final results will be the best output of n_init consecutive runs in terms of inertia.
            max_iter: Maximum number of iterations of the k-means algorithm for a single run.
            random_state: Determines random number generation for centroid initialization. Use an int to make the randomness deterministic.
        
        Returns:
            Cluster labels array
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return np.array([])
        
        self._log(f"Performing KMeans Clustering (n_clusters={n_clusters}, n_init={n_init}, max_iter={max_iter}) for phase: {phase if phase else 'all'} - n_structures: {len(structures)}")
        feature_matrix = self.feature_matrix()
        model = KMeans(n_clusters=n_clusters,
                       n_init=n_init,
                       max_iter=max_iter,
                       random_state=random_state)
        labels = model.fit_predict(feature_matrix)
        self.labels = labels  # Cache labels for later use

        return labels


    # ================================ Obtain the Cluster Representatives ================================
    def get_cluster_representatives(self, method: str = "centroid", labels: Optional[np.ndarray] = None) -> List[StructureData]:
        """  
        Get the representative structures for each cluster

        Method == "centroid": For each cluster computes the centroid in the normalized feature space and selects the structure closest to the centroid (Euclidean distance)

        Args:
            labels: Cluster^er labels array. If None, uses cached labels from the last clustering run. If no cached labels, raises  ValueError
        
        Returns:
            List of StructureData objects representing the cluster centers
        """
        if labels is None:
            labels = self.labels
        if labels is None:
            raise ValueError("No cluster labels available. Run a clustering method first")

        feature_mat = self.feature_matrix(normalize=True)
        unique_labels = np.unique(labels)
        representatives = []
        for label in unique_labels:
            if label == -1:
                # Skip noise points form DBSCAN/HDBSCAN
                continue
            mask = labels == label
            cluster_features  = feature_mat[mask]
            cluster_indices = np.where(mask)[0]

            centroid = cluster_features.mean(axis=0)
            distances = np.linalg.norm(cluster_features - centroid, axis=1)
            closest_idx = cluster_indices[np.argmin(distances)]
            rep = self.structures[closest_idx]
            self._log(f"Cluster {label}: representative index {closest_idx}, distance to centroid: {distances.min():.3f}, E={rep.energy:.6f} Ha")
            representatives.append(rep)
        
        self._log(f"Total clusters (excluding noise): {len(unique_labels) - (1 if -1 in unique_labels else 0)}. Representatives obtained: {len(representatives)}")
        return representatives

    

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
        


    def tsne(self, n_components: int = 2, perplexity: float = 30.0, n_iter: int = 1000, random_state: int = 42, phase: Optional[str] = None) -> np.ndarray:
        """   
        Performs t-SNE dimensionality reduction on the feature matrix.

        Args:
            n_components: Number of t-SNE dimensions (typically 2)
            perplexity: Perplexity parameter for t-SNE
            n_iter: Number of iterations for optimization
            random_state: Random seed for reproducibility
            phase: Specified phase to analyze. If None, analyzes all structures.

        Returns:
            t-SNE embedding array of shape (n_structures, n_components)
        """
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return np.array([])
        
        feature_mat = self.feature_matrix()
        # Adjust perplexity if needed: must be less than n_samples
        effective_perplexity = min(perplexity, len(structures) - 1)
        tsne_model = TSNE(n_components=n_components, perplexity=effective_perplexity,
                          n_iter=n_iter, random_state=random_state)
        tsne_result = tsne_model.fit_transform(feature_mat)
        return tsne_result

    def umap(self, n_components: int = 2, n_neighbors: int = 15, min_dist: float = 0.1, random_state: int = 42, phase: Optional[str] = None) -> np.ndarray:
        """   
        Performs UMAP dimensionality reduction on the feature matrix.

        Args:
            n_components: Number of UMAP dimensions (typically 2)
            n_neighbors: Number of neighbors for UMAP
            min_dist: Minimum distance parameter for UMAP
            random_state: Random seed for reproducibility
            phase: Specified phase to analyze. If None, analyzes all structures.

        Returns:
            UMAP embedding array of shape (n_structures, n_components)
        """
        if UMAP is None:
            raise ImportError("UMAP is not installed. Install with `pip install umap-learn`.")
        
        if phase:
            structures = self.phases[phase]
        else:
            structures = self.structures
        if not structures:
            self._log(f"No structures found for phase: {phase}", level="warning")
            return np.array([])
        
        feature_mat = self.feature_matrix()
        # Adjust n_neighbors if needed
        effective_neighbors = min(n_neighbors, len(structures) - 1)
        umap_model = UMAP(n_components=n_components, n_neighbors=effective_neighbors,
                          min_dist=min_dist, random_state=random_state)
        umap_result = umap_model.fit_transform(feature_mat)
        return umap_result

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


    def plot_pca_clustered(self, 
                           n_clusters: Optional[int] = None,
                           phase: Optional[str] = None, 
                           cluster_method: str = "agglomerative",
                           **cluster_kwargs):
        """   
        Plots the PCA projection colored by cluster labels.

        Args:
            n_clusters: Number of clusters to form
            phase: Specified phase to analyze. If None, analyzes all structures.
            cluster_method: Clustering method ("agglomerative", "kmeans")
            **cluster_kwargs: Passed to the clustering method
        """
        pca_result, explained_variance = self.pca(phase=phase)


        if cluster_method in ["dbscan", "hdbscan"]:
            labels = self.cluster(method=cluster_method, phase=phase, **cluster_kwargs)
            title_suffix = f"clusters={len(set(labels) - {-1})}, noise={np.sum(labels == -1)}"
        else:
            if n_clusters is None:
                raise ValueError("n_clusters must be specified for agglomerative or kmeans clustering")
            labels = self.cluster(method=cluster_method, n_clusters=n_clusters, phase=phase, **cluster_kwargs)
            title_suffix = f"k={n_clusters}"




        plt.figure(figsize=(10, 6))
        scatter = plt.scatter(pca_result[:, 0], pca_result[:, 1], c=labels, cmap='tab10', alpha=0.7)
        plt.xlabel(f'PC1 ({explained_variance[0]*100:.1f}%)')
        plt.ylabel(f'PC2 ({explained_variance[1]*100:.1f}%)')
        plt.title(f'PCA + {cluster_method.capitalize()} ({title_suffix}) - Phase: {phase if phase else "All"}')
        plt.colorbar(scatter, label='Cluster Label')
        plt.grid(True, alpha=0.3)
        plt.savefig(f"figures/pca_{cluster_method}_{phase if phase else 'all'}.png", dpi=200, bbox_inches='tight')
        plt.close()

    def plot_tsne_clustered(self,
                            n_clusters: Optional[int] = None,
                            phase: Optional[str] = None,
                            cluster_method: str = "agglomerative",
                            n_iter: int = 1000,
                            random_state: int = 42,
                            colored_by: str = "cluster",
                            **cluster_kwargs) -> None:
        """
        Plots t-SNE projection colored by cluster labels or energy.
        4 subplots with perplexities [20, 40, 60, 80].
        """
        if cluster_method in ["dbscan", "hdbscan"]:
            labels = self.cluster(method=cluster_method, phase=phase, **cluster_kwargs)
            title_suffix = f"clusters={len(set(labels) - {-1})}, noise={np.sum(labels == -1)}"
        else:
            if n_clusters is None:
                raise ValueError("n_clusters must be specified for agglomerative or kmeans clustering")
            labels = self.cluster(method=cluster_method, n_clusters=n_clusters, phase=phase, **cluster_kwargs)
            title_suffix = f"k={n_clusters}"

        perplexities = [20, 40, 60, 80]
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        for ax, perplexity in zip(axes.flatten(), perplexities):
            tsne_result = self.tsne(n_components=2, perplexity=perplexity, n_iter=n_iter, random_state=random_state, phase=phase)
            if colored_by == "cluster":
                scatter = ax.scatter(tsne_result[:, 0], tsne_result[:, 1], c=labels, cmap='tab10', alpha=0.7)
                plt.colorbar(scatter, ax=ax, label='Cluster Label')
            elif colored_by == "energy":
                energies = np.array([s.energy for s in (self.phases[phase] if phase else self.structures)])
                scatter = ax.scatter(tsne_result[:, 0], tsne_result[:, 1], c=energies, cmap='viridis', alpha=0.7)
                plt.colorbar(scatter, ax=ax, label='Energy')
            ax.set_title(f't-SNE (perp={perplexity}) - {colored_by}', fontsize=12)
            ax.set_xlabel('t-SNE Dim 1')
            ax.set_ylabel('t-SNE Dim 2')
            ax.grid(True, alpha=0.3)
        plt.suptitle(f"t-SNE + {cluster_method.capitalize()} ({title_suffix}) - Phase: {phase if phase else 'All'}", fontsize=14, fontweight='bold')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(f"figures/tsne_{cluster_method}_{phase if phase else 'all'}_{colored_by}.png", dpi=200, bbox_inches='tight')
        plt.close()

    def plot_umap_clustered(self, 
                            n_clusters: Optional[int] = None,
                            phase: Optional[str] = None,
                            cluster_method: str = "agglomerative",
                            random_state: int = 42,
                            colored_by: str = "cluster",
                            **cluster_kwargs) -> None:
        """
        Plots UMAP projection colored by cluster labels or energy.
        4 subplots with n_neighbors [10, 15, 30, 50].
        """
        if UMAP is None:
            raise ImportError("UMAP is not installed.")
        
        if cluster_method in ["dbscan", "hdbscan"]:
            labels = self.cluster(method=cluster_method, phase=phase, **cluster_kwargs)
            title_suffix = f"clusters={len(set(labels) - {-1})}, noise={np.sum(labels == -1)}"
        else:
            if n_clusters is None:
                raise ValueError("n_clusters must be specified for agglomerative or kmeans clustering")
            labels = self.cluster(method=cluster_method, n_clusters=n_clusters, phase=phase, **cluster_kwargs)
            title_suffix = f"k={n_clusters}"
         

     
        n_neighbors_list = [10, 15, 30, 50]
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        for ax, n_neighbors in zip(axes.flatten(), n_neighbors_list):
            umap_result = self.umap(n_components=2, n_neighbors=n_neighbors, random_state=random_state, phase=phase)
            if colored_by == "cluster":
                scatter = ax.scatter(umap_result[:, 0], umap_result[:, 1], c=labels, cmap='tab10', alpha=0.7)
                plt.colorbar(scatter, ax=ax, label='Cluster Label')
            elif colored_by == "energy":
                energies = np.array([s.energy for s in (self.phases[phase] if phase else self.structures)])
                scatter = ax.scatter(umap_result[:, 0], umap_result[:, 1], c=energies, cmap='viridis', alpha=0.7)
                plt.colorbar(scatter, ax=ax, label='Energy')
            ax.set_title(f'UMAP (n_neighbors={n_neighbors}) - {colored_by}', fontsize=12)
            ax.set_xlabel('UMAP Dim 1')
            ax.set_ylabel('UMAP Dim 2')
            ax.grid(True, alpha=0.3)
        plt.suptitle(f"UMAP + {cluster_method.capitalize()} ({title_suffix}) - Phase: {phase if phase else 'All'}", fontsize=14, fontweight='bold')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(f"figures/umap_{cluster_method}_{phase if phase else 'all'}_{colored_by}.png", dpi=200, bbox_inches='tight')
        plt.close()

    def log_representative_features(
            self,
            representatives: List[StructureData],
            labels: Optional[np.ndarray] = None
        ) -> None:
        """ 
        Log the feature matrix rows for the given representative structures

        Args:
            representatives: List of StructureData objects representing the cluster centers
            labels: Cluster labels array to identify which representative belongs to which cluster. If None, uses cached labels.
        """
        if not self.logger:
            return
        self._ensure_feature_matrix()
        raw = self._feature_matrix_raw
        normalized = self._feature_matrix_normalized

        # Find the indices of the representatives
        rep_indices = []
        for rep in representatives:
            try:
                idx = self.structures.index(rep)
                rep_indices.append(idx)
            except ValueError:
                self._log(f"Representative structure not found in structures list: {rep}", level="warning")
                continue
        if not rep_indices:
            self._log("No valid representative indices found for logging features.", level="warning")
            return
    
        if labels is not None:
            unique_non_noise = [l for l in np.unique(labels) if l != -1]
            row_labels = [f"Cluster {labels[idx]}" for idx in rep_indices]
        else:
            row_labels = [f"Rep {i}" for i in range(len(rep_indices))]
        
        col_labels = [
            "Delta_E", "RMSD", "Rg",
            "Rot_A", "Rot_B", "Rot_C",
            "N_Hbond", "Avg_HB_Ang", "Avg_DA_Dist",
        ]
        n_interatomic = raw.shape[1] - 13
        col_labels += [f"IAD_{j}" for j in range(n_interatomic)]
        col_labels += ["Asphericity", "Acylindricity", "Kappa^2", "Dipole_Mag"]

        rep_matrix_raw = raw[rep_indices]
        rep_matrix_norm = normalized[rep_indices]

        self.logger.log_matrix(
            matrix=rep_matrix_raw,
            row_labels=row_labels,
            col_labels=col_labels,
            title="Representative Structures (Raw)",
            fmt=".4f",
            max_col_width=12
        )
        self.logger.log_matrix(
            matrix=rep_matrix_norm,
            row_labels=row_labels,
            col_labels=col_labels,
            title="Representative Structures (Normalized)",
            fmt=".4f",
            max_col_width=12
        )

