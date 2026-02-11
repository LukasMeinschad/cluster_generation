"""  
Module for various clustering methods and analysis of the sampled molecular configurations
"""
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional, Dict, Any, Union
from dataclasses import dataclass, field
from collections import defaultdict
import seaborn as sns
from sklearn.cluster import DBSCAN, KMeans, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA


from molecule_class import Molecule
from transformations import GeometryOps


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
    def __init__(self, name: str = "BHMC_Analysis"):
        self.name = name 
        self.structures: List[StructureData] = []
        self.phases: Dict[str, List[StructureData]] = defaultdict(list)

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
            structures = self.phase[phase]
        else:
            structures = self.structures
        n = len(structures)
        rmsd_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                rmsd = self._calculate_rmsd(structures[i].molecule.coordinates, structures[j].molecule.coordinates)
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd  # Symmetric matrix
        return rmsd_matrix


    @staticmethod 
    def _calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
        """   
        Calculate RMSD between two sets of coordinates
        """
        # Center the coordinates
        coords1_centered = coords1 - np.mean(coords1, axis=0)
        coords2_centered = coords2 - np.mean(coords2, axis=0)
        diff = coords1_centered - coords2_centered
        rsmd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        return rsmd

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
        I = GeometryOps.inertia_tenstor(coords_centered, masses)
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

        
    


    # ====================== Clustering Methods =======================

    def feature_matrix(self, normalize = True) -> np.ndarray:
        """   
        Current Feature Matrix wifth Columns

        1. Delta E to lowest energy structure
        2. RMSD to lowest energy structure
        3. Radius of Gyration
        4. Rotational Constants (A,B,C)
        """
        delta_e = np.array([s.energy - self.get_lowest_energy_structure().energy for s in self.structures])
        rmsd_values = np.array([self._calculate_rmsd(s.molecule.coordinates, self.get_lowest_energy_structure().molecule.coordinates) for s in self.structures])
        rg_values = np.array([self.radius_of_gyration(s.molecule.coordinates, s.molecule.masses) for s in self.structures])
        rotational_constants = np.array([self.determine_rotational_constants(s.molecule.coordinates, s.molecule.masses) for s in self.structures])
        feature_matrix = np.column_stack((delta_e, rmsd_values, rg_values, rotational_constants))

        # Normalize using StandardScaler
        scaler = StandardScaler()
        if normalize:
            feature_matrix = scaler.fit_transform(feature_matrix)


        return feature_matrix
    
    def AgglomerativeClustering(self, n_clusters: int = 5, phase: Optional[str] = None, linkage: str = "ward") -> np.ndarray:
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
        clustering_model = AgglomerativeClustering(n_clusters=n_clusters, linkage=linkage)
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
        labels = self.AgglomerativeClustering(n_clusters=n_clusters, phase=phase, linkage=linkage)
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
            
     
    






class MolecularCluster:
    """  
    A class to perform clustering analysis on sampled molecular configurations
    """

    def __init__(
            self,
            sampled_molecules: List,
            energies: List[float],
            reference_molecule: Optional["Molecule"] = None,
            logger = None, 
            sampling_region: Optional[Dict[str, Any]] = None
        ):
        """
        Initializes the MolecularCLuster with the sampled molecules and their energies

        Args:
            sampled_molecules (List): List of sampled Molecule objects
            energies (List[float]): List of energies corresponding to the sampled molecules
            logger: Logger object
        """
        # Core Data 
        self.sampled_molecules = sampled_molecules 
        self.energies = np.array(energies)
        self.reference_molecule = reference_molecule

        # Clustering Attributes
        self.feature_matrix = None 
        self.feature_names = None
        self.labels = None 
        self.cluster_centers = None 
        self.cluster_info = None 
        self.scaler = StandardScaler()

        # Optional components 
        self.logger = logger
        self.sampling_region = sampling_region if sampling_region is not None else {}

        # Computed properteis
        self.stats = self._compute_statistics()

        # Hydrogen Bond Analysis (Delayed initialization)
        self.hbond_statistics = None
        self.hbond_data = None

        # Symmetry Analysis Attributes
        self.total_symmetry_elements = None


    def _compute_statistics(self) -> Dict[str, float]:
        """  
        Computes the basic statistics of the molecule energies
        """ 
        max_energy = np.max(self.energies)
        min_energy = np.min(self.energies)
        mean_energy = np.mean(self.energies)
        std_energy = np.std(self.energies)
        median_energy = np.median(self.energies)

        # Compute the IQR
        q75, q25 = np.percentile(self.energies, [75 ,25])
        iqr_energy = q75 - q25
        # Compute Mean of q75 and q25
        mid_iqr_energy = (q75 + q25) / 2
        
        statistics = {
            "max_energy": max_energy,
            "min_energy": min_energy,
            "mean_energy": mean_energy,
            "std_energy": std_energy,
            "median_energy": median_energy,
            "iqr_energy": iqr_energy,
            "mid_iqr_energy": mid_iqr_energy
        }

        if self.logger:
            self.logger.write_cluster_statistics(statistics)

        return statistics
    
    # =================== Symmetry Analysis =======================

    @staticmethod 
    def _analyze_single_molecule_symmetry(mol, tolerance: float = 1e-3) -> int: 
        """   
        Helper function to analyze a single molecule's symmetry
        Must be static to work with multiprocessing
        """
        Sym_Analyzer = SymmetryAnalyzer(molecule=mol)
        symmetry_info = Sym_Analyzer.analyze_full_symmetry(tolerance=tolerance)
        return symmetry_info.get('total_symmetry_elements', 0)

    def get_total_symmetry_elements_samples(
            self,
            tolerance: float = 1e-3,
            n_processes: Optional[int] = None,
            use_multiprocessing: bool = True
        ) -> List[int]:
        """   
        Computes the total number of symmetry elements for each sampled molecule 

        Args:
            tolerance (float): Tolerance for symmetry analysis
            n_processes (Optional[int]): Number of processes to use for multiprocessing
            use_multiprocessing (bool): Whether to use multiprocessing
        """
        start_time = time.time()
        if use_multiprocessing and len(self.sampled_molecules) > 1:
            # Use multiprocessing for large datasets 
            if n_processes is None:
                n_processes = max(1, cpu_count() - 1)
            
            # Create partial function with fixed tolerance
            analyze_func = partial(self._analyze_single_molecule_symmetry, tolerance=tolerance)

            # Pool 
            with Pool(processes=n_processes) as pool:
                total_symmetry_elements = pool.map(analyze_func, self.sampled_molecules)
        else:
            # Sequential processing
            total_symmetry_elements = [] 
            for mol in self.sampled_molecules:
                Sym_Analyzer = SymmetryAnalyzer(molecule=mol)
                symmetry_info = Sym_Analyzer.analyze_full_symmetry(tolerance=tolerance)
                total_symmetry_elements.append(symmetry_info.get('total_symmetry_elements', 0))
        self.total_symmetry_elements = total_symmetry_elements
        end_time = time.time()
        print(f"Symmetry elements analysis completed in {end_time - start_time:.2f} seconds.")

    
    # =================== Hydrogen Bond Analysis ===================== 
    
    @staticmethod
    def _analyze_single_molecule(mol, angle_threshold: float = 150.0) -> List:
        """   
        Helper function to analyze the single molecule's Hydrogen bond
        Must be static to work with multiprocessing
        """
        return mol.get_valid_hbond_configurations(angle_threshold=angle_threshold)
    
    def analyze_hydrogen_bonds(
            self, 
            angle_threshold: float = 150.0,
            n_processes: Optional[int] = None,
            use_multiprocessing: bool = True
        ) -> Tuple[Dict[str, float], List]:
        """  
        Analyzes the hydrogen bond configurations across all sampled molecules

        Args:
            angle_threshold (float): Angle threshold for valid H-bond configurations
            n_processes (Optional[int]): Number of processes to use for multiprocessing
            use_multiprocessing (bool): Whether to use multiprocessing
        Returns:
            Tuple containing H-bond statistics and detailed H-bond data
        """
        time_start = time.time()
        
        if use_multiprocessing and len(self.sampled_molecules) > 1: 
            # Use multiprocessing for large datasets 
            if n_processes is None:
                n_processes = max(1, cpu_count() - 1)
            
            # Create partial function with fixed angle_threshold
            analyze_func = partial(self._analyze_single_molecule, angle_threshold=angle_threshold)

            # Use multiprocessing Pool
            with Pool(processes=n_processes) as pool:
                hbond_data = pool.map(analyze_func, self.sampled_molecules)
        else:
            hbond_data = [
                mol.get_valid_hbond_configurations(angle_threshold=angle_threshold) 
                for mol in self.sampled_molecules
            ]
        
        # Compute statistics
        hbond_counts = [len(hbonds) for hbonds in hbond_data]
        max_hbonds = np.max(hbond_counts)
        min_hbonds = np.min(hbond_counts)
        mean_hbonds = np.mean(hbond_counts)
        std_hbonds = np.std(hbond_counts)
        median_hbonds = np.median(hbond_counts)
        hbond_statistics = {
            "max_hbonds": max_hbonds,
            "min_hbonds": min_hbonds,
            "mean_hbonds": mean_hbonds,
            "std_hbonds": std_hbonds,
            "median_hbonds": median_hbonds
        }
        self.hbond_statistics = hbond_statistics
        self.hbond_data = hbond_data
        time_end = time.time()
        print(f"H-bond analysis completed in {time_end - time_start:.2f} seconds.")
        if self.logger:
            self.logger.write_hbond_analysis(hbond_statistics)
        return hbond_statistics, hbond_data


    # =================== Feature Calculation ======================
    """  
    In this section we set up a feature matrix for the later machine learning

    1. RMSD to reference molecule
    """
    def compute_feature_matrix(
            self, 
            include_energy: bool = True,
            include_rmsd: bool = True,
            include_hbonds: bool = True,
            inculde_radius_of_gyration: bool = True,
            include_symmetry: bool = True,
            normalize: bool = True
        ) -> np.ndarray:
        """
        Computes the feature matrix for clustering analysis

        Args:
            include_energy (bool): Whether to include energy as a feature
            include_rmsd (bool): Whether to include RMSD to reference molecule as a feature
            include_hbonds (bool): Whether to include number of H-bonds as a feature
            include_symmetry (bool): Whether to include total symmetry elements as a feature
        Returns:
            np.ndarray: Feature matrix
        """
        features_list = []
        features_names = []

        n_samples = len(self.sampled_molecules)

        # 1. Energy Feature
        if include_energy:
            energy_feature = self.energies.reshape(-1,1)  # Transforms into column vector
            features_list.append(energy_feature)
            features_names.append("energy")
            
        if include_rmsd and self.reference_molecule is not None: 
            rmsd_values = self.compute_rmsd_to_reference() 
            features_list.append(rmsd_values.reshape(-1,1))
            features_names.append("rmsd")

        if include_hbonds and self.hbond_statistics is not None:
            hbond_counts = np.array([len(hbonds) for hbonds in self.hbond_data]).reshape(-1,1)
            features_list.append(hbond_counts)
            features_names.append("hbond_counts")

        if inculde_radius_of_gyration:
            rg_values = self.compute_radius_of_gyration()
            features_list.append(rg_values.reshape(-1,1))
            features_names.append("radius_of_gyration")

        

        if not features_list:
            raise ValueError("No features selected for feature matrix computation.")
        
        feature_matrix = np.hstack(features_list)
        print(f"Total feature matrix shape: {feature_matrix.shape}")
        print(f"Included features: {features_names}")
        if normalize:
            feature_matrix = self.scaler.fit_transform(feature_matrix)
        
        self.feature_matrix = feature_matrix
        self.feature_names = features_names

        

        
    
    
    def compute_rmsd_to_reference(self) -> np.ndarray:
        """  
        Computes the RMSD of all sampled molecules to the reference molecule
        """
        if self.reference_molecule is None:
            raise ValueError("Reference molecule not set for RMSD calculation.")
        
        rmsd_values = []
        ref_coords = self.reference_molecule.coordinates
        
        for mol in self.sampled_molecules:
            mol_coords = mol.coordinates
            if mol_coords.shape != ref_coords.shape:
                raise ValueError("Molecule and reference molecule have different number of atoms.")
            
            # Calculate RMSD
            diff = mol_coords - ref_coords
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
            rmsd_values.append(rmsd)
        
        return np.array(rmsd_values)
    
    def compute_radius_of_gyration(self) -> np.ndarray:
        """
        Computes the radius of gyration (R_g) which is the root-mean-square distance of each part of the molecule from
        its center of mass.

        R_g = sqrt( (1/N) * sum((r_i - r_cm)^2) )

        where r_i is the position vector of atom i, r_cm is the center of mass position vector, and N is the number of atoms.
        """
        rg_values = []
        for mol in self.sampled_molecules:
            coords = mol.coordinates
            masses = mol.masses
            com = np.sum(coords.T * masses, axis=1) / np.sum(masses)
            N = coords.shape[0]
            diff = coords - com
            rg = np.sqrt(np.sum(masses * np.sum(diff**2, axis=1)) / np.sum(masses))
            rg_values.append(rg)
        return np.array(rg_values) 
    
    # =================== USR Descriptor Calculation ======================
    def compute_usr_descriptors(self, use_multiprocessing: bool = True, n_processes: Optional[int] = None) -> np.ndarray:
        """   
        Computes USR descriptors for all sampled molecules 

        Args:
            use_multiprocessing (bool): Whether to use multiprocessing
            n_processes (Optional[int]): Number of processes to use for multiprocessing

        Returns:
            np.ndarray: USR descriptor matrix
        """
        start_time = time.time()

        if use_multiprocessing and len(self.sampled_molecules) > 1: 
            if n_processes is None: 
                n_processes = max(1, cpu_count() -1) 
            
            with Pool(processes = n_processes) as pool:
                usr_descriptors = pool.map(self._compute_single_usr, self.sampled_molecules)
        else:
            # Sequential processing
            usr_descriptors = [self._compute_single_usr(mol) for mol in self.sampled_molecules]

        end_time = time.time()
        print(f"USR descriptor computation completed in {end_time - start_time:.2f} seconds.")
        return np.array(usr_descriptors)
    
    @staticmethod 
    def _compute_single_usr(molecule) -> np.ndarray:
        """    
        Static method to compute the USR descriptor for a single molecule 
        Required form multiprocessing        
        
        Args:
            molecule: Molecule object
        Returns:
            np.ndarray: USR descriptor
        """
        usr_calc = USRDescriptor()
        return usr_calc.compute_usr_descriptor(molecule)
    
    def compute_usr_similarity_matrix(self, usr_descriptors: np.ndarray = None) -> np.ndarray:
        """   
        Computes a pairwise similarity matrix based on the USR descritpros
        Uses the Manhatten Distance as in the original USR paper 

        Args:
            usr_descriptors (np.ndarray): Precomputed USR descriptors. If None, computes them. 
        Returns:
            np.ndarray: Pairwise similarity matrix
        """
        if usr_descriptors is None:
            usr_descriptors = self.compute_usr_descriptors()

        n_molecules = len(usr_descriptors)
        similarity_matrix = np.zeros((n_molecules, n_molecules))

        for i in range(n_molecules):
            for j in range(i + 1, n_molecules):
                # Use Manhattan Distance
                distance = np.sum(np.abs(usr_descriptors[i] - usr_descriptors[j]))
                similarity_matrix[i, j] = distance
                similarity_matrix[j, i] = distance  # Symmetric matrix
        return similarity_matrix
    
    


    # ==================== Clustering Methods ======================



    
    # ================= Plotting ====================

    def plot_energy_distribution(self, bins: int = 100):
        """  
        Plots the energy distribution of the sampled molecules

        
        Makes one plot of the total energy distribution and subsequently a 4 subplot of the energy in the 
        25% quantile ranges.
        """ 
        # Overall Energy Distribution
        plt.figure(figsize=(10, 6))
        plt.hist(self.energies, bins=bins, color='blue', alpha=0.7)
        plt.title('Energy Distribution of Sampled Molecules')
        plt.xlabel('Energy')
        plt.ylabel('Frequency')
        plt.grid(True)
        plt.savefig("figures/energy_distribution_overall.png")
        plt.close()
        # Make 4 subplots of the energy in the following ranges
        # 0 - 25% quantile
        # 25% - 50% quantile
        # 50% - 75% quantile
        # 75% - 100% quantile
        q25 = np.percentile(self.energies, 25)
        q50 = np.percentile(self.energies, 50)
        q75 = np.percentile(self.energies, 75)
        ranges = [(self.energies <= q25, f'0 - 25% Quantile (<= {q25:.2f})'),
                  ((self.energies > q25) & (self.energies <= q50), f'25% - 50% Quantile ({q25:.2f} - {q50:.2f})'),
                  ((self.energies > q50) & (self.energies <= q75), f'50% - 75% Quantile ({q50:.2f} - {q75:.2f})'),
                  (self.energies > q75, f'75% - 100% Quantile (> {q75:.2f})')]
        plt.figure(figsize=(12, 10))
        for i, (condition, title) in enumerate(ranges, 1):
            plt.subplot(2, 2, i)
            plt.hist(self.energies[condition], bins=bins, color='green', alpha=0.7)
            plt.title(f'Energy Distribution: {title}')
            plt.xlabel('Energy')
            plt.ylabel('Frequency')
            plt.grid(True)
        plt.tight_layout()
        plt.savefig("figures/energy_distribution_quantiles.png")
        plt.close()

    def plot_hydrogen_bond_statistics(self):
        """
        Plots a histogram of the hydrogen bond statistics 
        """
        if self.hbond_statistics is None:
            raise ValueError("Hydrogen bond statistics not computed. Please run analyze_hydrogen_bonds() first.")
        
        # Extract data
        hbond_counts = [len(hbonds) for hbonds in self.hbond_data]
        
        plt.figure(figsize=(10, 6))
        plt.hist(hbond_counts, bins=range(0, max(hbond_counts)+2), color='purple', alpha=0.7, align='left')
        plt.title('Hydrogen Bond Count Distribution')
        plt.xlabel('Number of Hydrogen Bonds')
        plt.ylabel('Frequency')
        plt.grid(True)
        plt.savefig("figures/hydrogen_bond_distribution.png")
        plt.close()

    def plot_symmetry_elements_distribution(self):
        """
        Plots a histogram of the total symmetry elements across sampled molecules
        """
        if self.total_symmetry_elements is None:
            raise ValueError("Total symmetry elements not computed. Please run get_total_symmetry_elements_samples() first.")
        
        plt.figure(figsize=(10, 6))
        plt.hist(self.total_symmetry_elements, bins=range(0, max(self.total_symmetry_elements)+2), color='orange', alpha=0.7, align='left')
        plt.title('Total Symmetry Elements Distribution')
        plt.xlabel('Number of Symmetry Elements')
        plt.ylabel('Frequency')
        plt.grid(True)
        plt.savefig("figures/symmetry_elements_distribution.png")
        plt.close()

    def plot_usr_similarity_heatmap(self, usr_descriptors: np.ndarray = None):
        """   
        Plots a heatmap of the USR similarity matrix

        Args:
            usr_descriptors (np.ndarray): Precomputed USR descriptors. If None, computes them. 
        """
        similarity_matrix = self.compute_usr_similarity_matrix(usr_descriptors=usr_descriptors)

        plt.figure(figsize=(10, 8))
        plt.imshow(similarity_matrix, cmap='hot', interpolation='nearest')
        plt.colorbar(label='USR Similarity (Manhattan Distance)')
        plt.title('USR Similarity Heatmap')
        plt.xlabel('Molecule Index')
        plt.ylabel('Molecule Index')
        plt.savefig("figures/usr_similarity_heatmap.png")
        plt.close()