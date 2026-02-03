import numpy as np
from sklearn.cluster import DBSCAN, KMeans, AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler   
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional, Dict, Any
from multiprocessing import Pool, cpu_count
from functools import partial
import time

from symmetry import SymmetryAnalyzer

class USRDescriptor:
    """  
    Class to compute the Ultra-fast Shape Recognition (USR) descriptor for a molecule

    Paper : https://doi.org/10.1098/rspa.2007.1823
    """
    def __init__(self):
        pass 

    def calculate_centroid(self,molecule) -> np.ndarray:
        """  
        Calculates the geometric center of the molecule

        Args:
            molecule: Molecule object with coordinates attribute

        Returns:
            np.ndarray: Centroid coordinates
        """
        return np.mean(molecule.coordinates, axis=0)
    
    def euclidean_distance(self, coord1: np.ndarray, coord2: np.ndarray) -> float:
        """   
        Computes the Euclidean Distance between two 3D coords

        Args:
            coord1 (np.ndarray): First coordinate
            coord2 (np.ndarray): Second coordinate

        Returns:
            float: Euclidean distance
        """
        return np.sqrt(np.sum((coord1 - coord2) ** 2))
    
    def find_closest_and_furthest_to_centroid(self, molecule) -> Tuple[int, int]:
        """  
        Finds the atoms closest and furthest from the centroid.

        Args:
            molecule: Molecule object with coordinates attribute
        Returns:
            Tuple[int, int]: Indices of closest and furthest atoms
        """
        centroid = self.calculate_centroid(molecule)
        distances = []
        for i, coord in enumerate(molecule.coordinates):
            dist = self.euclidean_distance(coord, centroid)
            distances.append((i, dist))
        distances.sort(key=lambda x: x[1])

        # Find closest and furthest
        closest_index = distances[0][0]
        furthest_index = distances[-1][0]
        return closest_index, furthest_index
    
    def find_furthest_from_atom(self, molecule, atom_idx: int) -> int:
        """
        Find the atom furthest from a specified atom

        Args:
            molecule: Molecule object with coordinates attribute
            atom_idx (int): Index of the reference atom
        
        Returns:
            int: Index of the furthest atom
        """
        ref_coords = molecule.coordinates[atom_idx]
        max_distance = 0
        furthest_index = -1
        for i, coord in enumerate(molecule.coordinates):
            dist = self.euclidean_distance(coord, ref_coords)
            if dist > max_distance:
                max_distance = dist
                furthest_index = i
        return furthest_index
    
    def calculate_moments_to_point(self, molecule, reference_point: np.ndarray) -> Tuple[float, float, float]:
        """   
        Calculates the first three moments (mean, variance and skewness) of distances to a reference point

        Args:
            molecule: Molecule object with coordinates attribute
            reference_point (np.ndarray): Reference point coordinates

        Returns:
            Tuple[float, float, float]: First three moments 
        """
        distances = []
        for coord in molecule.coordinates:
            dist = self.euclidean_distance(coord, reference_point)
            distances.append(dist) 

        distances = np.array(distances) # Convert to np array 
        n_atoms = len(distances)

        # First Moment (Mean)
        mean = np.mean(distances)


        # Second Moment (Variance)
        variance = np.sum((distances - mean) ** 2) / n_atoms 

        # Third Moment (Skewness)
        if variance > 0:
            numerator = np.sum((distances - mean) ** 3)
            denominator = n_atoms * (variance ** 1.5) # Take 1.5 power because variance is squared
            skewness = numerator / denominator
        else:
            skewness = 0.0

        return mean, variance, skewness

    def compute_usr_descriptor(self, molecule) -> np.ndarray:
        """   
        Computes the 12-dimensional USR descriptor for the molecule

        This descriptor contains:
        - 3 Moments from centroid 
        - 3 Moments from closest atom to centroid
        - 3 Moments from furthest atom to centroid
        - 3 Moments from furthest atom from the furthest atom to centroid
        """ 
        # Calculate the reference points 
        centroid = self.calculate_centroid(molecule)
        closest_idx, furthest_idx = self.find_closest_and_furthest_to_centroid(molecule)
        furthest_from_furthest_idx = self.find_furthest_from_atom(molecule, furthest_idx)

        # Calculate the moments from each reference point 
        moments_centroid = self.calculate_moments_to_point(molecule, centroid)
        moments_closest = self.calculate_moments_to_point(molecule, molecule.coordinates[closest_idx])
        moments_furthest = self.calculate_moments_to_point(molecule, molecule.coordinates[furthest_idx])
        moments_furthest_from_furthest = self.calculate_moments_to_point(molecule, molecule.coordinates[furthest_from_furthest_idx])

        # Combine all moments into a single descriptor
        usr_descriptor = np.array([
            *moments_centroid,
            *moments_closest,
            *moments_furthest,
            *moments_furthest_from_furthest
        ])  

        return usr_descriptor


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