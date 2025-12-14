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
