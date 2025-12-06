import numpy as np
from sklearn.cluster import DBSCAN, KMeans, AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler   
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional, Dict, Any
import time



class MolecularCluster:
    """ 
    A class to perform various clustering algorithms on the sampled molecular configurations
    """
    def __init__(self, sampled_molecules: List, energies: List[float], reference_molecule: Optional["Molecule"], logger = None, sampling_region: Optional[Dict[str, Any]] = None):
        """
        Initializes the MolecularCluster with sampled molecules and their respective energies
        
        Args:
            sampled_molecules (List): List of sampled Molecule objects
            energies (List[float]): List of energies corresponding to the sampled molecules
            logger: Logger object
        """

        self.sampled_molecules = sampled_molecules
        self.energies = energies
        self.feature_matrix = None
        self.labels = None

        self.reference_molecule = reference_molecule

        self.cluster_centers = None 
        self.cluster_info = None
        self.scaler = StandardScaler()
        self.logger = logger

        # Optional for plotting are Sampling region parameters
        self.sampling_region = sampling_region if sampling_region is not None else {}
        
        # Sampling Stats
        self.stats = self.basic_statistics()


        # Hydrogen bond analysis attributed 
        self.sampled_molecules_valid_hbond = []

    def basic_statistics(self) -> Dict[str, Any]:
        """   
        Computes basic statistics of the molecular cluster
        """
        total_molecules = len(self.sampled_molecules)
        avg_energy = np.mean(self.energies)
        min_energy = np.min(self.energies)
        max_energy = np.max(self.energies)
        median_energy = np.median(self.energies)
        std_energy = np.std(self.energies)

        # calculate 25th and 75th percentiles
        percentile_25 = np.percentile(self.energies, 25)
        percentile_75 = np.percentile(self.energies, 75)

        # Calculate interquartile range (IQR)
        iqr = percentile_75 - percentile_25



        stats = {
            "total_molecules": total_molecules,
            "average_energy": avg_energy,
            "minimum_energy": min_energy,
            "maximum_energy": max_energy,
            "median_energy": median_energy,
            "std_energy": std_energy,
            "25_perc_energy": percentile_25,
            "75_perc_energy": percentile_75,
            "iqr_energy": iqr
        }
        if self.logger:
            msg = "Molecular Cluster Basic Statistics:\n"
            for key, value in stats.items():
                msg += f"{key.replace('_', ' ').title()}: {value}\n"
            self.logger.write_message_block(msg)





    def plot_energies_sampling_region(self, diff_to_min: bool = False, log_scale: bool = False,
                                       top_percent: Optional[float] = None,
                                       diff_to_mean: bool = False):
        """   
        Function to plot the energies of the sampled molecules in relation to their sampling region
        """
        if not self.sampling_region:
            raise ValueError("Sampling region information not available.")

        if diff_to_min:
            min_energy = min(self.energies)
            energies = [energy - min_energy for energy in self.energies]
        else:
            energies = self.energies

        if diff_to_mean:
            mean_energy = np.mean(self.energies)
            energies = [energy - mean_energy for energy in self.energies]


        if top_percent is not None:
            if not (0 < top_percent < 100):
                raise ValueError("top_percent must be between 0 and 100")
            n_mols_to_plot = int(len(energies) * (top_percent / 100))
            sorted_indices = np.argsort(energies)
            selected_indices = sorted_indices[:n_mols_to_plot]
            energies = [energies[i] for i in selected_indices]
        
        if self.sampling_region["form"] == "rectangle":
            # For this we plot a 2D scatter plot with the two direction vectors as axes
            dir_vector1 = self.sampling_region["parameters"]["dir_vector1"]
            dir_vector2 = self.sampling_region["parameters"]["dir_vector2"]
            length1 = self.sampling_region["parameters"]["length1"]
            length2 = self.sampling_region["parameters"]["length2"]
            center_point = self.sampling_region["parameters"]["center_point"]

            # We calculate the center of each molecule and the difference to the center point
            coords_list = []
            for mol in self.sampled_molecules:
                coords = mol.coordinates
                mol_center = np.mean(coords, axis=0)
                diff_vector = mol_center - center_point
                coord1 = np.dot(diff_vector, dir_vector1) / np.linalg.norm(dir_vector1)
                coord2 = np.dot(diff_vector, dir_vector2) / np.linalg.norm(dir_vector2)
                coords_list.append((coord1, coord2))
            coords_array = np.array(coords_list)

            if top_percent is not None:
                coords_array = coords_array[selected_indices]

            plt.figure(figsize=(8, 6))
            scatter = plt.scatter(coords_array[:,0], coords_array[:,1], c=energies, cmap='viridis', alpha=0.7)
            if diff_to_min:
                plt.colorbar(scatter, label=f'Energy Difference to Minimum {min_energy:.4f} (Hartree)')
            elif diff_to_mean:
                plt.colorbar(scatter, label=f'Energy Difference to Mean {mean_energy:.4f} (Hartree)')
            else:
                plt.colorbar(scatter, label='Energy (Hartree)')
            plt.title('Energies of Sampled Molecules in Sampling Rectangle')
            plt.xlabel('Projection onto Direction Vector 1')
            plt.ylabel('Projection onto Direction Vector 2')
            
            
            if log_scale:
                plt.yscale('log')
            
            plt.grid(True)
            plt.savefig("energies_sampling_rectangle.png")

    def analyze_h_bond_configurations(self, plot_info: bool = False):
        """
        Function to analyze possible hydrogen bond configurations in the sampled molecules

        For each sampled molecule, the covalent and hydrogen bonds are identified. Then, the hydrogen bond configurations and their angles are calculated and stored.
        """
        # For all the sampled molecules run the get_bonds function
        start_time = time.time()
        covalent_bonds_list = []
        hydrogen_bonds_list = []
        for mol in self.sampled_molecules:
            cov_bonds, h_bonds = mol.get_bonds()
            covalent_bonds_list.append(cov_bonds)
            hydrogen_bonds_list.append(h_bonds) 

        # Get hydrogen bond configurations of all molecules
        # Get angles as well
        hbond_configurations_all = []
        hbond_configurations_angles_all = []
        valid_hbond_configurations_all = []
        valid_hbond_configurations_angles = []
        for mol in self.sampled_molecules:
            hbond_configs = mol.hbond_configurations()
            hbond_configurations_all.append(hbond_configs)
            angles = mol.get_angles_of_hbond_configurations()
            hbond_configurations_angles_all.append(angles)
            valid_configs = mol.get_valid_hydrogen_bond_configurations(angle_threshold=150.0)
            valid_hbond_configurations_all.append(valid_configs)
            valid_angles = mol._get_angle_of_valid_hbond_configurations()
            valid_hbond_configurations_angles.append(valid_angles)


        # Obtain mols with valid hydrogen bond configurations
        filtered_mols = []
        for mol, valid_configs in zip(self.sampled_molecules, valid_hbond_configurations_all):
            if len(valid_configs) > 0:
                filtered_mols.append(mol)
        
        self.sampled_molecules_valid_hbond = filtered_mols
        


        end_time = time.time()
        elapsed_time = end_time - start_time
        


        if self.logger:
            # Make a message of basic statistics
            total_mols = len(self.sampled_molecules)
            avg_hbonds = np.mean([len(hbonds) for hbonds in hydrogen_bonds_list])
            msg = f"Analyzed {total_mols} sampled molecules for hydrogen bond configurations.\n"
            msg += f"Average number of hydrogen bonds per molecule: {avg_hbonds:.2f}\n"
            # Find maximum number of hydrogen bonds in any molecule
            max_hbonds = max([len(hbonds) for hbonds in hydrogen_bonds_list])
            msg += f"Maximum number of hydrogen bonds in a molecule: {max_hbonds}\n"
            msg += f"Time taken for analysis: {elapsed_time:.2f} seconds\n"
            self.logger.write_message_block(msg)


        if plot_info:
            # Plot a histogramm of the number of hydrogen bonds per molecule
            num_hbonds = [len(hbonds) for hbonds in hydrogen_bonds_list]
            plt.figure(figsize=(8, 6))
            plt.hist(num_hbonds, bins=range(0, max(num_hbonds)+2), alpha=0.7, color='blue', edgecolor='black')
            plt.title('Distribution of Hydrogen Bonds per Sampled Molecule')
            plt.xlabel('Number of Hydrogen Bonds')
            plt.ylabel('Frequency')
            plt.grid(axis='y', alpha=0.75)
            plt.savefig("hydrogen_bond_distribution.png")
            plt.close()

            # Plot histogramm of total hydrogen bond configurations angles
            all_angles = [angle for angles in hbond_configurations_angles_all for angle in angles]
            plt.figure(figsize=(8, 6))
            plt.hist(all_angles, bins=30, alpha=0.7, color='green', edgecolor='black')
            plt.title('Distribution of Hydrogen Bond Configuration Angles')
            plt.xlabel('Angle (degrees)')
            plt.ylabel('Frequency')
            plt.grid(axis='y', alpha=0.75)
            plt.savefig("hydrogen_bond_configuration_angles.png")
            plt.close()

            # Plot histogramm of valid hydrogen bond configurations angles
            valid_hbond_configurations_angles = [angle for angles in valid_hbond_configurations_angles for angle in angles]

            plt.figure(figsize=(8, 6))
            plt.hist(valid_hbond_configurations_angles, bins=30, alpha=0.7, color='orange', edgecolor='black')
            plt.title('Distribution of Valid Hydrogen Bond Configuration Angles')
            plt.xlabel('Angle (degrees)')
            plt.ylabel('Frequency')
            plt.grid(axis='y', alpha=0.75)
            plt.savefig("valid_hydrogen_bond_configuration_angles.png")
            plt.close()

    # TODO refactor this plotting function to make it better visually appealing
    def plot_energy_distribution(self, bins: int = 50, density: bool = True, top_percent: Optional[float] = None, diff_to_min: bool = False):
        """
        Function that plots the energy distribution of the sampled molecules
        
        Args:
            bins (int): Number of bins for the histogram
            density (bool): Whether to normalize the histogram
            top_percent (float, optional): If provided, only plots the top percentage of lowest energy molecules
            diff_to_min (bool): If True, plots the energy difference to the minimum energy instead of absolute energies
        """

        if top_percent is not None:
            if not (0 < top_percent < 100):
                raise ValueError("plot_top_percent must be between 0 and 100")
            n_mols_to_plot = int(len(self.energies) * (top_percent / 100))
            sorted_indices = np.argsort(self.energies)
            selected_indices = sorted_indices[:n_mols_to_plot]
            energies_to_plot = [self.energies[i] for i in selected_indices]
        else:
            energies_to_plot = self.energies

        if diff_to_min:
            min_energy = min(energies_to_plot)
            energies_to_plot = [energy - min_energy for energy in energies_to_plot]
        


        plt.figure(figsize=(8, 6))
        plt.hist(energies_to_plot, bins=bins, alpha=0.7, color='cyan', edgecolor='black', density=density)
        plt.title('Energy Distribution of Sampled Molecules')
        plt.xlabel('Energy / (Hartree)')
        plt.ylabel('Frequency')
        plt.grid(axis='y', alpha=0.75)
        if diff_to_min:
            plt.xlabel(f'Energy Difference to Minimum {min_energy:.2f} / (Hartree)')
        plt.tight_layout()
        plt.savefig("energy_distribution_sampled_molecules.png")
        plt.close() 


    def plot_energy_distribution_of_valid_hbond_molecules(self):
        """  
        Function that plots the energy distribution of the sampled molecules that have valid hydrogen bond configurations
        """            
        if not self.sampled_molecules_valid_hbond:
            raise ValueError("No valid hydrogen bond configurations analyzed yet. Please run analyze_h_bond_configurations first.")
        
        # Obtain energies of the valid hbond molecules
        valid_energies = []
        for mol in self.sampled_molecules_valid_hbond:
            idx = self.sampled_molecules.index(mol)
            valid_energies.append(self.energies[idx])
        
        plt.figure(figsize=(8, 6))
        plt.hist(valid_energies, bins=30, alpha=0.7, color='purple', edgecolor='black')
        plt.title('Energy Distribution of Sampled Molecules with Valid Hydrogen Bond Configurations')
        plt.xlabel('Energy (Hartree)')
        plt.ylabel('Frequency')
        plt.grid(axis='y', alpha=0.75)
        plt.savefig("energy_distribution_valid_hbond_molecules.png")
        plt.close()


    def select_n_lowest_energy(self, n: int =5) -> Tuple[List, List[float]]: 
        """ 
        Selects the n lowest energy configurations from the sampled molecules
        """
        if n > len(self.sampled_molecules):
            raise ValueError("n is larger than the number of sampled molecules")
            n = len(self.sampled_molecules)
        
        energy_indices = np.argsort(self.energies)[:n]
        lowest_energy_mols = [self.sampled_molecules[i] for i in energy_indices]
        lowest_energies = [self.energies[i] for i in energy_indices]
        return lowest_energy_mols, lowest_energies

    def calculate_distance_features(self,
                                    submol_atom_labels: List[List[str]],
                                    center_type: str = "centroid",
                                    include_energy: bool = True, 
                                    normalize_features: bool = True):
        """
        Calculates feature matrix based on distances between submolecules

        Args:
            submol_atom_labels (List[List[str]]): List of lists containing atom labels for each submolecule
            center_type (str): Type of center to calculate ('centroid' or 'com')
            include_energy (bool): Whether to include energy as a feature
            normalize_features (bool): Whether to normalize features using StandardScaler

        Motivation of Scaling --> Range of distance and energy can differ significantly, which may bias clustering algorithms
        """
        features = []
        for i, mol in enumerate(self.sampled_molecules):
            coords = mol.coordinates 
            atom_labels = mol.atom_labels
            atom_masses = mol.masses
            mol_features = []

            # Calculate distances between all pairs of submolecules
            n_submols = len(submol_atom_labels)
            for j in range(n_submols):
                for k in range(j + 1, n_submols):
                    # Obtain the indices for the two submolecules
                    indices_submol_j = [idx for idx, label in enumerate(atom_labels) if label in submol_atom_labels[j]]
                    indices_submol_k = [idx for idx, label in enumerate(atom_labels) if label in submol_atom_labels[k]]

                    if indices_submol_j == [] or indices_submol_k == []:
                        raise ValueError(f"Submolecule atom labels not found in molecule: {submol_atom_labels[j]} or {submol_atom_labels[k]}")
                    
                    if center_type == "com":
                        # Obtain the masses for the atoms in the submolecules
                        masses_j = [atom_masses[idx] for idx in indices_submol_j]
                        masses_k = [atom_masses[idx] for idx in indices_submol_k]
                        center_j = self._calculate_center(coords, indices_submol_j, center_type, masses_j)
                        center_k = self._calculate_center(coords, indices_submol_k, center_type, masses_k)
                    else:
                        center_j = self._calculate_center(coords, indices_submol_j, center_type)
                        center_k = self._calculate_center(coords, indices_submol_k, center_type)

                    distance = np.linalg.norm(center_j - center_k)
                    mol_features.append(distance)

            if include_energy:
                mol_features.append(self.energies[i])
            features.append(mol_features)

        self.feature_matrix = np.array(features)                                   
        if normalize_features:
            self.feature_matrix = self.scaler.fit_transform(self.feature_matrix)
         

    def _calculate_center(self, coords: np.ndarray, indices: List[int], center_type: str, masses: Optional[List[float]] = None) -> np.ndarray:
        """ 
        Calculates the center (centroid or COM) of a set of atoms given their indices

        Args:
            coords (np.ndarray): Coordinates of all atoms in the molecule
            indices (List[int]): Indices of the atoms forming the submolecule
            center_type (str): Type of center to calculate ('centroid' or 'com')
            masses (List[float], optional): Masses of the atoms, required for COM calculation
        """
        if center_type == "centroid":
            center = np.mean(coords[indices], axis=0)
        elif center_type == "com":
            if masses is None:
                raise ValueError("Masses must be provided for COM calculation")
            total_mass = sum(masses)
            center = sum(coords[idx] * mass for idx, mass in zip(indices, masses)) / total_mass
        else:
            raise ValueError("Invalid center_type. Choose 'centroid' or 'com'.")
        return center


    def calculate_rmsd_features(self, reference_molecule = None, include_energy: bool = True, normalize_features: bool = True):
        """ 
        Calculates a feature matrix based on RMSD to a reference molecule

        If no reference molecule is provided the lowest energy molecule is used as reference

        Args:
            reference_molecule: The molecule to calculate RMSD against
            include_energy (bool): Whether to include energy as a feature
            normalize_features (bool): Whether to normalize features using StandardScaler
        
        The rmst is given by:
            RMSD = sqrt( (1/N) * sum_i^N ||r_i - r_i_ref||^2)
        """
        if reference_molecule is None:
            reference_molecule = self.sampled_molecules[np.argmin(self.energies)]

        ref_coords = reference_molecule.coordinates
        features = []
        for i, mol in enumerate(self.sampled_molecules):
            coords = mol.coordinates 
            diff = coords - ref_coords
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
            mol_features = [rmsd]
            if include_energy:
                mol_features.append(self.energies[i])
            features.append(mol_features)
        self.feature_matrix = np.array(features)

        if include_energy:
            self.feature_matrix = np.array(features)

        if normalize_features:
            self.feature_matrix = self.scaler.fit_transform(self.feature_matrix)

    def calculate_pairwise_distance_features(self, include_energy: bool = True, normalize_features: bool = True):
        """
        Calculates features based on pairwise atomic distances
        """
        features = []
        for mol in self.sampled_molecules:
            coords = mol.coordinates
            # Calculate pairwise distance matrix
            dist_matrix = np.linalg.norm(coords[:, None, :] - coords[None,: ,:], axis=-1)
            # Extract upper triangle distances to avoid redundancy
            mol_features = []

            upper_triangle = dist_matrix[np.triu_indices(len(coords), k=1)]

            mol_features.extend([
                np.mean(upper_triangle), # Average distance
                np.std(upper_triangle),  # Standard deviation of distances
                np.min(upper_triangle),  # Minimum distance
                np.max(upper_triangle)   # Maximum distance
            ])

            if include_energy:
                mol_features.append(self.energies[i])
            features.append(mol_features)
        self.feature_matrix = np.array(features)

        if normalize_features:
            self.feature_matrix = self.scaler.fit_transform(self.feature_matrix)

    def cluster_kmeans(self, n_cluster: int=4) -> Dict[str, Any]:
        """ 
        Performs KMeans Clustering on the feature matrix
        Args:
            n_cluster (int): Number of clusters for KMeans
        """
        if self.feature_matrix is None:
            raise ValueError("Feature matrix not calculated. Please run calculate_distance_features or calculate_rmsd_features first.")
        
        kmeans = KMeans(n_clusters=n_cluster, random_state=42)
        self.labels = kmeans.fit_predict(self.feature_matrix)
        self.cluster_centers = kmeans.cluster_centers_

        silhouette_avg = silhouette_score(self.feature_matrix, self.labels)

        self.cluster_info = {
            "method": "KMeans",
            "n_clusters": n_cluster,
            "silhouette_score": silhouette_avg
        }
        return self.cluster_info
    
    def cluster_kmeans_silhouette_analysis(self, max_clusters: int=10, plot_analysis: bool = False) -> List[Dict[str, Any]]:
        """ 
        Performs KMeans Clustering for a range of cluster numbers and evaluates this using silhouette score
        
        Args:
            max_clusters (int): Maximum number of clusters to evaluate
        """
        if self.feature_matrix is None:
            raise ValueError("Feature matrix not calculated. Please run calculate_distance_features or calculate_rmsd_features first.")
        
        cluster_analysis = []
        for n_cluster in range(2, max_clusters + 1):
            kmeans = KMeans(n_clusters=n_cluster, random_state=42)
            labels = kmeans.fit_predict(self.feature_matrix)
            silhouette_avg = silhouette_score(self.feature_matrix, labels)
            cluster_info = {
                "n_clusters": n_cluster,
                "silhouette_score": silhouette_avg
            }
            cluster_analysis.append(cluster_info)
        
        # Obtain the best clustering based on silhouette score
        best_clustering = max(cluster_analysis, key=lambda x: x["silhouette_score"])
        self.labels = KMeans(n_clusters=best_clustering["n_clusters"], random_state=42).fit_predict(self.feature_matrix)
        self.cluster_centers = KMeans(n_clusters=best_clustering["n_clusters"], random_state=42).fit(self.feature_matrix).cluster_centers_

        print(f"Best clustering: {best_clustering['n_clusters']} clusters with silhouette score {best_clustering['silhouette_score']:.4f}")
        
        if plot_analysis:
            n_clusters_list = [info["n_clusters"] for info in cluster_analysis]
            silhouette_scores = [info["silhouette_score"] for info in cluster_analysis]

            plt.figure(figsize=(8, 6))
            plt.plot(n_clusters_list, silhouette_scores, marker='o')
            plt.title('KMeans Silhouette Analysis')
            plt.xlabel('Number of Clusters')
            plt.ylabel('Silhouette Score')
            plt.xticks(n_clusters_list)
            plt.grid(True)
            plt.savefig("kmeans_silhouette_analysis.png")
            plt.close()
        
        return cluster_analysis


    def obtain_cluster_representatives(self, method: str ="closest_to_center") -> List:
        """ 
        Obtains the representative molecule for each cluster based on the specified method
        """ 
        if self.labels is None:
            raise ValueError("Clustering not performed yet. Please run a clustering method first.")
        
        representatives = []
        n_clusters = len(set(self.labels)) # Number of unique clusters
        for cluster_id in range(n_clusters):
            cluster_indices = np.where(self.labels == cluster_id)[0]
            cluster_molecules = [self.sampled_molecules[i] for i in cluster_indices]
            cluster_features = self.feature_matrix[cluster_indices]

            if method == "closest_to_center":
                center = self.cluster_centers[cluster_id]
                distances_to_center = np.linalg.norm(cluster_features - center, axis=1)
                closest_index = cluster_indices[np.argmin(distances_to_center)]
                representatives.append(self.sampled_molecules[closest_index])
            if method == "lowest_energy":
                cluster_energies = [self.energies[i] for i in cluster_indices]
                lowest_energy_index = cluster_indices[np.argmin(cluster_energies)]
                representatives.append(self.sampled_molecules[lowest_energy_index])


        return representatives
    
    def cluster_agglomerative(self,n_clusters: int=4, metric: str = "euclidean", linkage: str = "ward") -> Dict[str, Any]:
        """ 
        Performs Agglomerative Clustering on the feature matrix
        
        Args:
            n_clusters (int): Number of clusters for Agglomerative Clustering
            metric (str): Distance metric to use
            linkage (str): Linkage criterion to use

        + Ward linkage minimizes the variance of the clusters being merged.
        + Average uses the average of the distances of each observation of the two sets
        + Complete or "maximum" linkage uses the maximum distances between all observations of the two sets.
        + "single" uses the minmum distances between all observations of the two sets
        """
        if self.feature_matrix is None:
            raise ValueError("Feature matrix not calculated. Please run calculate_distance_features or calculate_rmsd_features first.")
        
        agglomerative = AgglomerativeClustering(n_clusters=n_clusters, affinity=metric, linkage=linkage)
        self.labels = agglomerative.fit_predict(self.feature_matrix)

        self.cluster_info = {
            "method": "Agglomerative",
            "n_clusters": n_clusters,
            "metric": metric,
            "linkage": linkage
        }
        return self.cluster_info
    
    def plot_dendrogram(self, method: str = "ward"):
        """ 
        Plots the dendrogram for the hierachical cluster method
        
        Args:
            method (str): Linkage method to use for dendrogram
        """
        if self.feature_matrix is None:
            raise ValueError("Feature matrix not calculated. Please run calculate_distance_features or calculate_rmsd_features first.")
        
       
        linked = linkage(self.feature_matrix, method=method)

        plt.figure(figsize=(10, 7))
        dendrogram(linked,
                   orientation='top',
                   distance_sort='descending',
                   show_leaf_counts=True)
        plt.title('Hierarchical Clustering Dendrogram')
        plt.xlabel('Sample Index')
        plt.ylabel('Distance')
        plt.savefig("dendrogram.png")
        plt.close()


    def visalize_cluster_2d(self):
        """ 
        Visualizes the clustering results in 2D using the first two features
        """
        if self.feature_matrix is None or self.labels is None:
            raise ValueError("Feature matrix or labels not available. Please run clustering first.")
        
        plt.figure(figsize=(8, 6))
        scatter = plt.scatter(self.feature_matrix[:, 0], self.feature_matrix[:, 1], c=self.labels, cmap='viridis', alpha=0.7)
        plt.title('Molecular Clustering Results')
        plt.xlabel('Feature 1')
        plt.ylabel('Feature 2')
        plt.colorbar(scatter, label='Cluster Label')
        plt.grid(True)
        plt.legend(*scatter.legend_elements(), title="Clusters")
        plt.savefig("clustering_results_2d.png")
        plt.close()