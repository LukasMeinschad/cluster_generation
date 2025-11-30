import numpy as np
from sklearn.cluster import DBSCAN, KMeans, AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler   
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional, Dict, Any



class MolecularCluster:
    """ 
    A class to perform various clustering algorithms on the sampled molecular configurations
    """
    def __init__(self, sampled_molecules: List, energies: List[float]):
        """
        Initializes the MolecularCluster with sampled molecules and their respective energies
        
        Args:
            sampled_molecules (List): List of sampled Molecule objects
            energies (List[float]): List of energies corresponding to the sampled molecules
        """

        self.sampled_molecules = sampled_molecules
        self.energies = energies
        self.feature_matrix = None
        self.labels = None
        self.cluster_centers = None 
        self.cluster_info = None
        self.scaler = StandardScaler()

        

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