import numpy as np
from sklearn.cluster import DBSCAN, KMeans

class MolecularCluster:
    """ 
    A class to perform various clustering algorithms on the sampled molecular configurations
    """
    def __init__(self, sampled_molecules, energies):
        """
        Initializes the MolecularCluster with sampled molecules and their respective energies
        """
        self.sampled_molecules = sampled_molecules
        self.energies = energies
        self.feature_matrix = None
        self.labels = None
        self.cluster_centers = None

    def calculate_com_features(self, submol_atom_labels, center_type="centroid"):
        """ 
        Calculate a feature matrix based on the COM/geometric distances and energies
        Uses the submol_atom_labels to identify the submolecule within each sampled molecule

        At the moment this function works only for two submolecules in the system
        """
        features = []
        submol_atom_labels = [set(labels) for labels in submol_atom_labels]
        print(submol_atom_labels)
        for mol in self.sampled_molecules:
            coords = mol.coordinates
            atom_labels = mol.atom_labels

            # Identify the indices for the submolecule
            indices_submol1 = [i for i, label in enumerate(atom_labels) if label in submol_atom_labels[0]]
            indices_submol2 = [i for i, label in enumerate(atom_labels) if label in submol_atom_labels[1]]

            # TODO Implement com calculation here as well, do this by using the molecule class --> get_com function
            if center_type == "centroid":
                center1 = np.mean(coords[indices_submol1], axis=0)
                center2 = np.mean(coords[indices_submol2], axis=0)
            else:
                raise NotImplementedError("Only centroid center_type is implemented")
            
            distance = np.linalg.norm(center1 - center2)
            # Create a feature vector: [distance, energy]
            feature_vector = [distance, self.energies[self.sampled_molecules.index(mol)]]
            features.append(feature_vector)
        self.feature_matrix = np.array(features)
        
    def cluster_kmeans(self,n_clusters = 4):
        """ 
        Perform KMeans clustering on the feature matrix
        
        Args:
            n_clusters: Number of clusters for KMeans
        """
        if self.feature_matrix is None:
            raise ValueError("Feature matrix not calculated. Run calculate_com_features() first.")
        
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        self.labels = kmeans.fit_predict(self.feature_matrix)
        self.cluster_centers = kmeans.cluster_centers_
        print(f"KMeans clustering completed with {n_clusters} clusters.")

    def visualize_cluster_com(self):
        """
        Visualizes the clusters based on the COM distance and energy features
        """
        if self.feature_matrix is None or self.labels is None:
            raise ValueError("Feature matrix or labels not available. Ensure clustering has been performed.")
        
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8,6))
        scatter = plt.scatter(self.feature_matrix[:,0], self.feature_matrix[:,1], c=self.labels, cmap='viridis', alpha=0.7)
        plt.xlabel("COM Distance (Å)")
        plt.ylabel("Energy (Hartree)")
        plt.title("KMeans Clustering of Molecular Configurations")
        plt.colorbar(scatter, label='Cluster Label')
        plt.grid()
        plt.savefig("kmeans_cluster_com.png")
        plt.show()
        