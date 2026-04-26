"""
Integration of the DScribe Package used as a featurizer for subsequent machine learning.
DScribe provides methods to transform atomic structures into fixed-size numeric vectors. These vectords are
built in a way that tehay efficiently summarize the contents of the input structures. Such a transformation is very useful for machine learning applications
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Optional, Any
from dataclasses import dataclass

# Import Molecule class 
from modules.molecule_class import Molecule


# Import DScribe descriptors
from dscribe.descriptors import SOAP, ACSF, MBTR

@dataclass
class FeaturizerConfig:
    descriptor_type: str = "SOAP"
    soap_rcut: float = 6 # Cutoff for local region
    soap_nmax: int = 12 # Number of radial basis functions
    soap_lmax: int = 8 # Maximum degree of spherical harmonics
    soap_sigma: float = 0.5 # Std of the Gaussians used to expand the atomic density
    soap_rbf: str = "gto" # Type of radial basis functions (e.g. "gto", "polynomial")
    soap_average: str = "inner" # Averaging method for SOAP features ("inner", "outer", "off")





class Featurizer:
    local_descriptors = ["ACE", "ACSF", "SOAP"]
    global_descriptors = ["MBTR"]
    def __init__(self, config: FeaturizerConfig):
        """ 
        Initialize the featurizer with a specific DScribe descriptor object
        """
        self.config = config

    def _initialize_descriptor(self, species: List[str], **kwargs) -> Any:
        """  
        Internal method to initialize the DScribe descriptor based on the provided chemical species

        We group into local and global descriptors --> For this see paper
        https://www.nature.com/articles/s41524-022-00721-x
        """ 
        # get the descripor type from the config
        desc_type = self.config.descriptor_type.upper()
        if desc_type == "SOAP":
            return SOAP(
                species=species,
                r_cut=self.config.soap_rcut,
                n_max=self.config.soap_nmax,
                l_max=self.config.soap_lmax,
                sigma=self.config.soap_sigma,
                rbf=self.config.soap_rbf,
                average=self.config.soap_average
            )
        elif desc_type == "ACSF":
            # Initialize ACSF descriptor with appropriate parameters
            pass
        elif desc_type == "MBTR":
            return MBTR(
                geometry=self.config.mbtr_geometry,
                grid=self.config.mbtr_grid,
                weighting=self.config.mtbr_weighting,
                species=species
            )
        else:
            raise ValueError(f"Unsupported descriptor type: {desc_type}")
        
    def compute_soap_features(self, molecules: List[Molecule]) -> np.ndarray:
        """ 
        Compute SOAP features for a list of Molecule objects
        """
        if self.config.descriptor_type.upper() != "SOAP":
            raise ValueError("Descriptor type must be 'SOAP' to compute SOAP features.")
        
        # Extract unique chemical species from the molecules
        species = set()
        for mol in molecules:
            # transform to ASE atoms
            atoms = mol.to_ase_atoms()
            species.update(atoms.get_chemical_symbols())

        species = sorted(species)  # Sort species for consistent ordering
        # Initialize the SOAP descriptor
        soap_descriptor = self._initialize_descriptor(species=species)

        # Compute SOAP features for each molecule and aggregate into a feature matrix
        feature_matrix = []
        for mol in molecules:
            features = soap_descriptor.create(mol.to_ase_atoms(), n_jobs=-1)
            feature_matrix.append(features)
        
        return np.array(feature_matrix)

    def plot_soap_feature_histograms_average(self, feature_mat: np.ndarray, n_features_to_plot: int = 5):
        """ 
        Plot histograms of the average SOAP feature values across the dataset for a specified number of features
        """
        # Choose n random feature dimensions to plot
        n_total_features = feature_mat.shape[1]
        feature_indices = np.random.choice(n_total_features, size=n_features_to_plot, replace=False)
        plt.figure(figsize=(15, 5))
        for i, idx in enumerate(feature_indices):
            plt.subplot(1, n_features_to_plot, i + 1)
            plt.hist(feature_mat[:, idx], bins=30, alpha=0.7)
            plt.title(f"Average SOAP Feature {idx}")
            plt.xlabel("Feature Value")
            plt.ylabel("Frequency")
        plt.tight_layout()
        plt.savefig("figures/average_soap_feature_histograms.png")
        plt.close()
        
        







if __name__ == "__main__":
    # Import test molecules from file
    test_molecules = np.load("initial_molecules.npy", allow_pickle=True)[0]
    print(test_molecules[1])
    # Initialize featurizer
    config = FeaturizerConfig(descriptor_type="SOAP")
    featurizer = Featurizer(config)
    # Compute features for test molecules
    features = featurizer.compute_soap_features(test_molecules)
    print("Feature matrix shape:", features.shape)

    # plot histograms of some feature dimensions 
    featurizer.plot_soap_feature_histograms_average(features, n_features_to_plot=3)

