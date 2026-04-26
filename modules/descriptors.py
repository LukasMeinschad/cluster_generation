import numpy as np 
from numba import njit, prange
from multiprocessing import Pool, cpu_count
from typing import List, Tuple

from modules.geometry import GeometryOps
from modules.molecule_class import Molecule


# ====================== Numba Core Functions =====================

@njit(fastmath=True)
def _fast_gyration_tensor_features(coords: np.ndarray, masses: np.ndarray) -> Tuple[float, float, float]:
    """ 
    Computes the Gyration Tensor and Asphericity, Acylindricity and Kappa² for single molecules
    """
    total_mass = np.sum(masses) 
    com = np.zeros(3)
    for i in range(len(masses)):
        com += masses[i] * coords[i]
    com /= total_mass

    # Build Gyration Tensor
    S = np.zeros((3, 3))
    for i in range(len(masses)):
        diff = coords[i] - com 
        for j in range(3):
            for k in range(3):
                S[j,k] += masses[i] * diff[j] * diff[k]
    S /= total_mass

    # Eigenvalues
    eigenvalues = np.linalg.eigvalsh(S)
    eigenvalue = np.sort(eigenvalues)  # Ensure sorted order
    l1, l2, l3 = eigenvalues[0], eigenvalues[1], eigenvalues[2]

    asphericity = l3 - 0.5 * (l1 + l2) 
    acylindricity = l2 - l1
    trace = l1 + l2 + l3

    kappa_squared = 0.0
    if trace > 1e-12:
        kappa_squared = (3.0/2.0) * (l1**2 + l2**2 + l3**2) / (trace**2) - 0.5
    

    return asphericity, acylindricity, kappa_squared

@njit(fastmath=True)
def _fast_radius_of_gyration(coords: np.ndarray, masses: np.ndarray) -> float:
    total_mass = np.sum(masses)
    com = np.zeros(3)
    for i in range(len(masses)):
        com += masses[i] * coords[i]
    com /= total_mass

    rg_sq = 0.0
    for i in range(len(masses)):
        diff = coords[i] - com 
        dist_sq = diff[0]**2 + diff[1]**2 + diff[2]**2
        rg_sq += masses[i] * dist_sq

    return np.sqrt(rg_sq / total_mass)

def worker_compute_features(args) -> np.ndarray:
    """  
    Multiprocessing worker function.
    Given a single structure's data, it computes all fast metrics and returns a feature vector.
    """
    coords, masses, submol_indices = args

    # 1. Radius of Gyration
    rg = _fast_radius_of_gyration(coords, masses)

    # 2. Gyration Tensor Features
    asphericity, acylindricity, kappa_squared = _fast_gyration_tensor_features(coords, masses)

    # 3. Intermolecular distance (using center of mass of first 2 submolecules)
    inter_dist = 0.0
    if submol_indices and len(submol_indices) >= 2:
        idx1, idx2 = submol_indices[0], submol_indices[1]
        com1 = GeometryOps.center_of_mass(coords[idx1], masses[idx1])
        com2 = GeometryOps.center_of_mass(coords[idx2], masses[idx2])
        inter_dist = np.linalg.norm(com1 - com2)

    return np.array([rg, inter_dist, asphericity, acylindricity, kappa_squared])


# ================= Parallel Feature Extractor 

class FeatureExtractor:
    """  
    Handels the parallel extraction of physical based features from structures
    """
    def __init__(self, n_processes: int = None):
        self.n_processes = n_processes or max(1, cpu_count() - 1)

    def extract_fast_features(self, structures: List, submolecule_indices: List[List[int]]) -> np.ndarray:
        """ 
        Extracts fast features for a list of structures in parallel.
        Returns a 2D array where each row corresponds to the feature vector of a structure.
        """
        worker_args = [
            (s.molecule.coordinates, s.molecule.masses, submolecule_indices)
            for s in structures
        ]

        with Pool(self.n_processes) as pool:
            features = pool.map(worker_compute_features, worker_args)

        return np.vstack(features)