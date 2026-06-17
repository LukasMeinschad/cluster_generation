"""  
Structure Comparison Module
Goal of this module is to internally keep track of the obtained structures compare them to each other
in terms of their structural similarity and evaluate if they propose a local minima on the PES by analysis of their Hessian
"""

import numpy as np
from typing import Optional, Tuple, List, Dict, Any
import os    
from modules.molecule_class import Molecule
from modules.geometry import GeometryOps
from modules.calculator import EnergyEvaluator

from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from dataclasses import dataclass

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist
from dscribe.descriptors import CoulombMatrix

@dataclass
class StructureAnalysisConfig:
    """  
    Configuration for the StructureAnalysis class
    """
    calculator_backend: str = "xtb"
    calculator_qm_method: str = "mp2"
    calculator_qm_basis: str = "cc-pvdz"
    calculator_xtb_method: str = "GFN2-xTB"
    calculator_gpaw_mode: str = "lcao"
    calculator_gpaw_basis: str = "dzp"
    calculator_gpaw_xc: str = "B3LYP"

    # rmsd threshold for considering two structures as similar
    rmsd_threshold: float = 0.5  # in Angstrom
    sparsity_zero_threshold: float = 1e-1 # Threshold to consider a element in the hessian as zero
    
    


class StructureAnalysis:
    """ 
    Class for Analysis of the obtained structures during the BHMC algorithm
    Tracks both initial and optimized molecular geometries
    """
    def __init__(self, logger, mols: Optional[List[Molecule]] = None, config: Optional[StructureAnalysisConfig] = None):
        self.logger = logger
        self.geometry_ops = GeometryOps()

        # State Arrays to prevent tracking issues
        self.initial_mols: List[Molecule] = mols if mols is not None else []
        self.optimized_mols: Optional[List[Molecule]] = None
        self.optimized_energies: Optional[List[float]] = None

        self.config = config if config is not None else StructureAnalysisConfig()
        self._setup_calculator()

    def _setup_calculator(self):
        """  
        Sets up the calculator based on the configuration
        """
        self.calculator = EnergyEvaluator(
            backend=self.config.calculator_backend,
            qm_method=self.config.calculator_qm_method,
            qm_basis=self.config.calculator_qm_basis,
            xtb_method=self.config.calculator_xtb_method,
            gpaw_mode=self.config.calculator_gpaw_mode,
            gpaw_basis=self.config.calculator_gpaw_basis,
            gpaw_xc=self.config.calculator_gpaw_xc,
        )




    def _log(self, message: str):
        if self.logger:
            self.logger.info(message)
        else:
            # Fallback to printing if no logger is provided
            print(message)

    def _get_target_mols(self, use_optimized: bool) -> List[Molecule]:
        """  
        Helper to resolve targed coordinates savely
        """
        if use_optimized:
            if self.optimized_mols is None:
                self._log("No optimized molecules available, falling back to initial molecules.")
                return self.initial_mols
            return self.optimized_mols
        else:
            return self.initial_mols

    def load_mols_from_xyz(self, trajectory_path: str) -> list[Molecule]:
        """  
        Loads a given xyz trajectory file and returns a list of Molecule objects
        """ 
        with open(trajectory_path, 'r') as f:
            content = f.read()

        mols = [Molecule.from_xyz(frame) for frame in self._split_xyz_frames(content)]
        self._log(f"Loaded {len(mols)} molecules from {trajectory_path}")
        self.initial_mols = mols
        self.optimized_mols = None # Reset optimized mols when loading new structures
        self.optimized_energies = None
        return mols


    def from_mols(self, mols: List[Molecule]):
        """
        Initializes the StructureAnalysis object with a list of Molecule objects
        """
        self.initial_mols = mols
        self.optimized_mols = None
        self.optimized_energies = None
        self._log(f"Initialized StructureAnalysis with {len(mols)} structures")

    def compute_pairwise_rmsd(self, use_optimized: bool = False) -> np.ndarray:
        """ 
        Computes the pairwise RMSD between all structures.
        """
        mols = self._get_target_mols(use_optimized=use_optimized)
        n = len(mols)
        rmsd_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                rmsd = self.geometry_ops.compute_optimal_correspondence_rmsd(mols[i].coordinates, mols[j].coordinates)
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd
        
        self._log(f"Computed pairwise RMSD matrix (use_optimized={use_optimized})")
        return rmsd_matrix

    def get_unique_structures(self, use_optimized: bool = False) -> Tuple[List[int], List[Molecule]]:
        """
        Filters out structurally redundant geometries based on the configured RMSD threshold
        """
        self._log(f"Identifying unique structures using RMSD threshold of {self.config.rmsd_threshold} Å (use_optimized={use_optimized})")
        mols = self._get_target_mols(use_optimized=use_optimized)
        rmsd_matrix = self.compute_pairwise_rmsd(use_optimized=use_optimized)

        unique_indices = []
        for i in range(len(mols)):
            if all(rmsd_matrix[i, j] >= self.config.rmsd_threshold for j in unique_indices):
                unique_indices.append(i)
        unique_mols = [mols[i] for i in unique_indices]
        self._log(f"Identified {len(unique_mols)} unique structures (use_optimized={use_optimized})")
        return unique_indices, unique_mols

    def compute_coulomb_matrix_eigenspectra(self, use_optimized: bool = False) -> np.ndarray:
        """ 
        Computes rotation-invariant representation eigenspectra of the Coulomb matrix for all structures.
        """
        mols = self._get_target_mols(use_optimized=use_optimized)

        if not mols:
            self._log("No molecules available to compute Coulomb matrix eigenspectra.")
            return np.array([])
        max_atoms = max(mol.coordinates.shape[0] for mol in mols)
        cm = CoulombMatrix(n_atoms_max=max_atoms, permutation="eigenspectrum")
        ase_mols = [mol.to_ase_atoms() for mol in mols]
        return cm.create(ase_mols)
    
    @staticmethod
    def _compare_eigenspectra(eig1: np.ndarray, eig2: np.ndarray, metric: str = "euclidean") -> float:
        if metric == "euclidean":
            return np.linalg.norm(eig1 - eig2)
        if metric == "manhattan":
            return np.sum(np.abs(eig1 - eig2))
        if metric == "cosine":
            return 1 - np.dot(eig1, eig2) / (np.linalg.norm(eig1) * np.linalg.norm(eig2))
        if metric == "maximum":
            return np.max(np.abs(eig1 - eig2))
        raise ValueError(f"Unsupported metric: {metric}")
    
    def compute_eigenspectra_distance_matrix(self, metric: str = "euclidean", use_optimized: bool = False) -> np.ndarray:
        """  
        Computes a distance matrix between the Coulomb matrix eigenspectra of the structures using the specified metric.
        """
        eigenspectra = self.compute_coulomb_matrix_eigenspectra(use_optimized=use_optimized)
        n = eigenspectra.shape[0]
        distance_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i + 1, n):
                d = self._compare_eigenspectra(eigenspectra[i], eigenspectra[j], metric=metric)
                distance_matrix[i, j] = d
                distance_matrix[j, i] = d
        return distance_matrix

    def _save_plot(self, save_path: str):
        """  
        Helper to ensure path exists and save cleanly, further free plot memory
        """
        parent = os.path.dirname(save_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()


    def plot_eigenspectra_distance_heatmap(self, metric: str = "euclidean", use_optimized: bool = False, save_path: str = "figures/eigenspectra_distance_heatmap.png"):
        """ 
        Plots a distance map evaluation derived from Coloumb matrix eigenspectra
        """ 
        mols = self._get_target_mols(use_optimized=use_optimized)
        distance_matrix = self.compute_eigenspectra_distance_matrix(metric=metric, use_optimized=use_optimized)

        plt.figure(figsize=(10, 8))
        sns.heatmap(distance_matrix, cmap="viridis", xticklabels=np.arange(len(mols)), yticklabels=np.arange(len(mols)), annot=True, fmt=".2f")
        plt.title(f"Coulomb Matrix Eigenspectra Distance Heatmap ({metric} metric)")
        plt.xlabel("Structure Index")
        plt.ylabel("Structure Index")
        self._save_plot(save_path)

    @staticmethod  
    def compute_distance_matrix_flattened(mol) -> np.ndarray:
        """ 
        Computes a sorted 1D array of all unique pairise interatomic distances
        """
        return np.sort(pdist(mol.coordinates)) 

    @staticmethod
    def compare_distance_matrices(dist1: np.ndarray, dist2: np.ndarray, metric: str = "euclidean") -> float:
        if metric == "euclidean":
            return np.linalg.norm(dist1 - dist2)
        if metric == "manhattan":
            return np.sum(np.abs(dist1 - dist2))
        if metric == "cosine":
            return 1 - np.dot(dist1, dist2) / (np.linalg.norm(dist1) * np.linalg.norm(dist2))
        if metric == "maximum":
            return np.max(np.abs(dist1 - dist2))
        raise ValueError(f"Unsupported metric: {metric}")
    
    def plot_distance_matrix_heatmap(self, metric: str = "euclidean", use_optimized: bool = False,  save_path: str = "figures/distance_matrix_heatmap.png"):
        """  
        Plots a heatmap of the distance matrix between the flattened pairwise interatomic distances of the structures using the specified metric.
        """
        mols = self._get_target_mols(use_optimized=use_optimized)
        n = len(mols)
        distance_matrix = np.zeros((n,n))

        # Precompute all flattened distance matrices
        flat_distances = [self.compute_distance_matrix_flattened(mol) for mol in mols]

        for i in range(n):
            for j in range(i + 1, n):
                d = self.compare_distance_matrices(flat_distances[i], flat_distances[j], metric=metric)
                distance_matrix[i, j] = d
                distance_matrix[j, i] = d

        plt.figure(figsize=(10, 8))
        sns.heatmap(distance_matrix, cmap="viridis", xticklabels=np.arange(len(mols)), yticklabels=np.arange(len(mols)), annot=True, fmt=".2f")
        plt.title(f"Interatomic Distance Matrix Heatmap ({metric} metric)")
        plt.xlabel("Structure Index")
        plt.ylabel("Structure Index")
        self._save_plot(save_path) 

    def plot_pairwise_rmsd_heatmap(self, use_optimized: bool = False, save_path: str = "figures/pairwise_rmsd_heatmap.png"):
        """  
        Plots a heatmap of the pairwise RMSD matrix
        """
        mols = self._get_target_mols(use_optimized=use_optimized)
        rmsd_matrix = self.compute_pairwise_rmsd(use_optimized=use_optimized)

        plt.figure(figsize=(10, 8))
        sns.heatmap(rmsd_matrix, cmap="viridis", xticklabels=np.arange(len(mols)), yticklabels=np.arange(len(mols)), annot=True, fmt=".2f")
        plt.title(f"Pairwise RMSD Heatmap (use_optimized={use_optimized})")
        plt.xlabel("Structure Index")
        plt.ylabel("Structure Index")
        self._save_plot(save_path)

    def optimize_geometries(self, n_workers: int = 4) -> List[Molecule]:
        """
        Optimizes the geometries of all structures in parallel using multiprocessing
        """
        optimized_mols_tracker: Dict[int, Molecule] = {}
        optimized_energies_tracker: Dict[int, float] = {}
        ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as executor:
            futures = {
                executor.submit(self._optimize_geometry_worker, mol, self.config): idx
                for idx, mol in enumerate(self.initial_mols)
            }
            for future in as_completed(futures):
                idx = futures[future]
                try:
                    optimized_mols_tracker[idx], optimized_energies_tracker[idx] = future.result()
                except Exception as e:
                    self._log(f"Error optimizing molecule at index {idx}: {e}")

        # Re-Sort to maintain original order
        self.optimized_mols = [optimized_mols_tracker[i] for i in range(len(self.initial_mols)) if i in optimized_mols_tracker]
        self.optimized_energies = [optimized_energies_tracker[i] for i in range(len(self.initial_mols)) if i in optimized_energies_tracker]
        self._log(f"Completed geometry optimization for {len(self.optimized_mols)} structures")
        return self.optimized_mols



    @staticmethod
    def _optimize_geometry_worker(mol: Molecule, calculator_config: StructureAnalysisConfig) -> Tuple[Molecule, float]:
        calculator = EnergyEvaluator(
            backend=calculator_config.calculator_backend,
            qm_method=calculator_config.calculator_qm_method,
            qm_basis=calculator_config.calculator_qm_basis,
            xtb_method=calculator_config.calculator_xtb_method,
            gpaw_mode=calculator_config.calculator_gpaw_mode,
            gpaw_basis=calculator_config.calculator_gpaw_basis,
            gpaw_xc=calculator_config.calculator_gpaw_xc,
        )
        return calculator.optimize_geometry(mol, optimizer="LBFGS")

    def analyze_hessians(self, use_optimized: bool = True, n_workers: int = 4) -> List[Any]:
        """Computes tracking matrices for stationary point criteria evaluation."""
        mols = self._get_target_mols(use_optimized)
        results, errors = self.calculator.compute_hessians_parallel(mols, n_workers)
        
        for mol, e in errors:
            self._log(f"Error analyzing Hessian for molecule {mol.name}: {e}")
            
        self._log(f"Completed Hessian analysis for all structures (use_optimized={use_optimized})")
        for mol, analysis in results:
            self._log(f"Molecule: {mol.name} order_stationary_point={analysis['order_stationary_point']} global_minimum={analysis['global_minimum']}")
            self._log(f"Frequencies (cm^-1): {analysis['frequencies']}")
            
        self.hessian_analysis_results = results
        return results

    def analyze_hessian_sparsity(self, use_optimized: bool = True, n_workers: int = 4, save_path: str = "figures/hessian_sparsity.png"):
        """Determines the element distribution matrix densities inside computed Hessians."""
        if not hasattr(self, 'hessian_analysis_results'):
            self.analyze_hessians(use_optimized=use_optimized, n_workers=n_workers)
        
        sparsity_values = []
        for mol, analysis in self.hessian_analysis_results:
            hessian = analysis["hessian"]
            n_elements = hessian.size
            n_zero_elements = np.sum(np.abs(hessian) < self.config.sparsity_zero_threshold)
            sparsity = n_zero_elements / n_elements
            sparsity_values.append(sparsity)
            
        plt.figure(figsize=(8, 6))
        sns.histplot(sparsity_values, bins=20, kde=True)
        plt.title("Sparsity of Hessian Matrices")
        plt.xlabel("Sparsity (Fraction of Near-Zero Elements)")
        plt.ylabel("Frequency")
        self._save_plot(save_path)
    
    def plot_hessian_matrix(self, save_path: str = "figures/hessian_matrix.png"):
        """Plots structural stiffness profiles side by side over an axes arrangement grid."""
        n = len(self.hessian_analysis_results)
        if n == 0:
            self._log("No Hessian analysis results to plot.")
            return
        
        n_cols = 4
        n_rows = (n + n_cols - 1) // n_cols
        sns.set_theme(style="white")
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
        axes_flat = axes.flatten() if n_rows > 1 else [axes]

        for i, (mol, analysis) in enumerate(self.hessian_analysis_results):
            ax = axes_flat[i]
            sns.heatmap(analysis["hessian"], ax=ax, cmap="viridis", square=True, cbar=True, xticklabels=False, yticklabels=False)
            ax.set_title(f"Structure {i} Hessian", fontsize=14)

        for j in range(i + 1, n_rows * n_cols):
            fig.delaxes(axes_flat[j])

        self._save_plot(save_path)

    def plot_rmsd_comparison(self, save_path: str = "figures/rmsd_comparison.png"):
        """Generates cross-examination heatmaps contrasting structural profiles before and after optimization."""
        if self.optimized_mols is None:
            self._log("No optimized structures found. Run optimize_geometries() first.")
            return

        def _get_matrix(mols: List[Molecule]) -> np.ndarray:
            size = len(mols)
            mat = np.zeros((size, size))
            for idx_i in range(size):
                for idx_j in range(idx_i + 1, size):
                    val = self.geometry_ops.compute_optimal_correspondence_rmsd(mols[idx_i].coordinates, mols[idx_j].coordinates)
                    mat[idx_i, idx_j] = val
                    mat[idx_j, idx_i] = val
            return mat

        rmsd_before = _get_matrix(self.initial_mols)
        rmsd_after = _get_matrix(self.optimized_mols)

        fig, (ax_before, ax_after) = plt.subplots(1, 2, figsize=(16, 6))
        
        sns.heatmap(rmsd_before, ax=ax_before, cmap="viridis", annot=True, fmt=".1f", annot_kws={"size": 8},
                    xticklabels=np.arange(len(self.initial_mols)), yticklabels=np.arange(len(self.initial_mols)))
        ax_before.set_title("Pairwise RMSD — Before Optimization")
        
        sns.heatmap(rmsd_after, ax=ax_after, cmap="viridis", annot=True, fmt=".1f", annot_kws={"size": 8},
                    xticklabels=np.arange(len(self.optimized_mols)), yticklabels=np.arange(len(self.optimized_mols)))
        ax_after.set_title("Pairwise RMSD — After Optimization")

        self._save_plot(save_path)
        self._log(f"Saved RMSD comparison plot to {save_path}")

    def plot_optimized_energy_distribution(self, save_path: str = "figures/optimized_energy_distribution.png"):
        """Plots a barplot of the optimized energies across configurations.""" 
        if self.optimized_energies is None:
            self._log("No optimized energies found. Run optimize_geometries() first.")
            return

        # Compute differences from lowest energy for better visualization
        min_energy = min(self.optimized_energies)
        energy_differences = [e - min_energy for e in self.optimized_energies]
        
        # sort the structures by energy for better visualization
        sorted_indices = np.argsort(energy_differences)
        energy_differences = [energy_differences[i] for i in sorted_indices]


        plt.figure(figsize=(10, 6))
        x = np.arange(len(energy_differences))
        sns.barplot(x=x, y=energy_differences, hue=x, palette="viridis", legend=False)
        plt.xticks(np.arange(len(energy_differences)), [f"{i}" for i in sorted_indices])
        plt.title("Optimized Energy Distribution Across Configurations")
        plt.xlabel("Structure Index")
        plt.ylabel("Optimized Energy (Hartree)")
        self._save_plot(save_path)
        self._log(f"Saved optimized energy distribution plot to {save_path}")

    @staticmethod
    def _split_xyz_frames(text: str) -> List[str]:
        lines = text.splitlines()
        frames = []
        i = 0
        while i < len(lines):
            if not lines[i].strip():
                i += 1
                continue
            try:
                n_atoms = int(lines[i].strip())
            except ValueError:
                i += 1
                continue
            frame_lines = lines[i : i + n_atoms + 2]
            frames.append("\n".join(frame_lines))
            i += n_atoms + 2
        return frames

if __name__ == "__main__":
    """
    Run small self test to ensure that the StructureAnalysis class is working correctly
    """
    xyz_traj_path = "/media/storage_6/lme/master_thesis/cluster_generation/trajectories/initial_candidates.xyz"
    struc_config = StructureAnalysisConfig(rmsd_threshold = 0.5)
    structure_analysis = StructureAnalysis(logger=None, config=struc_config)
    structure_analysis.load_mols_from_xyz(xyz_traj_path)

    # Analyize raw input candidates
    structure_analysis.plot_pairwise_rmsd_heatmap(save_path="figures/pairwise_rmsd_heatmap_initial.png")

    # Run processing piplein
    optimized_structures = structure_analysis.optimize_geometries(n_workers=4)

    # Plot tracking verification
    structure_analysis.plot_distance_matrix_heatmap(metric="euclidean", use_optimized=True, save_path="figures/distance_matrix_heatmap_optimized.png")
