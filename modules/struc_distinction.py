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

# Import concurrent futures
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# import dataclass
from dataclasses import dataclass

# Plotting 
import matplotlib.pyplot as plt
import seaborn as sns

@dataclass
class StructureAnalysisConfig:
    """  
    Configuration for the StructureAnalysis class
    """
    calculator_backend: str = "psi4"
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
    """
    def __init__(self, logger, mols: Optional[List[Molecule]] = None, config: Optional[StructureAnalysisConfig] = None):
        self.logger = logger
        self.geometry_ops = GeometryOps()
        self.mols = mols
        self.mols_before_opt: Optional[List[Molecule]] = None
        self.config = config if config is not None else StructureAnalysisConfig()

        # Setup the calculator based on the configuration
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
        self.logger.info(message)
    
    def load_mols_from_xyz(self, trajectory_path: str) -> list[Molecule]:
        """  
        Loads a given xyz trajectory file and returns a list of Molecule objects

        for this we open up the xyz file and read the coordinates and energies for each frame, then we create a Molecule object for each frame and store it in a list
        """
        with open(trajectory_path, 'r') as f:
            content = f.read()

        mols = []
        for frame in self._split_xyz_frames(content):
            mols.append(Molecule.from_xyz(frame))
        if self.logger:
            # Log the number of loaded structures
            self._log(f"Loaded {len(mols)} structures from {trajectory_path}")
        # store the mols in the class variable
        self.mols = mols
        return mols
    
    def from_mols(self, mols: List[Molecule]):
        """  
        Initializes the StructureAnalysis object with a list of Molecule objects
        """
        self.mols = mols
        self._log(f"Initialized StructureAnalysis with {len(mols)} structures")

    def compute_pairwise_rmsd(self) -> np.ndarray:
        """ 
        Computes the pairwise RMSD between all structures.
        Use the geometryops class to find the optimal correspondence between the atoms and then compute the RMSD based on the aligned coordinates 
        """
        n = len(self.mols)
        rmsd_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                mol_i = self.mols[i]
                mol_j = self.mols[j]
                rmsd = self.geometry_ops.compute_optimal_correspondence_rmsd(mol_i.coordinates, mol_j.coordinates)
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd  # Symmetric matrix
        self._log("Computed pairwise RMSD matrix")
        # store rmsd matris
        self.rmsd_matrix = rmsd_matrix

    def get_unique_structures(self) -> Tuple[List[int], List[Molecule]]:
        """
        Returns unique structures from the current self.mols without modifying state.
        Computes a fresh pairwise RMSD matrix on self.mols and filters by rmsd_threshold.
        Returns (unique_indices, unique_mols).
        """
        n = len(self.mols)
        rmsd_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                r = self.geometry_ops.compute_optimal_correspondence_rmsd(
                    self.mols[i].coordinates, self.mols[j].coordinates
                )
                rmsd_matrix[i, j] = r
                rmsd_matrix[j, i] = r

        unique_indices = []
        for i in range(n):
            if all(rmsd_matrix[i, j] >= self.config.rmsd_threshold for j in range(i)):
                unique_indices.append(i)

        unique_mols = [self.mols[i] for i in unique_indices]
        self._log(
            f"Found {len(unique_indices)} unique structures out of {n} "
            f"(threshold={self.config.rmsd_threshold} Å)"
        )
        return unique_indices, unique_mols



    def plot_pairwise_rmsd_heatmap(self, save_path: str = "figures/pairwise_rmsd_heatmap.png"):
        """  
        Plots a heatmap of the pairwise RMSD matrix
        """
        import matplotlib.pyplot as plt
        import seaborn as sns

        if not hasattr(self, 'rmsd_matrix'):
            self._log("Pairwise RMSD matrix not computed yet. Computing now...")
            self.compute_pairwise_rmsd()

        plt.figure(figsize=(10, 8))
        # Annote x and y with indices
        sns.heatmap(self.rmsd_matrix, cmap="viridis", xticklabels=np.arange(len(self.mols)), yticklabels=np.arange(len(self.mols)), annot=True, fmt=".2f")
        plt.title("Pairwise RMSD Heatmap")
        plt.xlabel("Structure Index")
        plt.ylabel("Structure Index")
        plt.tight_layout()
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300)
        plt.close()

    def optimize_geometries(self, n_workers: int = 4):
        """
        Optimizes the geometries of all structures in parallel using multiprocessing
        """
        self.mols_before_opt = list(self.mols)
        optimized_mols = []
        ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as executor:
            futures = {executor.submit(self._optimize_geometry_worker, mol, self.config): mol for mol in self.mols}
            for future in as_completed(futures):
                mol = futures[future]
                try:
                    optimized_mol = future.result()
                    optimized_mols.append(optimized_mol)
                except Exception as e:
                    msg = f"Error optimizing geometry for molecule: {e}"
                    if self.logger:
                        self._log(msg)
                    else:
                        print(msg)
        if self.logger:
            self._log("Completed geometry optimization for all structures")
        # Update the mols with the optimized geometries
        self.mols = optimized_mols
        return optimized_mols

    @staticmethod
    def _optimize_geometry_worker(mol: Molecule, calculator_config: StructureAnalysisConfig) -> Molecule:
        """ 
        Worker function to carry out a geometry optimization
        TODO: Problem is when we switch between methods for example go from xtb to psi4 we need to also optimize our
              structures with the new method to get a correct hessian analysis, but this is very time consuming, so maybe we can do a single point hessian analysis with the new method and then only optimize the structures that are close to a local minima based on the hessian analysis
        """
        calculator = EnergyEvaluator(
            backend=calculator_config.calculator_backend,
            qm_method=calculator_config.calculator_qm_method,
            qm_basis=calculator_config.calculator_qm_basis,
            xtb_method=calculator_config.calculator_xtb_method,
            gpaw_mode=calculator_config.calculator_gpaw_mode,
            gpaw_basis=calculator_config.calculator_gpaw_basis,
            gpaw_xc=calculator_config.calculator_gpaw_xc,
        )
        optimized_mol = calculator.optimize_geometry(mol, optimizer="LBFGS")
        return optimized_mol
        

    def analyze_hessians(self, n_workers: int = 4):
        """
        Computes the Hessian for each structure in parallel and analyzes eigenvalues to
        determine if each structure is a local minimum.
        """
        results, errors = self.calculator.compute_hessians_parallel(self.mols, n_workers)
        for mol, e in errors:
            msg = f"Error analyzing Hessian for molecule: {e}"
            if self.logger:
                self._log(msg)
            else:
                print(msg)
        if self.logger:
            self._log("Completed Hessian analysis for all structures")
            for mol, analysis in results:
                self._log(f"Molecule: {mol.name} order_stationary_point={analysis['order_stationary_point']} global_minimum={analysis['global_minimum']}")
                self._log(f"Molecule: {mol.name} Frequencies (cm^-1): {analysis['frequencies']}")
        self.hessian_analysis_results = results
        return results

    def analyze_hessian_sparsity(self, n_workers: int = 4, save_path: str = "figures/hessian_sparsity.png"):
        """  
        Determines the sparsity of the hessian matrices of the analyzed structures and plots a histogram
        """
        if not hasattr(self, 'hessian_analysis_results'):
            self._log("Hessian analysis results not computed yet. Computing now...")
            self.analyze_hessians(n_workers=n_workers)
        
        sparsity_values = []
        for mol, analysis in self.hessian_analysis_results:
            hessian = analysis["hessian"]
            n_elements = hessian.size
            n_zero_elements = np.sum(np.abs(hessian) < self.config.sparsity_zero_threshold)
            sparsity = n_zero_elements / n_elements
            sparsity_values.append(sparsity)
        self._log("Computed sparsity values for Hessian matrices")
        self._log(f"Sparsity values: {sparsity_values}")
        plt.figure(figsize=(8, 6))
        sns.histplot(sparsity_values, bins=20, kde=True)
        plt.title("Sparsity of Hessian Matrices")
        plt.xlabel("Sparsity (Fraction of Near-Zero Elements)")
        plt.ylabel("Frequency")
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()
    
    def plot_hessian_matrix(self, save_path: str = "figures/hessian_matrix.png"):
        """ 
        Helper function to plot all the hessian matrices of the analyzed structures in a grid
        """
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
            hessian = analysis["hessian"]

            sns.heatmap(
                hessian,
                ax=ax,
                cmap = "viridis",
                square=True,
                cbar=True,
                xticklabels=False,
                yticklabels=False,
            )
            ax.set_title(f"Structure {i} Hessian", fontsize=14)

        # remove empty subplots
        for j in range(i + 1, n_rows * n_cols):
            fig.delaxes(axes_flat[j])

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")   
        plt.close()

    

                

    def plot_rmsd_comparison(self, save_path: str = "figures/rmsd_comparison.png"):
        """
        Plots pairwise RMSD heatmaps before and after geometry optimization side by side.
        Requires optimize_geometries() to have been called first.
        """
        if self.mols_before_opt is None:
            self._log("No pre-optimization structures stored. Run optimize_geometries() first.")
            return

        def _pairwise_rmsd(mols: List[Molecule]) -> np.ndarray:
            n = len(mols)
            mat = np.zeros((n, n))
            for i in range(n):
                for j in range(i + 1, n):
                    r = self.geometry_ops.compute_optimal_correspondence_rmsd(
                        mols[i].coordinates, mols[j].coordinates
                    )
                    mat[i, j] = r
                    mat[j, i] = r
            return mat

        rmsd_before = _pairwise_rmsd(self.mols_before_opt)
        rmsd_after  = _pairwise_rmsd(self.mols)

        # format and log the RMSD matrices
        self._log("Pairwise RMSD values before optimization:")
        for i in range(len(self.mols_before_opt)):
            row = "  ".join(f"{rmsd_before[i, j]:.2f}" for j in range(len(self.mols_before_opt)))
            self._log(f"Structure {i}: {row}")
        self._log("Pairwise RMSD values after optimization:")
        for i in range(len(self.mols)):
            row = "  ".join(f"{rmsd_after[i, j]:.2f}" for j in range(len(self.mols)))
            self._log(f"Structure {i}: {row}")

        n_before = len(self.mols_before_opt)
        n_after  = len(self.mols)

        fig, (ax_before, ax_after) = plt.subplots(1, 2, figsize=(16, 6))
        sns.heatmap(
            rmsd_before, ax=ax_before, cmap="viridis", annot=True, fmt=".1f",
            annot_kws={"size": 8},
            xticklabels=np.arange(n_before), yticklabels=np.arange(n_before),
        )
        ax_before.set_title("Pairwise RMSD — Before Optimization")
        ax_before.set_xlabel("Structure Index")
        ax_before.set_ylabel("Structure Index")

        sns.heatmap(
            rmsd_after, ax=ax_after, cmap="viridis", annot=True, fmt=".1f",
            annot_kws={"size": 8},
            xticklabels=np.arange(n_after), yticklabels=np.arange(n_after),
        )
        ax_after.set_title("Pairwise RMSD — After Optimization")
        ax_after.set_xlabel("Structure Index")
        ax_after.set_ylabel("Structure Index")

        plt.tight_layout()
        parent = os.path.dirname(save_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
        self._log(f"Saved RMSD comparison plot to {save_path}")

    @staticmethod
    def _split_xyz_frames(text: str) -> list[str]:
        """Split a concatenated multi-frame XYZ string into individual single-frame strings."""
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
    xyz_traj_path = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/water_dimer_local_minima/h2o_2_traj.xyz"
    #xyz_traj_path = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/water_dimer_local_minima/h2o_2_globalmin.xyz"
    struc_config = StructureAnalysisConfig(
        calculator_backend="psi4",
        calculator_qm_method="mp2",
        calculator_qm_basis="cc-pvdz",
        calculator_xtb_method="GFN2-xTB",
        calculator_gpaw_mode="lcao",
        calculator_gpaw_basis="dzp",
        calculator_gpaw_xc="B3LYP",
        rmsd_threshold=0.5,
    )
    structure_analysis = StructureAnalysis(logger=None, config=struc_config)
    mols = structure_analysis.load_mols_from_xyz(xyz_traj_path)
    print(f"Loaded {len(mols)} molecules from trajectory")
    # Optimize
    optimized_mols = structure_analysis.optimize_geometries(n_workers=4)
    # Analyze Hessians
    results = structure_analysis.analyze_hessians(n_workers=4)

    for mol, analysis in results:
        print(f"Molecule: {mol.name} order_stationary_point={analysis['order_stationary_point']} global_minimum={analysis['global_minimum']}")

    structure_analysis.plot_hessian_matrix(save_path="figures/hessian_matrices.png")