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

@dataclass
class StructureAnalysisConfig:
    """  
    Configuration for the StructureAnalysis class
    """
    calculator_backend: str = "psi"
    calculator_qm_method: str = "mp2"
    calculator_qm_basis: str = "cc-pvdz"
    calculator_xtb_method: str = "GFN2-xTB"
    calculator_gpaw_mode: str = "lcao"
    calculator_gpaw_basis: str = "dzp"
    calculator_gpaw_xc: str = "B3LYP"

    # rmsd threshold for considering two structures as similar
    rmsd_threshold: float = 0.5  # in Angstrom
    
    


class StructureAnalysis:
    """ 
    Class for Analysis of the obtained structures during the BHMC algorithm
    """
    def __init__(self, logger, mols: Optional[List[Molecule]] = None, config: Optional[StructureAnalysisConfig] = None):
        self.logger = logger
        self.geometry_ops = GeometryOps()  
        self.mols = mols    
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
        sns.heatmap(self.rmsd_matrix, cmap="viridis", xticklabels=False, yticklabels=False, annot=True, fmt=".2f")
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
                    if self.logger:
                        self._log(f"Error optimizing geometry for molecule: {e}")
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
        optimized_mol = calculator.optimize_geometry(mol)
        return optimized_mol
        

    def analyze_hessians(self, n_workers: int = 4):
        """  
        Computes the Hessian for each structure and analyzes its eigenvalues to determine if the structure is a local minimum
        """
        results = []
        ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as executor:
            futures = {executor.submit(self._analyze_hessian_worker, mol, self.config): mol for mol in self.mols}
            for future in as_completed(futures):
                mol = futures[future]
                try:
                    result = future.result()
                    results.append((mol, result))
                except Exception as e:
                    if self.logger:
                        self._log(f"Error analyzing Hessian for molecule: {e}")
        if self.logger:
            self._log("Completed Hessian analysis for all structures")
        self.hessian_analysis_results = results
        # Log the Hessian analysis results
        if self.logger:
            for mol, analysis in results:
                self._log(f"Molecule: {mol.name} Ha: {analysis['n_negative']} negative, {analysis['n_zero']} zero, {analysis['n_positive']} positive eigenvalues")
                # Log the frequencies
                self._log(f"Molecule: {mol.name} Frequencies (cm^-1): {analysis['frequencies']}")

        return results

    @staticmethod
    def _analyze_hessian_worker(mol: Molecule, calculator_config: StructureAnalysisConfig):
        """ 
        Worker function for multiprocessing to compute the Hessian and analyze its eigenvalues for a given molecule
        """  
        
        # Setup the Calculator 
        calculator = EnergyEvaluator(
            backend=calculator_config.calculator_backend,
            qm_method=calculator_config.calculator_qm_method,
            qm_basis=calculator_config.calculator_qm_basis,
            xtb_method=calculator_config.calculator_xtb_method,
            gpaw_mode=calculator_config.calculator_gpaw_mode,
            gpaw_basis=calculator_config.calculator_gpaw_basis,
            gpaw_xc=calculator_config.calculator_gpaw_xc,
        )
        # Compute the hessian
        hessian, frequencies, order_saddle_point, local_minimum = calculator.compute_hessian(mol)

        # Return the order of the saddle point wether it is a minimum and the frequencies
        return {
            "order_saddle_point": order_saddle_point,
            "local_minimum": local_minimum,
            "frequencies": frequencies,
        } 

    

                

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
    #xyz_traj_path = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/water_dimer_local_minima/h2o_2_traj.xyz"
    xyz_traj_path = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/water_dimer_local_minima/h2o_2_globalmin.xyz"
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
    # Optimize
    optimized_mols = structure_analysis.optimize_geometries(n_workers=1)
    # Analyze Hessians
    results = structure_analysis.analyze_hessians(n_workers=1)

    for mol, analysis in results:
        print(f"Molecule: {mol.name} Order of saddle point: {analysis['order_saddle_point']} Local minimum: {analysis['local_minimum']} Frequencies (cm^-1): {analysis['frequencies']}")  
