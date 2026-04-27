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
from modules.calculator import Calculator

# import dataclass
from dataclasses import dataclass

@dataclass
class StructureAnalysisConfig:
    """  
    Configuration for the StructureAnalysis class
    """
    calculator_backend: str = "xtb"
    calculator_qm_method: str = "hf"
    calculator_qm_basis: str = "sto-3g"
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
    def __init__(self, logger, mols: Optional[List[Molecule]] = None):
        self.logger = logger
        self.geometry_ops = GeometryOps()  
        self.mols = mols    


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

        self._log(f"Loaded {len(mols)} structures from {trajectory_path}")
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
