"""
Calculator wrapper for energy evaluation. This is ASE-based and uses the ASE calculator interface
"""
try:
    from modules.molecule_class import Molecule
except ImportError:
    from molecule_class import Molecule

# ====== ASE imports ======
from ase import Atoms

# XTB and EMT calculators
from xtb.ase.calculator import XTB
from ase.calculators.emt import EMT

# Import ase io for write
from ase.io import write

# Import read for Hessian matrix analysis
from ase.io import read

# Import GPAW for DFT calculations
from gpaw import GPAW

# PSI4 calculator
from ase.calculators.psi4 import Psi4

# Import ASE optimizers
from ase.optimize import BFGS
from ase.optimize import LBFGS

# Import ASE Vibrational Analysis
from ase.vibrations import Vibrations

# Multiprocessing for parallel Hessian evaluation
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# Seaborn and mpl for plotting stuff
import matplotlib.pyplot as plt
import seaborn as sns

from typing import Optional, List, Tuple
import numpy as np
import os
import tempfile
from pathlib import Path

class EnergyEvaluator:
    """  
    Flexible energy evaluator that supports multiple computational methods
    """
    # Conversion factors
    EV_TO_HARTREE = 0.0367493  # 1 eV in Hartree
    E_ANG_TO_DEBYE = 4.80320427  # 1 eÅ in Debye

    def __init__(self, backend: str = "psi4",
                       qm_method: str = "hf",
                       qm_basis: str = "6-31g",
                       xtb_method: str = "GFN2-xTB",
                       gpaw_mode: str = "lcao",
                       gpaw_basis: str = "dzp",
                       gpaw_xc: str = "B3LYP"):
        """  
        Initialize the energy evaluator with the specified backend and method parameters
        """
        self.backend = backend.lower()
        self.qm_method = qm_method
        self.qm_basis = qm_basis
        self.xtb_method = xtb_method
        self.gpaw_mode = gpaw_mode
        self.gpaw_basis = gpaw_basis
        self.gpaw_xc = gpaw_xc
        self.calculator = None

        if self.backend == "xtb":
            self.calculator = XTB(method=self.xtb_method)
        elif self.backend == "ase_emt":
            self.calculator = EMT()
        elif self.backend == "psi4":
            self.calculator = Psi4(method=self.qm_method, basis=self.qm_basis)
        elif self.backend == "gpaw":
            self.calculator = GPAW(mode=gpaw_mode, basis=gpaw_basis, xc=gpaw_xc)
        else:
            raise ValueError(f"Unsupported backend: {self.backend}. Supported backends are 'psi4', 'xtb', and 'ase_emt'.")
        
    def _ensure_psi_scratch(self):
        """ 
        Validates or sets the PSI_SCRATCH environment variable
        """
        scratch_path = os.environ.get("PSI_SCRATCH")

        if not scratch_path:
            # Option A: create a default local temp directoy
            default_scratch = os.path.join(os.getcwd(), ".psi_scratch")
            os.makedirs(default_scratch, exist_ok=True)
            os.environ["PSI_SCRATCH"] = default_scratch
            print(f"Warning: PSI_SCRATCH not set. Using default scratch directory at {default_scratch}")
        else:
            # Validate that the existing path is writable
            path = Path(scratch_path)
            if not path.exists() or not os.access(path, os.W_OK):
                raise PermissionError(f"PSI_SCRATCH path '{scratch_path}' is not writable. Please set PSI_SCRATCH to a valid directory.")

    def evaluate(self, molecule: Molecule):
        """ 
        Evaluate the energy of the molecule
        """
        if self.backend == "psi4":
            self._ensure_psi_scratch()

        # Convert the molecule to a ASE Atoms object
        molecule_ase = molecule.to_ase_atoms()
        # Set the calculator and calculate energy
        molecule_ase.set_calculator(self.calculator)
        energy = molecule_ase.get_potential_energy()  # Energy in eV for ASE
        energy_hartree = energy * self.EV_TO_HARTREE  # Convert to Hartree
        return energy_hartree
    
    def evaluate_forces(self, molecule: Molecule):
        """ 
        Evaluate the forces (gradients) on the molecule
        """
        if self.backend == "psi4":
            self._ensure_psi_scratch()

        atoms = molecule.to_ase_atoms()
        atoms.set_calculator(self.calculator)
        forces = atoms.get_forces()  # Forces in eV/Å for ASE
        forces_hartree_per_angstrom = forces * self.EV_TO_HARTREE  # Convert to Hartree/Å
        return forces_hartree_per_angstrom
    
    def optimize_geometry(self, molecule: Molecule,
                            optimizer: str = "BFGS",
                            trajectory_fp: Optional[str] = None) -> Tuple[Molecule, float]:
        """  
        Helper function to optimize the geometry of a molecule using the specified optimizer
        Returns the optimized molecule and its converged energy in Hartree.
        """
        if self.backend == "psi4":
            self._ensure_psi_scratch()
        atoms = molecule.to_ase_atoms()
        atoms.set_calculator(self.calculator)

        def write_xyz():
            if trajectory_fp is not None:
                write(trajectory_fp, atoms, append=True)
            
        if trajectory_fp is not None:
            with open(trajectory_fp, "w") as f:
                f.write("")  # Clear the file contents

        # Select the optimizer
        if optimizer.lower() == "bfgs":
            opt = BFGS(atoms)
        elif optimizer.lower() == "lbfgs":
            opt = LBFGS(atoms)
        else:
            raise ValueError(f"Unsupported optimizer: {optimizer}. Supported optimizers are 'BFGS' and 'LBFGS'.")
        
        opt.attach(write_xyz, interval=1)
        opt.run(fmax=0.001)         
        energy_hartree = atoms.get_potential_energy() * self.EV_TO_HARTREE
        optimized_molecule = Molecule.from_ase_atoms(atoms)
        return optimized_molecule, energy_hartree

    def compute_hessian(self, molecule, delta=0.01):
        """
        Evaluate vibrational properties and return the 2D Hessian.
        Strips rotational/translational modes before classifying eigenvalues.
        Returns (hessian, frequencies, order_stationary_point, global_minimum).
        """
        atoms = molecule.to_ase_atoms()

        dof_total = 3 * len(atoms)
        dof_rot_trans = 6  # 3 translations + 3 rotations for a non-linear molecule

        with tempfile.TemporaryDirectory() as tmpdir:
            if self.backend == "psi4":
                # Ensure psi4 scratch is set for the subprocess
                self._ensure_psi_scratch()
                # Create a new calculator instance for the subprocess to avoid issues with shared state
                calc = Psi4(method=self.qm_method, basis=self.qm_basis)
            else:
                calc = self.calculator
            atoms.set_calculator(calc)
            vib_name = os.path.join(tmpdir, "vib")
            vib = Vibrations(atoms, name=vib_name, delta=delta)
            vib.run()
            vibrations = vib.get_vibrations()
            frequencies = vibrations.get_frequencies()
            hessian = vibrations.get_hessian_2d()


        eigenvalues, _ = np.linalg.eigh(hessian)
        # Sort eigenvalues from negative to positive
        eigenvalues_sorted = np.sort(eigenvalues)
        # The six lowest eigenvalues correspond to the 6 rotational + translational modes
        vibrational_eigenvalues = eigenvalues_sorted[dof_rot_trans:] 
        threshold = 1e-2 # Threshold to classify zero eigenvalues
        order_stationary_point = np.sum(np.abs(vibrational_eigenvalues) < threshold)
        global_minimum = order_stationary_point == 0
        
        return hessian, frequencies, order_stationary_point, global_minimum

        

    def compute_hessians_parallel(self, molecules: List[Molecule], n_workers: int = 4):
        """
        Compute Hessians for a list of molecules in parallel using multiprocessing.
        Returns (results, errors) where results is a list of (mol, analysis_dict) tuples
        and errors is a list of (mol, exception) tuples for any failed molecules.
        """
        calculator_kwargs = {
            "backend": self.backend,
            "qm_method": self.qm_method,
            "qm_basis": self.qm_basis,
            "xtb_method": self.xtb_method,
            "gpaw_mode": self.gpaw_mode,
            "gpaw_basis": self.gpaw_basis,
            "gpaw_xc": self.gpaw_xc,
        }
        results = []
        errors = []
        ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as executor:
            futures = {
                executor.submit(EnergyEvaluator._hessian_worker, mol, calculator_kwargs): mol
                for mol in molecules
            }
            for future in as_completed(futures):
                mol = futures[future]
                try:
                    results.append((mol, future.result()))
                except Exception as e:
                    errors.append((mol, e))
        return results, errors

    @staticmethod
    def _hessian_worker(mol: Molecule, calculator_kwargs: dict) -> dict:
        """Worker function for parallel Hessian computation (runs in subprocess)."""
        calc = EnergyEvaluator(**calculator_kwargs)
        hessian, frequencies, order_stationary_point, global_minimum = calc.compute_hessian(mol)
        return {
            "order_stationary_point": order_stationary_point,
            "global_minimum": global_minimum,
            "frequencies": frequencies,
            "hessian": hessian,
        }
    
    def plot_hessian_matrix_trajectory(self, trajectory_fp: str, save_path: str, n_workers: int = 4):
        """  
        Computes and Plots the Hessian matrices for each frame with a shared global colorbar. 
        """
        traj = read(trajectory_fp, index=":")
        molecules = [Molecule.from_ase_atoms(frame) for frame in traj]
        results, errors = self.compute_hessians_parallel(molecules, n_workers=n_workers)

        if errors:
            print(f"Warning: {len(errors)} molecules failed. Check logs.")

        n_frames = len(results)
        if n_frames == 0:
            return results, errors

        # 1. Find global min/max for a consistent color scale
        all_values = [analysis["hessian"] for _, analysis in results]
        vmin = min(h.min() for h in all_values)
        vmax = max(h.max() for h in all_values)

        n_cols = min(4, n_frames)
        n_rows = (n_frames + n_cols - 1) // n_cols

        # Increase figsize slightly to accommodate the global colorbar on the right
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols + 1, 4 * n_rows))
        axes_flat = np.array(axes).flatten()

        for i, (mol, analysis) in enumerate(results):
            ax = axes_flat[i]
            # 2. Use vmin/vmax so the colors are comparable across frames
            sns.heatmap(
                analysis["hessian"], 
                cmap="viridis", 
                cbar=False, 
                ax=ax, 
                vmin=vmin, 
                vmax=vmax,
                square=True
            )
            ax.set_title(f"Frame {i}")
            ax.set_xlabel("Coord Index")
            ax.set_ylabel("Coord Index")

        # Hide unused axes
        for j in range(n_frames, len(axes_flat)):
            axes_flat[j].set_visible(False)

        # 3. Create the global colorbar
        # We create a scalar mappable for the colorbar using the same cmap and limits
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([]) # Empty array for the mappable

        # Add colorbar to the figure, specifying the axes it should reference
        # 'ax' can be the flattened list; 'shrink' helps it fit vertically
        cbar = fig.colorbar(sm, ax=axes_flat, orientation='vertical', fraction=0.02, pad=0.04)
        cbar.set_label('Hessian Value Units', rotation=270, labelpad=15)

        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        # Note: tight_layout sometimes struggles with global colorbars; 
        # if it looks cramped, try fig.subplots_adjust(right=0.85)
        fig.tight_layout()
        fig.savefig(save_path, dpi=100, bbox_inches="tight")
        plt.close(fig)

        return results, errors

    


def run_self_tests():
    """ 
    Run the self tests of the EnergyEvaluator class to ensure that all backends are functioning correctly
    """ 
    h2o_dimer_path = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/water_clusters/h2o_dimer.xyz"
    #h2o_dimer_path = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/old/h2o_2_nonopt.xyz"
    #h2o_dimer_path = "/media/storage_6/lme/master_thesis/molpro_calcs/h2o_dimer/cyclic_c2/h2o_2_cyclic_c2.xyz"
    #h2o_dimer_path = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/water_clusters/h2o_timer.xyz"
    #h2o_dimer_path = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/h2o.xyz"
    h2o_xyz = []
    with open(h2o_dimer_path, "r") as f:
        h2o_xyz = f.read()
    molecule = Molecule.from_xyz(h2o_xyz)
    print("Testing PSI4 backend...")
    psi4_evaluator = EnergyEvaluator(backend="psi4", qm_method="mp2", qm_basis="cc-pvdz")
    # Optimize and calculate hessian
    optimized_molecule_psi4, energy_psi4 = psi4_evaluator.optimize_geometry(molecule, optimizer="LBFGS", trajectory_fp="test_molecules/h2o_opt.xyz")
    hessian_psi4, frequencies_psi4, order_stationary_point, global_minimum = psi4_evaluator.compute_hessian(optimized_molecule_psi4)

    print("Optimized energy (Hartree):", energy_psi4) 


    print("hessian shape:", hessian_psi4.shape)
    print("hessian:", hessian_psi4)
    
    print("frequencies shape:", frequencies_psi4.shape)
    print("PSI4 tests passed!")
    print("Vibrational frequencies (cm^-1):", frequencies_psi4)
    print("Order of stationary point:", order_stationary_point)
    print("Global minimum:", global_minimum)

    # Read in the trajectory and plot Hessian matrices for each frame
    print("Testing Hessian matrix plotting for trajectory...")
    psi4_evaluator.plot_hessian_matrix_trajectory(trajectory_fp="test_molecules/h2o_opt.xyz", save_path="figures/h2o_hessian_matrices.png", n_workers=20)




    

if __name__ == "__main__":
    run_self_tests()

    
    
