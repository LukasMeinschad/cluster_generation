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

# Import GPAW for DFT calculations
from gpaw import GPAW

# PSI4 calculator
from ase.calculators.psi4 import Psi4

# Import ASE optimizers
from ase.optimize import BFGS
from ase.optimize import LBFGS

# Import ASE Vibrational Analysis
from ase.vibrations import Vibrations


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
                            write_trajectory: bool = False) -> Molecule:
        """  
        Helper function to optimize the geometry of a molecule using the specified optimizer
        """
        if self.backend == "psi4":
            self._ensure_psi_scratch()
        atoms = molecule.to_ase_atoms()
        atoms.set_calculator(self.calculator)

        

        # Select the optimizer
        if optimizer.lower() == "bfgs":
            if write_trajectory:
                opt = BFGS(atoms, trajectory="optimization.traj")
            else:
                opt = BFGS(atoms)
        elif optimizer.lower() == "lbfgs":
            if write_trajectory:
                opt = LBFGS(atoms, trajectory="optimization.traj")
            else:
                opt = LBFGS(atoms)
        else:
            raise ValueError(f"Unsupported optimizer: {optimizer}. Supported optimizers are 'BFGS' and 'LBFGS'.")
        opt.run(fmax=0.001)  # Tight convergence required for reliable Hessian/frequencies
        optimized_molecule = Molecule.from_ase_atoms(atoms)
        return optimized_molecule

    def compute_hessian(self, molecule, delta=0.01):
        """ 
        Evaluate the vibrational properties and the return the 2D Hessian Matrix 
        """ 
        if self.backend == "psi4":
            self._ensure_psi_scratch()

        atoms = molecule.to_ase_atoms()
        dof = len(atoms) * 3
        n_vib = dof - 6  # Number of vibrational modes (3N-6 for non-linear molecules)
        atoms.set_calculator(self.calculator)
        vib = Vibrations(atoms, delta=delta)
        vib.clean()  # Remove stale cache files from previous runs
        vib.run()
        vibrations = vib.get_vibrations() # Returns the VibrationsData Object
        # Extract the Vibrational Frequencies
        frequencies = vibrations.get_frequencies() # Returns a 1D array of frequencies in cm^-1
        # Imaginary frequencies → negative by convention (saddle-point modes)
        hessian = vibrations.get_hessian_2d() # Returns the 2D Hessian matrix
        # Determine Eigenvalues, Eigenvectors, and the order of the saddle point
        eigenvalues, eigenvectors = np.linalg.eigh(hessian)

        # Set negative eigenvalues and near zero eigenvalues to zero for stability
        threshold = 1e-4  # Threshold for zero eigenvalues
        eigenvalues[np.abs(eigenvalues) < threshold] = 0.0 
        eigenvalues[eigenvalues < 0] = 0.0
        print("Eigenvalues of the Hessian:", eigenvalues)
        # Check if there are 3n-6 positive eigenvalues if not the difference is the order of the saddle point
        n_positive = np.sum(eigenvalues > 0)
        order_saddle_point = dof - 6 - n_positive
        local_minimum = order_saddle_point == 0

        

        return hessian, frequencies, order_saddle_point, local_minimum
    
    
    
    
        
 

    

def run_self_tests():
    """ 
    Run the self tests of the EnergyEvaluator class to ensure that all backends are functioning correctly
    """ 
    h2o_dimer_path = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/water_clusters/h2o_dimer.xyz"
    h2o_xyz = []
    with open(h2o_dimer_path, "r") as f:
        h2o_xyz = f.read()
    molecule = Molecule.from_xyz(h2o_xyz)
    print("Testing PSI4 backend...")
    psi4_evaluator = EnergyEvaluator(backend="psi4", qm_method="hf", qm_basis="cc-pvdz")
    # Optimize and calculate hessian
    optimized_molecule_psi4 = psi4_evaluator.optimize_geometry(molecule, optimizer="BFGS")
    hessian_psi4, frequencies_psi4, order_saddle_point, local_minimum = psi4_evaluator.compute_hessian(optimized_molecule_psi4)

    # for the water dimer we expect 3N = 18 mode
    # 3N-6 = 12 vibrational modes and 6 zero-frequency modes (3 translations + 3 rotations)

    print("hessian shape:", hessian_psi4.shape)
    print("frequencies shape:", frequencies_psi4.shape)
    print("PSI4 tests passed!")
    # Print the frequencies
    print("Vibrational frequencies (cm^-1):", frequencies_psi4)
    print("Order of saddle point:", order_saddle_point)
    print("Is local minimum:", local_minimum)
    



    

if __name__ == "__main__":
    run_self_tests()

    
    
