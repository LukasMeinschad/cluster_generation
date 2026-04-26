import numpy as np
from typing import Tuple, List, Optional

from molecule_class import Molecule
from cluster_generation.modules.class_test_results.psi4_interface import direct_energy

# Import ASE based calculators
try:
    from ase import Atoms
    from ase.calculators.emt import EMT
except ImportError:
    print("ASE not found. ASE-based energy evaluation will be unavailable.")

# Import XTB calculator
try:
    from xtb.ase.calculator import XTB
except ImportError:
    print("XTB not found. XTB-based energy evaluation will be unavailable.")





class EnergyEvaluator:
    """
    Flexible energy evaluator that supports multiple computational methods

    Backends:
        - Psi4: Quantum chemical calculations (e.g., HF, DFT)
        - ASE: Semi-empirical methods (e.g., xTB) if ASE is available
    """

    # Conversion Factors
    EV_TO_HARTREE = 0.0367493  # 1 eV in Hartree
    E_ANG_TO_DEBYE = 4.80320427  # 1 eÅ in Debye


    
    def __init__(self, backend: str = "psi4", method: str = "hf", basis: str = "sto-3g", xtb_method: str = "GFN2-xTB"):
        self.backend = backend.lower()
        self.method = method
        self.basis = basis
        self.xtb_method = xtb_method

        # Initialize the arse calculator 
        self._calculator = None 
        if self.backend == "xtb":
            self._calculator = XTB(method=self.xtb_method)
        elif self.backend == "ase_emt":
            self._calculator = EMT()
    
    def evaluate(self, molecule: Molecule) -> Tuple[Optional[float], Optional[List[float]], Optional[float]]:
        """ 
        Evaluate energy and the dipole moment.

        Returns:
            Tuple of (energy in Hartree, dipole vector [dx,dy,dz] in Debye, dipole magnitude in Debye)
            or (None, None, None, None) if evaluation fails.
        """
        if self.backend == "psi4":
            return self._evaluate_psi4(molecule)
        elif self.backend in ["xtb", "ase_emt"]:
            return self._evaluate_ase(molecule)
        else:
            print(f"Unknown backend requested: {self.backend}")
            return None, None, None, None 

    def _evaluate_psi4(self, molecule: Molecule) -> Tuple[Optional[float], Optional[List[float]], Optional[float]]:
        """  
        Evaluate the Energy based on Psi4
        """ 
        try:
            return direct_energy(molecule, method=self.method, basis=self.basis)
        except Exception as e:
            print(f"Psi4 energy evaluation failed: {e}")
            return None, None, None, None
        
    def _evaluate_ase(self, molecule: Molecule) -> Tuple[Optional[float], Optional[List[float]], Optional[float]]:
        """ 
        Energy evaluation using ASE-based calculators like xTB or EMT
        """
        try:
            symbols = molecule.get_symbols()
            if symbols is None:
                print("Failed to extract symbols from molecule for ASE evaluation.")
                return None, None, None, None
            
            atoms = Atoms(symbols=symbols, positions=molecule.coordinates)
            atoms.calc = self._calculator

            # Calculate the Energy and convert to Hartree
            energy_ev = atoms.get_potential_energy()
            energy_hartree = energy_ev * self.EV_TO_HARTREE

            # Try to obtain the dipole moment
            try:
                dipole_moments = atoms.get_dipole_moment()
                # Convert from eÅ to Debye
                dipole_moments_debye = [d * self.E_ANG_TO_DEBYE for d in dipole_moments]
                dipole_magnitude = np.linalg.norm(dipole_moments)

                            

            except Exception as e:
                print(f"ASE dipole moment evaluation failed: {e}")
                dipole_moments = None
                dipole_magnitude = None 

            # Try to obtain partial charges --> Not the same thing as mulliken!
            try:
                mulliken_charges = atoms.get_charges()

            except Exception as e:
                print(f"ASE Mulliken charge evaluation failed: {e}")
                mulliken_charges = None


           
            return energy_hartree, dipole_moments_debye, dipole_magnitude, mulliken_charges


        except Exception as e:
            print(f"ASE energy evaluation failed: {e}")
            return None, None, None, None

    
