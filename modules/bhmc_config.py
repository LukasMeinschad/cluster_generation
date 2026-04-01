import time
from dataclasses import dataclass


@dataclass
class BHMCConfig:
    """Configuration for Basin Hopping Monte Carlo sampling."""
    temperature: float = 400.0  # Temperature in Kelvin
    method: str = "hf"          # Quantum chemistry method
    basis: str = "sto-3g"       # Basis set
    backend: str = "psi4"       # Backend for energy evaluation 
    xtb_method: str = "GFN2-xTB" # Method for XTB if backend is set to xtb
    


    verbose: bool = False       # Enable debug output in workers
    adaptive_operators: bool = True  # Use adaptive scaling for operators


    # In-worker adaptive box control (Phase) A
    adaptive_box: bool = True
    box_update_interval: int = 10  # Update every N steps
    box_target_acceptance: float = 0.6 # Desired acceptance rate 
    box_acceptance_window: float = 0.05 # no update of box rate inside [box_target_acceptance - window, box_target_acceptance + window]
    box_growth_kp: float = 0.6 # propotional gain
    box_growth_max: float = 1.15 # cap per update
    box_max_scale: float = 2 # Relative to initial box size
    box_stable_windows: int = 3 # stop growing after this many stable windows
    
    # Adaptive Temperature Control in Phase B
    max_temperature: float = 50 



