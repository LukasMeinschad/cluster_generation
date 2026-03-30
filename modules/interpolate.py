""" 
Module for construction of an interpolation path between two distinct minimum configurations
"""
from typing import List, Tuple, Optional
import numpy as np

from geometry import GeometryOps
from molecule_class import Molecule

if __name__ == "__main__":

    h2o_2_global_min = """
    6
    CCSD(T)-F12B/CVTZ-F12  ENERGY=-152.87631795
    H           0.0000000000       -0.0706396306       -0.5331292670
    H           0.0000000000       -0.8138599744       -1.8622054507
    H           0.7602291948        0.4099801780        1.7599189902
    H          -0.7602291948        0.4099801780        1.7599189902
    O           0.0000000000        0.0669822407       -1.4877213295
    O           0.0000000000       -0.0629156802        1.4168673828
    """
    h2o_2_local_min = """ 
    6
    CCSD(T)-F12B/CVTZ-F12  ENERGY=-152.87516258
    O           0.0628156806       -0.0251981904       -1.4028788915
    O          -0.0628156806        0.0251981904        1.4028788915
    H          -0.4574953073        0.4039197362       -2.0827543373
    H          -0.5778157355       -0.3632857456       -0.7713911917
    H           0.4574953073       -0.4039197362        2.0827543373
    H           0.5778157355        0.3632857456        0.7713911917
    """
    mol1 = Molecule.from_xyz(h2o_2_global_min)
    mol2 = Molecule.from_xyz(h2o_2_local_min)
    print("Molecules initialized")
    