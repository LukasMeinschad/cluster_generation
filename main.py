import matplotlib.pyplot as plt
import numpy as np
from mendeleev.fetch import fetch_table
import itertools
import networkx as nx


import sys
from pathlib import Path

module_dir = Path(__file__).parent / "modules"
sys.path.insert(0, str(module_dir))

# Import modules
from molecule_class import Molecule
from transformations import Transformation

        


if __name__ == "__main__":
    with open('/media/storage_6/lme/master_thesis/initial_tests/h2o_2/h2o_2.xyz', 'r') as file:
        xyz_content = file.read()
   

    molecule = Molecule.from_xyz(xyz_content)
    cov_bonds, hydrogen_bonds = molecule.degree_of_covalence()
    submolecules = molecule.find_submolecules(cov_bonds)
    transformation = Transformation()
    transformation.com(molecule)
