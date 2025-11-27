import matplotlib.pyplot as plt
import numpy as np
from mendeleev.fetch import fetch_table
import itertools
import networkx as nx
import multiprocessing as mp

import sys
from pathlib import Path

module_dir = Path(__file__).parent / "modules"
sys.path.insert(0, str(module_dir))

# Import modules
from molecule_class import Molecule
from transformations import Transformation
from plotting import Plotting
from psi4_interface import Psi4Calculator
from cluster import MolecularCluster
from modules.args import get_args
from modules.ConfigSampler import ConfigSampler
from modules.logger import Logger 


if __name__ == "__main__":
    
    # Set multiprocessing start method
    # Start method is the way to start the child process in python
    # 'spawn' creates a new process
    # 'fork' would copy the parent process
    mp.set_start_method('spawn', force=True)

    # Parse the Arguments
    args = get_args()

    with open(args.i[0], "r") as file:
        xyz_content = file.read()
    

    logger = Logger(out_file="cluster_generation.log", mode="w")
    logger.write_header()

    molecule = Molecule.from_xyz(xyz_content)
    covalent_bonds, hydrogen_bonds = molecule.get_bonds()
    submolecules = molecule.fragment_by_connectivity()
    Transformation = Transformation()
    ref_frame, center_point = Transformation.set_reference_frame(molecule=molecule, method="com",plot_frame=False) 

    Sampler = ConfigSampler(reference_frame=ref_frame, center_point=center_point, molecule=molecule)
    h_bond_unit_vector = Transformation.axis_along_bond(molecule=molecule,atom_label1=hydrogen_bonds[0][0], atom_label2=hydrogen_bonds[0][1])
    print("H-bond unit vector:", h_bond_unit_vector)
    
    sampled_mols  = Sampler.sample_submol_cone(submolecule=submolecules[1], apex=molecule.get_coords_by_label(hydrogen_bonds[0][0])[0],
                                               axis = h_bond_unit_vector, base_radius=2, height=3, num_points=100)
    