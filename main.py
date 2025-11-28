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
import time
import timeit


if __name__ == "__main__":
    time_start = time.time()

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
    # Switch to append mode for further logging
    logger = Logger(out_file="cluster_generation.log", mode="a")

    molecule = Molecule.from_xyz(xyz_content)
    covalent_bonds, hydrogen_bonds = molecule.get_bonds()
    submolecules = molecule.fragment_by_connectivity()

    logger.write_molecule_info(molecule)

    Transformation = Transformation()
    ref_frame, center_point = Transformation.set_reference_frame_submolecule(submolecule=submolecules[0],parent_molecule=molecule,method="com")



    Sampler = ConfigSampler(reference_frame=ref_frame, center_point=center_point, molecule=molecule) 
    
    sampled_mols = Sampler.sample_submol_sphere(submolecule=submolecules[1],radius = 3, num_points=10, rotation=True,rotational_grid_dim=4,plot_sampling=False)
    time_end = time.time() 

    Calculator = Psi4Calculator(ls_of_molecules=sampled_mols, method="MP2")
    results = Calculator.batch_single_point_calc(parallel=True,n_processes=8)
    logger.write_scf_batch_result(results)


    Cluster_obj = MolecularCluster(sampled_molecules=sampled_mols, energies=[energy for _, energy in results])
    lowest_energy_mols, lowest_energies = Cluster_obj.select_n_lowest_energy(n=2)

    #Perform geomtry optimization with the 5 lowest energy mols
    Calculator = Psi4Calculator(ls_of_molecules=lowest_energy_mols, method="MP2")
    opt_results = Calculator.batch_geometry_optimization(parallel=True,n_processes=8)
    logger.write_optimization_xyz_batch(opt_results)
    
    