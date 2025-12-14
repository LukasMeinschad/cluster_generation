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
from psi4_interface import Psi4Calculator, Psi4Config
from cluster import MolecularCluster
from symmetry import SymmetryAnalyzer
from graph import MolecularGraph


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
    
    # Initialize Logger
    logger = Logger(out_file="cluster_gen.out", mode= "w")
    logger.write_header()
    # Swith to append mode for further logging
    logger = Logger(out_file="cluster_gen.out", mode="a")

    molecule = Molecule.from_xyz(xyz_content)

    
    if args.test:
        # Run through all tests and check which one was selected
        for test in args.test:
            if test == "method_basis_combinations":
                # Initialize Psi4 Calculator 
                config = Psi4Config(method = "hf", basis_set="cc-pvdz", memory="1 GB", num_threads=1)
                calc = Psi4Calculator(config=config, verbose=False)
                results = calc.determine_method_and_basis_set_combinations(molecule=molecule)
                logger.write_method_basis_combinations(results)

    #calc.test_basis_set_convergence(molecule=molecule, method="CCSD(T)")

    # Test computation of frequency 
    config = Psi4Config(method = "hf", basis_set="cc-pvdz", memory="2 GB", num_threads=1)
    calc = Psi4Calculator(config=config, verbose=True)

    

    molecule.compute_bonds()

    SymmetryAnalyzer = SymmetryAnalyzer(molecule=molecule)
    SymmetryAnalyzer.check_linearity()
    SymmetryAnalyzer.check_inversion_center()
    SymmetryAnalyzer.check_planarity()
    SymmetryAnalyzer.find_symmetry_planes(tolerance=1e-4)    
    SymmetryAnalyzer.find_all_rotation_axes()
    # Log Molecule Info
    logger.write_molecule_info(molecule)

    # Get Submolecules
    submolecules = molecule.fragment_by_connectivity()
    logger.write_submolecule_info(submolecules)
        
    # Find possible H-bond configurations 
    configurations = molecule.find_hbond_configurations()
    valid_configurations = molecule.get_valid_hbond_configurations(angle_threshold= 150.0)
    logger.write_hbond_configurations(configurations)

    coords_sanity = submolecules[1].coordinates.copy()

    # Set Reference Frame to submolecule 0
    Transformation = Transformation()
    ref_frame = Transformation.set_reference_frame_submolecule(submolecule=submolecules[0])
    logger.write_reference_frame_info(ref_frame)

    submolecules = molecule.fragment_by_connectivity()

    
        

    # Initialize Sampler
    Sampler = ConfigSampler(reference_frame=ref_frame)
    #Sampler.test_sampling_speed(submolecule=submolecules[1], method="sphere", rotation=False)
    sampled_mols = Sampler.sample_submol_sphere(submolecule=submolecules[1],
                                                center = ref_frame.origin, 
                                                radius=5,  
                                                num_points=2, 
                                                rotation=True, 
                                                rotation_grid_dim=4, 
                                                rotation_method="spherical")
    
    #Sampler.test_convergence_metrics(submolecule=submolecules[1],
    #                                       center=ref_frame.origin,
    #                                       radius=5,
    #                                       method="sphere",
    #                                       rotation=False)
    
    calc = Psi4Calculator(config=config, verbose=False)
    results = calc.batch_single_point_energy(molecules=sampled_mols, parallel = True, n_processes=30)
    
    energies = [energy for _, energy in results]
    


    Sampler.calculate_sampling_statistics(sampled_mols, submolecule=submolecules[1])
    logger.write_sampling_statistics(Sampler)
    logger.write_trajectory_sampling(sampled_mols)

    Cluster = MolecularCluster(sampled_molecules = sampled_mols,
                               energies = energies,
                                reference_molecule = submolecules[1],
                                logger=logger)
    Cluster.analyze_hydrogen_bonds()
    Cluster.plot_hydrogen_bond_statistics()
    Cluster.plot_energy_distribution()


