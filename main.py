import matplotlib.pyplot as plt
import numpy as np
from mendeleev.fetch import fetch_table
import itertools
import networkx as nx
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore", category=UserWarning)  # Supress library warnings for cleaner output

import sys
from pathlib import Path

module_dir = Path(__file__).parent / "modules"
sys.path.insert(0, str(module_dir))

# Import modules
from molecule_class import Molecule
from transformations import Transformation
from plotting import Plotting
from psi4_interface import Psi4Calculator, Config  # Changed from Psi4Config to Config
from cluster import MolecularCluster, BHMCAnalyzer
from symmetry import SymmetryAnalyzer
from graph import MolecularGraph
from coord_projector import CoordinateProjector
from bhmc import MultiPhaseBHMC, BHMCConfig
import bhmc as bhmc_module


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
                config = Config(method = "hf", basis="cc-pvdz")
                calc = Psi4Calculator(config=config, verbose=False)
                results = calc.determine_method_and_basis_set_combinations(molecule=molecule)
                logger.write_method_basis_combinations(results)

    #calc.test_basis_set_convergence(molecule=molecule, method="CCSD(T)")

    # Test computation of frequency 
    config = Config(method = "hf", basis="cc-pvdz")
    calc = Psi4Calculator(config=config, verbose=True)


    molecule.compute_bonds()
    Transformation = Transformation()

    submolecules = molecule.fragment_by_connectivity()
    # Get submolecule indices
    submol_indices = [submol.get_index_in_parent() for submol in submolecules]

    # Set up BHCM Config
    bhmc_config = BHMCConfig(temperature=300.0, max_steps=10, step_size=0.2,
                             method="hf", basis="cc-pvdz")
    # Initialize BHMC Sampler
    bhmc_sampler = MultiPhaseBHMC(config=bhmc_config) 
    phase_a_candidates = bhmc_sampler.run_phase_a(initial_molecule=molecule, submolecule_indices=submol_indices, n_structures_per_worker=10, n_processes=5)
    # phase a = (structure, energy)
    # obtain all structures
    phase_a_structures = [structure for structure, energy in phase_a_candidates]
    logger.write_trajectory(phase_a_structures) 

    representatives = bhmc_sampler.analyse_phase_a_results(phase_a_candidates) 
     
  
#    # Test computation of frequency 
    config = Config(method = "hf", basis="cc-pvdz")  # Changed basis_set to basis
    calc = Psi4Calculator(config=config, verbose=True)
    optimization_results = calc.batch_optimize_parallel_unordered(representatives, n_processes=20)
    optimized_mols = [result.molecule for result in optimization_results if result.success]
    logger.write_trajectory(optimized_mols, filename="optimized_representatives.xyz")

#
#    optimization_results = calc.batch_optimize_parallel_unordered(phase_a_structures,n_processes=10)
#
#    # Extract successful results
#    optimized_structures = [
#        (result.molecule, result.energy) 
#        for result in optimization_results 
#        if result.success
#    ]
#    optimized_mols = [mol for mol, energy in optimized_structures]
#    logger.write_trajectory(optimized_mols, filename="optimized_structures.xyz")
