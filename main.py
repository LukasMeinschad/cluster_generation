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
from symmetry import MoleculeSymmetry
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

    config = Psi4Config(method = "hf", basis_set="cc-pvdz", memory="1 GB", num_threads=1)

    calc = Psi4Calculator(config=config, verbose=True)

    calc.test_basis_set_convergence(molecule=molecule, method="ccsd(t)")




    

    molecule.compute_bonds()

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
                                                num_points=100, 
                                                rotation=True, 
                                                rotation_grid_dim=4, 
                                                rotation_method="spherical")
    
    #Sampler.test_convergence_metrics(submolecule=submolecules[1],
    #                                       center=ref_frame.origin,
    #                                       radius=5,
    #                                       method="sphere",
    #                                       rotation=False)
    
    calc = Psi4Calculator(config=config, verbose=True)
    #results = calc.batch_single_point_energy(molecules=sampled_mols, parallel = True, n_processes=10)

    #print(calc.get_results_summary())


    Sampler.calculate_sampling_statistics(sampled_mols, submolecule=submolecules[1])
    logger.write_sampling_statistics(Sampler)
    logger.write_trajectory_sampling(sampled_mols)



  #  Calculator = Psi4Calculator(ls_of_molecules=sampled_mols, method="HF")
  #  results = Calculator.batch_single_point_calc(parallel=True,n_processes=20)
  #  logger.write_scf_batch_result(results)
#
  #  sucessful_mols, energies = zip(*results)
  #  Cluster_obj = MolecularCluster(sampled_molecules=sucessful_mols, energies=np.array(energies), logger=logger, reference_molecule=molecule, sampling_region=Sampler.sampling_region)
  #  Cluster_obj.analyze_h_bond_configurations(plot_info=True)
  #  Cluster_obj.plot_energy_distribution(bins=50, top_percent=20, diff_to_min=True)
  #  Cluster_obj.plot_energy_distribution_of_valid_hbond_molecules()
  #  #Cluster_obj.plot_energies_sampling_region(diff_to_mean=True)

#    Cluster_obj.calculate_rmsd_features(include_energy=True, normalize_features=True)
#    cluster_info = Cluster_obj.cluster_kmeans_silhouette_analysis(max_clusters=10, plot_analysis=True)
#
#    Cluster_obj.visalize_cluster_2d() 
#
#    cluster_reps = Cluster_obj.obtain_cluster_representatives(method="lowest_energy")
#
#    
#
#
#    #Perform geomtry optimization with the 5 lowest energy mols
#    Calculator = Psi4Calculator(ls_of_molecules=cluster_reps, method="MP2")
#    opt_results = Calculator.batch_geometry_optimization(parallel=True,n_processes=10)
#    logger.write_optimization_xyz_batch(opt_results)
    
    