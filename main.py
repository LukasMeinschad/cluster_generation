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

    # Set Reference Frame to submolecule 0
    Transformation = Transformation()
    ref_frame = Transformation.set_reference_frame(molecule=submolecules[0], method="com")
    logger.write_reference_frame_info(ref_frame)


    Sampler = ConfigSampler(reference_frame=ref_frame, center_point=center_point, molecule=molecule, logger=logger) 
#
#    # Test sampling along a plane
#    #Sampler.sample_plane_along_vector(center_point=[0,0,0], direction_vector=[1,1,1], normal_vector=[-1,1,0], s_length=5, t_params=[-5,5], num_points=200, plot_sampling=True)
#
#
    sampled_mols = Sampler.sample_submol_sphere(submolecule=submolecules[1],radius = 5, num_points=100, rotation=True,rotational_grid_dim=6,plot_sampling=False)
#    sampled_mols = Sampler.sample_submol_cone(submolecule=submolecules[1],apex=cone_apex,axis=cone_axis,height=8,base_radius=4,num_points=300,rotation=True,rotational_grid_dim=6,plot_sampling=False)
    #sampled_mols = Sampler.sample_submol_rectangle(submolecule=submolecules[1], center_point=center_point, dir_vector1=ref_frame[:,0], dir_vector2=ref_frame[:,1], length1=5, length2=5, num_points=20, rotation=True, rotational_grid_dim=4, plot_sampling=False)
    
    logger.write_trajectory_sampling(sampled_mols)

#    
    time_end = time.time() 
#
    Calculator = Psi4Calculator(ls_of_molecules=sampled_mols, method="HF")
    results = Calculator.batch_single_point_calc(parallel=True,n_processes=20)
    logger.write_scf_batch_result(results)
#
    sucessful_mols, energies = zip(*results)
    Cluster_obj = MolecularCluster(sampled_molecules=sucessful_mols, energies=np.array(energies), logger=logger, reference_molecule=molecule, sampling_region=Sampler.sampling_region)
    Cluster_obj.analyze_h_bond_configurations(plot_info=True)
    Cluster_obj.plot_energy_distribution(bins=50, top_percent=20, diff_to_min=True)
    Cluster_obj.plot_energy_distribution_of_valid_hbond_molecules()
    #Cluster_obj.plot_energies_sampling_region(diff_to_mean=True)

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
    
    