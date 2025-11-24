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
from transformations import Transformation, ConfigSampler
from plotting import Plotting
from psi4_interface import Psi4Calculator
from cluster import MolecularCluster



        


if __name__ == "__main__":
    
    # Set multiprocessing start method
    # Start method is the way to start the child process in python
    # 'spawn' creates a new process
    # 'fork' would copy the parent process
    mp.set_start_method('spawn', force=True)
    
    with open('/media/storage_6/lme/master_thesis/initial_tests/h2o_2/h2o_2.xyz', 'r') as file:
        xyz_content = file.read()


   

    molecule = Molecule.from_xyz(xyz_content)
    cov_bonds, hydrogen_bonds = molecule.degree_of_covalence()
    submolecules = molecule.find_submolecules(cov_bonds)
    submolecules_labels = [submol.atom_labels for submol in submolecules]
    transformation = Transformation() 
    
    ref_frame, center_point = transformation.set_reference_frame_submolecule(submolecules = submolecules, submol_index=0, molecule=molecule, method="com", print_info=False, plot_frame=True)

    calculator = Psi4Calculator(molecule=molecule, method="scf", basis_set="sto-3g", memory="400 MB", num_threads=1)
    calculator.single_point_calc()

    config_sampler = ConfigSampler(molecule, ref_frame, center_point)
    
    sampling_cube = config_sampler.create_cube_volume(plot_cube=True, size=5)

    sampled_mols = config_sampler.place_submolecule_uniformly(submolecules, submol_trans_id=1,submol_fixed_id=0, num_samples=100, volume_type="cube", plot_samples=True)

    calculator = Psi4Calculator(ls_of_molecules=sampled_mols)
    results = calculator.batch_single_point_calc(parallel=True,n_processes=8)

    # Initialize Cluster Analysis
    energies = [energy for mol, energy in results]
    cluster_analysis = MolecularCluster(sampled_molecules=sampled_mols, energies=energies)
    cluster_analysis.calculate_com_features(submol_atom_labels=submolecules_labels)
    cluster_analysis.cluster_kmeans(n_clusters=3)
    cluster_analysis.visualize_cluster_com()
