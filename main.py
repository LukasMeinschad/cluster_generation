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
    
    #Plotting().plot_molecule_3d(molecule, cov_bonds=[(bond[0], bond[1]) for bond in cov_bonds], hydrogen_bonds=[(bond[0], bond[1]) for bond in hydrogen_bonds])
    transformation = Transformation() 
    
    ref_frame, center_point = transformation.set_reference_frame_submolecule(submolecules = submolecules, submol_index=0, molecule=molecule, method="com", print_info=False, plot_frame=True)


    config_sampler = ConfigSampler(molecule, ref_frame, center_point)
    
    sampling_cube = config_sampler.create_cube_volume(plot_cube=True)

    sampled_mols = config_sampler.place_submolecule_uniformly(submolecules, submol_trans_id=1,submol_fixed_id=0, num_samples=10, volume_type="cube", plot_samples=True)

    calculator = Psi4Calculator(ls_of_molecules=sampled_mols)

    calculator.batch_single_point_calc()
    calculator.batch_single_point_calc_parallel(n_processes=5)
    
    #calculator.build_geometry_string()
    #calculator.geometry_optimization()
