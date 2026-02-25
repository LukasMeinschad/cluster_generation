import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore", category=UserWarning)  # Supress library warnings for cleaner output

import sys
from pathlib import Path
import time

module_dir = Path(__file__).parent / "modules"
sys.path.insert(0, str(module_dir))

# Module imports
from molecule_class import Molecule
from transformations import Transformation
from psi4_interface import Psi4Calculator, Config
from cluster import BHMCAnalyzer
from bhmc import MultiPhaseBHMC, BHMCConfig, benchmark_temperature_acceptance
from init import ClusterInitializer, InitializationConfig
from args import get_args
from logger import Logger






if __name__ == "__main__":
    time_start = time.time()
    mp.set_start_method('spawn', force=True)

    # Parse the Arguments
    args = get_args()

    with open(args.i[0], "r") as file:
        xyz_content = file.read()
    
    logger = Logger(out_file="cluster_gen.out", mode="w")
    logger.write_header()

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


    init_config = InitializationConfig(
        method="mp2",
        basis="cc-pvdz",
        box_type="sphere",
        box_scale_factor = 5,
        min_distance=1.8,
        optimize_submolecules=True,
        verbose=True
    )
    initializer = ClusterInitializer(config=init_config)
    initial_molecule, submol_indices, simulation_box = initializer.initialize_from_xyz(args.i[0])


 #   # Run benchmark
 #   results = benchmark_temperature_acceptance(
 #       initial_molecule=initial_molecule,
 #       submolecule_indices=submol_indices,
 #       simulation_box=simulation_box,
 #       method="mp2",
 #       basis="cc-pvdz",
 #       n_steps=100,
 #       n_trials=5,
 #       save_plot=True,
 #       plot_filename="figures/temperature_acceptance.png"
 #   )

    # Set up BHMC Config
    bhmc_config = BHMCConfig(
        temperature=500,
        method="mp2",
        basis="cc-pvdz"
    )
    
    # Initialize BHMC Sampler
    bhmc_sampler = MultiPhaseBHMC(config=bhmc_config, simulation_box=simulation_box)

    phase_a_candidates = bhmc_sampler.run_phase_a(
        initial_molecule=initial_molecule, 
        submolecule_indices=submol_indices,
        n_structures_per_worker=500,
        n_processes=20
    )


    # Obtain all structures
    phase_a_structures = [structure for structure, energy in phase_a_candidates]
    logger.write_trajectory(phase_a_structures) 

    # Analyze results and get cluster representatives
    representatives = bhmc_sampler.analyse_phase_a_results(
        phase_a_candidates, 
        submolecule_indices=submol_indices,
        n_clusters=10,
        simulation_box=simulation_box
    ) 

    logger.write_trajectory(representatives, filename="representatives.xyz") 


    config = Config(method = "mp2", basis="cc-pvdz")  # Changed basis_set to basis
    calc = Psi4Calculator(config=config, verbose=False)
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
