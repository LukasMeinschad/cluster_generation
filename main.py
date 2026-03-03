import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

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
    
    # ── Main summary logger (cluster_gen.out) ───────────────────
    logger = Logger(name="main", log_file="cluster_gen.out", file_mode="w")
    logger.write_program_header()  # Don't reassign — returns None

    molecule = Molecule.from_xyz(xyz_content)

    if args.test:
        for test in args.test:
            if test == "method_basis_combinations":
                config = Config(method="hf", basis="cc-pvdz")
                calc = Psi4Calculator(config=config, verbose=False)
                results = calc.determine_method_and_basis_set_combinations(molecule=molecule)

    # ── Initialization ──────────────────────────────────────────
    # Appends to same summary log
    init_logger = Logger(name="init", log_file="cluster_gen.out", file_mode="a")

    init_config = InitializationConfig(
        method="hf",
        basis="cc-pvdz",
        box_type="sphere",
        box_scale_factor=5,
        min_distance=1.8,
        optimize_submolecules=True,
        verbose=False
    )
    initializer = ClusterInitializer(config=init_config, logger=init_logger)
    initial_molecule, submol_indices, simulation_box = initializer.initialize_from_xyz(args.i[0])

    # ── BHMC Phase A ────────────────────────────────────────────
    # Summary logger → cluster_gen.out
    # Worker detail  → bhmc_workers.log (separate file)
    bhmc_logger = Logger(name="bhmc", log_file="cluster_gen.out", file_mode="a")

    bhmc_config = BHMCConfig(
        temperature=800,
        method="hf",
        basis="cc-pvdz",
        verbose=False,
        adaptive_operators=True
    )

    bhmc_sampler = MultiPhaseBHMC(
        config=bhmc_config,
        simulation_box=simulation_box,
        logger=bhmc_logger,
        worker_log_file="bhmc_workers.log"  # Worker detail goes here
    )

    phase_a_candidates = bhmc_sampler.run_phase_a(
        initial_molecule=initial_molecule,
        submolecule_indices=submol_indices,
        n_structures_per_worker=28,
        n_processes=28
    )

    # Obtain all structures
    phase_a_structures = [structure for structure, energy in phase_a_candidates]

    # ── Analysis ────────────────────────────────────────────────
    logger_analysis = Logger(name="analysis", log_file="cluster_gen.out", file_mode="a")
    representatives = bhmc_sampler.analyse_phase_a_results(
        phase_a_candidates,
        submolecule_indices=submol_indices,
        n_clusters=10,
        simulation_box=simulation_box, 
        logger = logger_analysis
    )

    # ── Local Optimization ──────────────────────────────────────
    logger_opt = Logger(name="optimization", log_file="cluster_gen.out", file_mode="a")
    logger_opt.header("Local Optimization")

    config = Config(method="mp2", basis="cc-pvdz")
    calc = Psi4Calculator(config=config, verbose=False)
    optimization_results = calc.batch_optimize_parallel_unordered(representatives, n_processes=20)
    optimized_mols = [result.molecule for result in optimization_results if result.success]
    optimized_energies = [result.energy for result in optimization_results if result.success]

    logger_opt.info(f"Optimized {len(optimized_mols)}/{len(representatives)} structures successfully")

    # --- Write Trajectories
    # Phase A candidates
    phase_a_energies = [energy for structure, energy in phase_a_candidates]
    logger_opt.write_xyz_trajectory(
        molecules=phase_a_structures,
        filepath="trajectories/phase_a_candidates.xyz",
        energies=phase_a_energies,
    )

    # Write optimizes representatives
    logger_opt.write_xyz_trajectory(
        molecules=optimized_mols,
        filepath="trajectories/optimized_representatives.xyz",
        energies=optimized_energies
    )

    


    # ── Timing ──────────────────────────────────────────────────
    elapsed = time.time() - time_start
    logger_opt.separator(char="=")
    logger_opt.info(f"Total runtime: {elapsed:.1f}s ({elapsed/60:.1f} min)")
    logger_opt.separator(char="=")
