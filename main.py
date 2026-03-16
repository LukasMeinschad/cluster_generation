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

    N_WORKERS = 40

    # Parse the Arguments
    args = get_args()

    with open(args.i[0], "r") as file:
        xyz_content = file.read()
    
    # ── Main summary logger (cluster_gen.out) ───────────────────
    logger = Logger(name="main", log_file="cluster_gen.out", file_mode="w")
    logger.write_program_header()

    molecule = Molecule.from_xyz(xyz_content)

    if args.test:
        for test in args.test:
            if test == "method_basis_combinations":
                config = Config(method="hf", basis="cc-pvdz")
                calc = Psi4Calculator(config=config, verbose=False)
                results = calc.determine_method_and_basis_set_combinations(molecule=molecule)

    # ── Initialization ──────────────────────────────────────────
    init_logger = Logger(name="init", log_file="cluster_gen.out", file_mode="a")

    init_config = InitializationConfig(
        method="hf",
        basis="cc-pvdz",
        box_type="sphere",
        box_scale_factor=2.0,
        min_distance=1.5,
        optimize_submolecules=True,
        verbose=False
    )
    initializer = ClusterInitializer(config=init_config, logger=init_logger)
    initial_molecules, submol_indices, simulation_box = initializer.initialize_from_xyz(
        args.i[0], 
        n_configurations=N_WORKERS
    )

    # ── BHMC Phase A: Global Exploration ────────────────────────
    bhmc_logger = Logger(name="bhmc", log_file="cluster_gen.out", file_mode="a")

    bhmc_config = BHMCConfig(
        temperature=800,
        method="hf",
        basis="cc-pvdz",
        verbose=False,
        adaptive_operators=True,
        adaptive_box = True,
        box_update_interval=5, # Update box every 5 steps
        box_target_acceptance =0.6,
        box_acceptance_window = 0.05,
        box_growth_max = 1.15,
        box_max_scale = 4.0,
        box_stable_windows = 3
    )

    bhmc_sampler = MultiPhaseBHMC(
        config=bhmc_config,
        simulation_box=simulation_box,
        logger=bhmc_logger,
        worker_log_file="bhmc_workers.log"
    )

    phase_a_candidates = bhmc_sampler.run_phase_a(
        initial_molecules=initial_molecules,
        submolecule_indices=submol_indices,
        n_structures_per_worker=20,
        n_processes=N_WORKERS
    )

    

    # ── Analysis & Clustering ───────────────────────────────────
    logger_analysis = Logger(name="analysis", log_file="cluster_gen.out", file_mode="a")
    representatives = bhmc_sampler.analyse_phase_a_results(
        phase_a_candidates,
        submolecule_indices=submol_indices,
        cluster_method ="dbscan",
        eps=1.8,
        min_samples=5,
        simulation_box=simulation_box, 
        logger=logger_analysis
    )

    # Kabsch Interpolation
    interpolated_candidates = bhmc_sampler.generate_interpolated_candidates(
        representatives=representatives,
        submolecule_indices=submol_indices,
        n_interpolations= 2  
        )

    # Extract molecules from interpolated candidates for Phase B
    interpolated_molecules = [candidate.molecule for candidate in interpolated_candidates]
    logger_analysis.write_xyz_trajectory(
        molecules=interpolated_molecules,
        filepath="trajectories/interpolated_candidates.xyz",
        energies=None
    )

    # ── BHMC Phase B: Local Refinement ──────────────────────────
    # Switch to lower temperature and disable adaptive scaling
    bhmc_sampler.config.temperature = 150
    bhmc_sampler.config.adaptive_operators = False

    phase_b_results = bhmc_sampler.run_phase_b(
        representatives=interpolated_candidates,
        submolecule_indices=submol_indices,
        n_steps_per_worker=50,
    )

    # Extract best structures from Phase B (rather than raw representatives)
    phase_b_molecules = []
    phase_b_energies = []
    for result in phase_b_results:
        best = result["best_structure"]
        if best is not None:
            phase_b_molecules.append(best[0])
            phase_b_energies.append(best[1])

    # ── Local Optimization ──────────────────────────────────────
    logger_opt = Logger(name="optimization", log_file="cluster_gen.out", file_mode="a")
    logger_opt.header("Local Optimization")

    config = Config(method="mp2", basis="cc-pvdz")
    calc = Psi4Calculator(config=config, verbose=False)
    optimization_results = calc.batch_optimize_parallel_unordered(phase_b_molecules, n_processes=20)
    optimized_mols = [result.molecule for result in optimization_results if result.success]
    optimized_energies = [result.energy for result in optimization_results if result.success]

    logger_opt.info(f"Optimized {len(optimized_mols)}/{len(phase_b_molecules)} structures successfully")

    # ── Write Trajectories ──────────────────────────────────────
    phase_a_structures = [candidate[0] for candidate in phase_a_candidates]
    phase_a_energies = [candidate[1] for candidate in phase_a_candidates] 

    logger_opt.write_xyz_trajectory(
        molecules=phase_a_structures,
        filepath="trajectories/phase_a_candidates.xyz",
        energies=phase_a_energies,
    )

    logger_opt.write_xyz_trajectory(
        molecules=[rep.molecule for rep in representatives],
        filepath="trajectories/representatives.xyz",
        energies=None
    )    

    

    logger_opt.write_xyz_trajectory(
        molecules=phase_b_molecules,
        filepath="trajectories/phase_b_best.xyz",
        energies=phase_b_energies
    )

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
