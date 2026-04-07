import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

import sys
from pathlib import Path
import time
import matplotlib.pyplot as plt

module_dir = Path(__file__).parent / "modules"
sys.path.insert(0, str(module_dir))

# Module imports
from molecule_class import Molecule
from transformations import Transformation
from psi4_interface import Psi4Calculator, Config
from cluster import BHMCAnalyzer
from bhmc import MultiPhaseBHMC, BHMCConfig
from init import ClusterInitializer, InitializationConfig
from args import get_args
from logger import Logger
from interpolate import Interpolator

# Random visualization delete sometimes
from visualizations import plot_contour_vdw_radii, plot_pairwise_distance_heatmap, plot_initial_molecules_centers


if __name__ == "__main__":
    time_start = time.time()
    mp.set_start_method('spawn', force=True)

    N_WORKERS = 5

    # Parse the Arguments
    args = get_args()

    with open(args.i[0], "r") as file:
        xyz_content = file.read()
    
    # ── Main summary logger (cluster_gen.out) ───────────────────
    logger = Logger(name="main", log_file="cluster_gen.out", file_mode="w")
    logger.write_program_header()

    molecule = Molecule.from_xyz(xyz_content)

    # Plot contour
    plot_contour_vdw_radii(molecule, save_path="figures/contour_vdw_radii.png")
    plot_pairwise_distance_heatmap(molecule, save_path="figures/pairwise_distance_heatmap.png")



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
        box_scale_factor=1.8,
        min_distance=1.5,
        optimize_submolecules=True,
        verbose=False
    )
    initializer = ClusterInitializer(config=init_config, logger=init_logger)
    initial_molecules, submol_indices, simulation_box = initializer.initialize_from_xyz(
        args.i[0], 
        n_configurations=N_WORKERS
    )

    plot_initial_molecules_centers(initial_molecules, submol_indices, simulation_box, save_path="figures/initial_molecules_centers.png")

    # ── BHMC Phase A: Global Exploration ────────────────────────
    bhmc_logger = Logger(name="bhmc", log_file="cluster_gen.out", file_mode="a")

    bhmc_config = BHMCConfig(
        temperature=5000,
        method="hf",
        basis="cc-pvdz",
        backend="xtb",
        verbose=False,
        adaptive_operators=True,
        adaptive_box = True,
        box_update_interval=5, # Update box every 5 steps
        box_target_acceptance =0.6,
        box_acceptance_window = 0.05,
        box_growth_max = 1.15,
        box_max_scale = 2,
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
        n_structures_per_worker=300,
        n_processes=N_WORKERS
    )

    # write phase a trajectory
    bhmc_logger.write_xyz_trajectory(
        molecules=[candidate[0] for candidate in phase_a_candidates],
        filepath="trajectories/phase_a_candidates.xyz",
        energies=[candidate[1] for candidate in phase_a_candidates] 
    )

    

    # ── Analysis & Clustering ───────────────────────────────────
    logger_analysis = Logger(name="analysis", log_file="cluster_gen.out", file_mode="a")
    representatives = bhmc_sampler.analyse_phase_a_results(
        phase_a_candidates,
        submolecule_indices=submol_indices,
        cluster_method ="kmeans",
        n_clusters=5,
        simulation_box=simulation_box, 
        logger=logger_analysis
    )
    if len(representatives) >= 2:
        logger_analysis.header("Crossover Interpolation of Phase A Representatives")
        interpolator = Interpolator()
        all_crossover_mols = []

        rep_mols = [rep.molecule for rep in representatives]
        for i in range(len(rep_mols)):
            for j in range(i + 1, len(rep_mols)):
                logger_analysis.info(f"Interpolating between Representative {i} and {j}...")
                crossover_mols = interpolator.crossover_interpolate(
                    mol_a=rep_mols[i],
                    mol_b=rep_mols[j],
                    submolecule_indices=submol_indices,
                    n_permutations=2
                )
                all_crossover_mols.extend(crossover_mols)
        
        logger_analysis.info(f"Generated {len(all_crossover_mols)} crossover structures from {len(representatives)} representatives")

        logger_analysis.write_xyz_trajectory(
            molecules=all_crossover_mols,
            filepath="trajectories/phase_a_crossover_structures.xyz",
            energies=None
        )


    # ── BHMC Phase B: Local Refinement ──────────────────────────
    # Switch to lower temperature and disable adaptive scaling
    bhmc_sampler.config.temperature = 150
    bhmc_sampler.config.adaptive_operators = False

    phase_b_results = bhmc_sampler.run_phase_b(
        representatives=representatives,
        submolecule_indices=submol_indices,
        n_steps_per_worker=50,
    )

    # Extract best structures from Phase B
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

    # Psi4 Config für Optimierung (und spätere Interpolations-Punkte)
    config = Config(method="mp2", basis="cc-pvdz", memory="2GB")
    calc = Psi4Calculator(config=config, verbose=False)
    
    optimization_results = calc.batch_optimize_parallel_unordered(phase_b_molecules, n_processes=N_WORKERS)
    
    # Erfolgreiche Optimierungen extrahieren und nach Energie sortieren (wichtig für die Interpolation)
    successful_opts = [res for res in optimization_results if res.success]
    successful_opts.sort(key=lambda x: x.energy)  # Index 0 ist jetzt das absolute globale Minimum!
    
    optimized_mols = [result.molecule for result in successful_opts]
    optimized_energies = [result.energy for result in successful_opts]

    logger_opt.info(f"Optimized {len(optimized_mols)}/{len(phase_b_molecules)} structures successfully")

    # ── Write Intermediate Trajectories ─────────────────────────
    logger_opt.write_xyz_trajectory(
        molecules=[candidate[0] for candidate in phase_a_candidates],
        filepath="trajectories/phase_a_candidates.xyz",
        energies=[candidate[1] for candidate in phase_a_candidates] 
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

    # ── Interpolation of Energy Paths ───────────────────────────
    if len(optimized_mols) >= 2:
        logger_opt.header("Rigid-Body Interpolation between Opt. Minima")
        interpolator = Interpolator()
        
        # Wir berechnen Pfade vom tiefsten Minimum (Index 0) zu den anderen n Minima
        n_paths = min(len(optimized_mols) - 1, 3) # Limitiert auf max. 3 Pfade, um Zeit zu sparen
        
        # Für die Energiebewertungen der Pfade reicht ggf. Hartree-Fock, um Zeit zu sparen.
        # Falls du MP2 willst, nutze einfach den bestehenden 'calc' weiter.
        eval_config = Config(method="hf", basis="cc-pvdz", memory="2GB")
        eval_calc = Psi4Calculator(config=eval_config, verbose=False)
        
        for i in range(1, n_paths + 1):
            logger_opt.info(f"Interpolating path from Global Min (0) to Local Min ({i})...")
            
            path_mols, path_energies = interpolator.slerp_interpolate(
                mol_a=optimized_mols[0],
                mol_b=optimized_mols[i],
                submolecule_indices=submol_indices,
                n_points=20,  # 12 Schritte pro Pfad
                calculator=eval_calc
            )
            

            logger_opt.write_xyz_trajectory(
                molecules=path_mols,
                filepath=f"trajectories/interpolated_path_0_to_{i}.xyz",
                energies=path_energies
            )

            # Plot the energy profiles
            valid_indices = [idx for idx, energy in enumerate(path_energies) if energy is not None]
            valid_energies = [energy for energy in path_energies if energy is not None]
            if valid_energies:
                plt.figure(figsize=(8, 5))
                plt.plot(valid_indices, valid_energies, marker='o', linestyle='-', color='#1f77b4')
                plt.xlabel("Interpolation Step")
                plt.ylabel("Energy (Hartree)")
                plt.title(f"Energy Profile: Global Min to Local Min {i}")
                plt.grid(True, linestyle="--", alpha=0.7)
                plt.savefig(f"trajectories/energy_profile_0_to_{i}.png")
                plt.close()
                logger_opt.info(f"Saved energy profile plot for path 0 to {i}")
            else:
                logger_opt.warning(f"No valid energies for path 0 to {i}, skipping plot.")


        
    # ── Timing ──────────────────────────────────────────────────
    elapsed = time.time() - time_start
    logger_opt.separator(char="=")
    logger_opt.info(f"Total runtime: {elapsed:.1f}s ({elapsed/60:.1f} min)")
    logger_opt.separator(char="=")
