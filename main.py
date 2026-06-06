import sys
from pathlib import Path
import numpy as np

module_dir = Path(__file__).parent / "modules"
sys.path.insert(0, str(module_dir))

from init import ClusterInitializer, InitializationConfig
from bhmc import BHMC
from bhmc_config import BHMCConfig
from args import get_args
from logger import Logger
from struc_distinction import StructureAnalysis, StructureAnalysisConfig
from molecule_class import Molecule
from calculator import EnergyEvaluator
from input_reader import read_ml_bhop_input, read_operator_distribution
import time
import warnings

# Silence Warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="networkx")
warnings.filterwarnings("ignore", category=UserWarning, module="pyfiglet")

if __name__ == "__main__":
    args = get_args()

    logger = Logger(name="cluster_gen", log_file="cluster_gen.out", file_mode="w")
    logger.program_header()

    logger.separator() 
    if args.s is not None:
        settings_path = args.s[0]
        logger.info(f"Reading settings from: {settings_path}")
        try:
            settings = read_ml_bhop_input(settings_path)
            logger.info("Settings successfully read:")
            for section, params in settings.items():
                logger.info(f"  [{section}]")
                for key, value in params.items():
                    logger.info(f"    {key} = {value}")
        except Exception as e:
            logger.error(f"Error reading settings file: {e}")
            sys.exit(1)
        logger.separator()

    if args.o is not None:
        operator_dist_path = args.o[0]
        logger.info(f"Reading operator distribution from: {operator_dist_path}")
        try:
            operator_distribution = read_operator_distribution(operator_dist_path)
            logger.info("Operator distribution successfully read:")
            for name, weight in operator_distribution:
                logger.info(f"  {name}: {weight}")
        except Exception as e:
            logger.error(f"Error reading operator distribution file: {e}")
            sys.exit(1)
        logger.separator()
    

    # ------------------------------------------------------------------
    # Step 1: Initial SP Estimation
    # ------------------------------------------------------------------
    # Open input xyz
    with open(args.i[0], "r") as f:
        content = f.read()
    # Parse Molecule
    mol = Molecule.from_xyz(content)
    
    # Obtain settings for initial SP calculation
    settings_isp = settings.get("InitialSPConfig", {})
    calc_sp = EnergyEvaluator(backend=settings_isp["backend"], xtb_method=settings_isp["xtb_method"])
    timer_start = time.time()
    sp_energy = calc_sp.evaluate(mol)
    timer_end = time.time()
    logger.info(f"Initial SP Energy: {sp_energy:.4f} eV (Calculated in {timer_end - timer_start:.2f} seconds)")


    # ------------------------------------------------------------------
    # Step 2: Initialization
    # ------------------------------------------------------------------

    # obtain settings for initialization
    settings_init = settings.get("InitializationConfig", {})

    # ── Initialization ────────────────────────────────────────────────────────
    init_config = InitializationConfig(
        backend = settings_init["backend"],
        qm_method = settings_init["qm_method"],
        qm_basis = settings_init["qm_basis"],
        xtb_method = settings_init["xtb_method"],
        gpaw_mode = settings_init["gpaw_mode"],
        gpaw_xc = settings_init["gpaw_xc"],
        classifier_backend = settings_init["classifier_backend"],
        classifier_kwargs = settings_init.get("classifier_kwargs", {}),
        box_type = settings_init["box_type"],
        box_scale_factor = settings_init["box_scale_factor"],
        eta_factor = settings_init["eta_factor"],
        min_distance = settings_init["min_distance"],
        max_placement_attempts = settings_init["max_placement_attempts"],
        optimize_submolecules = settings_init["optimize_submolecules"],
        optimize_parallel = settings_init["optimize_parallel"],
        optimize_cluster_representatives = settings_init["optimize_cluster_representatives"],
        verbose = settings_init["verbose"],
    )

    
    initializer = ClusterInitializer(config=init_config, logger=logger)
    
    # obtain settings for xyz initialization
    settings_xyz_init = settings.get("InitializeXYZ", {})
    
    initial_molecules, submol_indices, simulation_box, umap_model, svm_model = initializer.initialize_from_xyz(
        args.i[0],
        n_workers=settings_xyz_init["n_workers"],
        n_configurations=settings_xyz_init["n_configurations"],
        n_sampling_workers=settings_xyz_init["n_sampling_workers"],
        placing_method=settings_xyz_init["placing_method"],
        grid_spacing = settings_xyz_init["grid_spacing"],
        n_theta = settings_xyz_init["n_theta"],
        n_phi = settings_xyz_init["n_phi"],
        n_r = settings_xyz_init["n_r"],
    )
     
    logger.write_xyz_trajectory(
        molecules=initial_molecules,
        filepath="trajectories/initial_candidates.xyz",
        energies=None,
    )

    # ---------------------------------------------------------------------
    # Step 3: Structure Analysis and Optimization of Initial Candidates
    # ---------------------------------------------------------------------

    settings_struc_analysis = settings.get("StructureAnalysisInt", {})

    # Initialize the StructureAnalysis Class
    StructureAnalysisConfig = StructureAnalysisConfig(
        calculator_backend=settings_struc_analysis["calculator_backend"],
        calculator_qm_method=settings_struc_analysis["calculator_qm_method"],
        calculator_qm_basis=settings_struc_analysis["calculator_qm_basis"],
        calculator_xtb_method=settings_struc_analysis["calculator_xtb_method"],
        calculator_gpaw_mode=settings_struc_analysis["calculator_gpaw_mode"],
        calculator_gpaw_basis=settings_struc_analysis["calculator_gpaw_basis"],
        calculator_gpaw_xc=settings_struc_analysis["calculator_gpaw_xc"],
        rmsd_threshold=settings_struc_analysis["rmsd_threshold"],
    )

    logger.info("\n")
    logger.header("Structure Analysis and Optimization of Initial Candidates")
    time_struc_analysis = time.time() 
    structure_analysis = StructureAnalysis(logger=logger, mols= initial_molecules, config=StructureAnalysisConfig)
    # Optimize the initial structures and write to xyz
    optimized_mols = structure_analysis.optimize_geometries(n_workers=20)
    logger.write_xyz_trajectory(
        molecules=optimized_mols,
        filepath="trajectories/optimized_initial_candidates.xyz",
        energies=None,
    )
    structure_analysis.plot_rmsd_comparison(save_path="figures/rmsd_comparison.png")
    unique_indices, unique_mols = structure_analysis.get_unique_structures()
    logger.info(f"Number of unique structures after optimization: {len(unique_mols)}")

    logger.write_xyz_trajectory(
        molecules=unique_mols,
        filepath="trajectories/unique_optimized_initial_candidates.xyz",
        energies=None,
    )

    time_struc_analysis_end = time.time()
    elapsed = time_struc_analysis_end - time_struc_analysis
    logger.info(f"Structure analysis and optimization completed in {elapsed:.2f} seconds")

    # ---------------------------------------------------------------------
    # Step 4: BHMC Sampling
    # ---------------------------------------------------------------------
    


    timer_start_bhmc = time.time()
    
    # Obtain settings for BHMC sampling
    settings_bhmc = settings.get("BHMCConfig", {})

    # If n_workers > unique_mols randomly sample out of unique mols to full up n_workers, otherwise use all unique mols
    if settings_bhmc["n_workers"] > len(unique_mols):
        logger.info(f"n_workers ({settings_bhmc['n_workers']}) is greater than the number of unique molecules ({len(unique_mols)}). Randomly sampling with replacement to fill up workers.")
        sampled_indices = np.random.choice(len(unique_mols), size=settings_bhmc["n_workers"], replace=True)
        # obtain the corresponding molecules, submolecule indices are the same
        unique_mols = [unique_mols[i] for i in sampled_indices] 
    else:
        logger.info(f"Using all {len(unique_mols)} unique molecules for BHMC sampling.")


    bhmc_config = BHMCConfig(
        backend=settings_bhmc["backend"],
        qm_method = settings_bhmc["qm_method"],
        qm_basis = settings_bhmc["qm_basis"],
        xtb_method = settings_bhmc["xtb_method"],
        gpaw_mode = settings_bhmc["gpaw_mode"],
        gpaw_basis = settings_bhmc["gpaw_basis"],
        gpaw_xc = settings_bhmc["gpaw_xc"],
        temperature = settings_bhmc["temperature"],
        verbose = settings_bhmc["verbose"],
        adaptive_operators = settings_bhmc["adaptive_operators"],
        adaptive_box = settings_bhmc["adaptive_box"],
        box_update_interval = settings_bhmc["box_update_interval"],
        box_target_acceptance = settings_bhmc["box_target_acceptance"],
        box_acceptance_window = settings_bhmc["box_acceptance_window"],
        box_growth_kp = settings_bhmc["box_growth_kp"],
        box_growth_max = settings_bhmc["box_growth_max"],
        box_max_scale = settings_bhmc["box_max_scale"],
        box_stable_windows = settings_bhmc["box_stable_windows"],
        operators = operator_distribution,  # Use the operator distribution read from the file
    )


    bhmc_sampler = BHMC(
        config=bhmc_config,
        simulation_box=simulation_box,
        logger=logger,
        worker_log_file="bhmc_workers.log",
    )

    accepted_structures = bhmc_sampler.run(
        initial_molecules=unique_mols,
        submolecule_indices=submol_indices,
        n_steps_per_worker=settings_bhmc["n_steps_per_worker"],
        n_processes=len(unique_mols),
    )

    logger.info(f"BHMC complete: {len(accepted_structures)} structures accepted")

    logger.write_xyz_trajectory(
        molecules=[m for m, _ in accepted_structures],
        filepath="trajectories/bhmc_structures.xyz",
        energies=[e for _, e in accepted_structures],
    )

    bhmc_sampler.analyze_results(umap_model=umap_model, classifier=svm_model)




    
