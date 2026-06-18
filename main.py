import sys
import warnings
from pathlib import Path
import numpy as np

# Silence warnings before importing modules that trigger them at import time
# (logger.py imports pyfiglet/pkg_resources, molecule_class.py imports networkx)
warnings.filterwarnings("ignore", category=RuntimeWarning, module="networkx")
warnings.filterwarnings("ignore", category=UserWarning, module="pyfiglet")
warnings.filterwarnings("ignore", category=UserWarning, module="pkg_resources")

module_dir = Path(__file__).parent / "modules"
sys.path.insert(0, str(module_dir))

from init import ClusterInitializer, InitializationConfig
from bhmc import BHMC
from bhmc_config import BHMCConfig
from args import get_args
from logger import Logger
from struc_analysis import StructureAnalysis, StructureAnalysisConfig
from molecule_class import Molecule
from calculator import EnergyEvaluator
from input_reader import read_ml_bhop_input, read_operator_distribution
import time

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
            operator_dict = read_operator_distribution(operator_dist_path)
            logger.info("Operator distribution successfully read:")
            for phase, operators in operator_dict.items():
                logger.info(f"  [{phase}]")
                for op, weight in operators.items():
                    logger.info(f"    {op} = {weight}")
        except Exception as e:
            logger.error(f"Error reading operator distribution file: {e}")
            sys.exit(1)
        logger.separator()

    # Get the operators for training
    o_nl_train = operator_dict.get("NL-OpsTraining", {})
    o_l_train = operator_dict.get("L-OpsTraining", {})
    o_train = list({**o_nl_train, **o_l_train}.items())  # Combine into a list of (name, weight) tuples

    # Get the operators for sampling
    o_nl_sampling = operator_dict.get("NL-OPsMain", {})
    o_l_sampling = operator_dict.get("L-OPsMain", {})
    o_sampling = list({**o_nl_sampling, **o_l_sampling}.items()) # Combine into a list of (name, weight) tuples


    

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
        padding_factor = settings_init["padding_factor"],
        min_distance = settings_init["min_distance"],
        max_placement_attempts = settings_init["max_placement_attempts"],
        optimize_submolecules = settings_init["optimize_submolecules"],
        optimize_parallel = settings_init["optimize_parallel"],
        optimize_cluster_representatives = settings_init["optimize_cluster_representatives"],
        filter_com_outliers = settings_init["filter_com_outliers"],
        com_threshold = settings_init["com_threshold"],
        verbose = settings_init["verbose"],
    )

    
    initializer = ClusterInitializer(config=init_config, logger=logger)
    
    # obtain settings for xyz initialization
    settings_xyz_init = settings.get("InitializeXYZ", {})
    
    initial_molecules, submol_indices, simulation_box, init_clustering, feature_mat_raw = initializer.initialize_from_xyz(
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
    optimized_mols, optimized_energies = structure_analysis.optimize_geometries(n_workers=20)
    logger.write_xyz_trajectory(
        molecules=optimized_mols,
        filepath="trajectories/optimized_initial_candidates.xyz",
        energies=None,
    )
    # Plot RMSD comparison of initial vs optimized structures
    structure_analysis.plot_rmsd_comparison(save_path="figures/rmsd_comparison.png")
    # Plot distribution of optimized energies across configurations
    structure_analysis.plot_optimized_energy_distribution(save_path="figures/optimized_energy_distribution.png")
    # Compute Eigenspectra of Coloumb matrix
    structure_analysis.plot_eigenspectra_distance_heatmap(metric="euclidean", use_optimized=True, save_path="figures/eigenspectra_distance_heatmap.png")
    # Compute Distance Matrix of Coloumb Matrix
    structure_analysis.plot_distance_matrix_heatmap(metric="euclidean", use_optimized=True, save_path="figures/distance_matrix_heatmap.png")
    # Cluster the optimized structures and obtain unique representatives
    unique_indices, unique_mols = structure_analysis.get_unique_structures(use_optimized=True)
    unique_energies = [optimized_energies[i] for i in unique_indices]
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

    # 4.a --> Start of with a higher change of picking non-local operators
    # Then perform prediction and retraining of the operator distribution
    # This are 25 % of total steps
    n_training_steps = int(settings_bhmc["n_steps_per_worker"] * settings_bhmc["training_frac"])
    logger.info(f"Starting BHMC sampling with {n_training_steps} training steps per worker (total steps: {settings_bhmc['n_steps_per_worker']})")
    
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
        operators = o_train,
        clustering_method = settings_bhmc["clustering_method"],
        classifier_backend = settings_bhmc["classifier_backend"],
        classifier_kwargs = settings_bhmc.get("classifier_kwargs", None),
    )


    bhmc_sampler = BHMC(
        config=bhmc_config,
        simulation_box=simulation_box,
        logger=logger,
        worker_log_file="bhmc_workers.log",
    )
    # Run the training phase of BHMC sampling
    accepted_structures = bhmc_sampler.run(
        initial_molecules=unique_mols,
        submolecule_indices=submol_indices,
        n_steps_per_worker=n_training_steps,
        n_processes=len(unique_mols),
    )
    logger.info(f"BHMC complete: {len(accepted_structures)} structures accepted")
    logger.write_xyz_trajectory(
        molecules=[m for m, _ in accepted_structures],
        filepath="trajectories/bhmc_structures.xyz",
        energies=[e for _, e in accepted_structures],
    )
    updated_clustering, training_representatives = bhmc_sampler.analyze_training_results(
        reference_clustering=init_clustering,
        feature_mat_init=feature_mat_raw,
        mols_init=init_clustering.raw_molecules,
    )
    if updated_clustering is not None:
        init_clustering = updated_clustering

    timer_end_bhmc = time.time()
    elapsed_bhmc = timer_end_bhmc - timer_start_bhmc
    logger.info(f"BHMC training phase completed in {elapsed_bhmc:.2f} seconds")

    # 4.b -> Structure Analysis of the new representatives
    logger.header("Structure Analysis of Representatives from BHMC Training Phase")
    structure_analysis_bhmc = StructureAnalysis(logger=logger, mols=training_representatives, config=StructureAnalysisConfig)
    # Optimize the Representatives Compare RMSD and write to xyz
    optimized_representatives_training, optimized_energies_training = structure_analysis_bhmc.optimize_geometries(n_workers=20)
    logger.write_xyz_trajectory(
        molecules=optimized_representatives_training,
        filepath="trajectories/optimized_bhmc_representatives_training.xyz",
        energies=None,
    )
    structure_analysis_bhmc.plot_rmsd_comparison(save_path="figures/rmsd_comparison_bhmc_training.png")
    # Get new unique representatives from the training phase
    unique_indices_bhmc, unique_mols_bhmc = structure_analysis_bhmc.get_unique_structures(use_optimized=True)
    unique_energies_bhmc = [optimized_energies_training[i] for i in unique_indices_bhmc]
    logger.info(f"Number of unique structures after BHMC training phase: {len(unique_mols_bhmc)}")

    # Compare the init pool against the BHMC training representatives (energy + structure)
    pes_comparison_training = structure_analysis_bhmc.compare_pes(
        unique_mols, unique_energies,
        unique_mols_bhmc, unique_energies_bhmc,
        label_a="init", label_b="bhmc_training",
    )
    



    # 4.b -> Continue Sampling with an updated operator distribution
    n_sampling_steps = settings_bhmc["n_steps_per_worker"] - n_training_steps
    
    logger.info(f"Continuing BHMC sampling with {n_sampling_steps} sampling steps per worker (total steps: {settings_bhmc['n_steps_per_worker']})")
    
    # Setup new config
    bhmc_config_sampling = BHMCConfig(
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
        operators = o_sampling,
    )
    bhmc_sampler_sampling = BHMC(
        config=bhmc_config_sampling,
        simulation_box=simulation_box,
        logger=logger,
        worker_log_file="bhmc_sampling_workers.log",
    )
    accepted_structures_sampling = bhmc_sampler_sampling.run(
        initial_molecules=unique_mols,
        submolecule_indices=submol_indices,
        n_steps_per_worker=n_sampling_steps,
        n_processes=len(unique_mols),
    )
    logger.info(f"BHMC sampling complete: {len(accepted_structures_sampling)} structures accepted")
    logger.write_xyz_trajectory(
        molecules=[m for m, _ in accepted_structures_sampling],
        filepath="trajectories/bhmc_sampling_structures.xyz",
        energies=[e for _, e in accepted_structures_sampling],
    )

    # Decide which model to analyze the sampling run against: the retrained model
    # from the training phase if it found an improvement, otherwise the original
    # init-run model.
    if updated_clustering is not None:
        sampling_reference_clustering = updated_clustering
        logger.info("Using the retrained model from the training phase for sampling analysis.")
    else:
        sampling_reference_clustering = init_clustering
        logger.info("No improved model from the training phase — using the init-run model for sampling analysis.")

    novel_indices_sampling, mols_sampling, feature_mat_sampling, predicted_labels_sampling, predicted_probabilities_sampling = bhmc_sampler_sampling.analyze_sampling_results(
        reference_clustering=sampling_reference_clustering,
    )

    