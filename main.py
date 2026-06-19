import sys
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

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
from box import SimulationBox
from cluster import Clustering
from args import get_args
from logger import Logger
from struc_analysis import StructureAnalysis, StructureAnalysisConfig
from molecule_class import Molecule
from calculator import EnergyEvaluator
from input_reader import read_ml_bhop_input, read_operator_distribution


"""
Main.py Loop
This script condenses all of the necessary functions for the algorithm to run
In the following the respective methods for the __main__ are defined
"""



# ============================================================================
# Setup: CLI arguments, settings file, operator distribution
# ============================================================================

def require_input_files(args, logger: Logger) -> None:
    """Exit early with a clear message if a required CLI argument is missing.

    Without this, a missing --s/--o silently skips the block that defines
    `settings`/`operator_dict`, surfacing as a confusing NameError far later.
    """
    missing = [flag for flag, value in (("--i", args.i), ("--s", args.s), ("--o", args.o)) if value is None]
    if missing:
        logger.error(f"Missing required argument(s): {', '.join(missing)}. See --help.")
        sys.exit(1)


def load_settings(settings_path: str, logger: Logger) -> Dict[str, Any]:
    """Read and log the ML-BHOP settings file, exiting on failure."""
    logger.info(f"Reading settings from: {settings_path}")
    try:
        settings = read_ml_bhop_input(settings_path)
    except Exception as e:
        logger.error(f"Error reading settings file: {e}")
        sys.exit(1)

    logger.info("Settings successfully read:")
    for section, params in settings.items():
        logger.info(f"  [{section}]")
        for key, value in params.items():
            logger.info(f"    {key} = {value}")
    logger.separator()
    return settings


def load_operator_distribution(
    operator_dist_path: str, logger: Logger
) -> Tuple[List[Tuple[str, float]], List[Tuple[str, float]]]:
    """Read the operator distribution file and split it into training/sampling operator pools."""
    logger.info(f"Reading operator distribution from: {operator_dist_path}")
    try:
        operator_dict = read_operator_distribution(operator_dist_path)
    except Exception as e:
        logger.error(f"Error reading operator distribution file: {e}")
        sys.exit(1)

    logger.info("Operator distribution successfully read:")
    for phase, operators in operator_dict.items():
        logger.info(f"  [{phase}]")
        for op, weight in operators.items():
            logger.info(f"    {op} = {weight}")
    logger.separator()

    o_nl_train = operator_dict.get("NL-OpsTraining", {})
    o_l_train = operator_dict.get("L-OpsTraining", {})
    o_train = list({**o_nl_train, **o_l_train}.items())

    o_nl_sampling = operator_dict.get("NL-OPsMain", {})
    o_l_sampling = operator_dict.get("L-OPsMain", {})
    o_sampling = list({**o_nl_sampling, **o_l_sampling}.items())
    return o_train, o_sampling


# ============================================================================
# Step 1: Initial single-point energy estimate
# ============================================================================

def run_initial_single_point(xyz_path: str, settings_isp: Dict[str, Any], logger: Logger) -> float:
    """Quick single-point energy of the raw input structure, mainly as a sanity check."""
    try:
        with open(xyz_path, "r") as f:
            content = f.read()
    except OSError as e:
        logger.error(f"Could not read input XYZ file '{xyz_path}': {e}")
        sys.exit(1)
    mol = Molecule.from_xyz(content)

    calc_sp = EnergyEvaluator(backend=settings_isp["backend"], xtb_method=settings_isp["xtb_method"])
    start = time.time()
    sp_energy = calc_sp.evaluate(mol)
    logger.info(f"Initial SP Energy: {sp_energy:.4f} eV (Calculated in {time.time() - start:.2f} seconds)")
    return sp_energy


# ============================================================================
# Step 2: Cluster initialization (placement, prescreening, clustering)
# ============================================================================

def run_initialization(
    xyz_path: str,
    settings_init: Dict[str, Any],
    settings_xyz_init: Dict[str, Any],
    logger: Logger,
) -> Tuple[List[Molecule], List[List[int]], SimulationBox, Clustering, np.ndarray]:
    """Build the candidate cluster pool and reduce it to one representative per cluster."""
    init_config = InitializationConfig(
        backend=settings_init["backend"],
        qm_method=settings_init["qm_method"],
        qm_basis=settings_init["qm_basis"],
        xtb_method=settings_init["xtb_method"],
        gpaw_mode=settings_init["gpaw_mode"],
        gpaw_xc=settings_init["gpaw_xc"],
        classifier_backend=settings_init["classifier_backend"],
        classifier_kwargs=settings_init.get("classifier_kwargs", {}),
        box_type=settings_init["box_type"],
        box_scale_factor=settings_init["box_scale_factor"],
        eta_factor=settings_init["eta_factor"],
        padding_factor=settings_init["padding_factor"],
        min_distance=settings_init["min_distance"],
        max_placement_attempts=settings_init["max_placement_attempts"],
        optimize_submolecules=settings_init["optimize_submolecules"],
        optimize_parallel=settings_init["optimize_parallel"],
        optimize_cluster_representatives=settings_init["optimize_cluster_representatives"],
        filter_com_outliers=settings_init["filter_com_outliers"],
        com_threshold=settings_init["com_threshold"],
        verbose=settings_init["verbose"],
    )
    initializer = ClusterInitializer(config=init_config, logger=logger)

    initial_molecules, submol_indices, simulation_box, init_clustering, feature_mat_raw = initializer.initialize_from_xyz(
        xyz_path,
        n_workers=settings_xyz_init["n_workers"],
        n_configurations=settings_xyz_init["n_configurations"],
        n_sampling_workers=settings_xyz_init["n_sampling_workers"],
        placing_method=settings_xyz_init["placing_method"],
        grid_spacing=settings_xyz_init["grid_spacing"],
        n_theta=settings_xyz_init["n_theta"],
        n_phi=settings_xyz_init["n_phi"],
        n_r=settings_xyz_init["n_r"],
    )
    logger.write_xyz_trajectory(
        molecules=initial_molecules,
        filepath="trajectories/initial_candidates.xyz",
        energies=None,
    )
    return initial_molecules, submol_indices, simulation_box, init_clustering, feature_mat_raw


# ============================================================================
# Step 3: Structure analysis & optimization of the initial candidate pool
# ============================================================================

def build_structure_analysis_config(settings_struc_analysis: Dict[str, Any]) -> StructureAnalysisConfig:
    return StructureAnalysisConfig(
        calculator_backend=settings_struc_analysis["calculator_backend"],
        calculator_qm_method=settings_struc_analysis["calculator_qm_method"],
        calculator_qm_basis=settings_struc_analysis["calculator_qm_basis"],
        calculator_xtb_method=settings_struc_analysis["calculator_xtb_method"],
        calculator_gpaw_mode=settings_struc_analysis["calculator_gpaw_mode"],
        calculator_gpaw_basis=settings_struc_analysis["calculator_gpaw_basis"],
        calculator_gpaw_xc=settings_struc_analysis["calculator_gpaw_xc"],
        rmsd_threshold=settings_struc_analysis["rmsd_threshold"],
    )


def analyze_initial_candidates(
    initial_molecules: List[Molecule],
    struc_analysis_config: StructureAnalysisConfig,
    logger: Logger,
) -> Tuple[List[Molecule], List[float]]:
    """Optimize, plot, and de-duplicate the raw initial candidate pool."""
    logger.info("\n")
    logger.header("Structure Analysis and Optimization of Initial Candidates")
    start = time.time()

    structure_analysis = StructureAnalysis(logger=logger, mols=initial_molecules, config=struc_analysis_config)
    optimized_mols, optimized_energies = structure_analysis.optimize_geometries(n_workers=20)
    logger.write_xyz_trajectory(
        molecules=optimized_mols,
        filepath="trajectories/optimized_initial_candidates.xyz",
        energies=None,
    )

    structure_analysis.plot_rmsd_comparison(save_path="figures/rmsd_comparison.png")
    structure_analysis.plot_optimized_energy_distribution(save_path="figures/optimized_energy_distribution.png")
    structure_analysis.plot_eigenspectra_distance_heatmap(
        metric="euclidean", use_optimized=True, save_path="figures/eigenspectra_distance_heatmap.png"
    )
    structure_analysis.plot_distance_matrix_heatmap(
        metric="euclidean", use_optimized=True, save_path="figures/distance_matrix_heatmap.png"
    )

    unique_indices, unique_mols = structure_analysis.get_unique_structures(use_optimized=True)
    unique_energies = [optimized_energies[i] for i in unique_indices]
    logger.info(f"Number of unique structures after optimization: {len(unique_mols)}")

    logger.write_xyz_trajectory(
        molecules=unique_mols,
        filepath="trajectories/unique_optimized_initial_candidates.xyz",
        energies=None,
    )
    logger.info(f"Structure analysis and optimization completed in {time.time() - start:.2f} seconds")
    return unique_mols, unique_energies


# ============================================================================
# Step 4: BHMC sampling (training phase + main sampling phase)
# ============================================================================

def build_bhmc_config(settings_bhmc: Dict[str, Any], operators: List[Tuple[str, float]]) -> BHMCConfig:
    """Build a BHMCConfig from settings with the given operator pool swapped in.

    Shared by both the training-phase and main-sampling-phase runs so the two
    configs can never drift out of sync except for their operator weights.
    """
    return BHMCConfig(
        backend=settings_bhmc["backend"],
        qm_method=settings_bhmc["qm_method"],
        qm_basis=settings_bhmc["qm_basis"],
        xtb_method=settings_bhmc["xtb_method"],
        gpaw_mode=settings_bhmc["gpaw_mode"],
        gpaw_basis=settings_bhmc["gpaw_basis"],
        gpaw_xc=settings_bhmc["gpaw_xc"],
        temperature=settings_bhmc["temperature"],
        verbose=settings_bhmc["verbose"],
        adaptive_operators=settings_bhmc["adaptive_operators"],
        adaptive_box=settings_bhmc["adaptive_box"],
        box_update_interval=settings_bhmc["box_update_interval"],
        box_target_acceptance=settings_bhmc["box_target_acceptance"],
        box_acceptance_window=settings_bhmc["box_acceptance_window"],
        box_growth_kp=settings_bhmc["box_growth_kp"],
        box_growth_max=settings_bhmc["box_growth_max"],
        box_max_scale=settings_bhmc["box_max_scale"],
        box_stable_windows=settings_bhmc["box_stable_windows"],
        operators=operators,
        clustering_method=settings_bhmc["clustering_method"],
        classifier_backend=settings_bhmc["classifier_backend"],
        classifier_kwargs=settings_bhmc.get("classifier_kwargs", None),
    )


def select_bhmc_workers(
    unique_mols: List[Molecule], unique_energies: List[float], n_workers: int, logger: Logger
) -> Tuple[List[Molecule], List[float]]:
    """Resize the unique-structure pool to exactly n_workers, sampling with replacement if too few."""
    if n_workers > len(unique_mols):
        logger.info(
            f"n_workers ({n_workers}) is greater than the number of unique molecules "
            f"({len(unique_mols)}). Randomly sampling with replacement to fill up workers."
        )
        sampled_indices = np.random.choice(len(unique_mols), size=n_workers, replace=True)
        return [unique_mols[i] for i in sampled_indices], [unique_energies[i] for i in sampled_indices]
    logger.info(f"Using all {len(unique_mols)} unique molecules for BHMC sampling.")
    return unique_mols, unique_energies


def run_bhmc_training_phase(
    unique_mols: List[Molecule],
    submol_indices: List[List[int]],
    simulation_box: SimulationBox,
    init_clustering: Clustering,
    feature_mat_raw: np.ndarray,
    settings_bhmc: Dict[str, Any],
    o_train: List[Tuple[str, float]],
    n_training_steps: int,
    logger: Logger,
) -> Tuple[Optional[Clustering], Optional[List[Molecule]]]:
    """Run the exploratory (training) phase of BHMC and retrain the reference model on what it finds."""
    bhmc_config = build_bhmc_config(settings_bhmc, o_train)
    bhmc_sampler = BHMC(config=bhmc_config, simulation_box=simulation_box, logger=logger, worker_log_file="bhmc_workers.log")

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
    return updated_clustering, training_representatives


def analyze_training_phase_structures(
    training_representatives: Optional[List[Molecule]],
    struc_analysis_config: StructureAnalysisConfig,
    unique_mols: List[Molecule],
    unique_energies: List[float],
    logger: Logger,
) -> Tuple[Optional[Dict[str, Any]], Optional[List[Molecule]], Optional[List[float]]]:
    """Optimize and PES-compare the BHMC training-phase representatives against the init pool.

    training_representatives is None whenever analyze_training_results found no
    structures novel enough to justify retraining the reference model (or had
    nothing to analyze at all) — in that case there is nothing new to optimize
    or compare, so this step is skipped instead of running on an empty pool.

    Returns (comparison, unique_mols_bhmc, unique_energies_bhmc) — the latter two
    are the optimized, deduplicated training-phase representatives, forwarded so
    the final run-wide aggregation step doesn't need to redo this work.
    """
    logger.header("Structure Analysis of Representatives from BHMC Training Phase")

    if not training_representatives:
        logger.info(
            "No novel structures were found during the BHMC training phase — "
            "skipping optimization and PES comparison for this step."
        )
        return None, None, None

    structure_analysis_bhmc = StructureAnalysis(logger=logger, mols=training_representatives, config=struc_analysis_config)
    optimized_representatives, optimized_energies = structure_analysis_bhmc.optimize_geometries(n_workers=20)
    logger.write_xyz_trajectory(
        molecules=optimized_representatives,
        filepath="trajectories/optimized_bhmc_representatives_training.xyz",
        energies=None,
    )
    structure_analysis_bhmc.plot_rmsd_comparison(save_path="figures/rmsd_comparison_bhmc_training.png")

    unique_indices_bhmc, unique_mols_bhmc = structure_analysis_bhmc.get_unique_structures(use_optimized=True)
    unique_energies_bhmc = [optimized_energies[i] for i in unique_indices_bhmc]
    logger.info(f"Number of unique structures after BHMC training phase: {len(unique_mols_bhmc)}")

    comparison = structure_analysis_bhmc.compare_pes(
        unique_mols, unique_energies,
        unique_mols_bhmc, unique_energies_bhmc,
        label_a="init", label_b="bhmc_training",
    )
    return comparison, unique_mols_bhmc, unique_energies_bhmc


def run_bhmc_sampling_phase(
    unique_mols: List[Molecule],
    submol_indices: List[List[int]],
    simulation_box: SimulationBox,
    settings_bhmc: Dict[str, Any],
    o_sampling: List[Tuple[str, float]],
    n_sampling_steps: int,
    init_clustering: Clustering,
    updated_clustering: Optional[Clustering],
    logger: Logger,
):
    """Run the main BHMC sampling phase, against the retrained model if the training phase produced one."""
    logger.info(
        f"Continuing BHMC sampling with {n_sampling_steps} sampling steps per worker "
        f"(total steps: {settings_bhmc['n_steps_per_worker']})"
    )
    bhmc_config_sampling = build_bhmc_config(settings_bhmc, o_sampling)
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

    if updated_clustering is not None:
        sampling_reference_clustering = updated_clustering
        logger.info("Using the retrained model from the training phase for sampling analysis.")
    else:
        sampling_reference_clustering = init_clustering
        logger.info("No improved model from the training phase — using the init-run model for sampling analysis.")

    sampling_results = bhmc_sampler_sampling.analyze_sampling_results(
        reference_clustering=sampling_reference_clustering,
    )
    return accepted_structures_sampling, sampling_results


def analyze_final_results(
    struc_analysis_config: StructureAnalysisConfig,
    unique_mols: List[Molecule],
    unique_energies: List[float],
    training_mols: Optional[List[Molecule]],
    training_energies: Optional[List[float]],
    accepted_structures_sampling: List[Tuple[Molecule, float]],
    logger: Logger,
) -> Optional[Dict[str, Any]]:
    """Optimize the BHMC sampling output, combine it with the init and training pools,
    PES-compare everything new against the init baseline, then sort the full final
    pool by energy and write it to one xyz trajectory.

    training_mols/training_energies are None whenever the training phase found
    nothing novel (see analyze_training_phase_structures).
    """
    logger.header("Final Aggregation Across the Full Run")

    new_mols: List[Molecule] = list(training_mols) if training_mols else []
    new_energies: List[float] = list(training_energies) if training_energies else []

    if accepted_structures_sampling:
        sampling_analysis = StructureAnalysis(
            logger=logger,
            mols=[m for m, _ in accepted_structures_sampling],
            config=struc_analysis_config,
        )
        _, optimized_sampling_energies = sampling_analysis.optimize_geometries(n_workers=20)
        sampling_unique_indices, sampling_unique_mols = sampling_analysis.get_unique_structures(use_optimized=True)
        sampling_unique_energies = [optimized_sampling_energies[i] for i in sampling_unique_indices]
        logger.info(f"Number of unique structures after BHMC sampling phase: {len(sampling_unique_mols)}")
        new_mols += sampling_unique_mols
        new_energies += sampling_unique_energies

    final_comparison = None
    if new_mols:
        comparator = StructureAnalysis(logger=logger, config=struc_analysis_config)
        final_comparison = comparator.compare_pes(
            unique_mols, unique_energies,
            new_mols, new_energies,
            label_a="init", label_b="full_run",
        )
    else:
        logger.info("No new structures were obtained from the training or sampling phases.")

    final_pool_analysis = StructureAnalysis(
        logger=logger, mols=unique_mols + new_mols, config=struc_analysis_config
    )
    final_unique_indices, final_unique_mols = final_pool_analysis.get_unique_structures(use_optimized=False)
    final_unique_energies = [(unique_energies + new_energies)[i] for i in final_unique_indices]
    logger.info(f"Final unique structure pool across the whole run: {len(final_unique_mols)}")

    logger.write_xyz_trajectory(
        molecules=final_unique_mols,
        filepath="trajectories/final_results.xyz",
        energies=final_unique_energies,
        sort_by_energy=True,
    )
    return final_comparison



# ============================================================================
# Main Function 
# ============================================================================

def main() -> None:
    args = get_args()

    logger = Logger(name="cluster_gen", log_file=args.f[0], include_timestamp=False, include_level=False)
    logger.program_header()
    logger.separator()

    require_input_files(args, logger)

    settings = load_settings(args.s[0], logger)
    o_train, o_sampling = load_operator_distribution(args.o[0], logger)

    # ------------------------------------------------------------------
    # Step 1: Initial single-point estimation
    # ------------------------------------------------------------------
    settings_isp = settings.get("InitialSPConfig", {})
    run_initial_single_point(args.i[0], settings_isp, logger)

    # ------------------------------------------------------------------
    # Step 2: Initialization
    # ------------------------------------------------------------------
    settings_init = settings.get("InitializationConfig", {})
    settings_xyz_init = settings.get("InitializeXYZ", {})
    initial_molecules, submol_indices, simulation_box, init_clustering, feature_mat_raw = run_initialization(
        args.i[0], settings_init, settings_xyz_init, logger,
    )

    # ------------------------------------------------------------------
    # Step 3: Structure analysis & optimization of initial candidates
    # ------------------------------------------------------------------
    settings_struc_analysis = settings.get("StructureAnalysisInt", {})
    struc_analysis_config = build_structure_analysis_config(settings_struc_analysis)
    unique_mols, unique_energies = analyze_initial_candidates(initial_molecules, struc_analysis_config, logger)

    # ------------------------------------------------------------------
    # Step 4: BHMC sampling — training phase, then main sampling phase
    # ------------------------------------------------------------------
    settings_bhmc = settings.get("BHMCConfig", {})
    unique_mols, unique_energies = select_bhmc_workers(unique_mols, unique_energies, settings_bhmc["n_workers"], logger)

    # Start with a higher chance of picking non-local operators, then retrain the
    # operator distribution on what gets found. Training is `training_frac` of total steps.
    n_training_steps = int(settings_bhmc["n_steps_per_worker"] * settings_bhmc["training_frac"])
    n_sampling_steps = settings_bhmc["n_steps_per_worker"] - n_training_steps
    logger.info(
        f"Starting BHMC sampling with {n_training_steps} training steps per worker "
        f"(total steps: {settings_bhmc['n_steps_per_worker']})"
    )

    timer_start_bhmc = time.time()
    updated_clustering, training_representatives = run_bhmc_training_phase(
        unique_mols, submol_indices, simulation_box, init_clustering, feature_mat_raw,
        settings_bhmc, o_train, n_training_steps, logger,
    )
    if updated_clustering is not None:
        init_clustering = updated_clustering
    logger.info(f"BHMC training phase completed in {time.time() - timer_start_bhmc:.2f} seconds")

    _, training_mols_bhmc, training_energies_bhmc = analyze_training_phase_structures(
        training_representatives, struc_analysis_config, unique_mols, unique_energies, logger
    )

    accepted_structures_sampling, _ = run_bhmc_sampling_phase(
        unique_mols, submol_indices, simulation_box,
        settings_bhmc, o_sampling, n_sampling_steps,
        init_clustering, updated_clustering, logger,
    )

    analyze_final_results(
        struc_analysis_config, unique_mols, unique_energies,
        training_mols_bhmc, training_energies_bhmc,
        accepted_structures_sampling, logger,
    )


if __name__ == "__main__":
    main()
