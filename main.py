import sys
from pathlib import Path

module_dir = Path(__file__).parent / "modules"
sys.path.insert(0, str(module_dir))

from init import ClusterInitializer, InitializationConfig
from bhmc import MultiPhaseBHMC
from bhmc_config import BHMCConfig
from args import get_args
from logger import Logger
from struc_distinction import StructureAnalysis


if __name__ == "__main__":
    args = get_args()

    logger = Logger(name="cluster_gen", log_file="cluster_gen.out", file_mode="w")

    # ── Initialization ────────────────────────────────────────────────────────
    init_config = InitializationConfig(
        backend="xtb",
        box_scale_factor=1.5,
        xtb_method="GFN2-xTB",
        verbose=False,
    )

    initializer = ClusterInitializer(config=init_config, logger=logger)
    initial_molecules, submol_indices, simulation_box = initializer.initialize_from_xyz(
        args.i[0],
        n_workers=10,
        n_configurations=10000,
        placing_method="sobol",
        energy_backend="xtb",
        energy_xtb_method="GFN2-xTB",
    )

    logger.info(f"Obtained {len(initial_molecules)} initial structures")
    logger.write_xyz_trajectory(
        molecules=initial_molecules,
        filepath="trajectories/initial_candidates.xyz",
        energies=None,
    )

    # Initialize the StructureAnalysis Class
    structure_analysis = StructureAnalysis(logger=logger, mols=initial_molecules)
    structure_analysis.compute_pairwise_rmsd()
    structure_analysis.plot_pairwise_rmsd_heatmap(save_path="figures/initial_pairwise_rmsd_heatmap.png")

    # ── Phase A: Global PES Exploration ──────────────────────────────────────
    bhmc_config = BHMCConfig(
        backend="xtb",
        xtb_method="GFN2-xTB",
        qm_method="hf",
        qm_basis="sto-3g",
        temperature=400.0,
        verbose=True,
        adaptive_operators=True,
        adaptive_box=True,
        box_update_interval=10,
        box_target_acceptance=0.6,
        box_acceptance_window=0.05,
        box_growth_max=1.15,
        box_max_scale=4.0,
        box_stable_windows=3,
    )

    bhmc_sampler = MultiPhaseBHMC(
        config=bhmc_config,
        simulation_box=simulation_box,
        logger=logger,
        worker_log_file="bhmc_workers.log",
    )

    phase_a_structures = bhmc_sampler.run_phase_a(
        initial_molecules=initial_molecules,
        submolecule_indices=submol_indices,
        n_structures_per_worker=500,
        n_processes=len(initial_molecules),
    )

    logger.info(f"Phase A complete: {len(phase_a_structures)} structures accepted")

    # Save all Phase A structures
    phase_a_molecules = [mol for mol, _ in phase_a_structures]
    phase_a_energies = [e for _, e in phase_a_structures]
    logger.write_xyz_trajectory(
        molecules=phase_a_molecules,
        filepath="trajectories/phase_a_structures.xyz",
        energies=phase_a_energies,
    )

    
