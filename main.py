import sys
from pathlib import Path
import numpy as np

module_dir = Path(__file__).parent / "modules"
sys.path.insert(0, str(module_dir))

from init import ClusterInitializer, InitializationConfig
from bhmc import MultiPhaseBHMC
from bhmc_config import BHMCConfig
from args import get_args
from logger import Logger
from struc_distinction import StructureAnalysis, StructureAnalysisConfig


if __name__ == "__main__":
    args = get_args()

    logger = Logger(name="cluster_gen", log_file="cluster_gen.out", file_mode="w")

    # ── Initialization ────────────────────────────────────────────────────────
    init_config = InitializationConfig(
        backend="xtb",
        box_scale_factor=1.5,
        xtb_method="GFN2-xTB",
        optimize_cluster_representatives=False, 
        verbose=False,
    )

    initializer = ClusterInitializer(config=init_config, logger=logger)
    initial_molecules, submol_indices, simulation_box, umap_model, knn_model = initializer.initialize_from_xyz(
        args.i[0],
        n_workers=10,
        n_configurations=10,
        n_sampling_workers=10,
        placing_method="sobol",
        energy_backend="xtb",
        energy_xtb_method="GFN2-xTB",
    )
    print(f"Obtained {len(initial_molecules)} initial structures")


    logger.info(f"Obtained {len(initial_molecules)} initial structures")
    logger.write_xyz_trajectory(
        molecules=initial_molecules,
        filepath="trajectories/initial_candidates.xyz",
        energies=None,
    )

    # Initialize the StructureAnalysis Class

    StructureAnalysisConfig = StructureAnalysisConfig(
        calculator_backend="psi4",
        calculator_qm_method="mp2",
        calculator_qm_basis="cc-pvdz",
        calculator_xtb_method="GFN2-xTB",
        calculator_gpaw_mode="lcao",
        calculator_gpaw_basis="dzp",
        calculator_gpaw_xc="B3LYP",
        rmsd_threshold=0.5,
    )


   
   # structure_analysis = StructureAnalysis(logger=logger, mols=initial_molecules, config=StructureAnalysisConfig)
   # structure_analysis.compute_pairwise_rmsd()
   # structure_analysis.plot_pairwise_rmsd_heatmap(save_path="figures/initial_pairwise_rmsd_heatmap.png")
   # structure_analysis.analyze_hessians(n_workers=10)


    # ── Phase A: Global PES Exploration ──────────────────────────────────────
    bhmc_config = BHMCConfig(
        backend="xtb",
        xtb_method="GFN2-xTB",
        qm_method="hf",
        qm_basis="sto-3g",
        temperature=400.0,
        verbose=False,
        adaptive_operators=True,
        adaptive_box=True,
        box_update_interval=10,
        box_target_acceptance=0.6,
        box_acceptance_window=0.05,
        box_growth_max=1.15,
        box_max_scale=2.0,
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
        n_structures_per_worker=8000,
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

    # Run the phase A analysis
    bhmc_sampler.analyze_phase_a_results(umap_model=umap_model, knn_model=knn_model, phase_a_structures=phase_a_structures)




    
