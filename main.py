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
from molecule_class import Molecule
from calculator import EnergyEvaluator
import time
import warnings

# Silence Warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="networkx")
warnings.filterwarnings("ignore", category=UserWarning, module="pyfiglet")

if __name__ == "__main__":
    args = get_args()

    logger = Logger(name="cluster_gen", log_file="cluster_gen.out", file_mode="w")
    logger.program_header()




    # ------------------------------------------------------------------
    # Step 1: Initial SP Estimation
    # ------------------------------------------------------------------
    # Open input xyz
    with open(args.i[0], "r") as f:
        content = f.read()
    # Parse Molecule
    mol = Molecule.from_xyz(content)
    
    calc_sp = EnergyEvaluator(backend="xtb", xtb_method="GFN2-xTB")
    timer_start = time.time()
    sp_energy = calc_sp.evaluate(mol)
    timer_end = time.time()
    logger.info(f"Initial SP Energy: {sp_energy:.4f} eV (Calculated in {timer_end - timer_start:.2f} seconds)")


    # ------------------------------------------------------------------
    # Step 2: Initialization
    # ------------------------------------------------------------------

    # ── Initialization ────────────────────────────────────────────────────────
    init_config = InitializationConfig(
        backend="xtb",
        box_scale_factor=1.0,
        xtb_method="GFN2-xTB",
        optimize_cluster_representatives=False, 
        verbose=False,
    )

    initializer = ClusterInitializer(config=init_config, logger=logger)
    initial_molecules, submol_indices, simulation_box, umap_model, svm_model = initializer.initialize_from_xyz(
        args.i[0],
        n_workers=20,
        n_configurations=1000,
        n_sampling_workers=20,
        placing_method="sobol",
        energy_backend="xtb",
        energy_xtb_method="GFN2-xTB",
        classifier_backend="svm",
    )
     
    logger.write_xyz_trajectory(
        molecules=initial_molecules,
        filepath="trajectories/initial_candidates.xyz",
        energies=None,
    )

    # Initialize the StructureAnalysis Class

    StructureAnalysisConfig = StructureAnalysisConfig(
        calculator_backend="xtb",
        calculator_qm_method="mp2",
        calculator_qm_basis="cc-pvdz",
        calculator_xtb_method="GFN2-xTB",
        calculator_gpaw_mode="lcao",
        calculator_gpaw_basis="dzp",
        calculator_gpaw_xc="B3LYP",
        rmsd_threshold=0.5,
    )


   
    structure_analysis = StructureAnalysis(logger=logger, mols=initial_molecules, config=StructureAnalysisConfig)
    structure_analysis.plot_pairwise_rmsd_heatmap()
    
    
    # Optimize the initial structures and write to xyz
    optimized_mols = structure_analysis.optimize_geometries(n_workers=20)
    logger.write_xyz_trajectory(
        molecules=optimized_mols,
        filepath="trajectories/optimized_initial_candidates.xyz",
        energies=None,
    )

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
        box_max_scale=5.0,
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
        n_structures_per_worker=1000,
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
    bhmc_sampler.analyze_phase_a_results(umap_model=umap_model, classifier=svm_model, phase_a_structures=phase_a_structures)




    
