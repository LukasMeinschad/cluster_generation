"""Initialization of molecular clusters for BHMC sampling.

The module provides a complete initialization workflow:
1. Load structures from XYZ files
2. Fragment into connected submolecules
3. Optionally optimize each submolecule
4. Build a simulation box
5. Generate and rank candidate initial configurations
"""

import numpy as np
from typing import List, Tuple, Optional
from dataclasses import dataclass
import matplotlib.pyplot as plt


from modules.molecule_class import Molecule
from modules.box import SimulationBox
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# Import sobol sequence generator
from scipy.stats import qmc

# Logger import
from modules.logger import Logger

# Import geometry ops for rmsd
from modules.geometry import GeometryOps

# Import energy evaluator for prescreening
from modules.calculator import EnergyEvaluator

# Import the Featurizer
from modules.featurizer import Featurizer, FeaturizerConfig
from modules.cluster import Clustering



@dataclass
class InitializationConfig:
    """Configuration options for cluster initialization and prescreening."""
    
    # Configuration of energy evaluation for prescreening
    backend: str = "psi4" # QM based energy evaluation
    qm_method: str = "hf"              # QM method for optimization
    qm_basis: str = "sto-3g"           # Basis set for optimization
    xtb_method: str = "GFN2-xTB"       # Method for xTB energy evaluation (if applicable)
    gpaw_mode: str = "lcao"
    gpaw_basis: str = "dzp"
    gpaw_xc: str = "B3LYP"

    box_type: str = "sphere"        # "sphere" or "cube"
    box_scale_factor: float = 2.5   # Scaling factor for box size
    min_distance: float = 2.0       # Minimum distance between submolecules (Angstrom)
    max_placement_attempts: int = 1000  # Max attempts to place molecules
    optimize_submolecules: bool = True  # Whether to optimize submolecules 
    optimize_parallel: bool = True  # Wether ot optimize the submolecules in parallel
    verbose: bool = True



        
    def rank_candidates(
            self,
            molecules: List[Molecule],
            n_keep: int
        ) -> List[Tuple[float, Molecule]]:
        """  
        Evaluate all candidates and return lowest n_keep as (energy, molecule)
        """
        evaluator = EnergyEvaluator(
            backend=self.backend,
            qm_method=self.qm_method,
            qm_basis=self.qm_basis,
            xtb_method=self.xtb_method,
            gpaw_mode=self.gpaw_mode,
            gpaw_basis=self.gpaw_basis,
            gpaw_xc=self.gpaw_xc,
        )

        scored: List[Tuple[Optional[float], Molecule]] = []
        for mol in molecules:
            try:
                e = evaluator.evaluate(mol)
                scored.append((e, mol))
            except Exception:
                continue

        # Sort according to lowest energy and return top n_keep
        scored.sort(key=lambda x: x[0])
        return scored[:n_keep]

    


class ClusterInitializer:
    """Handles initialization of molecular clusters for BHMC sampling."""
    
    def __init__(self, config: InitializationConfig, logger: Optional[Logger] = None):
        """Create a cluster initializer with configuration and optional logger."""
        self.config = config
        self.logger = logger

        # Set the energy evaluator for prescreening 
        self.energy_evaluator = EnergyEvaluator(
            backend=self.config.backend,
            qm_method=self.config.qm_method,
            qm_basis=self.config.qm_basis,
            xtb_method=self.config.xtb_method,
            gpaw_mode=self.config.gpaw_mode,
            gpaw_basis=self.config.gpaw_basis,
            gpaw_xc=self.config.gpaw_xc   
        )

    def _log(self, msg: str, level: str = "info") -> None:
        """Log a message to stdout and to the configured logger."""
        if self.config.verbose:
            print(msg)
        if self.logger:
           getattr(self.logger, level)(msg) 

    def plot_energy_distribution_filtering(self, energies_before: List[float], energies_after: List[float], save_path: Optional[str] = None) -> None:
        """  
        Uses Seaborn to plot the energy distribution before and after filtering high-energy outliers. 
        """
        import seaborn as sns

        # side by side plot
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        sns.histplot(energies_before, bins=30, kde=False, ax=axes[0], color='blue')
        # KDE plot + Mean and STD lines
        sns.kdeplot(energies_before, ax=axes[0], color='blue', label='KDE')
        mean_before = np.mean(energies_before)
        std_before = np.std(energies_before)
        axes[0].axvline(mean_before, color='red', linestyle='--', label=f'Mean: {mean_before:.2f}')
        axes[0].axvline(mean_before + std_before, color='orange', linestyle='--', label=f'Mean + 1 STD: {mean_before + std_before:.2f}')
        axes[0].axvline(mean_before - std_before, color='orange', linestyle='--', label=f'Mean - 1 STD: {mean_before - std_before:.2f}')
        axes[0].set_title('Energy Distribution - Initialization (Before Z-Score Filtering)')
        axes[0].set_xlabel('Energy (Hartree)')
        axes[0].set_ylabel('Count')
        axes[0].legend()
        sns.histplot(energies_after, bins=30, kde=False, ax=axes[1], color='green')
        # KDE plot + Mean and STD lines
        sns.kdeplot(energies_after, ax=axes[1], color='green', label='KDE')
        mean_after = np.mean(energies_after)
        std_after = np.std(energies_after)
        axes[1].axvline(mean_after, color='red', linestyle='--', label=f'Mean: {mean_after:.2f}')
        axes[1].axvline(mean_after + std_after, color='orange', linestyle='--', label=f'Mean + 1 STD: {mean_after + std_after:.2f}')
        axes[1].axvline(mean_after - std_after, color='orange', linestyle='--', label=f'Mean - 1 STD: {mean_after - std_after:.2f}')
        axes[1].set_title('Energy Distribution - Initialization (After Z-Score Filtering)')
        axes[1].set_xlabel('Energy (Hartree)')
        axes[1].set_ylabel('Count')
        axes[1].legend()
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        plt.close()


    # Helper Method for RMSD Computation
    @staticmethod
    def _calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
        """Return RMSD after optimal atom correspondence matching."""
        coords_2_opt = GeometryOps.find_optimal_correspondence(coords1, coords2)
        diff = coords1 - coords_2_opt
        return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))
    
    def filter_high_energy_outliers(self, molecules: List[Molecule], energies: List[float], Z_threshold: float = 3.0) -> Tuple[List[Molecule], List[float]]:
        """  
        Peforms a filtering of unrealistic high-energy outliers based on the Z-score of the energy distribution.
        Args:
            molecules: List of candidate molecules
            energies: Corresponding list of energies for the candidates
            Z_threshold: Z-score threshold for filtering (default: 3.0)
        """
        energy_array = np.array(energies)
        mean_energy = np.mean(energy_array)
        std_energy = np.std(energy_array)
        z_scores = (energy_array - mean_energy) / std_energy
        filtered_molecules, filtered_energies = [], []
        # Filter only the high energy outliers (Z > Z_threshold)
        for mol, energy, z in zip(molecules, energies, z_scores):
            if z <= Z_threshold:
                filtered_molecules.append(mol)
                filtered_energies.append(energy)
        return filtered_molecules, filtered_energies
    
    

    def initialize_from_xyz(
        self,
        xyz_file: str,
        n_workers: int = 1,
        n_configurations: int = 10000,
        placing_method: str = "random",
        grid_spacing: Optional[str] = "sobol",
        n_theta: Optional[int] = None,
        n_phi: Optional[int] = None,
        n_r: Optional[int] = None,
        energy_backend: str= "xtb",
        energy_xtb_method: str = "GFN2-xTB",
    ) -> Tuple[List[Molecule], List[List[int]], SimulationBox]:
        """
        Run the full initialization workflow from an XYZ structure.

        Args:
            xyz_file: Path to XYZ file containing molecular structure
            n_workers: Number of configurations to return (e.g. one per BHMC worker)
            n_configurations: Number of candidates generated for prescreening
            placing_method: Placement method ("random", "sobol", or "grid")
            grid_spacing: Method to space points in the grid ("linear", "equal_volume_grid", or "sobol")
            n_theta: Number of angular theta grid points (grid mode)
            n_phi: Number of angular phi grid points (grid mode)
            n_r: Number of radial grid points (grid mode)
            energy_backend: Backend for prescreening energies ("xtb", "ase_emt", or "psi4")
            energy_xtb_method: xTB method if the backend uses xTB

        Returns:
            Tuple of (selected_molecules, submolecule_indices, simulation_box)
        """
        if self.logger:
            self.logger.header("Cluster Initialization")
            msg = (
                f"Configuration:\n"
                f"  QM Method (Optimization): {self.config.qm_method}\n"
                f"  Basis Set (Optimization): {self.config.qm_basis}\n"
                f"  Simulation Box Type: {self.config.box_type}\n"
                f"  Simulation Box Scale Factor: {self.config.box_scale_factor}\n"
                f"  Minimum Distance Between Submolecules: {self.config.min_distance} Angstrom\n"
                f"  Candidate configurations generated: {n_configurations}\n"
                f"  Workers (n_keep) for energy ranking: {n_workers}\n"
                f"  Energy Evaluation Backend: {energy_backend}\n"
                f"  xTB Method (if applicable): {energy_xtb_method}\n"
            )
            self.logger.info(msg)

        
        # Step 1: Load molecule
        molecule = self._load_molecule(xyz_file)
    
        # Step 2: Fragment into submolecules
        submolecules = self._fragment_molecule(molecule)
    
        # Step 2b: Store submolecule indices BEFORE optimization
        submol_indices = [submol.get_index_in_parent() for submol in submolecules]
    
        # Step 3: Optimize submolecules (optional)
        if self.config.optimize_submolecules:
            submolecules = self._optimize_submolecules(submolecules)
    
        # Step 4: Create simulation box
        simulation_box = self._create_simulation_box(submolecules)
    
        # Step 5: Generate random initial configuration
        initial_molecules = self._generate_random_configurations(
            submolecules = submolecules,
            simulation_box = simulation_box,
            n_configurations = n_configurations,
            placing_method = placing_method,
            grid_spacing = grid_spacing,
            n_theta = n_theta,
            n_phi = n_phi,
            n_r = n_r,
        )

        # Step 5b: Energy prescreening and ranking via the shared calculator wrapper
        self.energy_evaluator = EnergyEvaluator(
            backend=energy_backend,
            qm_method=self.config.qm_method,
            qm_basis=self.config.qm_basis,
            xtb_method=energy_xtb_method,
            gpaw_mode=self.config.gpaw_mode,
            gpaw_basis=self.config.gpaw_basis,
            gpaw_xc=self.config.gpaw_xc,
        )

        scored: List[Tuple[float, Molecule]] = []
        failed = 0
        for mol in initial_molecules:
            try:
                e = self.energy_evaluator.evaluate(mol)
                scored.append((e, mol))
            except Exception as exc:
                self._log(f"Energy evaluation failed for one candidate: {exc}", level="error")
                failed += 1
        self._log(f"\nEnergy evaluation complete: {len(scored)} successful, {failed} failed")
        
        if len(scored) == 0:
            raise RuntimeError("All energy evaluations failed. Cannot select initial configurations.")

        if len(scored) < n_workers:
            self._log(f"Warning: Only {len(scored)} valid configurations found, but n_workers={n_workers}. Returning all valid configurations.")
            n_workers = len(scored)
        
        scored.sort(key=lambda x: x[0])


        # Obtain unfiltered energies and mols
        mols_raw = [mol for _, mol in scored]
        energies_raw = [e for e, _ in scored]
        e_min, e_max = min(energies_raw), max(energies_raw)


        # Filter high energy outliers based on Z-score
        mols, energies = self.filter_high_energy_outliers(mols_raw, energies_raw, Z_threshold=3.0)

        # Plot energy distribution before and after filtering
        self.plot_energy_distribution_filtering(energies_before=energies_raw, energies_after=energies, save_path="figures/energy_distribution_filtering.png")


        # Log how many configurations were filtered out as high-energy outliers
        n_filtered = len(scored) - len(mols)
        self._log(f"Filtered out {n_filtered} high-energy outliers based on Z-score thresholding.")

        # Build feature matrix and cluster candidates with k-means
        featurizer_config = FeaturizerConfig(descriptor_type="SOAP")
        featurizer = Featurizer(featurizer_config)
        feature_mat = featurizer.build_feature_matrix(mols, energies=energies, include_hbonds=True)
        n_clusters = min(n_workers, len(mols))
        clustering = Clustering(feature_matrix=feature_mat, normalize=True)
        labels = clustering.kmeans_clustering(n_clusters=n_clusters)

        # Plot UMAP embedding (labels are stored on the clustering object after kmeans)
        embedding = clustering.umap(n_components=2)
        clustering.plot_embedding(embedding, title="Initialization UMAP", save_path="figures/initialization_umap.png")

        # Evaluate and log clustering quality metrics
        metrics = clustering.evaluate(labels)
        self._log(f"\nClustering evaluation metrics:")
        for metric_name, metric_value in metrics.items():
            self._log(f"  {metric_name}: {metric_value:.4f}")

        # Pick one representative per cluster → these become the MC starting structures
        rep_indices = clustering.get_cluster_representatives(labels=labels, method="centroid")
        selected_molecules = [mols[i] for i in rep_indices.values()]

        # Locally Optimize the selected representatives to ensure they are at least in a local minimum before starting BHMC
        for i, mol in enumerate(selected_molecules):
            self._log(f"Optimizing selected representative {i+1}/{len(selected_molecules)}...")
            try:
                optimized = self.energy_evaluator.optimize_geometry(mol, optimizer="BFGS", write_trajectory=False)
                selected_molecules[i] = optimized
                self._log(f" Optimization successful")
            except Exception as exc:
                self._log(f" Optimization failed: {exc}", level="error")
                # Keep original if optimization fails
        

        self._log(f"\nPre-screening completed: {len(scored)} valid configurations out of {n_configurations} generated.")
        self._log(f"Statistics of scored configurations:")
        self._log(f"  Energy range: {e_min:.6f} to {e_max:.6f} Hartree")
        self._log(f"  Energy distribution: mean {np.mean(energies):.6f}, median {np.median(energies):.6f}, std {np.std(energies):.6f} Hartree")
        self._log(f"  Clusters formed: {n_clusters}, Representatives selected: {len(selected_molecules)}")

        summary_lines = [
            f"{'='*60}",
            f"Initialization Complete!",
            f"Total atoms: {len(initial_molecules[0].coordinates)}",
            f"Submolecules: {len(submolecules)}",
            f"Simulation box type: {simulation_box.box_type}",
            f"Candidates generated: {n_configurations}, Valid: {len(scored)}, Failed: {failed}",
            f"Selected for workers: {len(selected_molecules)}, Energy range: {e_min:.6f} to {e_max:.6f} Hartree",
        ]
        if simulation_box.box_type == "sphere":
            summary_lines.append(f"Box radius: {simulation_box.radius:.2f} Angstrom")
            summary_lines.append(f"Box volume: {simulation_box.get_volume():.2f} Angstrom^3")
        else:
            summary_lines.append(f"Box dimensions: {simulation_box.box_dimensions}")
            summary_lines.append(f"Box volume: {simulation_box.get_volume():.2f} Angstrom^3")

        summary = "\n".join(summary_lines)
        self._log(summary)

        return selected_molecules, submol_indices, simulation_box

    def _generate_random_configurations(
            self,
            submolecules: List[Molecule],
            simulation_box: SimulationBox,
            n_configurations: int = 1,
            placing_method: str = "random",
            grid_spacing: Optional[str] = "sobol",
            n_theta: Optional[int] = None,
            n_phi: Optional[int] = None,
            n_r: Optional[int] = None,
        ) -> List[Molecule]:
        """Generate multiple initial configurations inside the simulation box.

        Args:
            submolecules: List of submolecules to place.
            simulation_box: Simulation box used for placement.
            n_configurations: Number of candidate configurations to generate.
            placing_method: Placement method ("random", "sobol", or "grid").
            grid_spacing: Grid spacing strategy for grid mode.
            n_theta: Number of angular theta grid points (grid mode).
            n_phi: Number of angular phi grid points (grid mode).
            n_r: Number of radial grid points (grid mode).

        Returns:
            List of generated candidate molecules.
        """
        if n_configurations <= 0:
            raise ValueError("n_configurations must be > 0")

        method = (placing_method or "random").lower().strip()
        if method not in {"random", "sobol", "grid"}:
            raise ValueError(f"Invalid placing method: {placing_method}. Choose from 'random', 'sobol', or 'grid'.")

        spacing = (grid_spacing or "sobol").lower().strip()

        if method == "grid":
            if spacing not in {"linear", "equal_volume_grid", "sobol"}:
                raise ValueError(f"Invalid grid spacing method: {grid_spacing}. Choose from 'linear', 'equal_volume_grid', or 'sobol'.")

            if n_theta is None or n_phi is None or n_r is None:
                n_submols = len(submolecules)
                n_per_submol = max(1, n_configurations // n_submols)
                n_theta = n_phi = n_r = int(np.ceil(n_per_submol ** (1/3)))
                self._log(f"Grid spacing selected. Setting n_theta = n_phi = n_r = {n_theta} to generate at least {n_configurations} configurations.")
            self._log(f"Generating grid-based configurations with spacing: {spacing}, n_theta: {n_theta}, n_phi: {n_phi}, n_r: {n_r}")
        else:
            spacing = n_theta = n_phi = n_r = None

        configurations: List[Molecule] = []
        for _ in range(n_configurations):
            mol = self._generate_configuration(
                submolecules=submolecules,
                simulation_box=simulation_box,
                placing_method=method,
                grid_spacing=spacing,
                n_theta=n_theta,
                n_phi=n_phi,
                n_r=n_r,
            )
            configurations.append(mol)
        return configurations
        


    def _load_molecule(self, xyz_file: str) -> Molecule:
        """Load molecule from XYZ file."""
        self._log(f"\n[1/5] Loading molecule from {xyz_file}") 

        
        with open(xyz_file, 'r') as f:
            xyz_content = f.read()
        
        molecule = Molecule.from_xyz(xyz_content)

        self._log(f" Loaded: {len(molecule.coordinates)} atoms")
        
        return molecule
    
    def _fragment_molecule(self, molecule: Molecule) -> List[Molecule]:
        """Fragment molecule into connected components."""
        self._log(f"\n[2/5] Fragmenting molecule into submolecules")
        

        molecule.compute_bonds()
        submolecules = molecule.fragment_by_connectivity()
        
        self._log(f" Found {len(submolecules)} submolecules:")
        for i, submol in enumerate(submolecules):
            self._log(f"  Submolecule {i+1}: {len(submol.coordinates)} atoms")
        

        
        return submolecules
    
    def _optimize_submolecules(self, submolecules: List[Molecule]) -> List[Molecule]:
        """Optimize all submolecules, optionally in parallel."""
        # If optimize parallel is enabled and we have more than 1 submol -> optimize in parallel
        if self.config.optimize_parallel and len(submolecules) > 1:
            return self._optimize_submolecules_parallel(submolecules)
        
        optimized = []
        for i, submol in enumerate(submolecules):
            self._log(f"  Optimizing submolecule {i+1} with {len(submol.coordinates)} atoms...")
            try:
                optimized_submol = self.energy_evaluator.optimize_geometry(submol, optimizer="BFGS", write_trajectory=False)
                optimized.append(optimized_submol)
                self._log(f"   Optimization successful")
            except Exception as exc:
                self._log(f"   Optimization failed: {exc}", level="error")
                optimized.append(submol)  # Fall back to original if optimization fails
        return optimized
        

    def _optimize_submolecules_parallel(self, submolecules: List[Molecule]) -> List[Molecule]:
        """Optimize each submolecule in parallel using a process pool."""

        n_workers = min(len(submolecules), multiprocessing.cpu_count())
        self._log(f"Workers: {n_workers}")

        optimized = [None] * len(submolecules)
    
        # Extract serializable data from submolecules
        submol_data = [
            {
                'atom_labels': submol.atom_labels,
                'coordinates': submol.coordinates,
                'index': i
            }
            for i, submol in enumerate(submolecules)
        ]
        ctx = multiprocessing.get_context("spawn")  # Use spawn to avoid issues with fork in some environments
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as executor:
            future_to_idx = {
                executor.submit(
                    self._optimize_single_static,
                    data['atom_labels'].tolist(),
                    np.asarray(data['coordinates'], dtype=np.float64),
                    self.config.backend,
                    self.config.qm_method,
                    self.config.qm_basis,
                    self.config.xtb_method,
                    self.config.gpaw_mode,
                    self.config.gpaw_basis,
                    self.config.gpaw_xc,
                    data['index']
                ): data['index']
                for data in submol_data
            }
            
            # Collect results
            for future in as_completed(future_to_idx):
                idx = future_to_idx[future]
                try:
                    result_mol = future.result()
                    optimized[idx] = result_mol
                    self._log(f"  Submolecule {idx+1} optimized successfully")
                except Exception as e:
                    self._log(f"  Submolecule {idx+1} optimization failed: {e}", level="error")
                    # Fall back to original submolecule
                    optimized[idx] = submolecules[idx]
    
        return optimized

    @staticmethod
    def _optimize_single_static(
        atom_labels: List[str],
        coordinates: np.ndarray,
        backend: str,
        method: str,
        basis: str,
        xtb_method: str,
        gpaw_mode: str,
        gpaw_basis: str,
        gpaw_xc: str,
        task_id: int
    ) -> Molecule:
        """Optimize one submolecule in an isolated worker process.
        
        Args:
            atom_labels: List of atom symbols
            coordinates: Atomic coordinates
            backend: Backend for calculator wrapper
            method: QM method
            basis: Basis set
        
        Returns:
            Optimized molecule.
        """
        import os
        import traceback 
        import tempfile
        from pathlib import Path 

        from modules.molecule_class import Molecule
        from modules.calculator import EnergyEvaluator

        # 1 Thread per worker to avoid oversubscription
        os.environ["OMP_NUM_THREADS"] = "1"
        os.environ["MKL_NUM_THREADS"] = "1"
        os.environ["OPENBLAS_NUM_THREADS"] = "1"

        # Unique scratch per task
        base = Path.cwd() / ".psi4_scratch"
        base.mkdir(parents=True, exist_ok=True)
        task_scratch = Path(tempfile.mkdtemp(prefix=f"task_{task_id}_", dir=str(base)))
        os.environ["PSI_SCRATCH"] = str(task_scratch)

        try:
            labels = list(atom_labels)
            coords = np.asarray(coordinates, dtype=np.float64).copy()

            submol = Molecule.from_labels_and_coords(labels, coords)
            calc = EnergyEvaluator(
                backend=backend,
                qm_method=method,
                qm_basis=basis,
                xtb_method=xtb_method,
                gpaw_mode=gpaw_mode,
                gpaw_basis=gpaw_basis,
                gpaw_xc=gpaw_xc,
            )
            optimized_submol = calc.optimize_geometry(submol, optimizer="BFGS", write_trajectory=False)
            return optimized_submol
        except Exception as exc:
            traceback.print_exc()
            raise RuntimeError(f"Optimization failed for task {task_id}: {exc}")
        finally:
            # Clean up scratch directory
            if task_scratch.exists():
                for item in task_scratch.iterdir():
                    if item.is_file():
                        item.unlink()
                    elif item.is_dir():
                        import shutil
                        shutil.rmtree(item)
                task_scratch.rmdir()
        
    
    
    def _create_simulation_box(self, submolecules: List[Molecule]) -> SimulationBox:
        """Create simulation box based on submolecules."""
        self._log(f"\n[4/5] Creating simulation box")
        
        # Collect all covalent radii
        all_cov_radii = []
        total_atoms = 0
        
        for submol in submolecules:
            all_cov_radii.extend(submol.covalent_radii)
            total_atoms += len(submol.coordinates)
        
        # Create box
        simulation_box = SimulationBox.from_covalent_radii(
            covalent_radii=all_cov_radii,
            n_atoms=total_atoms,
            box_type=self.config.box_type,
            scale_factor=self.config.box_scale_factor
        )
        self._log(f" Box type:  {simulation_box.box_type}")
        if simulation_box.box_type == "sphere":
            self._log(f" Radius: {simulation_box.radius:.2f} Angstrom")
            self._log(f" Volume: {simulation_box.get_volume():.2f} Angstrom^3")
        else:
            self._log(f" Dimensions: {simulation_box.box_dimensions}")
            self._log(f" Volume: {simulation_box.get_volume():.2f} Angstrom^3")
        
        return simulation_box
    
    def _generate_configuration(
        self, 
        submolecules: List[Molecule], 
        simulation_box: SimulationBox,
        placing_method: str = "random",
        grid_spacing: Optional[str] = None,
        n_theta: Optional[int] = None,
        n_phi: Optional[int] = None,
        n_r: Optional[int] = None,
    ) -> Molecule:
        """
        Generate random initial configuration inside simulation box.
        
        Args:
            submolecules: List of submolecules to place
            simulation_box: Simulation box to place molecules in
            placing_method: Method to place molecules:
                - "random" for random placement with distance checks
                - "sobol" for low-discrepancy sequence placement
                - "grid" for systematic grid-based placement (requires additional parameters)
            
                
        """
        
        if len(submolecules) == 0:
            raise ValueError("No submolecules to place")
        
        method = (placing_method or "random").lower().strip(    )
        
        if method not in {"random", "grid", "sobol"}:
            raise ValueError(f"Invalid placing method: {placing_method}. "
                             "Choose from 'random', 'grid', or 'sobol'.")

        sobol_engine = None
        sobol_engine_rotation = None
        if method == "sobol":
            sobol_engine = qmc.Sobol(d=3, scramble=True)
            sobol_engine_rotation = qmc.Sobol(d=3, scramble=True)
        
        partition_points = None
        partition_perm = None
        partition_ptr = None
        if method == "grid":
            if simulation_box.box_type != "sphere":
                raise ValueError("Grid-based placement is currently only implemented for spherical boxes.")
            
            spacing = (grid_spacing or "sobol").lower().strip()
            if spacing not in {"linear", "equal_volume_grid", "sobol"}:
                raise ValueError(f"Invalid grid spacing method: {grid_spacing}. Choose from 'linear', 'equal_volume_grid', or 'sobol'.")
            
            # Choose n_theta, n_phi, n_r based on config or defaults
            n_theta = 10 if n_theta is None else int(n_theta)
            n_phi = 10 if n_phi is None else int(n_phi)
            n_r = 5 if n_r is None else int(n_r)
            if n_theta <= 0 or n_phi <= 0 or n_r <= 0:
                raise ValueError("n_theta, n_phi, and n_r must be positive integers for grid placement.")
            
            partition_points = self.generate_partition_points(
                center=simulation_box.center,
                radius=simulation_box.radius,
                n_partitions=len(submolecules),
                n_theta=n_theta,
                n_phi=n_phi,
                n_r=n_r,
                spacing=spacing,
                sobol_scramble=True,
                sobol_seed=None
            )
            # Randomize traversal per partition ot avoid starting always at same points
            partition_perm = [np.random.permutation(len(points)) for points in partition_points]
            # Pointers to track next point in each partition
            partition_ptr = [0] * len(submolecules)

        placed_coords = np.empty((0, 3))
        placed_atoms = []

        for i, submol in enumerate(submolecules):
            if self.config.verbose:
                print(f" Placing submolecule {i+1}/{len(submolecules)} with {len(submol.coordinates)} atoms...")
            
            # Get submolecule COM
            submol_com = np.average(submol.coordinates, axis=0, weights=submol.masses)
            centered_coords = submol.coordinates - submol_com
            
            placed = False
            for _ in range(self.config.max_placement_attempts):
                
                if method == "random":
                    position = simulation_box.get_random_position(n_points=1)
                elif method == "sobol":
                    u = sobol_engine.random(n=1)[0]  # Get one point

                    if simulation_box.box_type == "cube":
                        half_dims = simulation_box.box_dimensions / 2.0
                        position = simulation_box.center + (2.0 * u - 1.0) * half_dims

                    elif simulation_box.box_type == "sphere":
                        # Uniform in volume:
                        # r = R * u1^(1/3), cos(theta) = 2*u2 - 1, phi = 2*pi*u3
                        r = simulation_box.radius * (u[0] ** (1.0/3.0))
                        cos_theta = 2.0 * u[1] - 1.0
                        sin_theta = np.sqrt(max(0.0, 1.0 - cos_theta * cos_theta)) 
                        phi = 2.0 * np.pi * u[2]

                        direction = np.array([
                            sin_theta * np.cos(phi),
                            sin_theta * np.sin(phi),
                            cos_theta
                        ])  
                        position = simulation_box.center + r * direction
                    else:
                        raise ValueError(f"Unsupported box type for Sobol placement: {simulation_box.box_type}")
                
                # Grid based mode
                else:

                    # Get next point from the current partition's list
                    perm = partition_perm[i]
                    ptr = partition_ptr[i]
                    if ptr >= len(perm):
                        break

                    p_idx = perm[ptr]
                    partition_ptr[i] = ptr + 1
                    position = partition_points[i][p_idx]

                # Sample Random Rotation
                if method == "sobol":
                    rotation_matrix = self._sobol_rotation_matrix(sobol_engine_rotation)
                else:
                    rotation_matrix = self._random_rotation_matrix()
            
                rotated_coords = (rotation_matrix @ centered_coords.T).T
                new_coords = rotated_coords + position
                
                # Check distance constraints
                if i == 0 or self._check_min_distance(new_coords, placed_coords):
                    placed_coords = np.vstack([placed_coords, new_coords])                     
                    placed_atoms.extend(submol.atom_labels.tolist())
                    placed = True
                    break
            
            if not placed:
                raise RuntimeError(
                    f"Could not place submolecule {i+1} after "
                    f"{self.config.max_placement_attempts} attempts. "
                    "Try increasing box size or reducing min_distance."
                )
        
        initial_molecule = Molecule.from_labels_and_coords(
            atom_labels=placed_atoms,
            coordinates=placed_coords.copy(),  
        )
        
        return initial_molecule
    
    def _check_min_distance(
        self, 
        new_coords: np.ndarray, 
        existing_coords: np.ndarray 
    ) -> bool:
        """Check if new coordinates satisfy minimum distance constraint."""
        if len(existing_coords) == 0:
            return True
        
        # existing_coords is already a numpy array, no conversion needed
        diff = new_coords[:, np.newaxis, :] - existing_coords[np.newaxis, :, :]
        distances = np.linalg.norm(diff, axis=2)
        return np.all(distances >= self.config.min_distance)



    @staticmethod
    def _random_rotation_matrix() -> np.ndarray:
        """Generate a random 3D rotation matrix from a unit quaternion."""
        u = np.random.rand(3)
        q = np.array([
            np.sqrt(1 - u[0]) * np.sin(2 * np.pi * u[1]),
            np.sqrt(1 - u[0]) * np.cos(2 * np.pi * u[1]),
            np.sqrt(u[0]) * np.sin(2 * np.pi * u[2]),
            np.sqrt(u[0]) * np.cos(2 * np.pi * u[2])
        ])
        q0, q1, q2, q3 = q 
        return np.array([
            [1 - 2*(q2**2 + q3**2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
            [2*(q1*q2 + q0*q3), 1 - 2*(q1**2 + q3**2), 2*(q2*q3 - q0*q1)],
            [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1**2 + q2**2)]
        ])
    
    def _sobol_rotation_matrix(self, sobol_engine) -> np.ndarray:
        """Generate a quasi-random 3D rotation matrix from Sobol samples."""
        u = sobol_engine.random(n=1)[0]

        # Shoemake's algorithm like above
        q = np.array([
            np.sqrt(1 - u[0]) * np.sin(2 * np.pi * u[1]),
            np.sqrt(1 - u[0]) * np.cos(2 * np.pi * u[1]),
            np.sqrt(u[0]) * np.sin(2 * np.pi * u[2]),
            np.sqrt(u[0]) * np.cos(2 * np.pi * u[2])
        ])
        q0, q1, q2, q3 = q
        return np.array([
            [1 - 2*(q2**2 + q3**2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
            [2*(q1*q2 + q0*q3), 1 - 2*(q1**2 + q3**2), 2*(q2*q3 - q0*q1)],
            [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1**2 + q2**2)]
        ])

    @staticmethod
    def partition_sphere(center, radius, n):
        """Split azimuth angle [0, 2*pi) into `n` equal partitions.

        Args:
            center: Sphere center (unused, kept for interface compatibility).
            radius: Sphere radius (unused, kept for interface compatibility).
            n: Number of partitions.

        Returns:
            List of (theta_start, theta_end) tuples.
        """
        if n <= 0:
            raise ValueError("n must be a positive integer")

        partitions = []

        # Calculate the angle step
        angle_step = (2*np.pi) / n

        for i in range(n):
            start_angle = i * angle_step
            end_angle = (i + 1) * angle_step
            partitions.append((start_angle, end_angle))
        return partitions            

    @staticmethod
    def generate_partition_points(center,
                                  radius,
                                  n_partitions,
                                  n_theta = 10,
                                  n_phi = 10,
                                  n_r = 5,
                                  spacing: Optional[str] = "linear",
                                  sobol_scramble: bool = True,
                                  sobol_seed: Optional[int] = None):
        """Generate candidate points for each partition of a spherical box.

        Args:
            center: Center of the sphere (3D coordinates)
            radius: Radius of the sphere
            n_partitions: Number of partitions (submolecules)
            n_theta: Number of points in theta direction per partition
            n_phi: Number of points in phi direction
            n_r: Number of points in radial direction

        Spacing methods:
            - "linear": using a linear grid in (r, phi, theta)
            - "equal_volume_grid" a deterministic grid closer to uniform 3D density
            - "sobol" using quasi-random sampling uniformly in volume per partition        
        """
        center = np.asarray(center, dtype=np.float64)
        spacing = (spacing or "linear").lower().strip()

        if radius <= 0:
            raise ValueError("Radius must be positive")
        if n_partitions <= 0:
            raise ValueError("n_partitions must be > 0") 

        slices = ClusterInitializer.partition_sphere(center, radius, n_partitions) 
        all_partition_points = []

        if spacing not in {"linear", "equal_volume_grid", "sobol"}:
            raise ValueError(f"Unknown spacing method: {spacing}. Choose from 'linear', 'equal_volume_grid', or 'sobol'.")
    
        for p_idx, (theta_start, theta_end) in enumerate(slices):
            if spacing == "sobol":
                # Keep the point count consistent with linear grid
                n_points = max(1, n_theta * n_phi * n_r)
                
                # Generate engine per partition
                seed = None if sobol_seed is None else sobol_seed + p_idx
                engine = qmc.Sobol(d=3, scramble=sobol_scramble, seed=seed)
                u = engine.random(n=n_points)

                # Uniform sample in the wedge volume:
                # r^3 ~ U(0, R^3), cosphi = U(-1,1) theta ~ U(theta0, theta1)
                r = radius * np.cbrt(u[:, 0])  
                cos_phi = 2.0 * u[:, 1] - 1.0
                sin_phi = np.sqrt(np.clip(1.0 - cos_phi * cos_phi, 0.0, 1.0))
                theta = theta_start + (theta_end - theta_start) * u[:, 2]

                # Transform to cartesian
                x = center[0] + r * sin_phi * np.cos(theta) 
                y = center[1] + r * sin_phi * np.sin(theta)
                z = center[2] + r * cos_phi
                points = np.column_stack((x, y, z))
                all_partition_points.append(points)
                continue
            
            theta_values = np.linspace(theta_start, theta_end, n_theta, endpoint=False)

            # grid based methods
            if spacing == "linear":
                phi_values = np.linspace(0, np.pi, n_phi)
                r_values = np.linspace(0, radius, n_r)
            elif spacing == "equal_volume_grid":
                # Midpoints to avoid piling exactly at poles and center
                j = np.arange(n_phi)
                k = np.arange(n_r)

                cos_phi_values = -1.0 + 2.0 * (j + 0.5) / n_phi
                phi_values = np.arccos(np.clip(cos_phi_values, -1.0, 1.0))

                # r^3 uniformly
                r_values = radius * np.cbrt((k + 0.5) / n_r)

            partition_points = []
            for r in r_values:
                for phi in phi_values:
                    sin_phi = np.sin(phi)
                    cos_phi = np.cos(phi)
                    for theta in theta_values:
                        x = center[0] + r * sin_phi * np.cos(theta)
                        y = center[1] + r * sin_phi * np.sin(theta)
                        z = center[2] + r * cos_phi
                        partition_points.append([x, y, z])
            all_partition_points.append(np.array(partition_points))
        return all_partition_points
    

            
    



    
# =============================== Debugging and Testing ===============================

def test_initializer(xyz_file: str,
                     method: str = "hf",
                     basis: str = "cc-pvdz",
                     optimize: bool = True,
                     box_type: str = "sphere",
                     box_scale: float = 2.5,
                     min_distance: float = 2.0,
                     placing_method: str = "random",
                     n_configurations: int = 1,
                     save_output: bool = True):
    """Run an end-to-end initialization sanity check and print diagnostics.
    
    Args:
        xyz_file: Path to XYZ file
        method: QM method for optimization
        basis: Basis set for optimization
        optimize: Whether to optimize submolecules
        box_type: Type of simulation box ("sphere" or "cube")
        box_scale: Scaling factor for box size
        min_distance: Minimum distance between submolecules (Angstrom)
        placing_method: Placement method ("random", "sobol", or "grid")
        n_configurations: Number of candidate configurations to generate
        save_output: Whether to write generated configurations to disk
    """
    print(f"\n{'='*80}")
    print(f"Testing ClusterInitializer with {xyz_file}")
    print(f"{'='*80}\n")
    print(f"Method: {method}, Basis: {basis}, Optimize: {optimize}")
    print(f"Box type: {box_type}, Box scale: {box_scale}, Min distance: {min_distance} Angstrom")
    print(f"Placement method: {placing_method}")

    config = InitializationConfig(
        qm_method=method,
        qm_basis=basis,
        optimize_submolecules=optimize,
        box_type=box_type,
        box_scale_factor=box_scale,
        min_distance=min_distance,
        verbose=True
    )

    initializer = ClusterInitializer(config=config)
    try:
        initial_molecules, submol_indices, simulation_box = initializer.initialize_from_xyz(
            xyz_file,
            n_configurations=n_configurations,
            placing_method=placing_method,
        )
        
        # Use first configuration for diagnostics
        initial_molecule = initial_molecules[0]
        
        print(f"\nInitialization successful!")
        print(f"Generated {len(initial_molecules)} initial configuration(s)")
        print(f"Molecule Information:")
        print(f"  Total atoms: {len(initial_molecule.coordinates)}")
        print(f" Total mass: {np.sum(initial_molecule.masses):.2f} amu")
        print(f" Center of mass: {np.average(initial_molecule.coordinates, axis=0, weights=initial_molecule.masses)}")
        print(f"\n Submolecule Information:")
        print(f"  Number of submolecules: {len(submol_indices)}")
        for i, indices in enumerate(submol_indices):
            print(f"    Submolecule {i+1}: {len(indices)} atoms")
        print(f"\n Simulation Box Information:")
        print(f"  Box type: {simulation_box.box_type}")
        if simulation_box.box_type == "sphere":
            print(f"  Box radius: {simulation_box.radius:.2f} Angstrom")
            print(f"  Box volume: {simulation_box.get_volume():.2f} Angstrom^3")
        else:
            print(f"  Box dimensions: {simulation_box.box_dimensions}")
            print(f"  Box volume: {simulation_box.get_volume():.2f} Angstrom^3")
        all_inside = simulation_box.is_inside(initial_molecule.coordinates)
        print(f"  All atoms inside box: {'Yes' if all_inside else 'No'}")

        # Calculate pairwise distances 
        print(f"\n Distance Check:")
        min_dist = np.inf 
        max_dist = 0.0 
        coords = initial_molecule.coordinates
        for i in range(len(coords)):
            for j in range(i+1, len(coords)):
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < min_dist:
                    min_dist = dist
                if dist > max_dist:
                    max_dist = dist
        print(f"  Minimum interatomic distance: {min_dist:.2f} Angstrom")
        print(f"  Maximum interatomic distance: {max_dist:.2f} Angstrom")

        # Calculate center of mass distance between submolecules
        print(f"\n Submolecule Distance Check:")
        if len(submol_indices) > 1:
            coms = []
            for indices in submol_indices:
                sub_coords = initial_molecule.coordinates[indices]
                sub_masses = initial_molecule.masses[indices]
                com = np.average(sub_coords, axis=0, weights=sub_masses)
                coms.append(com)
            for i in range(len(coms)):
                for j in range(i+1, len(coms)):
                    dist = np.linalg.norm(coms[i] - coms[j])
                    print(f"  Distance between submolecule {i+1} and {j+1}: {dist:.2f} Angstrom")
        else:
            print(f"  Only one submolecule, skipping distance check.")
        # Save output files
        if save_output:
            # Write all initial configuration into 1 xyz file
            with open("initial_configurations.xyz", "w") as f:
                for idx, mol in enumerate(initial_molecules):
                    f.write(f"{len(mol.coordinates)}\n")
                    f.write(f"Configuration {idx+1}\n")
                    for label, coord in zip(mol.atom_labels, mol.coordinates):
                        f.write(f"{label} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")
            print(f"\nInitial configurations saved to initial_configurations.xyz") 
        
    except Exception as e:
        print(f"\nInitialization failed with error: {e}")

def test_submolecule_partition_mapping(
        initial_molecules: List[Molecule],
        submol_indices: List[List[int]],
        simulation_box: SimulationBox,
        atol: float = 1e-10,
        verbose: bool = True
        ) -> bool:
        """Verify that each submolecule COM is in its assigned sphere partition.

        Returns:
            True if the mapping is correct, False otherwise
        """
        if simulation_box.box_type != "sphere":
            raise ValueError("Submolecule partition mapping test is only applicable for spherical boxes.")
        
        if len(submol_indices) == 0:
            raise ValueError("No submolecules to test")
        
        center = np.asarray(simulation_box.center, dtype=np.float64)
        n_submols = len(submol_indices)
        partitions = ClusterInitializer.partition_sphere(center, simulation_box.radius, n_submols)

        failures = []

        for cfg_idx, mol in enumerate(initial_molecules):
            for i, atom_idx in enumerate(submol_indices):
                idx = np.asarray(atom_idx, dtype=int)
                sub_coords = mol.coordinates[idx]
                submol_masses = mol.masses[idx]

                # Calculate com
                com = np.average(sub_coords, axis=0, weights=submol_masses)
                rel = com - center

                # Azimuth angle in [0, 2pi)]
                theta = np.arctan2(rel[1], rel[0])
                if theta < 0:
                    theta += 2 * np.pi

                start_angle, end_angle = partitions[i]

                # Partitions are contigous [start, end) we allow tiny tolerance
                # for the last partition include the right edge ~2*pi
                if i == n_submols - 1:
                    ok = (theta >= start_angle - atol) and (theta <= end_angle + atol)
                else:
                    ok = (theta >= start_angle - atol) and (theta < end_angle + atol)

                if not ok:
                    failures.append((cfg_idx, i, theta, start_angle, end_angle))
                
        if verbose:
            if not failures:
                print(f"\nSubmolecule partition mapping test passed for all {len(initial_molecules)} configurations.")
            else:
                print(f"\nSubmolecule partition mapping test failed for {len(failures)} cases:")
                for cfg_idx, sub_idx, theta, start, end in failures:
                    print(f"  Configuration {cfg_idx+1}, Submolecule {sub_idx+1}: theta={theta:.4f} not in [{start:.4f}, {end:.4f})")
        return len(failures) == 0

if __name__ == "__main__":
    """  
    Run test when module is executed directly
    """
    xyz_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/co2_h2o.xyz"
    


    init_config = InitializationConfig(
        backend="xtb",
        qm_method="hf",
        qm_basis="cc-pvdz",
        optimize_submolecules=True,
        optimize_parallel=True,
        box_type="sphere",
        box_scale_factor=2.5,
        min_distance=2.0,
        verbose=True,
        max_placement_attempts=1000
    ) 
    initializer = ClusterInitializer(config=init_config)
    initial_molecules, submol_indices, simulation_box = initializer.initialize_from_xyz(
        xyz_file=xyz_file,
        n_workers=4,
        n_configurations=5000,
        placing_method="sobol",
    )
    

