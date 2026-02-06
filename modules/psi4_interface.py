import numpy as np 
from multiprocessing import Pool, cpu_count
from typing import List, Optional, Tuple, Dict, Any, Callable
from dataclasses import dataclass, field
from enum import Enum 
import tempfile 
import json 
import pickle
import os 
import time 
import shutil 
from tqdm import tqdm
from pathlib import Path 
from contextlib import contextmanager
import matplotlib.pyplot as plt
from molecule_class import Molecule


# =========================== DATA CLASSES ========================== 

# Enums here are used to define constant named values for calculation types 
class CalculationType(Enum):
    """  
    Different Types of Quantum Chemistry Calculations
    """
    SINGLE_POINT = "single_point"
    GEOMETRY_OPTIMIZATION = "geometry_optimization"
    FREQUENCY = "frequency"
    OPT_FREQ = "opt_freq"
    PROPERTY = "property"

@dataclass 
class CalculationResult:
    """  
    Writes the Result from A Quantum Chemistry Calculation
    """
    molecule: Molecule
    calculation_type: CalculationType
    sucess: bool 
    energy: Optional[float] = None 
    frequencies: Optional[List[float]] = None
    hessian: Optional[np.ndarray] = None 
    error: Optional[str] = None
    timestamp: str = field(default_factory = lambda: time.strftime("%Y%m%d_%H%M%S")) # Auto-generate the time-stamp
    method: Optional[str] = None 
    basis_set: Optional[str] = None 
    wall_time: float = 0.0 

    def to_dict(self) -> Dict[str, Any]:
        """   
        Convert to a dictionary for easy serialization
        """
        return {
            "molecule_name": self.molecule.name,
            "calculation_type": self.calculation_type.value,
            "sucess": self.sucess,
            "energy": self.energy,
            "frequencies": self.frequencies,
            "error": self.error,
            "timestamp": self.timestamp,
            "method": self.method,
            "basis_set": self.basis_set,
            "wall_time": self.wall_time
        }

@dataclass 
class Psi4Config:
    """
    Configuration settings for Psi4 calculations
    """
    method: str = "hf"
    basis_set: str = "sto-3g"
    memory: str = "500 MB"
    num_threads: int = 1
    max_iterations: int = 100 
    convergence: float = 1e-6

    def validate(self) -> None:
        """  
        Validate the configuration settings
        """
        if self.num_threads < 1:
            raise ValueError("Number of threads must be at least 1")
        if not self.memory.endswith(("MB", "GB")):
            raise ValueError("Memory must be specified in MB or GB")

# ============================ Helper Classes ===========================

class GeometryBuilder:
    """   
    Builds the geometry string for Psi4 from a Molecule object
    """

    @staticmethod
    def build_psi4_geometry(molecule: Molecule) -> str:
        """  
        Constructs the geometry string for Psi4 input
        """
        geom = f"{molecule.charge} {molecule.spin_mult}\n"
        for atom, (x,y,z) in zip(molecule.atom_labels, molecule.coordinates):
            geom += f" {GeometryBuilder.remove_digits(atom)} {x:.8f} {y:.8f} {z:.8f}\n"
        return geom
    
    @staticmethod
    def remove_digits(text: str) -> str:
        """Remove all digits from a string"""
        return ''.join(char for char in text if not char.isdigit())
    
class ScratchManager:
    """   
    Classs that manages the temporary scratch directories for Psi4 calculations
    """

    @staticmethod
    @contextmanager
    def temporary_scratch(prefix: str ="psi4_scratch_"):
        """   
        Context manager for temporary scratch directory 
        """
        scratch_dir = None
        try:
            scratch_dir = tempfile.mkdtemp(prefix=prefix)
            os.environ["PSI4_SCRATCH"] = scratch_dir
            os.environ["PSI_SCRATCH"] = scratch_dir
            yield scratch_dir
        finally:
            if scratch_dir and os.path.exists(scratch_dir):
                shutil.rmtree(scratch_dir, ignore_errors=True)

class ResultsManager:
    """   
    Manages the Calculation Results storage and Retrival
    """
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results: List[CalculationResult] = []

    def add_result(self, result: CalculationResult) -> None:
        """   
        Adds a CalculationResult to the results list
        """
        self.results.append(result)

    def get_sucessful_results(self) -> List[CalculationResult]:
        """   
        Returns a list of sucessful CalculationResults
        """
        return [res for res in self.results if res.sucess]
    
    def get_failed_results(self) -> List[CalculationResult]:
        """   
        Returns a list of failed CalculationResults
        """
        return [res for res in self.results if not res.sucess]
    
    def get_energies(self) -> List[Tuple[Molecule, float]]:
        """
        Get (molecule, energy) tuples for scuessful calculations
        """
        return [(res.molecule, res.energy) for res in self.results if res.sucess and res.energy is not None]

    def summary(self) -> Dict[str, Any]:
        """  
        Get Summary Statistics
        """
        total = len(self.results)
        sucessful = len(self.get_sucessful_results())
        failed = len(self.get_failed_results())

        energies = [res.energy for res in self.get_sucessful_results() if res.energy is not None]

        return {
            "total": total, 
            "sucessful": sucessful,
            "failed": failed,
            "sucess_rate": sucessful / total if total > 0 else 0.0,
            "min_energy": min(energies) if energies else None,
            "max_energy": max(energies) if energies else None,
            "mean_energy": np.mean(energies) if energies else None,
            "std_energy": np.std(energies) if energies else None
        }
    
# ============================ Worker Functions ===========================

class Psi4Worker:
    """  
    Static Psi4 Worker functions for parallel execution 
    """
    @staticmethod 
    def single_point_worker(args: Tuple) -> CalculationResult:
        """   
        Worker function for single point calculations

        Args:
            args: Tuple containing (Molecule, config_dict)

        Returns:
            CalculationResult object
        """
        molecule, config_dict = args
        config = Psi4Config(**config_dict)

        scratch_dir = tempfile.mkdtemp(prefix="psi4_scratch_")


        with ScratchManager.temporary_scratch():
            
            # Set the scratch directory for psi4
            os.environ["PSI4_SCRATCH"] = scratch_dir
            os.environ["PSI_SCRATCH"] = scratch_dir
            

            import psi4 
            psi4.core.set_output_file(os.devnull) # Supress output file 

            start_time = time.time()

            try:
                # Build geometry
                geom_str = GeometryBuilder.build_psi4_geometry(molecule)
                mol = psi4.geometry(geom_str)

                # Calculate energy
                method_str = f"{config.method}/{config.basis_set}"
                energy= psi4.energy(method_str, molecule=mol)

                wall_time = time.time() - start_time

                return CalculationResult(
                    molecule=molecule,
                    calculation_type = CalculationType.SINGLE_POINT,
                    sucess = True,
                    energy = energy,
                    method = config.method,
                    basis_set = config.basis_set,
                    wall_time = wall_time
                    )
            except Exception as e:
                wall_time = time.time() - start_time 
                return CalculationResult(
                    molecule=molecule,
                    calculation_type = CalculationType.SINGLE_POINT,
                    sucess = False,
                    error = str(e),
                    method = config.method,
                    basis_set = config.basis_set,
                    wall_time = wall_time
                )
            finally:
                psi4.core.clean() # Clean up psi4 core after worker

    @staticmethod 
    def gradient_worker(args: Tuple):
        """   
        Worker function for gradient calculations

        Args:
            args: Tuple containing (Molecule, config_dict)
        """
        molecule, config_dict = args
        config = Psi4Config(**config_dict)
        with ScratchManager.temporary_scratch():
            import psi4 
            psi4.core.set_output_file(os.devnull) 

            start_time = time.time()
            try:
                # Build geometry
                geom_str = GeometryBuilder.build_psi4_geometry(molecule)
                mol = psi4.geometry(geom_str)

                # Calculate gradient
                method_str = f"{config.method}/{config.basis_set}"
                gradient = psi4.gradient(method_str, molecule=mol)

                wall_time = time.time() - start_time

                #TODO clean this up
                return (molecule, np.array(gradient), True, None)
            except Exception as e:
                wall_time = time.time() - start_time 
                # TODO clean this up
                return (molecule, None, False, str(e))
            
            finally:
                psi4.core.clean() # Clean up psi4 core after worker


    @staticmethod
    def hessian_worker(args: Tuple):
        """ 
        Worker functio nfor the Hessian calculation

        Args:
            args: Tuple containing (Molecule, config_dict)
        """
        molecule, config_dict = args 
        config = Psi4Config(**config_dict)

        with ScratchManager.temporary_scratch():
            import psi4 
            psi4.core.set_output_file(os.devnull)

            start_time = time.time()

            try:
                geom_str = GeometryBuilder.build_psi4_geometry(molecule)
                mol = psi4.geometry(geom_str)

                # Calculate Hessian
                method_str = f"{config.method}/{config.basis_set}"
                hessian = psi4.hessian(method_str, molecule=mol)

                # convert to numpy array
                hessian_array = np.array(hessian)

                wall_time = time.time() - start_time
                return (molecule, hessian_array, True, None)
            except Exception as e:
                wall_time = time.time() - start_time
                return (molecule, None, False, str(e))
            finally:
                psi4.core.clean() # Clean up psi4 core after worker 

            


    @staticmethod
    def geometry_optimization_worker(args: Tuple) -> CalculationResult:
        """   
        Worker function for geometry optimization calculations

        Args:
            args: Tuple containing (Molecule, config_dict)
        Returns:
            CalculationResult object
        """
        molecule, config_dict = args
        config = Psi4Config(**config_dict)
        # create scratch for worker
        scratch_dir = tempfile.mkdtemp(prefix=f"psi4_scratch_{os.getpid()}_")

        try:
            # set enviroment in this scratch
            os.environ["PSI4_SCRATCH"] = scratch_dir
            os.environ["PSI_SCRATCH"] = scratch_dir

            import psi4 
            psi4.core.set_output_file(os.devnull) # Supress output file 

            start_time = time.time() 
            
            try:
                # Build geometry 
                geom_str = GeometryBuilder.build_psi4_geometry(molecule)
                mol = psi4.geometry(geom_str)

                # Optimize
                method_str = f"{config.method}/{config.basis_set}"
                energy = psi4.optimize(method_str, molecule=mol)

                opt_geom = np.array(mol.geometry())
                bohr_to_angstrom = 0.529177
                molecule.coordinates = opt_geom * bohr_to_angstrom

                wall_time = time.time() - start_time

                return CalculationResult(
                    molecule =molecule,
                    calculation_type = CalculationType.GEOMETRY_OPTIMIZATION,
                    sucess = True,
                    energy = energy,
                    method = config.method,
                    basis_set = config.basis_set,
                    wall_time = wall_time
                )
            except Exception as e:
                wall_time = time.time() - start_time
                return CalculationResult(
                    molecule =molecule,
                    calculation_type = CalculationType.GEOMETRY_OPTIMIZATION,
                    sucess = False,
                    error = str(e),
                    method = config.method,
                    basis_set = config.basis_set,
                    wall_time = wall_time
                )

        finally:
            if os.path.exists(scratch_dir):
                shutil.rmtree(scratch_dir, ignore_errors=True) # Clean up scratch directory after worker
            try:
                import psi4 
                psi4.core.clean() # Clean up psi4 core after worker
            except:
                pass


    def optimization_and_frequency_worker(args: Tuple) -> CalculationResult:
        """  
        Worker Function for combined optimization and frequency calculation
        Includes a verification if a true minimum is found by checking for imaginary frequencies
        
        1. Perform geometry optimization using psi4.optimize()
        2. Perform frequency calculation using psi4.frequency()
        3. Check if the structure is a true minimum (no imaginary frequencies)
        """ 
        molecule, config_dict = args
        config = Psi4Config(**config_dict)

        with ScratchManager.temporary_scratch():
            import psi4 
            psi4.core.set_output_file(os.devnull) # Supress output file 

            start_time = time.time() 

            try:
                # Build the geometry
                geom_str = GeometryBuilder.build_psi4_geometry(molecule)
                mol = psi4.geometry(geom_str)

                # Store original coordinates for RMSD calculation
                original_coords = molecule.coordinates.copy()  

                # Optimize 
                method_str = f"{config.method}/{config.basis_set}"
                energy = psi4.optimize(method_str, molecule=mol)

                # Update molecule coordinates
                opt_geom = np.array(mol.geometry())
                bohr_to_angstrom = 0.529177
                molecule.coordinates = opt_geom * bohr_to_angstrom


                # Frequency calculation
                scf_e, scf_wfn = psi4.frequency(method_str, molecule=mol, return_wfn=True)

                # Extract frequencies
                frequencies = scf_wfn.frequencies().to_array()

                # Check for imaginary frequencies 
                imaginary_freqs = [freq for freq in frequencies if freq < 0.0]
                is_minimum = len(imaginary_freqs) == 0
                num_imaginary = len(imaginary_freqs)
                lowest_freq = min(imaginary_freqs) if imaginary_freqs else None

                wall_time = time.time() - start_time

                result = CalculationResult(
                    molecule = molecule,
                    calculation_type = CalculationType.OPT_FREQ,
                    sucess = is_minimum,
                    energy = energy,
                    frequencies = scf_wfn.frequencies().to_array(),
                    hessian = scf_wfn.hessian().to_array(),
                    method = config.method,
                    basis_set = config.basis_set,
                    wall_time = wall_time
                )

                # TODO add minimum verification details to result
                return result
            except Exception as e:
                wall_time = time.time() - start_time
                return CalculationResult(
                    molecule = molecule,
                    calculation_type = CalculationType.OPT_FREQ,
                    sucess = False,
                    error = str(e),
                    method = config.method,
                    basis_set = config.basis_set,
                    wall_time = wall_time
                )
            finally:
                psi4.core.clean() # Clean up psi4 core after worker



        
# ============================ Main Calculator Class ===========================
class Psi4Calculator:
    """  
    Unified interface for Psi4 quantum chemistry calculations
    """
    def __init__(self,
                 config: Optional[Psi4Config] = None,
                 output_dir: str = "psi4_calculations",
                 save_results: bool = True,
                 verbose: bool = True):
        """   
        Initializes the Psi4Calculator object

        Args:
            config: Psi4Config object with calculation settings
            output_dir: Directory to save calculation outputs
            save_results: Whether to save results to files
            verbose: Whether to print progress messages
        """
        self.config = config or Psi4Config()
        # Validate config
        self.config.validate()

        # Output directory
        self.output_dir = Path(output_dir)
        self.save_results = save_results
        self.verbose = verbose

        # Results manager
        self.results_manager = ResultsManager(self.output_dir)

        # Configure Psi4
        self._configure_psi4()

    def _configure_psi4(self) -> None: 
        """  
        Configure global Psi4 settings
        """
        import psi4
        psi4.set_memory(self.config.memory)
        psi4.set_num_threads(self.config.num_threads)
        log_file = self.output_dir / "psi4_output.log"
        psi4.set_output_file(str(log_file), False)

    # ============================ Calculation Methods ===========================

    def single_point_energy(self, molecule: Molecule) -> Optional[float]:
        """  
        Calculate a single point enegry 

        Args:
            molecule: Molecule object for which to calculate the energy
        
        Returns:
            Energy in Hartree or None if failed
        """
        result = self._run_calculation(
            molecule = molecule,
            worker_func = Psi4Worker.single_point_worker
        )
        if self.save_results:
            self.results_manager.add_result(result)

        if self.verbose:
            self._print_result(result)

        return result.energy if result.sucess else None 
    
    def compute_gradient(self,molecule: Molecule) -> Optional[np.ndarray]:
        """   
        Computes the energy gradient for a molecule
        """
        config_dict = { 
            "method": self.config.method,
            "basis_set": self.config.basis_set,
            "memory": self.config.memory,
            "num_threads": self.config.num_threads,
        }
        molecule_copy, gradient, sucess, error = Psi4Worker.gradient_worker((molecule, config_dict))
        if self.verbose:
            if sucess:
                print(f"Gradient calculation for {molecule.name} succeeded.")
            else:
                print(f"Gradient calculation for {molecule.name} failed: {error}")
        if sucess:
            return gradient
        else:
            return None

    def compute_force(self,molecule) -> Optional[np.ndarray]:
        """   
        Computes the force of a molecule
        Force = - Gradient
        """
        gradient = self.compute_gradient(molecule)
        if gradient is not None:
            return -gradient
        else:
            return None

    def compute_hessian(self, molecule: Molecule) -> Optional[np.ndarray]:
        """  
        Computes the Hessian matrix of second derivatives for a molecule


        Returns the total non-mass-weighted electronic Hessian Matrix in Hartrees/Bohr^2
        """
        config_dict = {
            "method": self.config.method,
            "basis_set": self.config.basis_set,
            "memory": self.config.memory,
            "num_threads": self.config.num_threads,
        }
        molecule_copy, hessian, sucess, error = Psi4Worker.hessian_worker((molecule, config_dict))
        if self.verbose:
            if sucess:
                print(f"Hessian calculation for {molecule.name} succeeded.")
            else:
                print(f"Hessian calculation for {molecule.name} failed: {error}")
        if sucess:
            return hessian
        else:
            return None
    
    def geometry_optimization(self, molecule: Molecule) -> Optional[float]:
        """   
        Optimize molecular geometry

        Args:
            molecule: Molecule object to optimize
        Returns:
            Optimized energy in Hartree or None if failed
        """
        result = self._run_calculation(
            molecule = molecule,
            worker_func = Psi4Worker.geometry_optimization_worker 
        )

        if self.save_results:
            self.results_manager.add_result(result)
        if self.verbose:
            self._print_result(result)

        return result.energy if result.sucess else None
    
    def optimization_and_frequency_calculation(self,
                                    molecule: Molecule,
                                    freq_threshold: float = -5):
        """   
        Performs geometry optimization followed by a frequency calculation with verification of a true minimum 

        Args:
            molecule: Molecule object to optimize and analyze
            freq_threshold: Frequency threshold to consider imaginary frequencies
        Returns:
            Dictionary with Results 
        """
        config_dict = {
            "method": self.config.method,
            "basis_set": self.config.basis_set,
            "memory": self.config.memory,
            "num_threads": self.config.num_threads,
        }
        result = Psi4Worker.optimization_and_frequency_worker((molecule, config_dict))
        if self.save_results:
            self.results_manager.add_result(result)
        if self.verbose:
            self._print_result(result)
        return {
            "sucess": result.sucess,
            "energy": result.energy,
            "frequencies": result.frequencies,
            "hessian": result.hessian}
    

    
    # =========================== Batch Calculation Methods ==========================

    def batch_single_point_energy(self,
                                  molecules: List[Molecule],
                                  parallel: bool = False,
                                  n_processes: Optional[int] = None
                                  ) -> List[Tuple[Molecule, float]]:
        """
        Batch single point energy calculations

        Args:
            molecules: List of Molecule objects
            parallel: Whether to run calculations in parallel
            n_processes: Number of processes for parallel execution
        
        Returns:
            List of (Molecule, energy) tuples for successful calculations
        """
        return self._batch_calculation(
            molecules = molecules,
            worker_func=Psi4Worker.single_point_worker,
            parallel = parallel,
            n_processes = n_processes,
            calc_type="Single Point"
        )

    def batch_geometry_optimization(self,
                                    molecules: List[Molecule],
                                    parallel: bool = False,
                                    n_processes: Optional[int] = None
                                    ) -> List[Tuple[Molecule, float]]:
        """
        Batch geometry optimizations
        Args:
            molecules: List of Molecule objects
            parallel: Whether to run calculations in parallel
            n_processes: Number of processes for parallel execution
        Returns:
            List of (Molecule, energy) tuples for successful calculations
        """
        return self._batch_calculation(
            molecules = molecules,
            worker_func=Psi4Worker.geometry_optimization_worker,
            parallel = parallel,
            n_processes = n_processes,
            calc_type="Geometry Optimization"
        )
    
    # ============================ Internal Helper Methods ===========================
    def _run_calculation(self,
                         molecule: Molecule,
                         worker_func: Callable) -> CalculationResult:
        """  
        Run a single calculation using the specified worker function
        """
        config_dict = {
            "method": self.config.method,
            "basis_set": self.config.basis_set,
            "memory": self.config.memory,
            "num_threads": self.config.num_threads,
        }
        return worker_func((molecule, config_dict))
    
    def _batch_calculation(self,
                           molecules: List[Molecule],
                           worker_func: Callable,
                            parallel: bool,
                            n_processes: Optional[int],
                            calc_type: str) -> List[Tuple[Molecule, float]]:
        """  
        Generic batch calculation method

        Args:
            molecules: List of Molecule objects
            worker_func: Worker function to use for calculations
            parallel: Whether to run in parallel
            n_processes: Number of processes for parallel execution
            calc_type: String description of calculation type for logging

        Returns:
            List of (Molecule, energy) tuples for successful calculations    
        """
        if not molecules:
            raise ValueError("No molecules provided for batch calculation")
        
        start_time = time.time()

        if self.verbose:
            print(f"Starting batch {calc_type} for {len(molecules)} molecules...")
            print(f"Method: {self.config.method}, Basis Set: {self.config.basis_set}")
            print(f"Parallel: {parallel}, Processes: {n_processes if n_processes else 'Auto'}")

        config_dict = {
            "method": self.config.method,
            "basis_set": self.config.basis_set,
            "memory": self.config.memory,
            "num_threads": self.config.num_threads,
        }
        args = [(mol, config_dict) for mol in molecules]

        if parallel:
            results = self._execute_parallel(args, worker_func, n_processes)
        else:
            result = self._execute_serial(args, worker_func)

        # Store results
        if self.save_results:
            for res in results:
                self.results_manager.add_result(res)
        
        # Extract successful energies
        sucessful = [(r.molecule, r.energy) for r in results if r.sucess]

        elapsed_time = time.time() - start_time
        if self.verbose:
            print(f"Batch {calc_type} completed in {elapsed_time:.2f} seconds")
            summary = self.results_manager.summary()
            print(f"Total: {summary['total']}, Sucessful: {summary['sucessful']}, Failed: {summary['failed']}, Sucess Rate: {summary['sucess_rate']*100:.2f}%")

        return sucessful
    
    def _execute_serial(self,
                        args: List[Tuple],
                        worker_func: Callable) -> List[CalculationResult]:
        """  
        Execute calculations serially
        """

        results = []

        with tqdm(total=len(args),
                  desc="Processing molecules",
                  unit="mol",
                  bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]") as pbar:
            
            for arg in args:
                result = worker_func(arg)
                results.append(result)
                pbar.update(1)

                if self.verbose:
                    if result.sucess:
                        pbar.write(f"  ✓ {result.molecule.name}: {result.energy:.6f} Ha")
                    else:
                        pbar.write(f"  ✗ {result.molecule.name}: {result.error[:30]}")
        return results
    
    def _execute_parallel(self,
                          args: List[Tuple],
                          worker_func: Callable,
                            n_processes: Optional[int] = None,
                            timeout: int = 300) -> List[CalculationResult]:  # 5 min timeout
        """
        Execute calculations in parallel using multiprocessing
        """
        if n_processes is None:
            n_processes = max(1, cpu_count() - 2) # Leave some CPUs free

        results = []

        with Pool(processes=n_processes) as pool:
            with tqdm(total = len(args),
                      desc="Processing molecules",
                      unit="mol",
                      bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]") as pbar:

                # Use imap with timeout instead of imap_unordered
                async_results = [pool.apply_async(worker_func, (arg,)) for arg in args]
            
                for async_result in async_results:
                    try:
                        result = async_result.get(timeout=timeout)
                        results.append(result)
                        pbar.update(1)

                        if self.verbose:
                            if result.sucess:
                                pbar.write(f"  ✓ {result.molecule.name}: {result.energy:.6f} Ha")
                            else:
                                pbar.write(f"  ✗ {result.molecule.name}: {result.error[:30]}")
                    except TimeoutError:
                        pbar.write(f"  ✗ TIMEOUT after {timeout}s")
                        results.append(CalculationResult(
                            molecule=args[len(results)][0],
                            calculation_type=CalculationType.GEOMETRY_OPTIMIZATION,
                            sucess=False,
                            error=f"Timeout after {timeout} seconds"
                        ))
                        pbar.update(1)
                    except Exception as e:
                        pbar.write(f"  ✗ ERROR: {str(e)[:50]}")
                        pbar.update(1)
            
        return results

    def _print_result(self, result: CalculationResult) -> None:
        """   
        Print calculation result
        """
        calc_type = result.calculation_type.value.replace("_", " ").title()

        if result.sucess:
            print(f"{calc_type} for {result.molecule.name} succeeded: Energy = {result.energy} Hartree")
        else:
            print(f"{calc_type} for {result.molecule.name} failed: Error = {result.error}")

    # ============================= Results Managment Methods ============================

    def get_results_summary(self) -> Dict[str, Any]:
        """  
        Get summary statistics of all calculations
        """
        return self.results_manager.summary()
    
    def get_sucessful_results(self) -> List[CalculationResult]:
        """  
        Get list of sucessful CalculationResults
        """
        return self.results_manager.get_sucessful_results()
    
    def get_failed_results(self) -> List[CalculationResult]:
        """  
        Get list of failed CalculationResults
        """
        return self.results_manager.get_failed_results()
    
    # ============================== Test Methods =============================================

    def determine_method_and_basis_set_combinations(self,
                                                    molecule: Molecule,
                                                    methods: List[str] = ["hf", "b3lyp", "pbe0", "mp2", "ccsd", "ccsd(t)"],
                                                    basis_sets: List[str] = ["sto-3g", "3-21g", "6-31g", "6-311g(d,p)", "cc-pvdz", "cc-pvtz", "def2-svp", "def2-tzvp"]
                                                    ) -> List[Tuple[str, str, float]]:
        """"   
        Tests different method and basis set combinations for a given molecule
        """
        results = {}
        for method in methods:
            for basis in basis_sets:
                start_time = time.time()
                self.config.method = method 
                self.config.basis_set = basis
                energy = self.single_point_energy(molecule)
                end_time = time.time()
                elapsed = end_time - start_time

                if energy is not None:
                    results[(method, basis)] = (energy, elapsed)
                    if self.verbose:
                        print(f"Method: {method}, Basis Set: {basis}, Energy: {energy} Hartree, Time: {elapsed:.2f} seconds")
                else:
                    if self.verbose:
                        print(f"Method: {method}, Basis Set: {basis} failed.")
        
        # Make for each method a subplot showing basis set vs energy
        # We have 6 methods -> 2 rows, 3 columns
        fig, axs = plt.subplots(2, 3, figsize=(18, 10))
        axs = axs.flatten()
        for i, method in enumerate(methods):
            method_results = [(basis, results[(method, basis)][0]) for basis in basis_sets if (method, basis) in results]
            if not method_results:
                continue
            basis_labels = [b[0] for b in method_results]
            energies = [b[1] for b in method_results]
            axs[i].plot(basis_labels, energies, marker='o')
            axs[i].set_title(f'Method: {method.upper()}')
            axs[i].set_xlabel('Basis Set')
            axs[i].set_ylabel('Energy (Hartree)')
            axs[i].grid(True)
            axs[i].tick_params(axis='x', rotation=45)
        plt.tight_layout()
        plt.savefig("figures/method_basis_set_combinations.png")


        return results 

    def test_basis_set_convergence(self,
            molecule: Molecule,
            method="hf"):
            """
            Test the Basis set convergence for a given molecule and a given method:

            For this we use the following basis sets:

            + STO-3G
            + 3-21G
            + 6-311G(D,P)
            + cc-PVDZ
            + cc-PVTZ
            + cc-PVQZ                      
            """

            basis_sets = [
                "sto-3g",
                "3-21g",
                "6-311g(d,p)",
                "cc-pvdz",
                "cc-pvtz",
                "cc-pvqz"
            ]
            results = []
            times = []
            for basis in basis_sets:
                start_time = time.time()
                self.config.method = method
                self.config.basis_set = basis
                energy = self.single_point_energy(molecule)
                results.append((basis, energy))
                if self.verbose:
                    print(f"Basis Set: {basis}, Energy: {energy} Hartree")
                elapsed = time.time() - start_time
                times.append(elapsed)

            # Two subplots the call energy vs basis set and time vs basis set
            fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8,10))
            basis_labels = [b[0] for b in results]
            energies = [b[1] for b in results]
            ax1.plot(basis_labels, energies, marker='o')
            ax1.set_title(f'Basis Set Convergence using {method.upper()}')
            ax1.set_xlabel('Basis Set')
            ax1.set_ylabel('Energy (Hartree)')
            ax1.grid(True)
            times_minutes = [t/60.0 for t in times]
            ax2.plot(basis_labels, times_minutes, marker='o', color='orange')
            ax2.set_title(f'Calculation Time vs Basis Set using {method.upper()}')
            ax2.set_xlabel('Basis Set')
            ax2.set_ylabel('Time (Minutes)')
            ax2.grid(True)
            plt.tight_layout()
            plt.savefig("figures/basis_set_convergence.png")
            plt.close() 



    def test_threading_performance(self, 
            molecule: Molecule,
            method="hf",
            basis_set="cc-pvtz",
            max_threads: int = 10):
            """
            Small helper function to test speedup from threading
            """ 
            thread_counts = list(range(1, max_threads + 1))
            times = []
            for n_threads in thread_counts:
                self.config.method = method
                self.config.basis_set = basis_set 
                self.config.num_threads = n_threads
                start_time = time.time()
                energy = self.single_point_energy(molecule)
                elapsed = time.time() - start_time
                times.append(elapsed)
                if self.verbose:
                    print(f"Threads: {n_threads}, Energy: {energy} Hartree, Time: {elapsed:.2f} seconds")
            # Plotting the results
            plt.figure(figsize=(8,6))
            plt.plot(thread_counts, times, marker='o')
            plt.title(f'Threading Performance using {method.upper()} / {basis_set.upper()}')
            plt.xlabel('Number of Threads')
            plt.ylabel('Time (Seconds)')
            plt.grid(True)
            plt.savefig("figures/threading_performance.png")
            plt.close()