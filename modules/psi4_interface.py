"""
PSI4 Interface Module
"""

import os
import sys
import time
import tempfile
import shutil
import numpy as np
from typing import List, Optional, Tuple, Dict, Any
from dataclasses import dataclass, field
from enum import Enum
from multiprocessing import Process, Queue
from queue import Empty

from molecule_class import Molecule


# ==============================================================================
# DATA CLASSES - Simple containers for configuration and results
# ==============================================================================

class Status(Enum):
    """Calculation status codes"""
    SUCCESS = "success"
    FAILED = "failed"
    TIMEOUT = "timeout"
    NOT_CONVERGED = "not_converged"


@dataclass
class Config:
    """
    Configuration for Psi4 calculations.
    
    Attributes:
        method: QM method (e.g., 'hf', 'b3lyp', 'mp2')
        basis: Basis set (e.g., 'sto-3g', '6-31g*', 'cc-pvdz')
        charge: Molecular charge
        multiplicity: Spin multiplicity (1=singlet, 2=doublet, etc.)
        memory: Memory allocation (e.g., '2GB')
        max_iter: Maximum iterations for geometry optimization
        reference: Reference wavefunction ('rhf', 'uhf', 'rohf')
    """
    method: str = "hf"
    basis: str = "sto-3g"
    charge: int = 0
    multiplicity: int = 1
    memory: str = "2GB"
    max_iter: int = 50
    reference: str = "rhf"
    scf_max_iter: int = 50
    sp_timeout: int = 60 # Timeout for single point calculations
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for passing to worker processes"""
        return {
            'method': self.method,
            'basis': self.basis,
            'charge': self.charge,
            'multiplicity': self.multiplicity,
            'memory': self.memory,
            'max_iter': self.max_iter,
            'reference': self.reference, 
            'scf_max_iter': self.scf_max_iter,
            'sp_timeout': self.sp_timeout
        }


@dataclass
class Result:
    """
    Container for calculation results.
    
    Attributes:
        molecule: The molecule (with updated coordinates if optimization)
        energy: Final energy in Hartree (None if failed)
        status: Calculation status
        error_message: Error description if failed
        wall_time: Total calculation time in seconds
        iterations: Number of optimization iterations (if applicable)
    """
    molecule: Molecule
    energy: Optional[float] = None
    status: Status = Status.FAILED
    error_message: str = ""
    wall_time: float = 0.0
    iterations: int = 0
    
    @property
    def success(self) -> bool:
        """Check if calculation was successful"""
        return self.status == Status.SUCCESS
    
    def __repr__(self) -> str:
        if self.success:
            return f"Result(success=True, E={self.energy:.6f} Ha, time={self.wall_time:.1f}s)"
        return f"Result(success=False, status={self.status.value}, error='{self.error_message[:50]}')"


# ==============================================================================
# WORKER FUNCTIONS - Run in isolated subprocesses
# ==============================================================================
# These functions are intentionally simple and self-contained.
# They import psi4 inside the function to ensure clean process isolation.
# All state is passed via arguments and returned via the result.

def _remove_digits(text):
    """Utility function to remove digits from a string"""
    return ''.join([i for i in text if not i.isdigit()])

def _build_geometry_string(molecule: Molecule, charge: int, multiplicity: int) -> str:
    """
    Build Psi4 geometry string from Molecule object.
    
    Format:
        charge multiplicity
        symbol x y z
        ...
        units angstrom
        no_reorient
        no_com
    """
    lines = [f"{charge} {multiplicity}"]
    
    for symbol, coord in zip(molecule.atom_labels, molecule.coordinates):
        symbol_clean = _remove_digits(symbol)
        lines.append(f"{symbol_clean} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}")
    
    lines.extend(["units angstrom", "no_reorient", "no_com"])
    
    return "\n".join(lines)


def _single_point_worker(molecule: Molecule, config_dict: Dict) -> Dict:
    """
    Worker function for single point energy calculation.
    
    This runs in a separate process. It:
    1. Sets up Psi4 environment
    2. Builds the molecule
    3. Calculates the energy
    4. Returns results as a dictionary
    
    Args:
        molecule: Molecule object
        config_dict: Configuration dictionary
        
    Returns:
        Dictionary with keys: success, energy, error, wall_time
    """
    start_time = time.time()
    scratch_dir = tempfile.mkdtemp(prefix="psi4_sp_")
    
    result = {
        'success': False,
        'energy': None,
        'dipole_vec': None,
        'dipole_magnitude': None,
        'mulliken_charges': None,
        'error': '',
        'wall_time': 0.0,
        'coordinates': None
    }
    
    try:
        # Setup environment
        os.environ["PSI4_SCRATCH"] = scratch_dir
        os.environ["OMP_NUM_THREADS"] = "1"
        
        # Import psi4 here to ensure clean import in subprocess
        import psi4
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.set_output_file(os.devnull, False)
        psi4.set_memory(config_dict['memory'])
        psi4.set_num_threads(1)
        
        psi4.set_options({
            'reference': config_dict['reference'],
            'scf_type': 'df',
            'maxiter': config_dict.get('scf_max_iter', 50)
        })
        
        # Build geometry
        geom_str = _build_geometry_string(
            molecule, 
            config_dict['charge'], 
            config_dict['multiplicity']
        )
        psi4_mol = psi4.geometry(geom_str)
        
        # Calculate energy
        method_basis = f"{config_dict['method']}/{config_dict['basis']}"
        energy, wfn = psi4.energy(method_basis, molecule=psi4_mol, return_wfn=True)
        
        # Compute dipole moment
        psi4.oeprop(wfn, 'DIPOLE')
        vec = wfn.variable('CURRENT DIPOLE')
        vec = np.array(vec)
        magnitude = np.linalg.norm(vec)

        # Calculate Mulliken Atomic Charges 
        psi4.oeprop(wfn, 'MULLIKEN_CHARGES')
        mulliken = wfn.variable('MULLIKEN CHARGES')
        
        result['dipole_vec'] = vec.tolist()
        result['dipole_magnitude'] = magnitude
        result['mulliken_charges'] = mulliken

        result['success'] = True
        result['energy'] = float(energy)
        result['wall_time'] = time.time() - start_time

         
    except Exception as e:
        result['error'] = str(e)
        result['wall_time'] = time.time() - start_time
        
    finally:
        # Cleanup
        try:
            import psi4
            psi4.core.clean()
        except Exception:
            pass
        
        if os.path.exists(scratch_dir):
            shutil.rmtree(scratch_dir, ignore_errors=True)
    
    return result


def _geometry_opt_worker(molecule: Molecule, config_dict: Dict) -> Dict:
    """
    Worker function for geometry optimization.
    
    This runs in a separate process. It:
    1. Sets up Psi4 environment
    2. Builds the molecule
    3. Runs geometry optimization
    4. Extracts optimized coordinates
    5. Returns results as a dictionary
    
    Args:
        molecule: Molecule object
        config_dict: Configuration dictionary
        
    Returns:
        Dictionary with keys: success, energy, error, wall_time, coordinates, iterations
    """
    start_time = time.time()
    scratch_dir = tempfile.mkdtemp(prefix="psi4_opt_")
    
    result = {
        'success': False,
        'energy': None,
        'error': '',
        'wall_time': 0.0,
        'coordinates': None,
        'iterations': 0,
        'converged': False
    }
    
    try:
        # Setup environment
        os.environ["PSI4_SCRATCH"] = scratch_dir
        os.environ["OMP_NUM_THREADS"] = "1"
        
        # Import psi4 here to ensure clean import in subprocess
        import psi4
        psi4.core.set_output_file(os.devnull, False)
        psi4.set_memory(config_dict['memory'])
        psi4.set_num_threads(1)
        
        # Set options for optimization
        psi4.set_options({
            'reference': config_dict['reference'],
            'scf_type': 'df',
            'geom_maxiter': config_dict['max_iter'],
            'g_convergence': 'gau_loose',  # Looser convergence for speed
            'opt_coordinates': 'cartesian'
        })
        
        # Build geometry
        geom_str = _build_geometry_string(
            molecule, 
            config_dict['charge'], 
            config_dict['multiplicity']
        )
        psi4_mol = psi4.geometry(geom_str)
        
        # Run optimization
        method_basis = f"{config_dict['method']}/{config_dict['basis']}"
        
        try:
            energy = psi4.optimize(method_basis, molecule=psi4_mol)
            result['converged'] = True
        except psi4.OptimizationConvergenceError:
            # Get last energy even if not converged
            try:
                energy = psi4.core.variable('CURRENT ENERGY')
            except Exception:
                energy = None
            result['converged'] = False
            result['error'] = f"Did not converge in {config_dict['max_iter']} iterations"
        
        # Extract optimized coordinates
        opt_geom = np.array(psi4_mol.geometry())
        bohr_to_angstrom = 0.529177210903
        coordinates = opt_geom * bohr_to_angstrom
        
        # Get iteration count
        try:
            iterations = int(psi4.core.variable('OPTIMIZATION ITERATIONS'))
        except Exception:
            iterations = 0
        
        result['success'] = result['converged']
        result['energy'] = float(energy) if energy is not None else None
        result['coordinates'] = coordinates.tolist()
        result['iterations'] = iterations
        result['wall_time'] = time.time() - start_time
        
    except Exception as e:
        result['error'] = str(e)
        result['wall_time'] = time.time() - start_time
        
    finally:
        # Cleanup
        try:
            import psi4
            psi4.core.clean()
            psi4.core.clean_options()
        except Exception:
            pass
        
        if os.path.exists(scratch_dir):
            shutil.rmtree(scratch_dir, ignore_errors=True)
    
    return result


# ==============================================================================
# DIRECT FUNCTIONS - For use inside multiprocessing workers
# ==============================================================================

def direct_energy(molecule: Molecule, 
                  method: str = "hf", 
                  basis: str = "sto-3g",
                  charge: int = 0,
                  multiplicity: int = 1,
                  timeout: int = 60,
                  scf_max_iter: int = 100) -> Optional[float]:
    """
    Calculate energy directly WITHOUT spawning a subprocess.
    
    Use this when you are already inside a multiprocessing worker.
    
    Args:
        molecule: Molecule to calculate
        method: QM method
        basis: Basis set
        charge: Molecular charge
        multiplicity: Spin multiplicity
        timeout: Maximum time in seconds
        scf_max_iter: Maximum SCF iterations

    Returns:
        Tuple of (energy in Hartree, dipole vector [dx,dy,dz], dipole magnitude) and or None if failed
    """
    config_dict = {
        'method': method,
        'basis': basis,
        'charge': charge,
        'multiplicity': multiplicity,
        'memory': '2GB',
        'reference': 'rhf',
        'scf_max_iter': scf_max_iter,
    }
    
    try:
        result = _single_point_worker(molecule, config_dict)
        
        if result is None:
            print(f"direct_energy: _single_point_worker returned None")
            return None
            
        if result.get('success'):
            return result.get('energy'), result.get('dipole_vec'), result.get('dipole_magnitude'), result.get('mulliken_charges')
        else:
            # Zeige Fehler an für Debugging
            error = result.get('error', 'Unknown error')
            print(f"direct_energy failed: {error[:100]}")
            return None
            
    except Exception as e:
        # No printout here to avoid cluttering
        return None


def direct_optimize(molecule: Molecule,
                   method: str = "hf",
                   basis: str = "sto-3g", 
                   charge: int = 0,
                   multiplicity: int = 1,
                   max_iter: int = 50) -> Tuple[Optional[Molecule], Optional[float]]:
    """
    Optimize geometry directly WITHOUT spawning a subprocess.
    
    Use this when you are already inside a multiprocessing worker.
    
    Args:
        molecule: Molecule to optimize
        method: QM method
        basis: Basis set
        charge: Molecular charge
        multiplicity: Spin multiplicity
        max_iter: Maximum optimization iterations
        
    Returns:
        Tuple of (optimized_molecule, energy) or (None, None) if failed
    """
    config_dict = {
        'method': method,
        'basis': basis,
        'charge': charge,
        'multiplicity': multiplicity,
        'memory': '2GB',
        'reference': 'rhf',
        'max_iter': max_iter
    }
    
    result = _geometry_opt_worker(molecule, config_dict)
    
    if result.get('success') and result.get('coordinates') is not None:
        opt_mol = molecule.copy()
        opt_mol.coordinates = np.array(result['coordinates'])
        return opt_mol, result.get('energy')
    
    return None, None


# ==============================================================================
# PROCESS RUNNER - Handles process isolation and timeout
# ==============================================================================

def _worker_wrapper(func, args, queue):
    """
    Wrapper function for subprocess execution.
    """
    try:
        result = func(*args)
        queue.put(('success', result))
    except Exception as e:
        queue.put(('error', str(e)))


class ProcessRunner:
    """
    Runs a function in an isolated subprocess with timeout.
    """
    
    @staticmethod
    def run(worker_func, args: Tuple, timeout: int) -> Tuple[bool, Any]:
        """
        Run a worker function in a subprocess with timeout.
        
        Args:
            worker_func: The function to run (must be picklable)
            args: Arguments to pass to the function
            timeout: Maximum time in seconds
            
        Returns:
            Tuple of (completed, result)
            - completed: True if finished within timeout
            - result: Function result or None if timeout/error
        """
        queue = Queue()
        
        # Use module-level function instead of local lambda
        process = Process(
            target=_worker_wrapper,
            args=(worker_func, args, queue)
        )
        process.start()
        process.join(timeout=timeout)
        
        if process.is_alive():
            # Timeout - kill the process
            process.terminate()
            process.join(timeout=5)
            
            if process.is_alive():
                # Force kill if still alive
                process.kill()
                process.join()
            
            return False, None
        
        # Get result from queue
        try:
            status, result = queue.get_nowait()
            return status == 'success', result
        except Empty:
            return False, None


# ==============================================================================
# PARALLEL BATCH WORKER - For multiprocessing pool
# ==============================================================================

def _batch_optimize_worker(args: Tuple) -> Dict:
    """
    Worker for parallel batch optimization.
    Runs in a Pool process - does NOT spawn subprocesses.
    
    Args:
        args: Tuple of (molecule, config_dict, timeout)
        
    Returns:
        Dictionary with result data
    """
    molecule, config_dict = args
    
    # Run optimization directly (no subprocess)
    result = _geometry_opt_worker(molecule, config_dict)
    
    # Add molecule info to result for identification
    result['molecule_name'] = molecule.name
    result['original_molecule'] = molecule
    
    return result


def _batch_single_point_worker(args: Tuple) -> Dict:
    """
    Worker for parallel batch single point calculations.
    Runs in a Pool process - does NOT spawn subprocesses.
    """
    molecule, config_dict = args
    
    result = _single_point_worker(molecule, config_dict)
    result['molecule_name'] = molecule.name
    result['original_molecule'] = molecule
    
    return result


# ==============================================================================
# MAIN CALCULATOR CLASS - User-facing API
# ==============================================================================

class Psi4Calculator:
    """
    Main interface for Psi4 calculations.
    
    Supports both sequential and parallel batch processing.
    
    Example:
        config = Config(method='hf', basis='sto-3g')
        calc = Psi4Calculator(config)
        
        # Sequential (with timeout per molecule)
        results = calc.batch_optimize(molecules, timeout=300)
        
        # Parallel (faster, but no per-molecule timeout)
        results = calc.batch_optimize_parallel(molecules, n_processes=4)
    """
    
    def __init__(self, config: Optional[Config] = None, verbose: bool = True):
        """
        Initialize calculator.
        
        Args:
            config: Calculation configuration (default settings if None)
            verbose: Print progress and status messages
        """
        self.config = config or Config()
        self.verbose = verbose
        
        # Statistics tracking
        self._stats = {
            'total': 0,
            'success': 0,
            'failed': 0,
            'timeout': 0
        }
    
    # ==========================================================================
    # Single Calculations (with subprocess + timeout)
    # ==========================================================================
    
    def single_point(self, molecule: Molecule, timeout: int = 60) -> Result:
        """
        Calculate single point energy with timeout protection.
        Runs in isolated subprocess.
        """
        if self.verbose:
            print(f"Calculating energy for {molecule.name}...", end=" ", flush=True)
        
        completed, worker_result = ProcessRunner.run(
            _single_point_worker,
            args=(molecule, self.config.to_dict()),
            timeout=timeout
        )
        
        self._stats['total'] += 1
        
        if not completed:
            self._stats['timeout'] += 1
            if self.verbose:
                print(f"TIMEOUT ({timeout}s)")
            return Result(
                molecule=molecule,
                status=Status.TIMEOUT,
                error_message=f"Timeout after {timeout}s",
                wall_time=timeout
            )
        
        if worker_result is None or not worker_result.get('success', False):
            self._stats['failed'] += 1
            error_msg = worker_result.get('error', 'Unknown error') if worker_result else 'No result'
            if self.verbose:
                print(f"FAILED: {error_msg[:50]}")
            return Result(
                molecule=molecule,
                status=Status.FAILED,
                error_message=error_msg,
                wall_time=worker_result.get('wall_time', 0) if worker_result else 0
            )
        
        self._stats['success'] += 1
        if self.verbose:
            print(f"E = {worker_result['energy']:.6f} Ha ({worker_result['wall_time']:.1f}s)")
        
        return Result(
            molecule=molecule,
            energy=worker_result['energy'],
            status=Status.SUCCESS,
            wall_time=worker_result['wall_time']
        )
    
    def optimize(self, molecule: Molecule, timeout: int = 300) -> Result:
        """
        Optimize geometry with timeout protection.
        Runs in isolated subprocess.
        """
        if self.verbose:
            print(f"Optimizing {molecule.name}...", end=" ", flush=True)
        
        completed, worker_result = ProcessRunner.run(
            _geometry_opt_worker,
            args=(molecule, self.config.to_dict()),
            timeout=timeout
        )
        
        self._stats['total'] += 1
        
        if not completed:
            self._stats['timeout'] += 1
            if self.verbose:
                print(f"TIMEOUT ({timeout}s)")
            return Result(
                molecule=molecule,
                status=Status.TIMEOUT,
                error_message=f"Timeout after {timeout}s",
                wall_time=timeout
            )
        
        if worker_result is None:
            self._stats['failed'] += 1
            if self.verbose:
                print("FAILED: No result")
            return Result(
                molecule=molecule,
                status=Status.FAILED,
                error_message="No result returned"
            )
        
        if not worker_result.get('success', False):
            self._stats['failed'] += 1
            status = Status.NOT_CONVERGED if 'converge' in worker_result.get('error', '').lower() else Status.FAILED
            if self.verbose:
                print(f"FAILED: {worker_result.get('error', 'Unknown')[:50]}")
            return Result(
                molecule=molecule,
                energy=worker_result.get('energy'),
                status=status,
                error_message=worker_result.get('error', 'Unknown'),
                wall_time=worker_result.get('wall_time', 0),
                iterations=worker_result.get('iterations', 0)
            )
        
        self._stats['success'] += 1
        
        opt_molecule = molecule.copy()
        if worker_result.get('coordinates') is not None:
            opt_molecule.coordinates = np.array(worker_result['coordinates'])
        
        if self.verbose:
            print(f"E = {worker_result['energy']:.6f} Ha "
                  f"({worker_result['iterations']} iter, {worker_result['wall_time']:.1f}s)")
        
        return Result(
            molecule=opt_molecule,
            energy=worker_result['energy'],
            status=Status.SUCCESS,
            wall_time=worker_result['wall_time'],
            iterations=worker_result['iterations']
        )
    
    # ==========================================================================
    # Sequential Batch (with timeout per molecule)
    # ==========================================================================
    
    def batch_single_point(self, 
                          molecules: List[Molecule], 
                          timeout: int = 60) -> List[Result]:
        """
        Calculate energies sequentially with timeout per molecule.
        """
        if self.verbose:
            self._print_batch_header("Single Point", len(molecules), timeout)
        
        results = []
        for i, mol in enumerate(molecules):
            if self.verbose:
                print(f"[{i+1}/{len(molecules)}] ", end="")
            result = self.single_point(mol, timeout=timeout)
            results.append(result)
        
        if self.verbose:
            self._print_batch_summary(results)
        return results
    
    def batch_optimize(self, 
                      molecules: List[Molecule], 
                      timeout: int = 300) -> List[Result]:
        """
        Optimize geometries sequentially with timeout per molecule.
        """
        if self.verbose:
            self._print_batch_header("Geometry Optimization", len(molecules), timeout)
        
        results = []
        for i, mol in enumerate(molecules):
            if self.verbose:
                print(f"[{i+1}/{len(molecules)}] ", end="")
            result = self.optimize(mol, timeout=timeout)
            results.append(result)
        
        if self.verbose:
            self._print_batch_summary(results)
        return results
    
    # ==========================================================================
    # Parallel Batch (faster, no per-molecule timeout)
    # ==========================================================================
    
    def batch_single_point_parallel(self,
                                   molecules: List[Molecule],
                                   n_processes: Optional[int] = None,
                                   maxtasksperchild: int = 1) -> List[Result]:
        """
        Calculate energies in parallel using multiprocessing Pool.
        
        Note: No per-molecule timeout. Use maxtasksperchild=1 to prevent
        memory leaks and stuck workers.
        
        Args:
            molecules: List of molecules
            n_processes: Number of parallel processes (default: CPU count - 2)
            maxtasksperchild: Restart worker after N tasks (default: 1)
            
        Returns:
            List of Result objects
        """
        import multiprocessing as mp
        
        if n_processes is None:
            n_processes = max(1, mp.cpu_count() - 2)
        
        if self.verbose:
            print(f"\n{'='*60}")
            print(f"Parallel Single Point Calculations")
            print(f"Molecules: {len(molecules)}, Processes: {n_processes}")
            print(f"Method: {self.config.method}/{self.config.basis}")
            print(f"{'='*60}\n")
        
        # Prepare arguments
        config_dict = self.config.to_dict()
        args_list = [(mol, config_dict) for mol in molecules]
        
        # Run in parallel
        results = []
        with mp.Pool(processes=n_processes, maxtasksperchild=maxtasksperchild) as pool:
            for i, worker_result in enumerate(pool.imap(_batch_single_point_worker, args_list)):
                result = self._process_worker_result(worker_result, is_optimization=False)
                results.append(result)
                
                if self.verbose:
                    self._print_progress(i + 1, len(molecules), result)
        
        if self.verbose:
            self._print_batch_summary(results)
        
        return results
    
    def batch_optimize_parallel(self,
                               molecules: List[Molecule],
                               n_processes: Optional[int] = None,
                               maxtasksperchild: int = 1) -> List[Result]:
        """
        Optimize geometries in parallel using multiprocessing Pool.
        
        Note: No per-molecule timeout. Use maxtasksperchild=1 to prevent
        memory leaks and stuck workers.
        
        Args:
            molecules: List of molecules
            n_processes: Number of parallel processes (default: CPU count - 2)
            maxtasksperchild: Restart worker after N tasks (default: 1)
            
        Returns:
            List of Result objects
        """
        import multiprocessing as mp
        
        if n_processes is None:
            n_processes = max(1, mp.cpu_count() - 2)
        
        if self.verbose:
            print(f"\n{'='*60}")
            print(f"Parallel Geometry Optimization")
            print(f"Molecules: {len(molecules)}, Processes: {n_processes}")
            print(f"Method: {self.config.method}/{self.config.basis}")
            print(f"Max iterations: {self.config.max_iter}")
            print(f"{'='*60}\n")
        
        # Prepare arguments
        config_dict = self.config.to_dict()
        args_list = [(mol, config_dict) for mol in molecules]
        
        # Run in parallel with progress
        results = []
        with mp.Pool(processes=n_processes, maxtasksperchild=maxtasksperchild) as pool:
            for i, worker_result in enumerate(pool.imap(_batch_optimize_worker, args_list)):
                result = self._process_worker_result(worker_result, is_optimization=True)
                results.append(result)
                
                if self.verbose:
                    self._print_progress(i + 1, len(molecules), result)
        
        if self.verbose:
            self._print_batch_summary(results)
        
        return results
    
    def batch_optimize_parallel_unordered(self,
                                          molecules: List[Molecule],
                                          n_processes: Optional[int] = None,
                                          maxtasksperchild: int = 1) -> List[Result]:
        """
        Optimize geometries in parallel, returning results as they complete.
        
        Faster than batch_optimize_parallel because it doesn't wait for
        molecules to complete in order.
        
        Args:
            molecules: List of molecules
            n_processes: Number of parallel processes
            maxtasksperchild: Restart worker after N tasks
            
        Returns:
            List of Result objects (may be in different order than input)
        """
        import multiprocessing as mp
        
        if n_processes is None:
            n_processes = max(1, mp.cpu_count() - 2)
        
        if self.verbose:
            print(f"\n{'='*60}")
            print(f"Parallel Geometry Optimization (Unordered)")
            print(f"Molecules: {len(molecules)}, Processes: {n_processes}")
            print(f"Method: {self.config.method}/{self.config.basis}")
            print(f"{'='*60}\n")
        
        config_dict = self.config.to_dict()
        args_list = [(mol, config_dict) for mol in molecules]
        
        results = []
        completed = 0
        
        with mp.Pool(processes=n_processes, maxtasksperchild=maxtasksperchild) as pool:
            # imap_unordered returns results as soon as they're ready
            for worker_result in pool.imap_unordered(_batch_optimize_worker, args_list):
                result = self._process_worker_result(worker_result, is_optimization=True)
                results.append(result)
                completed += 1
                
                if self.verbose:
                    self._print_progress(completed, len(molecules), result)
        
        if self.verbose:
            self._print_batch_summary(results)
        
        return results
    
    # ==========================================================================
    # Helper Methods
    # ==========================================================================
    
    def _process_worker_result(self, worker_result: Dict, is_optimization: bool) -> Result:
        """Convert worker result dict to Result object"""
        original_mol = worker_result.get('original_molecule')
        
        self._stats['total'] += 1
        
        if not worker_result.get('success', False):
            self._stats['failed'] += 1
            error = worker_result.get('error', 'Unknown error')
            status = Status.NOT_CONVERGED if 'converge' in error.lower() else Status.FAILED
            
            return Result(
                molecule=original_mol,
                energy=worker_result.get('energy'),
                status=status,
                error_message=error,
                wall_time=worker_result.get('wall_time', 0),
                iterations=worker_result.get('iterations', 0)
            )
        
        self._stats['success'] += 1
        
        # For optimization, update coordinates
        if is_optimization and worker_result.get('coordinates') is not None:
            opt_mol = original_mol.copy()
            opt_mol.coordinates = np.array(worker_result['coordinates'])
        else:
            opt_mol = original_mol
        
        return Result(
            molecule=opt_mol,
            energy=worker_result['energy'],
            status=Status.SUCCESS,
            wall_time=worker_result.get('wall_time', 0),
            iterations=worker_result.get('iterations', 0)
        )
    
    def _print_batch_header(self, calc_type: str, n_molecules: int, timeout: int):
        """Print batch calculation header"""
        print(f"\n{'='*60}")
        print(f"Batch {calc_type}")
        print(f"Molecules: {n_molecules}, Timeout: {timeout}s each")
        print(f"Method: {self.config.method}/{self.config.basis}")
        print(f"{'='*60}\n")
    
    def _print_progress(self, current: int, total: int, result: Result):
        """Print progress for a single molecule"""
        if result.success:
            print(f"[{current}/{total}] ✓ {result.molecule.name}: "
                  f"E = {result.energy:.6f} Ha ({result.wall_time:.1f}s)")
        else:
            print(f"[{current}/{total}] ✗ {result.molecule.name}: "
                  f"{result.status.value} - {result.error_message[:40]}")
    
    def _print_batch_summary(self, results: List[Result]):
        """Print summary of batch results"""
        success = sum(1 for r in results if r.success)
        failed = sum(1 for r in results if r.status == Status.FAILED)
        timeout = sum(1 for r in results if r.status == Status.TIMEOUT)
        not_converged = sum(1 for r in results if r.status == Status.NOT_CONVERGED)
        
        total_time = sum(r.wall_time for r in results)
        
        print(f"\n{'='*60}")
        print(f"Batch Complete")
        print(f"  Success:       {success}/{len(results)}")
        print(f"  Failed:        {failed}")
        print(f"  Timeout:       {timeout}")
        print(f"  Not converged: {not_converged}")
        print(f"  Total time:    {total_time:.1f}s")
        
        # Energy statistics for successful calculations
        successful_energies = [r.energy for r in results if r.success and r.energy is not None]
        if successful_energies:
            print(f"\n  Energy range:  [{min(successful_energies):.6f}, {max(successful_energies):.6f}] Ha")
            print(f"  Mean energy:   {np.mean(successful_energies):.6f} Ha")
        
        print(f"{'='*60}\n")
    
    def get_stats(self) -> Dict[str, int]:
        """Get calculation statistics"""
        return self._stats.copy()
    
    def reset_stats(self):
        """Reset statistics counters"""
        for key in self._stats:
            self._stats[key] = 0


# ==============================================================================
# CONVENIENCE FUNCTIONS - For quick calculations
# ==============================================================================

def quick_energy(molecule: Molecule, method: str = "hf", basis: str = "sto-3g") -> Optional[float]:
    """
    Quick single point energy calculation.
    
    Args:
        molecule: Molecule to calculate
        method: QM method
        basis: Basis set
        
    Returns:
        Energy in Hartree or None if failed
    """
    config = Config(method=method, basis=basis)
    calc = Psi4Calculator(config, verbose=False)
    result = calc.single_point(molecule, timeout=120)
    return result.energy if result.success else None


def quick_optimize(molecule: Molecule, method: str = "hf", basis: str = "sto-3g") -> Optional[Molecule]:
    """
    Quick geometry optimization.
    
    Args:
        molecule: Molecule to optimize
        method: QM method
        basis: Basis set
        
    Returns:
        Optimized molecule or None if failed
    """
    config = Config(method=method, basis=basis)
    calc = Psi4Calculator(config, verbose=False)
    result = calc.optimize(molecule, timeout=300)
    return result.molecule if result.success else None


