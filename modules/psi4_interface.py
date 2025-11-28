import psi4 
import numpy as np
from multiprocessing import Pool, cpu_count
from typing import List, Tuple, Union, Optional, Dict, Any
import tempfile
import json
import pickle
import os
import time
from pathlib import Path
from molecule_class import Molecule


class Psi4Calculator:
    """
    A class to interface with Psi4 for quantum chemistry calculations.
    """
    def __init__(self, 
                 molecule: Optional[Molecule] = None,
                memory: str = "500 MB",
                num_threads: int =1,
                ls_of_molecules: Optional[List[Molecule]] = None,
                method: str ="ccsd(t)",
                basis_set: str ="cc-pvdz",
                output_dir: str = "psi4_calculations",
                save_results: bool = True):
        """ 
        Initializes a Psi4Calculator object with molecular structures and calculation setting

        Args:
            molecule (Molecule): A Molecule object for single calculations
            memory (str): Memory allocation for Psi4
            num_threads (int): Number of threads for Psi4
            ls_of_molecules (List[Molecule]): List of Molecule objects for batch calculations
            method (str): Quantum chemistry method to use
            basis_set (str): Basis set to use
            output_dir (str): Directory to save calculation outputs
            save_results (bool): Whether to save results to files
        """

        self.molecule = molecule
        self.ls_of_molecules = ls_of_molecules
        self.memory = memory
        self.num_threads = num_threads
        self.method = method
        self.basis_set = basis_set
        self.save_results = save_results
        self.output_dir = Path(output_dir)

        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
    

        # Configure Psi4 settings
        self._configure_psi4()

        # Results storage
        self.results = {}
        self.failed_calculations = []


    def _configure_psi4(self):
        """ 
        Configures the global Psi4 settings
        """
        psi4.set_memory(self.memory)
        psi4.set_num_threads(self.num_threads)
        psi4.set_output_file(str(self.output_dir / "psi4_output.log"), False)

    def build_geometry_string(self,molecule: Molecule) -> str:
        """ 
        Uses the Molecule object to build the geometry string for psi4
        """
        geom = f"{molecule.charge} {molecule.spin_mult}\n"
        for atom, (x,y,z) in zip(molecule.atom_labels, molecule.coordinates):
            geom += f" {self.remove_digits(atom)} {x:.8f} {y:.8f} {z:.8f}\n"
        return geom
    

    @staticmethod 
    def remove_digits(text: str) -> str:
        """Remove all digits from a string"""
        return ''.join(char for char in text if not char.isdigit())

    def _save_result(self, molecule: Molecule, result_type: str, result: Any, success: bool = True):
        """ 
        Saves the calculation results
        """
        if not self.save_results:
            return
        
        mol_name = molecule.name.replace(" ", "_")
        timestamp = time.strftime("%Y%m%d_%H%M%S")

        result_data = {
            "molecule_name": mol_name,
            "method": self.method,
            "basis_set": self.basis_set,
            "result_type": result_type,
            "success": success,
            "timestamp": timestamp,
            "result": result
        }

        key = f"{mol_name}_{result_type}_{timestamp}"
        self.results[key] = result_data

        if not sucess:
            self.failed_calculations.append(result_data)

    def single_point_calc(self,molecule: Optional[Molecule] = None) -> float:
        """ 
        Run a single point calculation with error handling

        Args:
            molecule: Molecule object for which to run the calculation. If None, uses self.molecule
        """
        target_molecule = molecule if molecule is not None else self.molecule
        if target_molecule is None:
            raise ValueError("No molecule provided for single point calculation")
        

        try:
            geom = self.build_geometry_string(target_molecule)
            mol = psi4.geometry(geom)
            method = f"{self.method}/{self.basis_set}"
            energy, wfn  = psi4.energy(method, molecule=mol, return_wfn=True)
        
            print(f"Single Point Energy for {target_molecule.name}: {energy} Hartree") 
            self.save_results(target_molecule, "single_point", energy, success=True)
            return energy
        except Exception as e:
            error_msg = f"Single point calculation falied for {target_molecule.name}: {e}"
            print(error_msg)
            self.save_results(target_molecule, "single_point", error_msg, success=False)

            return None
        
        

    def geometry_optimization(self, molecule: Optional[Molecule] = None) -> float:
        """
        Run a geometry optimization with error handling
        """
        target_molecule = molecule if molecule is not None else self.molecule
        if target_molecule is None:
            raise ValueError("No molecule provided for geometry optimization")
        
        try:
            geom = self.build_geometry_string(target_molecule)
            mol = psi4.geometry(geom)
            method = f"{self.method}/{self.basis_set}"
            optimized_energy = psi4.optimize(method, molecule=mol)
        
            # Update molecule coordinates after optimization
            optimized_coords = mol.save_string_xyz()
            new_molecule = Molecule.from_xyz(optimized_coords, name=target_molecule.name + "_optimized")
            target_molecule.coordinates = new_molecule.coordinates
            print(f"Optimized Energy for {target_molecule.name}: {optimized_energy} Hartree")

            # Save sucessful result
            self.save_results(target_molecule, "geometry_optimization", optimized_energy, success=True)
            return optimized_energy
        except Exception as e:
            error_msg = f"Geometry optimization failed for {target_molecule.name}: {e}"
            print(error_msg)
            self.save_results(target_molecule, "geometry_optimization", error_msg, success=False)
            return None
             
    
 
    def batch_single_point_calc(self, parallel: bool = False, n_processes: Optional[int] = None) -> List[Tuple[Molecule, float]]:
        """ 
        Run single point calculations for a list of Molecule objects

        If parallel is True, uses multiprocessing to speed up calculations

        """
        if self.ls_of_molecules is None:
            raise ValueError("No list of molecules provided for batch calculation")
        
        time_start = time.time()

        try:
            if parallel:
                results = self._batch_single_point_calc_parallel(n_processes=n_processes)
            else:
                results = self._single_point_calc_serial()
        except Exception as e:
            print(f"Batch calculation error: {e}")
            results = []

        time_end = time.time()
        print(f"Total time for batch calculations: {time_end - time_start:.2f} seconds")
        
        return results
    
    def _single_point_calc_serial(self) -> List[Tuple[Molecule, float]]:
        """ 
        Run a batch calculation serially 
        """
        energies = []
        for mol_obj in self.ls_of_molecules:
            try:
                geom = self.build_geometry_string(mol_obj)
                mol = psi4.geometry(geom)
                method = f"{self.method}/{self.basis_set}"
                energy = psi4.energy(method, molecule=mol)
                print(f"Molecule: {mol_obj.name}, Energy: {energy} Hartree")
                energies.append((mol_obj, energy))

                # Save successful result
                self._save_result(mol_obj, "single_point", energy, success=True)
            
            except Exception as e:
                error_msg = f"Single point calculation failed for {mol_obj.name}: {e}"
                print(error_msg)
                self._save_result(mol_obj, "single_point", error_msg, success=False)
        
        return energies

    def _batch_single_point_calc_parallel(self, n_processes: Optional[int] = None) -> List[Tuple[Molecule, float]]:
        """ 
        Run a batch calculation in parallel using multiprocessing
        """
        if n_processes is None:
            n_processes = cpu_count() - 3 # Leave some CPUs free

        args = [(m, self.method, self.basis_set) for m in self.ls_of_molecules]

        sucessful_results = []
        with Pool(processes=n_processes) as pool:
            for result in pool.imap_unordered(self._psi4_worker_single_point, args):
                if result["sucess"]:
                    sucessful_results.append((result["molecule"], result["energy"]))
                    print(f"Molecule: {result['molecule'].name}, Energy: {result['energy']} Hartree")
                else:
                    print(f"Calculation failed for {result['molecule'].name}: {result['error']}")
                    self.failed_calculations.append((result['molecule'], result['error']))

        return sucessful_results


    def batch_geometry_optimization(self,parallel: bool=False,n_processes: Optional[int] = None) -> List[Tuple[Molecule, float]]:
        """ 
        Run geometry optimizations for a list of Molecule objects
        
        If parallel is True, uses multiprocessing to speed up calculations
        """
        if self.ls_of_molecules is None:
            raise ValueError("No list of molecules provided for batch geometry optimization")
        
        time_start = time.time()

        try:
            if parallel:
                results = self._batch_geometry_optimization_parallel(n_processes=n_processes)
            else:
                results = self._batch_geometry_optimization_serial()
        except Exception as e:
            print(f"Error during batch geometry optimization: {e}")
            results = []
        


        time_end = time.time()
        print(f"Total time for batch geometry optimizations: {time_end - time_start:.2f} seconds")
        print(f"Number of successful optimizations: {len(results)} out of {len(self.ls_of_molecules)}")
        return results
    
    def _batch_geometry_optimization_serial(self) -> List[Tuple[Molecule, float]]:
        """
        Perform a batch geometry optimization serially
        """
        results = []
        for mol_obj in self.ls_of_molecules:
            try:
                geom = self.build_geometry_string(mol_obj)
                mol = psi4.geometry(geom)
                method = f"{self.method}/{self.basis_set}"
                optimized_energy = psi4.optimize(method, molecule=mol)
            
                # Update molecule coordinates after optimization
                optimized_coords = mol.save_string_xyz()
                new_molecule = Molecule.from_xyz(optimized_coords, name=mol_obj.name + "_optimized")
                mol_obj.coordinates = new_molecule.coordinates
                print(f"Molecule: {mol_obj.name}, Optimized Energy: {optimized_energy} Hartree")
                results.append((mol_obj, optimized_energy))

                # Save successful result
                self._save_result(mol_obj, "geometry_optimization", optimized_energy, success=True)
            
            except Exception as e:
                error_msg = f"Geometry optimization failed for {mol_obj.name}: {e}"
                print(error_msg)
                self._save_result(mol_obj, "geometry_optimization", error_msg, success=False)
        
        return results
    
    def _batch_geometry_optimization_parallel(self, n_processes: Optional[int]=None) -> List[Tuple[Molecule, float]]:
        """
        Performs a batch geometry optimization in parallel using multiprocessing
        """
        if n_processes is None:
            n_processes = cpu_count() - 3 # Leave some CPUs free

        args = [(m, self.method, self.basis_set) for m in self.ls_of_molecules]

        sucessful_results = []
        with Pool(processes=n_processes) as pool:
            for result in pool.imap_unordered(self._psi4_worker_geom_opt, args):
                if result["sucess"]:
                    sucessful_results.append((result["molecule"], result["energy"]))
                    print(f"Molecule: {result['molecule'].name}, Optimized Energy: {result['energy']} Hartree")
                else:
                    print(f"Geometry optimization failed for {result['molecule'].name}: {result['error']}")
                    self.failed_calculations.append((result['molecule'], result['error']))
        return sucessful_results

    @staticmethod
    def _psi4_worker_geom_opt(args) -> Dict[str, Any]:
        """ 
        Worker function for parallel execution of geometry optimizations
        
        Note each worker needs to reimport psi4 and set a temporary scratch
        """
        import tempfile
        import os
        import shutil

        mol_obj, method, basis_set = args

        # Set up scratch 
        scratch = tempfile.mkdtemp(prefix="psi4_scratch_")
        os.environ["PSI4_SCRATCH"] = scratch # Set scratch for this process
        os.environ["PSI_SCRATCH"] = scratch # Set scratch for this process

        import psi4
        psi4.core.set_output_file(os.devnull) # Supress output file

        try:
            geom = f"{mol_obj.charge} {mol_obj.spin_mult}\n"
            for atom, (x,y,z) in zip(mol_obj.atom_labels, mol_obj.coordinates):
                geom += f" {Psi4Calculator.remove_digits(atom)} {x:.8f} {y:.8f} {z:.8f}\n"
            mol = psi4.geometry(geom)
            full_method = f"{method}/{basis_set}"
            optimized_energy = psi4.optimize(full_method, molecule=mol)


            # Update molecule coordinates
            optimized_coords = mol.save_string_xyz()

            new_molecule = Molecule.from_xyz(optimized_coords, name=mol_obj.name + "_optimized")
            mol_obj.coordinates = new_molecule.coordinates

            return {
                "molecule": mol_obj,
                "energy": optimized_energy,
                "sucess": True,
                "error": None
            }
        
        except Exception as e:
            return {
                "molecule": mol_obj,
                "energy": None,
                "sucess": False,
                "error": str(e)
            }

        finally:
            # Clean up scratch
            import shutil
            shutil.rmtree(scratch, ignore_errors=True)
            psi4.core.clean() # Clean up psi4 core for this worker


    @staticmethod
    def _psi4_worker_single_point(args) -> Dict[str, Any]:
        """ 
        Worker function for parallel single point calculatins
        """
        import tempfile
        import os
        import shutil


        mol_obj, method, basis_set = args

        # Set up scratch 
        scratch = tempfile.mkdtemp(prefix="psi4_scratch_")
        os.environ["PSI4_SCRATCH"] = scratch # Set scratch for this process
        os.environ["PSI_SCRATCH"] = scratch # Set scratch for this process

        # Initialization of psi4 now
        import psi4
        psi4.core.set_output_file(os.devnull) # Supress output file


        try:
            geom = f"{mol_obj.charge} {mol_obj.spin_mult}\n"
            for atom, (x,y,z) in zip(mol_obj.atom_labels, mol_obj.coordinates):
                geom += f" {Psi4Calculator.remove_digits(atom)} {x:.8f} {y:.8f} {z:.8f}\n"
            mol = psi4.geometry(geom)
            full_method = f"{method}/{basis_set}"
            energy = psi4.energy(full_method, molecule=mol)

            return {
                "molecule": mol_obj,
                "energy": energy,
                "sucess": True,
                "error": None
            }
        
        except Exception as e:
            return {
                "molecule": mol_obj,
                "energy": None,
                "sucess": False,
                "error": str(e)
            }
        
        finally:
            # Clean up scratch
            import shutil
            shutil.rmtree(scratch, ignore_errors=True)
            psi4.core.clean() # Clean up psi4 core for this worker
        