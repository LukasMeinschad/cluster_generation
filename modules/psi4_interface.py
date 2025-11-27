import psi4 
import numpy as np
from multiprocessing import Pool, cpu_count
import tempfile
import os
import time
from molecule_class import Molecule


class Psi4Calculator:
    """
    A class to perform QM calculations using Psi4
    """
    def __init__(self, molecule = None, memory = "500 MB", num_threads=1,ls_of_molecules = None, method="ccsd(t)", basis_set="cc-pvdz"):
        """ 
        Initializes the Psi4Calculator with a Molecule object and calulation setting

        Args:
            molecule: Molecule object containing the molecular structure
            memory: Memory allocation for Psi4 calculations
            num_threads: Number of threads to use
            ls_of_molecules: List of Molecule objects for batch calculations
        """
        self.molecule = molecule
        self.ls_of_molecules = ls_of_molecules
        self.memory = memory
        self.num_threads = num_threads
        self.method = method
        self.basis_set = basis_set
        
        # Configure Psi4 settings 
        psi4.set_memory(self.memory)
        psi4.set_num_threads(self.num_threads)
        psi4.set_output_file("psi4_output.dat",False)


    def build_geometry_string(self,molecule):
        """ 
        Uses the Molecule object to build the geometry string for psi4
        """
        geom = f"{molecule.charge} {molecule.spin_mult}\n"
        for atom, (x,y,z) in zip(molecule.atom_labels, molecule.coordinates):
            geom += f" {self.remove_digits(atom)} {x:.8f} {y:.8f} {z:.8f}\n"
        return geom
    
    #TODO maybe remove this in the end
    @staticmethod 
    def remove_digits(text):
        """Remove all digits from a string"""
        return ''.join(char for char in text if not char.isdigit())



    def single_point_calc(self,molecule=None):
        """ 
        Run a single point calculation

        Args:
            molecule: Molecule object for which to run the calculation. If None, uses self.molecule
        """
        target_molecule = molecule if molecule is not None else self.molecule
        if target_molecule is None:
            raise ValueError("No molecule provided for single point calculation")
        
        geom = self.build_geometry_string(target_molecule)
        mol = psi4.geometry(geom)
        method = f"{self.method}/{self.basis_set}"
        energy, wfn  = psi4.energy(method, molecule=mol, return_wfn=True)
        
        print(f"Single Point Energy for {target_molecule.name}: {energy} Hartree") 
        return energy
        

    def geometry_optimization(self, molecule=None):
        """
        Run a geometry optimization
        """
        target_molecule = molecule if molecule is not None else self.molecule
        if target_molecule is None:
            raise ValueError("No molecule provided for geometry optimization")
        
        geom = self.build_geometry_string(target_molecule)
        mol = psi4.geometry(geom)
        method = f"{self.method}/{self.basis_set}"
        optimized_energy = psi4.optimize(method, molecule=mol)
        # Update molecule coordinates after optimization
        optimized_coords = mol.save_string_xyz()
        new_molecule = Molecule.from_xyz(optimized_coords, name=target_molecule.name + "_optimized")
        target_molecule.coordinates = new_molecule.coordinates
        print(f"Optimized Energy for {target_molecule.name}: {optimized_energy} Hartree")
        return optimized_energy

             
    
 
    def batch_single_point_calc(self, parallel=False, n_processes=None):
        """ 
        Run single point calculations for a list of Molecule objects

        If parallel is True, uses multiprocessing to speed up calculations

        """
        if self.ls_of_molecules is None:
            raise ValueError("No list of molecules provided for batch calculation")
        
        time_start = time.time()

        if parallel:
            results = self._batch_single_point_calc_parallel(n_processes=n_processes)
        else:
            results = self._single_point_calc_serial()

        time_end = time.time()
        print(f"Total time for batch calculations: {time_end - time_start:.2f} seconds")
        return results
    
    def _single_point_calc_serial(self):
        """ 
        Run a batch calculation serially 
        """
        energies = []
        for mol_obj in self.ls_of_molecules:
            geom = self.build_geometry_string(mol_obj)
            mol = psi4.geometry(geom)
            method = f"{self.method}/{self.basis_set}"
            energy = psi4.energy(method, molecule=mol)
            print(f"Molecule: {mol_obj.name}, Energy: {energy} Hartree")
            energies.append((mol_obj, energy))
        return energies

    def _batch_single_point_calc_parallel(self, n_processes=None):
        """ 
        Run a batch calculation in parallel using multiprocessing
        """
        if n_processes is None:
            n_processes = cpu_count() - 3 # Leave some CPUs free

        args = [(m, self.method, self.basis_set) for m in self.ls_of_molecules]

        with Pool(processes=n_processes) as pool:
            results = pool.map(self._psi4_worker, args)

        for mol, E in results:
            print(f"Molecule: {mol.name}, Energy: {E} Hartree")
        return results

    def batch_geometry_optimization(self,parallel=False,n_processes=None):
        """ 
        Run geometry optimizations for a list of Molecule objects
        
        If parallel is True, uses multiprocessing to speed up calculations
        """
        if self.ls_of_molecules is None:
            raise ValueError("No list of molecules provided for batch geometry optimization")
        
        time_start = time.time()

        if parallel:
            results = self._batch_geometry_optimization_parallel(n_processes=n_processes)
        else:
            results = self._geometry_optimization_serial()

        time_end = time.time()
        print(f"Total time for batch geometry optimizations: {time_end - time_start:.2f} seconds")
        return results
    
    def _geometry_optimization_serial(self):
        """
        Perform a batch geometry optimization serially
        """

        results = []
        for mol_obj in self.ls_of_molecules:
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
        return results
    
    def _batch_geometry_optimization_parallel(self, n_processes=None):
        """
        Performs a batch geometry optimization in parallel using multiprocessing
        """
        if n_processes is None:
            n_processes = cpu_count() - 3 # Leave some CPUs free

        args = [(m, self.method, self.basis_set) for m in self.ls_of_molecules]

        with Pool(processes=n_processes) as pool:
            results = pool.map(self._psi4_worker_geom_opt, args)

    
    @staticmethod
    def _psi4_worker_geom_opt(args):
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
            # Update molecule coordinates after optimization
            optimized_coords = mol.save_string_xyz()
            new_molecule = Molecule.from_xyz(optimized_coords, name=mol_obj.name + "_optimized")
            mol_obj.coordinates = new_molecule.coordinates
            return (mol_obj, optimized_energy)
        finally:
            # Clean up scratch
            import shutil
            shutil.rmtree(scratch, ignore_errors=True)
            psi4.core.clean() # Clean up psi4 core for this worker


    @staticmethod
    def _psi4_worker(args):
        """ 
        Worker function for parallel execution of single point batch calculations
        
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
            return (mol_obj, energy)
        finally:
            # Clean up scratch
            import shutil
            shutil.rmtree(scratch, ignore_errors=True)
            psi4.core.clean() # Clean up psi4 core for this worker
        