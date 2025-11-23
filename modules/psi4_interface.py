import psi4 
import numpy as np
from multiprocessing import Pool, cpu_count
import tempfile
import os
import time


#TODO maybe remove this in the end
def remove_digits(text):
    """Remove all digits from a string"""
    return ''.join(char for char in text if not char.isdigit())

class Psi4Calculator:
    """
    A class to perform QM calculations using Psi4
    """
    def __init__(self, molecule = None, memory = "500 MB", num_threads=1,ls_of_molecules = None):
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
        self.method = "ccsd(t)"
        self.basis_set = "cc-pvdz"
        self.geom = None # Here the Geometry Input for the Psi_4 calculation is written
        psi4.set_memory(self.memory)
        psi4.set_num_threads(self.num_threads)
        # Set the PSI4 Output Path
        psi4.set_output_file("psi4_output.dat",False)


    def build_geometry_string(self):
        """ 
        Uses the Molecule object to build the geometry string for psi4

        # Example how this geometry can look in Z-Matrix Coords
        h2o = psi4.geometry(
        O
        H 1 0.96
        H 1 0.96 2 104.5
        """
        geom = f"{self.molecule.charge} {self.molecule.spin_mult}\n"
        for atom, (x,y,z) in zip(self.molecule.atom_labels, self.molecule.coordinates):
            geom += f" {remove_digits(atom)} {x:.8f} {y:.8f} {z:.8f}\n"

        print(geom)
        self.geom = geom 
    
    def single_point_calc(self):
        """ 
        Run a single point calculation
        """
        geom = self.geom
        mol = psi4.geometry(geom)
        method = f"{self.method}/{self.basis_set}"
        energy = psi4.energy(method, molecule=mol)
        print(energy)

    def geometry_optimization(self):
        """
        Run a geometry optimization
        """
        geom = self.geom
        mol = psi4.geometry(geom) 
        method = f"{self.method}/{self.basis_set}"
        scf_e = psi4.optimize(method, molecule=mol) 
 
    def batch_single_point_calc(self):
        """ 
        Run single point calculations for a list of Molecule objects
        """
        # Start timer
        time_start = time.time()

        energies = []
        for mol_obj in self.ls_of_molecules:
            geom = f"{mol_obj.charge} {mol_obj.spin_mult}\n"
            for atom, (x,y,z) in zip(mol_obj.atom_labels, mol_obj.coordinates):
                geom += f" {remove_digits(atom)} {x:.8f} {y:.8f} {z:.8f}\n"
            mol = psi4.geometry(geom)
            method = f"{self.method}/{self.basis_set}"
            energy = psi4.energy(method, molecule=mol)
            energies.append(energy)
            print(f"Molecule: {mol_obj.name}, Energy: {energy}")

        time_end = time.time()
        print(f"Total time for batch calculations: {time_end - time_start:.2f} seconds")
        return energies
    
    @staticmethod
    def single_point_calc_worker(args):
        """ 
        Worker function for a parallel execution of single point calculations

        Its important to set the scratch file right here
        """
        mol_obj, method, basis_set = args
        scratch = tempfile.mkdtemp(prefix="psi4_scratch_")
        os.environ["PSI_SCRATCH"] = scratch 
        psi4.core.set_output_file(os.devnull)
        


        geom = f"{mol_obj.charge} {mol_obj.spin_mult}\n"
        for atom, (x,y,z) in zip(mol_obj.atom_labels, mol_obj.coordinates):
            geom += f" {remove_digits(atom)} {x:.8f} {y:.8f} {z:.8f}\n"
        mol = psi4.geometry(geom)
        
        energy = psi4.energy(f"{method}/{basis_set}", molecule=mol)

        return (mol_obj.name, energy) 
    
    @staticmethod
    def psi_worker(args):
        import psi4 
        import os, tempfile

        mol_obj, method, basis_set = args
        scratch = tempfile.mkdtemp(prefix="psi4_scratch_")
        os.environ["PSI_SCRATCH"] = scratch
        psi4.core.set_output_file(os.devnull)

        geom = f"{mol_obj.charge} {mol_obj.spin_mult}\n"
        for atom, (x,y,z) in zip(mol_obj.atom_labels, mol_obj.coordinates):
            geom += f" {remove_digits(atom)} {x:.8f} {y:.8f} {z:.8f}\n"
        mol = psi4.geometry(geom)
        E = psi4.energy(f"{method}/{basis_set}", molecule=mol)
        return (mol_obj.name, E)

    def batch_single_point_calc_parallel(self, n_processes=None):
        """ 
        Run single point calculations for a list of molecule objects in parallel
        
        Note that we need to create the right scratch directory for each process to avoid conflicts
        """
        time_start = time.time()

        if n_processes is None:
            n_processes = cpu_count()

        args = [(m, self.method, self.basis_set) for m in self.ls_of_molecules]

        with Pool(processes=n_processes) as pool:
            results = pool.map(Psi4Calculator.psi_worker, args)

        for name, E in results:
            print(f"Molecule: {name}, Energy: {E}")

        time_end = time.time()
        print(f"Total time for parallel batch calculations: {time_end - time_start:.2f} seconds")

    
