import os
from datetime import datetime

class Logger:
    """ 
    A simple Logger Class for writing a output file with timestamps
    """
    def __init__(self, out_file = "cluster_gen.out", mode="a"):
        """ 
        Initializes the Logger with the output file name and mode

        Args:
            out_file (str): Name of the output file
            mode (str): File mode, e.g., 'a' for append, 'w' for write
        """
        self.out_file = out_file
        self.mode = mode

    @staticmethod
    def remove_digits(text):
        """Remove all digits from a string"""
        return ''.join(char for char in text if not char.isdigit())
    
    def write_header(self):
        """ 
        Writes the general header to the output file
        """
        with open(self.out_file, self.mode) as f:
            f.write("="*50 + "\n")
            f.write(f"Cluster Generation Log - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("="*50 + "\n\n")

    def write_molecule_info(self,molecule):
        """ 
        Writes molecule information to the log file, i.e atom lables, coordinates, bonds, etc.
        """ 
        with open(self.out_file, self.mode) as f:
            f.write("-"*50 + "\n")
            f.write(f"Molecule Information - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("-"*50 + "\n")
            f.write(f"Number of Atoms: {len(molecule.atom_labels)}\n")
            f.write("Atom Labels and Coordinates:\n")
            for label, coord in zip(molecule.atom_labels, molecule.coordinates):
                f.write(f"{label}: {coord[0]:.6f}, {coord[1]:.6f}, {coord[2]:.6f}\n")
            covalent_bonds, hydrogen_bonds = molecule.get_bonds()
            f.write(f"Covalent Bonds ({len(covalent_bonds)}):\n")
            for bond in covalent_bonds:
                f.write(f"  {bond[0]} - {bond[1]} (Degree: {bond[2]:.3f})\n")
            f.write(f"Hydrogen Bonds ({len(hydrogen_bonds)}):\n")
            for bond in hydrogen_bonds:
                f.write(f"  {bond[0]} - {bond[1]} (Degree: {bond[2]:.3f})\n")
            f.write("\n")

    def write_scf_batch_result(self, results):
        """ 
        Writes the SCF batch calculation results to the log file
        """
        with open(self.out_file, self.mode) as f:
            f.write("-"*50 + "\n")
            f.write(f"SCF Batch Calculation Results - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("-"*50 + "\n")
            for idx, (sampled_mol ,energy) in enumerate(results): 
                f.write(f"Sampled Molecule {idx+1}: \n")
                f.write(f"SCF Energy: {energy:.6f} Hartree\n")
                f.write("\n")
    
    @staticmethod
    def write_trajectory_sampling(sampled_mols):
        """ 
        Writes a xyz trajectory file of the sampled molecules
        """
        with open("sampled_configurations_trj.xyz", "w") as trj_file:
            for i, mol in enumerate(sampled_mols):
                trj_file.write(f"{len(mol.atom_labels)}\n")
                trj_file.write(f"Sampled Configuration {i+1}\n")
                for label, coord in zip(mol.atom_labels, mol.coordinates):
                    trj_file.write(f"{Logger.remove_digits(label)} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")


    @staticmethod
    def write_optimization_xyz_batch(opt_results):
        """" 
        Makes a new directory 'optimized_geometries' and writes the optimized geometries in xyz files
        """
        output_dir = "optimized_geometries"
        os.makedirs(output_dir, exist_ok=True)
        
        for idx, (opt_mol, opt_energy) in enumerate(opt_results):
            file_path = os.path.join(output_dir, f"optimized_molecule_{idx+1}.xyz")
            with open(file_path, "w") as xyz_file:
                xyz_file.write(f"{len(opt_mol.atom_labels)}\n")
                xyz_file.write(f"Optimized Molecule {idx+1} Energy: {opt_energy:.6f} Hartree\n")
                for label, coord in zip(opt_mol.atom_labels, opt_mol.coordinates):
                    xyz_file.write(f"{Logger.remove_digits(label)} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")






