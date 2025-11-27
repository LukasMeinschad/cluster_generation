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
            f.write("Molecule Information:\n")
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







#TODO: I wrote this thing 4 times right now
def remove_digits(text):
    """Remove all digits from a string"""
    return ''.join(char for char in text if not char.isdigit())




def write_trajectory_sampling(sampled_mols):
    """
    Writes a Trajectory file of the sampled molecules.

    Args:
        sampled_mols (list): List of Molecule objects representing the sampled configurations.

    Basic Idea:

    - Open a _trj.xyz in write mode.
    - For each molecule in sampled mols:
        - Write the number of atoms
        - Write the comment line, e.g., "Sampled Configuration i"
        - for each atom in the molecule:
            - Write the atom label and coordinates in XYZ format.
    """
    with open("sampled_configurations_trj.xyz", "w") as trj_file:
        for i, mol in enumerate(sampled_mols):
            trj_file.write(f"{len(mol.atom_labels)}\n")
            trj_file.write(f"Sampled Configuration {i+1}\n")
            for label, coord in zip(mol.atom_labels, mol.coordinates):
                trj_file.write(f"{remove_digits(label)} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")

            


