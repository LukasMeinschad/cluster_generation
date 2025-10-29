import numpy as np
import itertools
from typing import Tuple, List
from mendeleev.fetch import fetch_table
import networkx as nx
from string import digits

class Molecule:
    

    # Elements table 
    ptable = fetch_table("elements")
    
    def __init__(self, name: str = "Unnamed Molecule"):
        """ 
        Initializes a Molecule object
        """
        self.name = name
        self.atom_labels = np.array([]) # empty array for element labels
        self.coordinates = np.empty((0,3), dtype=np.float64) # Empty 2D array
        self.masses = None 

    def get_atomic_mass(self,element_symbol):
        """ 
        Function to get the atomic_mass of a given element
        """
        try:
            row = self.ptable.loc[self.ptable["symbol"]==element_symbol]
            masses = row.get("atomic_weight")
            masses = masses.iloc[0]
            return masses
        except Exception as e:
            print(f"Warning: Could not get the mass for {element_symbol}, using Default = 1")
            return 1
    

    @classmethod
    def from_xyz(cls, xyz_content: str, name: str = None):
        """
        Create a Molecule from XYZ file content.
        
        Args:
            xyz_content (str): Content of the XYZ file as a string
            name (str, optional): Name for the molecule. If None, uses comment line from XYZ.
            
        Returns:
            Molecule: A new Molecule instance
        """
        lines = xyz_content.strip().split("\n")
        
        if len(lines) < 3:
            raise ValueError("Invalid XYZ file: insufficient lines")
        
        # Use provided name or comment line from XYZ
        if name is None:
            name = lines[1].strip() or "Parsed Molecule"
        
        molecule = cls(name)
        
        # Parse atoms and coordinates
        atom_labels = []
        coordinates_list = []
        
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                element = parts[0]
                coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                atom_labels.append(element)
                coordinates_list.append(coords)
      
        # Count Atoms and enumerate
        element_count = {}
        atom_labels_enum = []

        for element in atom_labels:
            if element not in element_count:
                element_count[element] = 0
            element_count[element] +=1
            atom_labels_enum.append(f"{element}{element_count[element]}")

        if atom_labels:
            molecule.add_atoms_batch(atom_labels_enum,np.array(coordinates_list))
        
        return molecule


    @classmethod
    def get_covalent_radius(cls,element):
        """ 
        Function that gives back the covalent radius for a given molecule using mendeleev 
        """
        ptable = fetch_table("elements")
        element = ''.join(filter(str.isalpha, element))
        row = ptable[ptable["symbol"]==element]
        if not row.empty:
            return row.iloc[0]["covalent_radius_pyykko"] / 100 # Angstrom
        else:
            raise ValueError(f"Element {element} not found in periodic table")


    def degree_of_covalence(self): 
        """ 
        Determines the degree of covalence of bonds and groups into covalent and hydrogen bonds
    
        Formula: exp(-[(r_ij)/(C_ij) - 1])
        where r_ij is the distance between atoms i and j 
        and C_ij is the sum of covalent radii of atoms i and j
        """
        cov_bonds = []
        hydrogen_bonds = []
    
        # Get all unique pairs of atom indices
        n_atoms = len(self.atom_labels)
        for i, j in itertools.combinations(range(n_atoms), 2):
            coords1 = self.coordinates[i]
            coords2 = self.coordinates[j]
            euclidean_distance = np.linalg.norm(coords1 - coords2)
            
            cov_radius1 = self.get_covalent_radius(self.atom_labels[i])
            cov_radius2 = self.get_covalent_radius(self.atom_labels[j])
            C_ij = cov_radius1 + cov_radius2
            
            degree = np.exp(-((euclidean_distance / C_ij) - 1))
            if degree >= 0.2 and degree < 0.7:
                hydrogen_bonds.append((self.atom_labels[i], self.atom_labels[j], degree))
            elif degree >= 0.7:
                cov_bonds.append((self.atom_labels[i], self.atom_labels[j], degree))

        # Filter hydrogen bonds to only consider O-H pairs
        hydrogen_bonds = [bond for bond in hydrogen_bonds 
                         if ('H' in bond[0] and 'O' in bond[1]) or ('H' in bond[1] and 'O' in bond[0])]
        
        return cov_bonds, hydrogen_bonds

    def find_submolecules(self, bonds):
        """ 
        Finds covalently connected components using networkx

        Args:
            bonds: List of bonds as tuples (atom_label1, atom_label2)

        Returns:
            List of Molecule object for each connected components
        """
        G = nx.Graph()

        for atom_label in self.atom_labels:
            G.add_node(atom_label)

        for bond in bonds:
            G.add_edge(bond[0], bond[1])

        # Find connected components 
        connected_components = list(nx.connected_components(G))

        # Create new molecule objects for each connected component
        submolecules = []

        for i,component in enumerate(connected_components):
            # Get indices of atoms in this components
            indices = [idx for idx, label in enumerate(self.atom_labels) if label in component]

            # Greate a new submolecule for this component
            submol = Molecule(name=f"{self.name}_fragment_{i+1}")
            submol.atom_labels = self.atom_labels[indices].copy()
            submol.coordinates = self.coordinates[indices].copy()
            submolecules.append(submol)

        return submolecules

    def add_atom(self, atom_label: str, coordinates):
        coords = np.array(coordinates, dtype=np.float64)
        if coords.shape != (3,):
            raise ValueError("Coordinates must be 3-dimensional (x, y, z)")
        self.atom_labels = np.append(self.atom_labels, atom_label)
        self.coordinates = np.vstack([self.coordinates, coords])
    
    def add_atoms_batch(self, atom_labels, coordinates,masses=None):
        if len(atom_labels) != coordinates.shape[0]:
            raise ValueError("Number of atom labels must match number of coordinate rows")
        if coordinates.shape[1] != 3:
            raise ValueError("Coordinates must have shape (n_atoms, 3)")
        self.atom_labels = np.append(self.atom_labels, atom_labels)
        self.coordinates = np.vstack([self.coordinates, coordinates])
        if masses is None:
            # Remove all the digits
            remove_digits = str.maketrans("","", digits)
            atom_labels_without_digits = [elem.translate(remove_digits) for elem in atom_labels]
            masses = [self.get_atomic_mass(elem) for elem in atom_labels_without_digits]
            self.masses = np.array(masses)
    
    
