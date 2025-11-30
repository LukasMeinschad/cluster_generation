import numpy as np
import itertools
from typing import Tuple, List, Optional # for type hints
from mendeleev.fetch import fetch_table
import networkx as nx
from string import digits
import itertools

class Molecule:
    """ 
    Base Class for representing a molecule
    """ 

    # Elements table 
    ptable = fetch_table("elements")
 

    # Internal list of hydrogen bond donors and acceptors
    hbond_donors = ["O", "N", "F"] 
    hbond_acceptors = ["H", "D"]
    
    def __init__(self, name: str = "Unnamed Molecule"):
        """ 
        Initializes a Molecule object
        """
        self.name = name
        self.atom_labels = np.array([],dtype=object) # empty array for element labels
        self.coordinates = np.empty((0,3), dtype=np.float64) # Empty 2D array
        self.masses = None 
        
        #TODO build something to automatically take care of spin and charge or ask for user input here
        self.charge = 0
        self.spin_mult = 1

        # Bonds
        self.covalent_bonds = []
        self.hydrogen_bonds = []

    @staticmethod 
    def _remove_digits_from_label(label: str) -> str:
        """ 
        Function to Remove digits from atom label to get element symbols
        """ 
        return ''.join(filter(str.isalpha, label))
    
    def get_atomic_mass(self,element_symbol: str): 
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
            return 1.0
    

    @classmethod
    def from_xyz(cls, xyz_content: str, name: Optional[str] = None) -> 'Molecule':
        """
        Create a Molecule from XYZ file content.
        
        Args:
            xyz_content (str): Content of the XYZ file as a string
            name (str, optional): Name for the molecule. If None, uses comment line from XYZ.
            
        Returns:
            Molecule: A new Molecule instance
        """
        lines = [line.strip() for line in xyz_content.strip().split("\n") if line.strip()]
        
        if len(lines) < 3:
            raise ValueError("Invalid XYZ file: insufficient lines")
        
        # Use provided name or comment line from XYZ
        if name is None:
            name = lines[1] or "Parsed Molecule"
        
        molecule = cls(name)
        
        # Parse atoms and coordinates
        atom_labels = []
        coordinates_list = []
        
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                element = parts[0]
                try:
                    coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                    atom_labels.append(element)
                    coordinates_list.append(coords)
                except ValueError as e:
                    raise ValueError(f"Invalid coordinate values in line: '{line}'") from e
                
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

    def get_coords_by_label(self,atom_label: str) -> np.ndarray:
        """ 
        Get coordinates of atoms by their label
        """
        indices = [i for i, label in enumerate(self.atom_labels) if label == atom_label]
        if not indices:
            raise ValueError(f"Atom label {atom_label} not found in molecule")
        return self.coordinates[indices]


    @classmethod
    def get_covalent_radius(cls,element_label: str) -> float:
        """ 
        Function that gives back the covalent radius for a given molecule using mendeleev 
        """
        element_symbol = cls._remove_digits_from_label(element_label)
        row = cls.ptable[cls.ptable["symbol"]==element_symbol]
        if not row.empty:
            radius = row.iloc[0]["covalent_radius_pyykko"]
            return radius / 100  # Convert pm to angstroms
        else:
            raise ValueError(f"Element {element} not found in periodic table")

    def get_bonds(self) -> Tuple[List[Tuple], List[Tuple]]:
        """ 
        Identify covalent and hydrogen bonds in the molecule
        """
        covalent_bonds = []
        hydrogen_bonds = []
        n_atoms = len(self.atom_labels)
        for i,j in itertools.combinations(range(n_atoms),2):
            distance = np.linalg.norm(self.coordinates[i] - self.coordinates[j])
            cov_radius_i = self.get_covalent_radius(self.atom_labels[i])
            cov_radius_j = self.get_covalent_radius(self.atom_labels[j])
            C_ij = cov_radius_i + cov_radius_j 
            degree = np.exp(-((distance / C_ij) - 1))

            if degree >= 0.7:
                covalent_bonds.append((self.atom_labels[i], self.atom_labels[j], degree))
            elif 0.3 <= degree < 0.7:
                # TODO hydrogen bonds only if H is involved
                labels = {self.atom_labels[i], self.atom_labels[j]}
                # Remove digits from labels
                labels = {self._remove_digits_from_label(label) for label in labels}
                if "H" in labels and "O" in labels:
                    hydrogen_bonds.append((self.atom_labels[i], self.atom_labels[j], degree))

        # Write Bonds to Molecule object
        self.covalent_bonds = covalent_bonds
        self.hydrogen_bonds = hydrogen_bonds
        return covalent_bonds, hydrogen_bonds
    
    def hbond_donor_configurations(self) -> List["SubMolecule"]:
        """
        Function that searches for hydrogen bond donors in the system and then gives back their neigbouring 
        atoms

        Example would be:
        
        + H-O-H
        + H-N-H-H
        + H-O-C
        """
        donor_indices = self._get_hbond_donors()
        configurations = []
        for donor_idx in donor_indices:
            # Find covalently bonded neighbours
            neighbours = []
            for bond in self.covalent_bonds:
                if self.atom_labels[donor_idx] in bond:
                    # Get the other atom in the bond
                    other_atom = bond[0] if bond[1] == self.atom_labels[donor_idx] else bond[1]
                    neighbours.append(other_atom)
            # Create SubMolecule for this configuration
            indices = [donor_idx] + [np.where(self.atom_labels == n)[0][0] for n in neighbours]
            submol = SubMolecule.from_parent_fragment(self, indices, fragment_index=donor_idx)
            configurations.append(submol)
        return configurations


    def _get_hbond_donors(self):
        """ 
        Helper function that searches for hydrogen bond donors in the molecule

        Returns:
            Indices of hydrogen bond donor atoms
        """
        donor_indices = []
        for idx, label in enumerate(self.atom_labels):
            element_symbol = self._remove_digits_from_label(label)
            if element_symbol in self.hbond_donors:
                donor_indices.append(idx)
        return donor_indices



    def fragment_by_connectivity(self) -> List["SubMolecule"]:
        """ 
        Fragments the molecule into connected components based on covalent bonds
        
        Returns:
            List of SubMolecule objects representing each connected component
        """
        # Check if bonds are already computed
        if not self.covalent_bonds:
            self.get_bonds()
        
        # Build Graph 
        G = nx.Graph() 
        G.add_nodes_from(self.atom_labels)
        for bond in self.covalent_bonds:
            # Add edges to Graph
            G.add_edge(bond[0], bond[1])

        # Find connected components
        connected_components = list(nx.connected_components(G))

        # Create the Submolecule objects
        submolecules = []
        for i,component in enumerate(connected_components):
            # Get indices of atoms in this components
            indices = [idx for idx, label in enumerate(self.atom_labels) if label in component]

            # Create SubMolecule from parent
            submol = SubMolecule.from_parent_fragment(self, indices, fragment_index=i+1)
            submolecules.append(submol)

        return submolecules


    def add_atom(self, atom_label: str, coordinates: List[float]) -> None:
        """ 
        Add a single atom to the molecule
        """
        coords = np.array(coordinates, dtype=np.float64)
        if coords.shape != (3,):
            raise ValueError("Coordinates must be 3-dimensional (x, y, z)")
        self.atom_labels = np.append(self.atom_labels, atom_label)
        self.coordinates = np.vstack([self.coordinates, coords])
        # Update masses
        element_symbol = self._remove_digits_from_label(atom_label)
        mass = self.get_atomic_mass(element_symbol)
        if self.masses is None:
            self.masses = np.array([mass])
        else:
            self.masses = np.append(self.masses, mass)
    
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

    def get_mass_by_label(self,atom_label: str) -> float:
        """ 
        Get atomic mass of atom by its label
        """
        index = np.where(self.atom_labels == atom_label)[0]
        if index.size == 0:
            raise ValueError(f"Atom label {atom_label} not found in molecule")
        return self.masses[index[0]]

    def parse_psi4_xyz(self, psi4_xyz: str) -> None:
        """ 
        Function to parse the Psi4 xyz format and initialize a new molecule object
        
        Coordinates are given in the following format:

        Charge SpinMultiplicity
        Element1   x1   y1   z2
        """
        lines = [line.strip() for line in psi4_xyz.strip().split("\n") if line.strip()]
        if len(lines) < 2:
            raise ValueError("Invalid Psi4 XYZ format: insufficient lines")
        
        # First line is charge and spin multiplicity
        charge_spin = lines[0].split()
        if len(charge_spin) != 2:
            raise ValueError("Invalid Psi4 XYZ format: first line must contain charge and spin multiplicity")
        self.charge = int(charge_spin[0])
        self.spin_mult = int(charge_spin[1])

        atom_labels = []
        coordinates_list = []
        for line in lines[1:]:
            parts = line.split()
            if len(parts) >= 4:
                element = parts[0]
                try:
                    coords = [float(parts[1]), float(parts[2]), float(parts[3])]
                    atom_labels.append(element)
                    coordinates_list.append(coords)
                except ValueError as e:
                    raise ValueError(f"Invalid coordinate values in line: '{line}'") from e
        # Count Atoms and enumerate
        element_count = {}
        atom_labels_enum = []
        for element in atom_labels:
            if element not in element_count:
                element_count[element] = 0
            element_count[element] +=1
            atom_labels_enum.append(f"{element}{element_count[element]}")
        if atom_labels:
            self.add_atoms_batch(atom_labels_enum,np.array(coordinates_list)) 
        return self 


class SubMolecule(Molecule):
    """   
    A subclass representing a submolecule within a larger molecule
    """
    def __init__(self,name: str = "Unnamed Submolecule", parent=None):
        super().__init__(name)
        self.parent = parent  # Reference to the parent Molecule object
        self.fragment_index = None  # Index of the fragment within the parent molecule

    @classmethod
    def from_parent_fragment(cls, parent_molecule: Molecule, atom_indices: List[int], fragment_index: int):
        """  
        Classmethod to generate a SubMolecule from a parent Molecule and atom indices
        """
        submolecule = cls(
            name=f"{parent_molecule.name}"
        )
        submolecule.atom_labels = parent_molecule.atom_labels[atom_indices].copy()
        submolecule.coordinates = parent_molecule.coordinates[atom_indices].copy()
        submolecule.masses = parent_molecule.masses[atom_indices].copy()
        submolecule.parent = parent_molecule
        submolecule.fragment_index = fragment_index
        return submolecule
    
    def get_atom_labels(self) -> List[str]:
        """ 
        Returns the list of atom labels in the submolecule
        """
        return self.atom_labels.tolist()
    
    def get_index_in_parent(self) -> List[int]:
        """ 
        Returns the indices of the submolecule's atom in the parent molecule
        """
        parent_labels = self.parent.atom_labels.tolist()
        indices = [parent_labels.index(label) for label in self.atom_labels]
        return indices
    
