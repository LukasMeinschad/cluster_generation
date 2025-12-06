import numpy as np
import itertools
from dataclasses import dataclass
from typing import Tuple, List, Optional # for type hints
from mendeleev.fetch import fetch_table
import networkx as nx
from string import digits


""" 
Note that we use dataclasses as a decorator to simplify classes that are mainly used to store data 
"""



@dataclass
class Bond:
    """   
    Represents a bond between two atoms
    """
    atom1: str 
    atom2: str 
    strength: float 

    def involves(self, atom_label: str) -> bool:
        """  
        Checks if a atom is involved in the bond
        """
        return atom_label in (self.atom1, self.atom2)
    
    def get_other_atom(self, atom_label: str) -> str: 
        """
        Get other atom in the bond 
        """
        if self.atom1 == atom_label:
            return self.atom2
        elif self.atom2 == atom_label:
            return self.atom1
        else:
            raise ValueError(f"Atom {atom_label} not involved in bond between {self.atom1} and {self.atom2}")

@dataclass 
class HBondConfiguration: 
    """  
    Represents a hydrogen bond configuration in the molecule: Donor-H...Acceptor
    """
    donor_idx: int 
    hydrogen_idx: int 
    acceptor_idx: int 
    angle: Optional[float] = None

    def is_valid(self, threshold: float = 150.0) -> bool:
        """  
        Check if configuration is valid based on an angle threshold
        """ 
        return self.angle is not None and self.angle >= threshold
    
class GeometryCalculator:
    """  
    Class that handles geometric calculations for molecules
    """
    @staticmethod 
    def calculate_distance(coords1: np.ndarray, coords2: np.ndarray) -> float: 
        """
        Function that calculates the Euclidean distance between two points 
        """ 
        return np.linalg.norm(coords1 - coords2)

    @staticmethod 
    def calculate_angle(coords1: np.ndarray, coords2: np.ndarray, coords3: np.ndarray) -> float: 
        """  
        Calculates the angle at coords2 formed by coords1 - coords2 - coords3
        """
        v1 = coords1 - coords2
        v2 = coords3 - coords2

        cos_theta = np.dot(v1,v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)) 
        # Clamp for numerical stability
        cos_theta = np.clip(cos_theta, -1.0, 1.0)  
        angle = np.arccos(cos_theta)
        return np.degrees(angle)

    @staticmethod
    def calculate_bond_strength(distance: float, cov_radius_sum: float) -> float: 
        """
        Calculates the bond strength based on the distance and sum of covalent radii
        """ 
        return np.exp(-((distance / cov_radius_sum) - 1))
    
class BondClassifier:
    """  
    Classifies bonds as covalent or hydrogen bonds
    """
    
    COVALENT_THRESHOLD = 0.7
    HYDROGEN_BOND_LOWER = 0.3
    HYDROGEN_BOND_UPPER = 0.7

    def __init__(self, molecule: " Molecule"):
        """  
        Initializes the BondClassifier with a molecule
        """
        self.molecule = molecule
        self.geometry = GeometryCalculator()

    def classify_bond(self, idx1: int, idx2: int) -> Optional[Bond]: 
        """  
        Classify the bond between two atoms
        Returns Bond object with type or None if no bond
        """
        label1 = self.molecule.atom_labels[idx1]
        label2 = self.molecule.atom_labels[idx2] 
        distance = self.geometry.calculate_distance(
            self.molecule.coordinates[idx1],
            self.molecule.coordinates[idx2]
        )
        radius_sum = (
            self.molecule.get_covalent_radius(label1) +
            self.molecule.get_covalent_radius(label2)
        )
        strength = self.geometry.calculate_bond_strength(distance, radius_sum)
        if strength >= self.COVALENT_THRESHOLD:
            return Bond(atom1=label1, atom2=label2, strength=strength)
        elif self.HYDROGEN_BOND_LOWER <= strength < self.HYDROGEN_BOND_UPPER:
            if self._is_hydrogen_bond(label1, label2):
                return Bond(atom1=label1, atom2=label2, strength=strength)

    def _is_hydrogen_bond(self, label1: str, label2: str) -> bool:
        """   
        Check if bond involves H and O/N/F for hydrogen bond
        """
        elements = {
            MMolecule._remove_digits_from_label(label1),
            MMolecule._remove_digits_from_label(label2)
        }
        return "H" in elements and any(elem in elements for elem in ["O", "N", "F"])
        
class HydrogenBondAnalyzer:
    """
    Analyzes hydrogen bond configurations in a molecule
    """ 
    def __init__(self,molecule: "Molecule"):
        """
        Initializes the HydrogenBondAnalyzer with a molecule
        """
        self.molecule = molecule
        self.geometry = GeometryCalculator()

    def find_configurations(self) -> List[HBondConfiguration]:
        """   
        Finds all possible H-Bond configurations: Donor-H...Acceptor
        """
        configurations = []
        donor_indices = self._get_donor_indices() 
        acceptor_indices = self._get_acceptor_indices()

    #TODO finish restructuring here

    def _get_donor_indices(self) -> List[int]: 
        """   
        Get the indices of potential H-bond donors
        """
        return [ 
            idx for idx, label in enumerate(self.molecule.atom_labels)
            if Molecule._remove_digits_from_label(label) in self.molecule.hbond_donors
        ]
    
    def _get_acceptor_indices(self) -> List[int]: 
        """   
        Get the indices of potential H-bond acceptors
        """
        return [ 
            idx for idx, label in enumerate(self.molecule.atom_labels)
            if Molecule._remove_digits_from_label(label) in self.molecule.hbond_accecptors
        ]

class Molecule:
    """ 
    Base Class for representing a molecule
    """ 

    # Elements table 
    ptable = fetch_table("elements")
 

    # Internal list of hydrogen bond donors and acceptors
    hbond_donors = ["O", "N", "F"] 
    hbond_accecptors = ["O", "N", "F"]
    
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

        # Hydrogen bond configurations
        self.hbond_configurations_list = []
        self.hbond_donor_configurations_list = []
        self.hbond_configurations_angles = []
        self.valid_hbond_configurations = []


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

    def get_coords_by_labels(self,atom_labels: List[str]) -> np.ndarray:
        """ 
        Get coordinates of multiple atoms by their labels
        """
        coords_list = []
        for label in atom_labels:
            coords = self.get_coords_by_label(label)
            coords_list.append(coords)
        return np.vstack(coords_list) 


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
    
    def hbond_configurations(self) -> List["SubMolecule"]:
        """
        Find possible HBond Configurations in a molecule

        For this we search for triples of atoms were:

        + Acceptor ----- H-Donor
        + Donor and H are covalently bonded
        """
        # First get all H-bond donors and acceptors
        donor_indices = self._get_hbond_donors()
        acceptor_indices = self._get_hbond_acceptors()
        configurations = []
        # Iterate over all Triples
        for donor_idx in donor_indices:
            for acceptor_idx in acceptor_indices:
                # Case that donor and acceptor are the same atom
                if donor_idx == acceptor_idx:
                    continue  
                # Obtain all covalent bonds of the donor
                cov_bonds = self._get_covalent_bonds_label(self.atom_labels[donor_idx])
                # Obtain all hydrogen bonds of the acceptor
                h_bonds = self._get_hydrogen_bonds_label(self.atom_labels[acceptor_idx])
                # Check if there is a common H atom
                for cov_bond in cov_bonds:
                    for h_bond in h_bonds:
                        # Find common atom
                        common_atoms = set(cov_bond).intersection(set(h_bond))
                        # Check that common atom is not the donor or acceptor itself
                        common_atoms.discard(self.atom_labels[donor_idx])
                        common_atoms.discard(self.atom_labels[acceptor_idx])
                        if common_atoms:
                            common_atom = common_atoms.pop()
                            # Create SubMolecule for this configuration in the form Donor-H  Acceptor
                            indices = [donor_idx] + [np.where(self.atom_labels == common_atom)[0][0]] + [acceptor_idx]
                            submol = SubMolecule.from_parent_fragment(self, indices, fragment_index=donor_idx)
                            configurations.append(submol)
        self.hbond_configurations_list = configurations
        return configurations  

    def get_angles_of_hbond_configurations(self) -> List[float]:
        """
        Calculates the angles in the hydrogen bond configurations of the molecule

        The hbond configurations is again defined as:

        + Acceptor ----- H-Donor

        And is invoked with the function hbond_configurations()
        """
        if not self.hbond_configurations_list:
            self.hbond_configurations_list = self.hbond_configurations()

        angles = []
        for config in self.hbond_configurations_list:
            if len(config.atom_labels) != 3:
                continue # Skip invalid configurations
            donor_label = config.atom_labels[0]
            hydrogen_label = config.atom_labels[1]
            acceptor_label = config.atom_labels[2]
            donor_coords = config.get_coords_by_label(donor_label)[0]
            hydrogen_coords = config.get_coords_by_label(hydrogen_label)[0]
            acceptor_coords = config.get_coords_by_label(acceptor_label)[0]
            
            v = donor_coords - hydrogen_coords
            w = acceptor_coords - hydrogen_coords
            cos_theta = np.dot(v, w) / (np.linalg.norm(v) * np.linalg.norm(w))
            angle = np.arccos(cos_theta) * (180.0 / np.pi)  # Convert to degrees
            angles.append(angle)
        
        # Add the angles to the molecule object
        self.hbond_configurations_angles = angles

        return angles
    
    def get_valid_hydrogen_bond_configurations(self, angle_threshold: float = 150.0) -> List["Submolecule"]:
        """
        Function that filters the hydrogen bond configurations based on an angle criterion 
        The angle is defined as the angle between the Donor-Hydrogen-Acceptor atoms
        If the angle is larger than the threshold, the configuration is considered valid
        """
        if not self.hbond_configurations_list:
            self.hbond_configurations_list = self.hbond_configurations()
        
        if not self.hbond_configurations_angles:
            self.hbond_configurations_angles = self.get_angles_of_hbond_configurations()
        
        valid_configurations = []
        for config, angle in zip(self.hbond_configurations_list, self.hbond_configurations_angles):
            if angle >= angle_threshold:
                valid_configurations.append(config)
        
        self.valid_hbond_configurations = valid_configurations
        return valid_configurations

    def _get_angle_of_valid_hbond_configurations(self) -> List[float]:
        """    
        Helper function that returns the angles of the valid hydrogen bond configurations
        """ 
        angles = []
        for config in self.valid_hbond_configurations:
            if len(config.atom_labels) != 3:
                continue # Skip invalid configurations
            donor_label = config.atom_labels[0]
            hydrogen_label = config.atom_labels[1]
            acceptor_label = config.atom_labels[2]
            donor_coords = config.get_coords_by_label(donor_label)[0]
            hydrogen_coords = config.get_coords_by_label(hydrogen_label)[0]
            acceptor_coords = config.get_coords_by_label(acceptor_label)[0]
            
            v = donor_coords - hydrogen_coords
            w = acceptor_coords - hydrogen_coords
            cos_theta = np.dot(v, w) / (np.linalg.norm(v) * np.linalg.norm(w))
            angle = np.arccos(cos_theta) * (180.0 / np.pi)  # Convert to degrees
            angles.append(angle)
        return angles
    
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

        self.hbond_donor_configurations_list = configurations
        return configurations


    def get_bond_vector(self, bond: Tuple[str, str]) -> np.ndarray:
        """ 
        Function to get the bond vector between a bond-tuple
        """
        coords1, coords2 = self._get_bond_coordinates(bond)
        bond_vector = coords2 - coords1
        return bond_vector

    def _get_covalent_bonds_label(self, label: str) -> List[Tuple[str, str]]:
        """
        Helper function that returns all covalent bonds where the given label is involved
        """
        bonds = []
        for bond in self.covalent_bonds:
            if label in bond:
                bonds.append((bond[0], bond[1]))
        return bonds
    
    def _get_hydrogen_bonds_label(self, label: str) -> List[Tuple[str, str]]:
        """
        Helper function that returns all hydrogen bond where the given label is involved
        """
        bonds = []
        for bond in self.hydrogen_bonds:
            if label in bond:
                bonds.append((bond[0], bond[1]))
        return bonds

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

    def _get_hbond_acceptors(self):
        """
        Helper function that searches for hydrogen bond acceptors in the molecule
        """
        acceptor_indices = []
        for idx, label in enumerate(self.atom_labels):
            element_symbol = self._remove_digits_from_label(label)
            if element_symbol in self.hbond_accecptors:
                acceptor_indices.append(idx)
        return acceptor_indices

    def _get_bond_coordinates(self, bond: Tuple[str, str]) -> Tuple[np.ndarray, np.ndarray]:
        """ 
        Helper function to get the coordinates of two atoms involved in a bond
        """
        atom1_label, atom2_label = bond[0], bond[1]
        coords1 = self.get_coords_by_label(atom1_label)
        coords2 = self.get_coords_by_label(atom2_label)
        return coords1, coords2
    



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
    

    
