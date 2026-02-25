import numpy as np
import itertools
from typing import Tuple, List, Optional, Dict
from dataclasses import dataclass
from mendeleev.fetch import fetch_table
import networkx as nx
from string import digits
import time 
import matplotlib.pyplot as plt

# Dataclass are used for simpler data structures


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
        Check if atom is involved in this bond
        """
        return atom_label in (self.atom1, self.atom2)
    
    def get_other_atom(self, atom_label: str) -> Optional[str]:
        """
        Get the other atom in the bond
        """
        if self.atom1 == atom_label:
            return self.atom2
        elif self.atom2 == atom_label:
            return self.atom1
        return None


@dataclass
class HBondConfiguration:
    """
    Represents a hydrogen bond configuration: Donor-H-Acceptor
    """
    donor_idx: int
    hydrogen_idx: int
    acceptor_idx: int
    angle: Optional[float] = None
    
    def is_valid(self, threshold: float = 150.0) -> bool:
        """
        Check if configuration meets angle threshold
        """
        return self.angle is not None and self.angle >= threshold


class GeometryCalculator:
    """
    Handles geometric calculations for molecules
    """
    
    @staticmethod
    def calculate_distance(coords1: np.ndarray, coords2: np.ndarray) -> float:
        """
        Calculate Euclidean distance between two points
        """
        return np.linalg.norm(coords1 - coords2)
    
    @staticmethod
    def calculate_angle(coords1: np.ndarray, coords2: np.ndarray, coords3: np.ndarray) -> float:
        """
        Calculate angle at coords2 formed by coords1-coords2-coords3
        Returns angle in degrees
        """
        v1 = coords1 - coords2
        v2 = coords3 - coords2
        
        cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        # Clamp to avoid numerical errors
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        angle_rad = np.arccos(cos_theta)
        return np.degrees(angle_rad)
    
    @staticmethod
    def calculate_bond_strength(distance: float, cov_radius_sum: float) -> float:
        """
        Calculate bond strength using exponential decay
        """
        return np.exp(-((distance / cov_radius_sum) - 1))


class BondClassifier:
    """
    Classifies bonds as covalent or hydrogen bonds
    """
    
    COVALENT_THRESHOLD = 0.7
    HYDROGEN_BOND_LOWER = 0.3
    HYDROGEN_BOND_UPPER = 0.7
    
    def __init__(self, molecule: 'Molecule'):
        self.molecule = molecule
        self.geometry = GeometryCalculator()
    
    def classify_bond(self, idx1: int, idx2: int) -> Optional[Bond]:
        """
        Classify bond between two atoms
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
            return Bond(label1, label2, strength)
        elif self.HYDROGEN_BOND_LOWER <= strength < self.HYDROGEN_BOND_UPPER:
            if self._is_hydrogen_bond(label1, label2):
                return Bond(label1, label2, strength)
        
        return None
    
    def _is_hydrogen_bond(self, label1: str, label2: str) -> bool:
        """Check if bond involves H and O (can be extended)"""
        elements = {
            Molecule._remove_digits_from_label(label1),
            Molecule._remove_digits_from_label(label2)
        }
        return "H" in elements and "O" in elements


class HydrogenBondAnalyzer:
    """Analyzes hydrogen bond configurations in molecules"""
    
    def __init__(self, molecule: 'Molecule'):
        self.molecule = molecule
        self.geometry = GeometryCalculator()
    
    def find_configurations(self) -> List[HBondConfiguration]:
        """Find all H-bond configurations: Donor-H-Acceptor"""
        configurations = []
        donor_indices = self._get_donor_indices()
        acceptor_indices = self._get_acceptor_indices()
        
        for donor_idx in donor_indices:
            configs = self._find_configs_for_donor(donor_idx, acceptor_indices)
            configurations.extend(configs)
        
        return configurations
    
    def _find_configs_for_donor(
        self, 
        donor_idx: int, 
        acceptor_indices: List[int]
    ) -> List[HBondConfiguration]:
        """Find all configurations for a specific donor"""
        configurations = []
        donor_label = self.molecule.atom_labels[donor_idx]
        donor_bonds = self.molecule.get_bonds_for_atom(donor_label, bond_type='covalent')
        
        for acceptor_idx in acceptor_indices:
            if donor_idx == acceptor_idx:
                continue
            
            acceptor_label = self.molecule.atom_labels[acceptor_idx]
            acceptor_bonds = self.molecule.get_bonds_for_atom(acceptor_label, bond_type='hydrogen')
            
            # Find common hydrogen atoms
            h_atoms = self._find_common_hydrogens(donor_bonds, acceptor_bonds, donor_label, acceptor_label)
            
            for h_label in h_atoms:
                h_idx = np.where(self.molecule.atom_labels == h_label)[0][0]
                config = HBondConfiguration(donor_idx, h_idx, acceptor_idx)
                config.angle = self._calculate_config_angle(config)
                configurations.append(config)
        
        return configurations
    
    def _find_common_hydrogens(
        self, 
        donor_bonds: List[Bond], 
        acceptor_bonds: List[Bond],
        donor_label: str,
        acceptor_label: str
    ) -> set:
        """Find hydrogen atoms common to donor and acceptor bonds"""
        donor_neighbors = {b.get_other_atom(donor_label) for b in donor_bonds}
        acceptor_neighbors = {b.get_other_atom(acceptor_label) for b in acceptor_bonds}
        
        common = donor_neighbors.intersection(acceptor_neighbors)
        common.discard(None)
        return common
    
    def _calculate_config_angle(self, config: HBondConfiguration) -> float:
        """Calculate angle for H-bond configuration"""
        coords_donor = self.molecule.coordinates[config.donor_idx]
        coords_h = self.molecule.coordinates[config.hydrogen_idx]
        coords_acceptor = self.molecule.coordinates[config.acceptor_idx]
        
        return self.geometry.calculate_angle(coords_donor, coords_h, coords_acceptor)
    
    def _get_donor_indices(self) -> List[int]:
        """Get indices of potential H-bond donors"""
        return [
            idx for idx, label in enumerate(self.molecule.atom_labels)
            if Molecule._remove_digits_from_label(label) in self.molecule.hbond_donors
        ]
    
    def _get_acceptor_indices(self) -> List[int]:
        """Get indices of potential H-bond acceptors"""
        return [
            idx for idx, label in enumerate(self.molecule.atom_labels)
            if Molecule._remove_digits_from_label(label) in self.molecule.hbond_acceptors
        ]


class Molecule:
    """Base class for representing a molecule"""
    
    _ptable = fetch_table("elements")
    # Pre-build lookup caches from the periodic table
    _mass_cache: Dict[str, float] = {}
    _covalent_radius_cache: Dict[str, float] = {}
    _vdw_radius_cache: Dict[str, float] = {}
    
    @classmethod
    def _init_caches(cls):
        """Build element lookup caches from periodic table (called once)."""
        if cls._mass_cache:  # already initialized
            return
        for _, row in cls._ptable.iterrows():
            sym = row["symbol"]
            cls._mass_cache[sym] = row["atomic_weight"]
            cov_r = row.get("covalent_radius_pyykko")
            cls._covalent_radius_cache[sym] = (cov_r / 100.0) if cov_r and not np.isnan(cov_r) else 0.0
            vdw_r = row.get("vdw_radius")
            cls._vdw_radius_cache[sym] = (vdw_r / 100.0) if vdw_r and not np.isnan(vdw_r) else 0.0
    
    hbond_donors = ["O", "N", "F"]
    hbond_acceptors = ["O", "N", "F"]
    
    def __init__(self, name: str = "Unnamed Molecule", logger: Optional[object] = None):
        """   
        Initializes the Molecule instance with name additionally the logger can be provided 
        """
        self.name = name
        self.atom_labels = np.array([], dtype=object)
        self.coordinates = np.empty((0, 3), dtype=np.float64)
        self.masses = np.array([], dtype=np.float64)
        self.vdw_radii = np.array([], dtype=np.float64)
        self.covalent_radii = np.array([], dtype=np.float64)

        # Initialize caches on first Molecule creation
        self._init_caches()


        self.charge = 0
        self.spin_mult = 1
        
        self._covalent_bonds: List[Bond] = []
        self._hydrogen_bonds: List[Bond] = []
        self._bonds_computed = False

        # Optional for Molecular Volume
        self.volume = None
    
    @staticmethod
    def from_labels_and_coords(
        atom_labels: List[str],
        coordinates: np.ndarray,
        name: Optional[str] = None
        ) -> 'Molecule':
        """Create Molecule from atom labels and coordinates"""
        molecule = Molecule(name or "Custom Molecule")
        molecule.add_atoms_batch(atom_labels, coordinates)
        return molecule

    # Static and class methods
    @staticmethod
    def _remove_digits_from_label(label: str) -> str:
        """Remove digits from atom label to get element symbol"""
        return ''.join(filter(str.isalpha, label))
    
    @classmethod
    def from_xyz(cls, xyz_content: str, name: Optional[str] = None) -> 'Molecule':
        """Create Molecule from XYZ file content"""
        lines = [line.strip() for line in xyz_content.strip().split("\n") if line.strip()]
        
        if len(lines) < 3:
            raise ValueError("Invalid XYZ file: insufficient lines")
        
        name = name or lines[1] or "Parsed Molecule"
        molecule = cls(name)
        
        atom_data = cls._parse_xyz_atoms(lines[2:])
        if atom_data:
            molecule.add_atoms_batch(atom_data['labels'], atom_data['coordinates'])
        
        return molecule
    
    @classmethod
    def _parse_xyz_atoms(cls, lines: List[str]) -> Dict:
        """Parse atom data from XYZ lines"""
        atoms = []
        coords = []
        
        for line in lines:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    atoms.append(parts[0])
                    coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
                except ValueError as e:
                    raise ValueError(f"Invalid coordinate values in line: '{line}'") from e
        
        if not atoms:
            return {}
        
        enumerated_labels = cls._enumerate_atom_labels(atoms)
        return {'labels': enumerated_labels, 'coordinates': np.array(coords)}
    
    @classmethod
    def _enumerate_atom_labels(cls, atom_labels: List[str]) -> List[str]:
        """Add numeric suffixes to atom labels"""
        element_count = {}
        enumerated = []
        
        for element in atom_labels:
            element_count[element] = element_count.get(element, 0) + 1
            enumerated.append(f"{element}{element_count[element]}")
        
        return enumerated
    
    @classmethod
    def get_covalent_radius(cls, element_label: str) -> float:
        """Get covalent radius for element in Angstroms (cached)."""
        element_symbol = cls._remove_digits_from_label(element_label)
        if element_symbol in cls._covalent_radius_cache:
            return cls._covalent_radius_cache[element_symbol]
        raise ValueError(f"Element {element_symbol} not found in periodic table")
    
    @classmethod 
    def get_vdw_radius(cls, element_label: str) -> float:
        """Get VDW radius for element in Angstroms (cached)."""
        element_symbol = cls._remove_digits_from_label(element_label)
        if element_symbol in cls._vdw_radius_cache:
            return cls._vdw_radius_cache[element_symbol]
        raise ValueError(f"Element {element_symbol} not found in periodic table")
    
    # Atom management
    def add_atom(self, atom_label: str, coordinates: List[float]) -> None:
        """Add single atom to molecule"""
        coords = np.array(coordinates, dtype=np.float64)
        if coords.shape != (3,):
            raise ValueError("Coordinates must be 3-dimensional (x, y, z)")
        
        self.atom_labels = np.append(self.atom_labels, atom_label)
        self.coordinates = np.vstack([self.coordinates, coords])
        
        element_symbol = self._remove_digits_from_label(atom_label)
        mass = self._get_atomic_mass(element_symbol)
        vdw_radius = self.get_vdw_radius(atom_label)
        self.vdw_radii = np.append(self.vdw_radii, vdw_radius)
        covalent_radius = self.get_covalent_radius(atom_label)
        self.covalent_radii = np.append(self.covalent_radii, covalent_radius)
        self.masses = np.append(self.masses, mass)
        
        self._bonds_computed = False
    
    def add_atoms_batch(self, atom_labels: List[str], coordinates: np.ndarray) -> None:
        """Add multiple atoms efficiently"""
        if len(atom_labels) != coordinates.shape[0]:
            raise ValueError("Number of labels must match number of coordinates")
        if coordinates.shape[1] != 3:
            raise ValueError("Coordinates must have shape (n_atoms, 3)")
        
        self.atom_labels = np.append(self.atom_labels, atom_labels)
        self.coordinates = np.vstack([self.coordinates, coordinates])
        
        masses = self._calculate_masses_batch(atom_labels)
        vdw_radii = np.array([self.get_vdw_radius(label) for label in atom_labels])
        covalent_radii = np.array([self.get_covalent_radius(label) for label in atom_labels])
        self.vdw_radii = np.append(self.vdw_radii, vdw_radii)
        self.covalent_radii = np.append(self.covalent_radii, covalent_radii)
        self.masses = np.append(self.masses, masses)
        
        self._bonds_computed = False
    
    def _calculate_masses_batch(self, atom_labels: List[str]) -> np.ndarray:
        """Calculate atomic masses for multiple atoms"""
        remove_digits = str.maketrans("", "", digits)
        elements = [label.translate(remove_digits) for label in atom_labels]
        return np.array([self._get_atomic_mass(elem) for elem in elements])
    
    def _get_atomic_mass(self, element_symbol: str) -> float:
        """Get atomic mass for element symbol (cached)."""
        if element_symbol in self._mass_cache:
            return self._mass_cache[element_symbol]
        print(f"Warning: Could not get mass for {element_symbol}, using default = 1")
        return 1.0
    
    # Coordinate access
    def get_coords_by_label(self, atom_label: str) -> np.ndarray:
        """Get coordinates of atom(s) by label"""
        indices = np.where(self.atom_labels == atom_label)[0]
        if indices.size == 0:
            raise ValueError(f"Atom label {atom_label} not found")
        return self.coordinates[indices]
    
    def get_coords_by_labels(self, atom_labels: List[str]) -> np.ndarray:
        """Get coordinates of multiple atoms"""
        coords_list = [self.get_coords_by_label(label) for label in atom_labels]
        return np.vstack(coords_list)
    
    def get_mass_by_label(self, atom_label: str) -> float:
        """Get atomic mass by atom label"""
        index = np.where(self.atom_labels == atom_label)[0]
        if index.size == 0:
            raise ValueError(f"Atom label {atom_label} not found")
        return self.masses[index[0]]
    
    # Bond analysis
    def compute_bonds(self) -> None:
        """Compute all bonds in molecule"""
        classifier = BondClassifier(self)
        covalent = []
        hydrogen = []
        
        n_atoms = len(self.atom_labels)
        for i, j in itertools.combinations(range(n_atoms), 2):
            bond = classifier.classify_bond(i, j)
            if bond:
                if bond.strength >= BondClassifier.COVALENT_THRESHOLD:
                    covalent.append(bond)
                else:
                    hydrogen.append(bond)
        
        self._covalent_bonds = covalent
        self._hydrogen_bonds = hydrogen
        self._bonds_computed = True
    
    def get_bonds(self, force_recompute: bool = False) -> Tuple[List[Bond], List[Bond]]:
        """Get covalent and hydrogen bonds"""
        if not self._bonds_computed or force_recompute:
            self.compute_bonds()
        return self._covalent_bonds, self._hydrogen_bonds
    
    def get_bonds_for_atom(self, atom_label: str, bond_type: str = 'covalent') -> List[Bond]:
        """Get all bonds involving specific atom"""
        if not self._bonds_computed:
            self.compute_bonds()
        
        bonds = self._covalent_bonds if bond_type == 'covalent' else self._hydrogen_bonds
        return [bond for bond in bonds if bond.involves(atom_label)]
    
    # Hydrogen bond analysis
    def find_hbond_configurations(self) -> List[HBondConfiguration]:
        """Find all hydrogen bond configurations"""
        if not self._bonds_computed:
            self.compute_bonds()
        
        analyzer = HydrogenBondAnalyzer(self)
        return analyzer.find_configurations()
    
    def get_valid_hbond_configurations(
        self, 
        angle_threshold: float = 150.0
        ) -> List[HBondConfiguration]:
        """
        Get hydrogen bond configurations meeting angle criterion
        """
        configs = self.find_hbond_configurations()
        return [c for c in configs if c.is_valid(angle_threshold)]
    
    def get_hbond_donor_environments(self) -> List["SubMolecule"]:
        """
        Get local environments around H-bond donors
        """
        if not self._bonds_computed:
            self.compute_bonds()
        
        donor_indices = [
            idx for idx, label in enumerate(self.atom_labels)
            if self._remove_digits_from_label(label) in self.hbond_donors
        ]
        
        submolecules = []
        for donor_idx in donor_indices:
            donor_label = self.atom_labels[donor_idx]
            bonds = self.get_bonds_for_atom(donor_label, 'covalent')
            
            neighbor_labels = [b.get_other_atom(donor_label) for b in bonds]
            neighbor_indices = [
                np.where(self.atom_labels == label)[0][0] 
                for label in neighbor_labels if label
            ]
            
            indices = [donor_idx] + neighbor_indices
            submol = SubMolecule.from_parent_fragment(self, indices, donor_idx)
            submolecules.append(submol)
        
        return submolecules

    def center_of_mass(self) -> np.ndarray:
        """Calculate center of mass of the molecule"""
        total_mass = np.sum(self.masses)
        com = np.sum(self.coordinates.T * self.masses, axis=1) / total_mass
        return com
    
    def geometric_center(self) -> np.ndarray:
        """Calculate geometric center of the molecule"""
        return np.mean(self.coordinates, axis=0)
    
    # Fragmentation
    def fragment_by_connectivity(self) -> List["SubMolecule"]:
        """Fragment molecule into connected components"""
        if not self._bonds_computed:
            self.compute_bonds()
        
        G = nx.Graph()
        G.add_nodes_from(self.atom_labels)
        
        for bond in self._covalent_bonds:
            G.add_edge(bond.atom1, bond.atom2)
        
        components = list(nx.connected_components(G))
        
        return [
            SubMolecule.from_parent_fragment(
                self,
                [idx for idx, label in enumerate(self.atom_labels) if label in comp],
                i + 1
            )
            for i, comp in enumerate(components)
        ]
    
    # Parsing
    def parse_psi4_xyz(self, psi4_xyz: str) -> 'Molecule':
        """
        Parse Psi4 XYZ format
        """
        lines = [line.strip() for line in psi4_xyz.strip().split("\n") if line.strip()]
        
        if len(lines) < 2:
            raise ValueError("Invalid Psi4 XYZ format")
        
        charge_spin = lines[0].split()
        if len(charge_spin) != 2:
            raise ValueError("First line must contain charge and spin multiplicity")
        
        self.charge = int(charge_spin[0])
        self.spin_mult = int(charge_spin[1])
        
        atom_data = self._parse_xyz_atoms(lines[1:])
        if atom_data:
            self.add_atoms_batch(atom_data['labels'], atom_data['coordinates'])
        
        return self

    def get_total_mass(self):
        """   
        Returns the total mass of the molecule
        """
        return np.sum(self.masses)
    
    def compute_volume_mc(self, n_points = 800000) -> float:
        """  
        Estimates the molecular volume using Monte Carlo Integration

        Here the Integral is estimated by I = V * (N_inside / N_total)
        where V is the volume of the bounding box, N_inside is the number of random points inside the molecule,
        """
        coords = self.coordinates 
        radii = self.vdw_radii 

        # Find bounding box with min/max coordinates, ad padding by max VDW radius
        min_coords = np.min(coords - radii[:, np.newaxis], axis=0)
        max_coords = np.max(coords + radii[:, np.newaxis], axis=0)
        box_dims = max_coords - min_coords
        box_volume = np.prod(box_dims)
        # Generate random points
        random_points = np.random.uniform(min_coords, max_coords, size=(n_points, 3))

        # Vectorize distance calculation
        coords_expanded = coords[:, np.newaxis, :]  # Shape (n_atoms, 1, 3)
        points_expanded = random_points[np.newaxis, :, :]  # Shape (1, n_points, 3)
        distances = np.linalg.norm(coords_expanded - points_expanded, axis=2)  # Shape (n_atoms, n_points)
        inside_mask = np.any(distances <= radii[:, np.newaxis], axis=0)  # Shape (n_points,)
        # Count points inside
        inside_count = np.sum(inside_mask)
        volume = (inside_count / n_points) * box_volume
        self.volume = volume
        return volume
    
    def get_coordinates_array(self) -> np.ndarray:
        """Get coordinates as numpy array"""
        return self.coordinates.copy()

    def get_atom_labels(self) -> List[str]:
        """Get list of atom labels"""
        return self.atom_labels.tolist() 
    
    def get_number_of_atoms(self) -> int:
        return len(self.coordinates)

    def copy(self) -> 'Molecule':
        """Create a lightweight copy of the molecule.
        
        Only copies mutable data (coordinates) as new arrays.
        Immutable/shared data (atom_labels, masses, radii) are shared references
        since they don't change during BHMC operations.
        """
        new_mol = Molecule.__new__(Molecule)
        new_mol.name = self.name
        new_mol.atom_labels = self.atom_labels  # shared (immutable during BHMC)
        new_mol.coordinates = self.coordinates.copy()  # deep copy - this is what changes
        new_mol.masses = self.masses  # shared
        new_mol.vdw_radii = self.vdw_radii  # shared
        new_mol.covalent_radii = self.covalent_radii  # shared
        new_mol.charge = self.charge
        new_mol.spin_mult = self.spin_mult
        new_mol._covalent_bonds = self._covalent_bonds
        new_mol._hydrogen_bonds = self._hydrogen_bonds
        new_mol._bonds_computed = self._bonds_computed
        new_mol.volume = self.volume
        return new_mol

    def deepcopy(self) -> 'Molecule':
        """Create a full deep copy of the molecule"""
        import copy
        return copy.deepcopy(self)
    
    def to_xyz_string(self) -> str:
        """Convert molecule to XYZ format string"""
        lines = [str(len(self.atom_labels)), self.name]
        for label, coord in zip(self.atom_labels, self.coordinates):
            lines.append(f"{self._remove_digits_from_label(label)} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}")
        return "\n".join(lines)

    def get_average_covalent_radii(self) -> List[float]:
        """   
        Calculate an Average Covalent Radius for each atom based on its covalent radius
        """
        # Count the different elemnt types
        element_counts = {}
        for label in self.atom_labels:
            element = self._remove_digits_from_label(label)
            element_counts[element] = element_counts.get(element, 0) + 1
        # Calculate average covalent radius for each element
        average_radii = {}
        for element, count in element_counts.items():
            radius = self.get_covalent_radius(element)
            average_radii[element] = radius / count
        # Map back to atom labels
        average_radii_list = []
        for label in self.atom_labels:
            element = self._remove_digits_from_label(label)
            average_radii_list.append(average_radii[element])
        return average_radii_list        


class SubMolecule(Molecule):
    """Represents a fragment of a larger molecule"""
    
    def __init__(self, name: str = "Unnamed Submolecule", parent: Optional[Molecule] = None):
        super().__init__(name)
        self.parent = parent
        self.fragment_index: Optional[int] = None
        self.volume: Optional[float] = None
    
    def __len__(self):
        """Return the Number of atoms"""
        return len(self.coordinates)

 

    @classmethod
    def from_parent_fragment(
        cls,
        parent_molecule: Molecule,
        atom_indices: List[int],
        fragment_index: int
    ) -> "SubMolecule":
        """
        Create SubMolecule from parent molecule and atom indices
        """
        submolecule = cls(name=f"{parent_molecule.name}_frag{fragment_index}")
        submolecule.atom_labels = parent_molecule.atom_labels[atom_indices].copy()
        submolecule.coordinates = parent_molecule.coordinates[atom_indices].copy()
        submolecule.masses = parent_molecule.masses[atom_indices].copy()
        submolecule.parent = parent_molecule
        submolecule.parent_indices = atom_indices.copy()
        submolecule.fragment_index = fragment_index
        submolecule.vdw_radii = parent_molecule.vdw_radii[atom_indices].copy()
        submolecule.covalent_radii = parent_molecule.covalent_radii[atom_indices].copy()
        return submolecule
    
    def get_atom_labels(self) -> List[str]:
        """Get list of atom labels"""
        return self.atom_labels.tolist()
    
    def get_index_in_parent(self) -> List[int]:
        """Get indices in parent molecule"""
        if self.parent is None:
            raise ValueError("SubMolecule has no parent")
        
        if self.parent_indices is not None:
            return self.parent_indices.copy()

        parent_labels = self.parent.atom_labels.tolist()
        indices = [parent_labels.index(label) for label in self.atom_labels]
        return indices 
    def compute_volume_mc(self, n_points = 800000) -> float:
        """   
        Estimate molecular volume using Monte Carlo Integration

        In the current implementation we use 800,000 random points for a good balance between accuracy and speed.
        """
        coords = self.coordinates 
        radii = self.vdw_radii

        # Find bounding box with min/max coordinates, ad padding by max VDW radius
        min_coords = np.min(coords - radii[:, np.newaxis], axis=0)
        max_coords = np.max(coords + radii[:, np.newaxis], axis=0)
        box_dims = max_coords - min_coords
        box_volume = np.prod(box_dims)

        # Generate random points 
        random_points = np.random.uniform(min_coords, max_coords, size=(n_points, 3))

        # Vectorizes distance calculation
        coords_expanded = coords[:, np.newaxis, :]  # Shape (n_atoms, 1, 3)
        points_expanded = random_points[np.newaxis, :, :]  # Shape (1, n_points, 3)

        distances = np.linalg.norm(coords_expanded - points_expanded, axis=2)  # Shape (n_atoms, n_points)

        inside_mask = np.any(distances <= radii[:, np.newaxis], axis=0)  # Shape (n_points,)

        # Count points inside
        inside_count = np.sum(inside_mask)
        volume = (inside_count / n_points) * box_volume

        self.volume = volume
        return volume
    
    # ===================== Submolecule Testing Functions ====================== 

    def check_convergence_and_speed_mc(self, n_points_list: List[int] = [1000, 5000, 10000, 50000, 100000, 200000, 500000, 800000]) -> Dict[int, float]:
        """   
        Test convergence and speed of Monte Carlo volume estimation
        """
        times = []
        volumes = []
        # Plot two subfigures: speed vs n_points and volume v.s n_points
        for n_points in n_points_list:
            start_time = time.time()
            volume = self.compute_volume_mc(n_points=n_points)
            end_time = time.time()
            elapsed = end_time - start_time

            times.append(elapsed)
            volumes.append(volume)
        fig = plt.figure(figsize=(12, 5))
        axs = fig.subplots(1, 2)
        axs[0].plot(n_points_list, times, marker='o')
        axs[0].set_xlabel("Number of Points")
        axs[0].set_ylabel("Computation Time (s)")
        axs[0].set_title("MC Volume Estimation Speed") 
        axs[0].grid(True)
        axs[1].plot(n_points_list, volumes, marker='o', color='orange')
        axs[1].set_xlabel("Number of Points")
        axs[1].set_ylabel("Estimated Volume (Å³)")
        axs[1].set_title("MC Volume Estimation Convergence")
        axs[1].grid(True)
        plt.tight_layout()
        plt.savefig("figures/mc_volume_convergence_speed.png")