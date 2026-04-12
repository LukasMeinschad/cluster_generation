import numpy as np
from typing import Tuple, List, Optional, Dict
from dataclasses import dataclass
from mendeleev.fetch import fetch_table
import networkx as nx
from string import digits
import time 
import matplotlib.pyplot as plt
from numba import njit, prange
# Dataclass are used for simpler data structures

from geometry import GeometryOps

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
    donor_acceptor_distance: Optional[float] = None 
    
    def is_valid(self, threshold: float = 150.0, max_distance: float = 3.5) -> bool:
        """
        Check if configuration meets angle threshold
        """
        angle_ok = self.angle is not None and self.angle >= threshold
        distance_ok = self.donor_acceptor_distance is not None and self.donor_acceptor_distance <= max_distance
        return angle_ok and distance_ok



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
    """Bond thresholds and numba helpers for pairwise bond classification."""
    
    COVALENT_THRESHOLD = 0.7
    HYDROGEN_BOND_LOWER = 0.3
    HYDROGEN_BOND_UPPER = 0.7
    
    @staticmethod
    @njit(fastmath=True)
    def _pairwise_bond_scan(
        coords: np.ndarray,
        cov_radii: np.ndarray,
        is_h: np.ndarray,
        is_o: np.ndarray,
        cov_thr: float,
        hb_lo: float,
        hb_hi: float
        ):
        """Numba kernel that classifies all atom pairs into bond buckets."""
        n = coords.shape[0]
        max_pairs = n *(n-1) // 2

        cov_i = np.empty(max_pairs, dtype=np.int64)
        cov_j = np.empty(max_pairs, dtype=np.int64)
        cov_s = np.empty(max_pairs, dtype=np.float64)
        n_cov = 0

        hb_i = np.empty(max_pairs, dtype=np.int64)
        hb_j = np.empty(max_pairs, dtype=np.int64)
        hb_s = np.empty(max_pairs, dtype=np.float64)
        n_hb = 0

        for i in range(n-1):
            xi, yi, zi = coords[i,0], coords[i,1], coords[i,2]
            ri = cov_radii[i]
            hi = is_h[i]
            oi = is_o[i]

            for j in range(i+1, n):
                dx = xi - coords[j,0]
                dy = yi - coords[j,1]
                dz = zi - coords[j,2]
                dist = np.sqrt(dx*dx + dy*dy + dz*dz)
                rsum = ri + cov_radii[j]
                if rsum <= 1e-14:
                    continue

                # calculate bond strength
                s = np.exp(-((dist/rsum) - 1))
                if s >= cov_thr:
                    cov_i[n_cov] = i
                    cov_j[n_cov] = j
                    cov_s[n_cov] = s
                    n_cov += 1
                elif hb_lo <= s < hb_hi:
                    # H-O only
                    if (hi and is_o[j]) or (oi and is_h[j]):
                        hb_i[n_hb] = i
                        hb_j[n_hb] = j
                        hb_s[n_hb] = s
                        n_hb += 1
    
        return (
            cov_i[:n_cov], cov_j[:n_cov], cov_s[:n_cov],
            hb_i[:n_hb], hb_j[:n_hb], hb_s[:n_hb]
        )

class HydrogenBondAnalyzer:
    """Analyzes hydrogen bond configurations in molecules"""
    
    def __init__(self, molecule: 'Molecule'):
        self.molecule = molecule
        self.geometry = GeometryCalculator()
    
    def find_configurations(self) -> List[HBondConfiguration]:
        """  
        Find all H-bond configurations of type Donor-H-Acceptor in molecule
        """
        # Ensure the bonds exist
        covalent_bonds, hydrogen_bonds = self.molecule.get_bonds()

        atom_labels = self.molecule.atom_labels
        coords = self.molecule.coordinates
        n_atoms = len(atom_labels)

        # One-time maps and masks
        label_to_idx = {label: i for i, label in enumerate(atom_labels)}
        elements = np.array(
            [Molecule._remove_digits_from_label(lbl) for lbl in atom_labels], dtype=object
        )
        donor_indices = np.where(np.isin(elements, self.molecule.hbond_donors))[0]
        acceptor_indices = np.where(np.isin(elements, self.molecule.hbond_acceptors))[0]
        hydrogen_indices = set(np.where(elements == "H")[0])

        # Build neighbor sets once (index based, no repeated lookups)
        cov_neighbors = [set() for _ in range(n_atoms)]
        for b in covalent_bonds:
            i = label_to_idx[b.atom1]
            j = label_to_idx[b.atom2]
            cov_neighbors[i].add(j)
            cov_neighbors[j].add(i)

        h_neighbors = [set() for _ in range(n_atoms)]
        for b in hydrogen_bonds:
            i = label_to_idx[b.atom1]
            j = label_to_idx[b.atom2]
            h_neighbors[i].add(j)
            h_neighbors[j].add(i)
        
        # Initialize configurations list
        configurations: List[HBondConfiguration] = []

        for donor_idx in donor_indices:
            donor_cov = cov_neighbors[donor_idx]
            if not donor_cov:
                continue

            donor_coords = coords[donor_idx]

            for acceptor_idx in acceptor_indices:
                if donor_idx == acceptor_idx:
                    continue

                # H must be covalently bound to donor and hydrogen bonded to acceptor
                common_h = donor_cov.intersection(h_neighbors[acceptor_idx]).intersection(hydrogen_indices)  
                if not common_h:
                    continue
                acceptor_coords = coords[acceptor_idx]
                da_dist = self.geometry.calculate_distance(donor_coords, acceptor_coords) 

                for h_idx in common_h:
                    h_coord = coords[h_idx]
                    angle = self.geometry.calculate_angle(donor_coords, h_coord, acceptor_coords) 

                    configurations.append(
                        HBondConfiguration(donor_idx, h_idx, acceptor_idx, angle=angle, donor_acceptor_distance=da_dist)
                    )  

        return configurations


    


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

    ELEMENT_COLORS = {
        'H': 'white', 'He': 'cyan',
        'C': 'gray', 'N': 'blue', 'O': 'red', 'F': 'green',
        'Ne': 'cyan', 'Na': 'purple', 'Mg': 'darkgreen',
        'P': 'orange', 'S': 'yellow', 'Cl': 'lime',
        'Ar': 'cyan', 'K': 'purple', 'Ca': 'darkgreen',
        'Br': 'darkred', 'I': 'darkviolet',
    }

    ELEMENT_SIZES = {
        'H': 50, 'He': 60,
        'C': 100, 'N': 110, 'O': 120, 'F': 100,
        'Ne': 80, 'Na': 140, 'Mg': 130,
        'P': 130, 'S': 140, 'Cl': 130,
        'Ar': 120, 'K': 160, 'Ca': 150,
        'Br': 140, 'I': 150
    }


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
    def get_element(self, atom_label: str) -> str:
        """  
        Get element symbol from atom label (strips digits)
        """
        return self._remove_digits_from_label(atom_label)  
    
    def get_atom_colors(self) -> list:
        """  
        Get CPK colors for all the atoms, based on element type
        """
        return [
            self.ELEMENT_COLORS.get(self._remove_digits_from_label(label), 'magenta')
            for label in self.atom_labels
        ]
    
    def get_atom_sizes(self) -> list:
        """  
        Get scatter plot sizes for all atoms, based on element type
        """
        return [
            self.ELEMENT_SIZES.get(self._remove_digits_from_label(label), 100)
            for label in self.atom_labels
        ]
    
    

    
    def get_mass_by_label(self, atom_label: str) -> float:
        """Get atomic mass by atom label"""
        index = np.where(self.atom_labels == atom_label)[0]
        if index.size == 0:
            raise ValueError(f"Atom label {atom_label} not found")
        return self.masses[index[0]]
    
    # Bond analysis
    def compute_bonds(self) -> None:
        """Compute covalent and hydrogen bonds from a single numba pair scan."""
        elements = np.array([self._remove_digits_from_label(label) for label in self.atom_labels], dtype=object)
        is_h = (elements == "H").astype(np.bool_)
        is_o = (elements == "O").astype(np.bool_)

        cov_i, cov_j, cov_s, hb_i, hb_j, hb_s = BondClassifier._pairwise_bond_scan(
            self.coordinates,
            self.covalent_radii,
            is_h,
            is_o,
            BondClassifier.COVALENT_THRESHOLD,
            BondClassifier.HYDROGEN_BOND_LOWER,
            BondClassifier.HYDROGEN_BOND_UPPER
        )
        self._covalent_bonds = [
            Bond(self.atom_labels[i], self.atom_labels[j], float(s))
            for i, j, s in zip(cov_i, cov_j, cov_s)
        
        ]
        self._hydrogen_bonds = [
            Bond(self.atom_labels[i], self.atom_labels[j], float(s))
            for i, j, s in zip(hb_i, hb_j, hb_s)
        ]
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
        angle_threshold: float = 150.0,
        max_distance: float = 3.5
        ) -> List[HBondConfiguration]:
        """
        Get hydrogen bond configurations meeting angle criterion
        """
        configs = self.find_hbond_configurations()
        return [c for c in configs if c.is_valid(angle_threshold, max_distance)]
    
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

    def compute_rmsd(self, other: 'Molecule', optimal_correspondence: bool = False) -> float:
        """ 
        Computes the RMSD between this molecule and another molecule.

        If optimal_correspondence is True, it will use the Hungarian algorithm to find the best atom mapping before calculating RMSD.
        """
        if len(self.coordinates) != len(other.coordinates):
            raise ValueError("Molecules must have the same number of atoms for RMSD calculation")
        if optimal_correspondence:
            coords2 = GeometryOps.find_optimal_correspondence(self.coordinates, other.coordinates)
        else:
            coords2 = other.coordinates
        diff = self.coordinates - coords2

        return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))



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

        
        # Count how many points are inside any atom's VDW sphere
        inside_count = self._count_points_inside_atoms(random_points, coords, radii)

        volume = (inside_count / n_points) * box_volume
        self.volume = volume
        return volume

    @staticmethod
    @njit(fastmath=True, parallel=True)
    def _count_points_inside_atoms(points: np.ndarray, coords: np.ndarray, radii: np.ndarray) -> int:
        """ 
        Numba-optimized function to count how many points are inside
        """
        inside = 0
        n_points = points.shape[0]
        n_atoms = coords.shape[0]
        for p in prange(n_points):
            px, py, pz = points[p, 0], points[p, 1], points[p, 2]
            hit = False
            for a in range(n_atoms):
                dx = px - coords[a, 0]
                dy = py - coords[a, 1]
                dz = pz - coords[a, 2]
                if dx*dx + dy*dy + dz*dz <= radii[a] * radii[a]:
                    hit = True
                    break
            if hit:
                inside += 1
        return inside

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
    
    def write_xyz(self, filename: str) -> None:
        """Write molecule to XYZ file"""
        with open(filename, "w") as f:
            f.write(self.to_xyz_string())

    def get_index_by_label(self, atom_label: str) -> int:
        """Get index of atom by label"""
        indices = np.where(self.atom_labels == atom_label)[0]
        if indices.size == 0:
            raise ValueError(f"Atom label {atom_label} not found")
        return int(indices[0]) 


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

    def get_symbols(self) -> List[str]:
        """Get list of element symbols for each atom"""
        return [self._remove_digits_from_label(label) for label in self.atom_labels]   



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

if __name__ == "__main__":
    h2o_2_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/h2o_2.xyz"
    with open(h2o_2_file, "r") as f:
        h2o_2_xyz = f.read()
    h2o_2 = Molecule.from_xyz(h2o_2_xyz, name="H2O Dimer")
    # compute and print bonds
    covalent_bonds, hydrogen_bonds = h2o_2.get_bonds()
    print("Covalent Bonds:")
    for bond in covalent_bonds:
        print(f"  {bond.atom1} - {bond.atom2}")
    print("Hydrogen Bonds:")
    for bond in hydrogen_bonds:
        print(f"  {bond.atom1} - {bond.atom2}") 

    # Test the hydrogen bond configurations
    time_start = time.time()
    hbond_configs = h2o_2.find_hbond_configurations()
    time_end = time.time()
    print(f"Found {len(hbond_configs)} H-bond configurations in {time_end - time_start:.4f} seconds")
    for config in hbond_configs:
        print(f"Donor: {h2o_2.atom_labels[config.donor_idx]}, Hydrogen: {h2o_2.atom_labels[config.hydrogen_idx]}, Acceptor: {h2o_2.atom_labels[config.acceptor_idx]}, Angle: {config.angle:.2f}, Distance: {config.donor_acceptor_distance:.2f}")

    # Test the compute bonds function for a bunch of alkanes
    ethane_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/alkanes_speed_test/ethan.xyz"
    propane_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/alkanes_speed_test/propan.xyz" 
    butan_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/alkanes_speed_test/butan.xyz"
    pentane_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/alkanes_speed_test/pentane.xyz"
    hexan_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/alkanes_speed_test/hexan.xyz"
    heptane_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/alkanes_speed_test/heptane.xyz"
    with open(ethane_file, "r") as f:
        ethane_xyz = f.read()
    with open(propane_file, "r") as f:
        propane_xyz = f.read()
    with open(butan_file, "r") as f:
        butan_xyz = f.read()
    with open(pentane_file, "r") as f:
        pentane_xyz = f.read()
    with open(hexan_file, "r") as f:
        hexan_xyz = f.read()
    with open(heptane_file, "r") as f:
        heptane_xyz = f.read()
    ethane = Molecule.from_xyz(ethane_xyz, name="Ethane")
    propane = Molecule.from_xyz(propane_xyz, name="Propane")
    butan = Molecule.from_xyz(butan_xyz, name="Butane")
    pentane = Molecule.from_xyz(pentane_xyz, name="Pentane")
    hexan = Molecule.from_xyz(hexan_xyz, name="Hexane")
    heptane = Molecule.from_xyz(heptane_xyz, name="Heptane")
    times = []
    # Compute bonds for each molecule and time it
    for mol in [ethane, propane, butan, pentane, hexan, heptane]:
        start_time = time.time()
        mol.compute_bonds()
        end_time = time.time()
        elapsed = end_time - start_time
        times.append(elapsed)
        print(f"{mol.name}: {elapsed:.4f} seconds")

    # Plot times vs number of atoms
    num_atoms = [len(mol.coordinates) for mol in [ethane, propane, butan, pentane, hexan, heptane]]
    plt.figure(figsize=(8, 5))
    plt.plot(num_atoms, times, marker='o')
    plt.xlabel("Number of Atoms")
    plt.ylabel("Bond Computation Time (s)")
    plt.title("Bond Computation Time vs Number of Atoms")
    plt.grid(True)
    plt.savefig("figures/bond_computation_time.png")