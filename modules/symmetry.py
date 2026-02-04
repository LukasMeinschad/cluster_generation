import numpy as np
from molecule_class import Molecule
from typing import List, Dict, Any, Tuple
from itertools import combinations
import time
import random 
import matplotlib.pyplot as plt

"""
Symmetry detection module based on:
https://doi.org/10.1016/j.cpc.2017.01.019
'SYVA: A program to analyze symmetry of molecules based on vector algebra'
Authors: László Gyevi-Nagy, Gyula Tasi
"""


class SymmetryAnalyzer:
    """
    A class to analyze the symmetry properties of a molecule.
    
    Attributes:
        molecule: Molecule object to analyze
        point_group: Detected point group (str)
        linear: Whether molecule is linear (bool)
        planar: Whether molecule is planar (bool)
        has_inversion: Whether molecule has inversion center (bool)
    """
    
    def __init__(self, molecule: Molecule):
        self.molecule = molecule
        self.point_group = None
        self.linear = False
        self.planar = False
        self.has_inversion = False

        # Store the symmetry info
        self.symmetry_info = {}


        self._coords = None
        self._com = None
        self._rel_coords = None
        self._atom_labels = None

    # ==================== Helper Methods ====================
    
    def _initialize_coordinates(self):
        """Cache coordinates and center of mass for performance."""
        if self._coords is None:
            self._coords = self.molecule.get_coordinates_array()
            self._com = self.molecule.center_of_mass()
            self._rel_coords = self._coords - self._com
            self._atom_labels = self.remove_digits(self.molecule.get_atom_labels())

    @staticmethod
    def remove_digits(symbols: List[str]) -> List[str]:
        """
        Removes digits from atom symbols.
        Example: ['C1', 'H2', 'O3'] -> ['C', 'H', 'O']
        """
        return [''.join(filter(str.isalpha, symbol)) for symbol in symbols]

    def check_structure_equivalence(
        self,
        coords1: np.ndarray,
        coords2: np.ndarray,
        labels1: List[str],
        labels2: List[str],
        tolerance: float = 1e-5
    ) -> bool:
        """
        Checks if two molecular structures are equivalent.
        
        Two structures are equivalent if each atom in structure 1 has a
        corresponding atom of the same element at the same position in structure 2.
        
        Args:
            coords1: Coordinates of first structure
            coords2: Coordinates of second structure
            labels1: Atom labels of first structure
            labels2: Atom labels of second structure
            tolerance: Distance tolerance for matching
            
        Returns:
            True if structures are equivalent
        """
        n_atoms = len(labels1)
        if n_atoms != len(labels2):
            return False

        matched = [False] * n_atoms

        for i in range(n_atoms):
            found_match = False
            for j in range(n_atoms):
                if matched[j]:
                    continue

                # Check same element type
                if labels1[i] != labels2[j]:
                    continue

                # Check position match
                distance = np.linalg.norm(coords1[i] - coords2[j])
                if distance < tolerance:
                    matched[j] = True
                    found_match = True
                    break

            if not found_match:
                return False

        return all(matched)

    # ==================== Basic Symmetry Checks ====================

    def check_linearity(self, tolerance: float = 1e-5) -> bool:
        """
        Checks if the molecule is linear.
        
        A molecule is linear when all position vectors (relative to COM)
        are parallel or anti-parallel:
        |(P_i - P_com) · (P_j - P_com)| / (||P_i|| * ||P_j||) = 1
        
        Args:
            tolerance: Numerical tolerance
            
        Returns:
            True if molecule is linear
        """
        self._initialize_coordinates()
        n_atoms = len(self._rel_coords)

        for i, j in combinations(range(n_atoms), 2):
            vec_i = self._rel_coords[i]
            vec_j = self._rel_coords[j]
            
            norm_i = np.linalg.norm(vec_i)
            norm_j = np.linalg.norm(vec_j)

            if norm_i < tolerance or norm_j < tolerance:
                continue

            cos_angle = np.abs(np.dot(vec_i, vec_j)) / (norm_i * norm_j)
            
            if not np.isclose(cos_angle, 1.0, atol=tolerance):
                self.linear = False
                return False

        self.linear = True
        return True

    def check_planarity(self, tolerance: float = 1e-5) -> bool:
        """
        Checks if all atoms lie in the same plane.
        
        Method:
        1. Find normal vector n from cross product of two non-parallel vectors
        2. Check if all position vectors are orthogonal to n
        
        Args:
            tolerance: Numerical tolerance
            
        Returns:
            True if molecule is planar
        """
        self._initialize_coordinates()
        n_atoms = len(self._rel_coords)

        if n_atoms < 3:
            self.planar = True
            return True

        # Find two non-parallel vectors
        normal_vector = None
        for i in range(n_atoms):
            if np.linalg.norm(self._rel_coords[i]) < tolerance:
                continue

            for j in range(i + 1, n_atoms):
                if np.linalg.norm(self._rel_coords[j]) < tolerance:
                    continue

                cross_prod = np.cross(self._rel_coords[i], self._rel_coords[j])
                cross_norm = np.linalg.norm(cross_prod)

                if cross_norm > tolerance:
                    normal_vector = cross_prod / cross_norm
                    break

            if normal_vector is not None:
                break

        if normal_vector is None:
            self.planar = True
            return True

        # Check orthogonality to normal
        for i in range(n_atoms):
            if np.linalg.norm(self._rel_coords[i]) < tolerance:
                continue

            dot_product = np.abs(np.dot(self._rel_coords[i], normal_vector))
            
            if dot_product > tolerance:
                self.planar = False
                return False

        self.planar = True
        return True

    def check_inversion_center(self, tolerance: float = 1e-5) -> bool:
        """
        Checks if molecule has an inversion center at the COM.
        
        An inversion center exists if for every atom at P_i, there exists
        an atom of the same element at -P_i (relative to COM).
        
        Args:
            tolerance: Distance tolerance
            
        Returns:
            True if inversion center exists
        """
        self._initialize_coordinates()
        n_atoms = len(self._atom_labels)
        matched = [False] * n_atoms

        for i in range(n_atoms):
            if matched[i]:
                continue

            # Atom at center of mass
            if np.linalg.norm(self._rel_coords[i]) < tolerance:
                matched[i] = True
                continue

            # Look for inverted partner
            inv_pos = -self._rel_coords[i]
            found_partner = False

            for j in range(n_atoms):
                if i == j or matched[j]:
                    continue

                if self._atom_labels[i] != self._atom_labels[j]:
                    continue

                distance = np.linalg.norm(self._rel_coords[j] - inv_pos)
                if distance < tolerance:
                    matched[i] = True
                    matched[j] = True
                    found_partner = True
                    break

            if not found_partner:
                self.has_inversion = False
                return False

        self.has_inversion = True
        return True

    def get_molecular_axis(self, tolerance: float = 1e-5) -> np.ndarray:
        """
        Returns the normalized direction vector of the molecular axis.
        Only valid for linear molecules.
        
        Returns:
            Normalized direction vector
            
        Raises:
            ValueError: If molecule is not linear
        """
        if not self.linear:
            raise ValueError("Molecule is not linear.")

        self._initialize_coordinates()

        for vec in self._rel_coords:
            norm = np.linalg.norm(vec)
            if norm > tolerance:
                return vec / norm

        raise ValueError("Cannot determine molecular axis.")

    def get_plane_normal(self, tolerance: float = 1e-5) -> np.ndarray:
        """
        Returns the normalized normal vector of the molecular plane.
        Only valid for planar molecules.
        
        Returns:
            Normalized normal vector
            
        Raises:
            ValueError: If molecule is not planar or normal cannot be determined
        """
        self._initialize_coordinates()
        n_atoms = len(self._rel_coords)

        for i in range(n_atoms):
            if np.linalg.norm(self._rel_coords[i]) < tolerance:
                continue

            for j in range(i + 1, n_atoms):
                if np.linalg.norm(self._rel_coords[j]) < tolerance:
                    continue

                cross_prod = np.cross(self._rel_coords[i], self._rel_coords[j])
                cross_norm = np.linalg.norm(cross_prod)

                if cross_norm > tolerance:
                    return cross_prod / cross_norm

        raise ValueError("Cannot determine plane normal.")

    # ==================== Transformation Operations ====================

    def reflect_point(
        self,
        point: np.ndarray,
        normal: np.ndarray,
        plane_point: np.ndarray = None
    ) -> np.ndarray:
        """
        Reflects a point across a plane.
        
        P' = P - 2 * proj_n(P - P0)
        where proj_n(v) = (v · n / ||n||²) * n
        
        Args:
            point: Point to reflect
            normal: Normal vector of plane
            plane_point: Point on plane (default: origin)
            
        Returns:
            Reflected point
        """
        if plane_point is None:
            plane_point = np.zeros(3)

        n = normal / np.linalg.norm(normal)
        v = point - plane_point
        proj = np.dot(v, n) * n

        return point - 2 * proj

    def rodrigues_rotation(
        self,
        point: np.ndarray,
        axis: np.ndarray,
        angle: float,
        axis_point: np.ndarray = None
    ) -> np.ndarray:
        """
        Rotates a point around an axis using Rodrigues' formula.
        
        P' = P + sin(φ) * (d × P) + (1 - cos(φ)) * (d × (d × P))
        
        Args:
            point: Point to rotate
            axis: Rotation axis (will be normalized)
            angle: Rotation angle in radians
            axis_point: Point on axis (default: origin)
            
        Returns:
            Rotated point
        """
        if axis_point is None:
            axis_point = np.zeros(3)

        d = axis / np.linalg.norm(axis)
        p = point - axis_point

        sin_phi = np.sin(angle)
        cos_phi = np.cos(angle)

        cross1 = np.cross(d, p)
        cross2 = np.cross(d, cross1)

        p_rotated = p + sin_phi * cross1 + (1 - cos_phi) * cross2

        return p_rotated + axis_point

    def improper_rotation(
        self,
        point: np.ndarray,
        axis: np.ndarray,
        order: int,
        axis_point: np.ndarray = None
    ) -> np.ndarray:
        """
        Performs improper rotation: rotation followed by reflection.
        
        S_n = C_n + σ_h
        
        Args:
            point: Point to transform
            axis: Rotation axis
            order: Order of rotation
            axis_point: Point on axis (default: origin)
            
        Returns:
            Transformed point
        """
        if axis_point is None:
            axis_point = np.zeros(3)

        angle = 2 * np.pi / order
        rotated = self.rodrigues_rotation(point, axis, angle, axis_point)
        reflected = self.reflect_point(rotated, axis, axis_point)

        return reflected

    # ==================== Symmetry Element Detection ====================

    def find_symmetry_planes(self, tolerance: float = 1e-5) -> List[Dict[str, Any]]:
        """
        Finds all symmetry planes (mirror planes).
        
        Returns:
            List of dicts with keys: 'normal', 'point', 'type', 'defining_atoms'
        """
        
        self._initialize_coordinates()
        n_atoms = len(self._atom_labels)
        symmetry_planes = []

        # 1. Check molecular plane for planar molecules
        if self.planar:
            try:
                plane_normal = self.get_plane_normal(tolerance)
                reflected_coords = np.array([
                    self.reflect_point(self._rel_coords[i], plane_normal)
                    for i in range(n_atoms)
                ])

                if self.check_structure_equivalence(
                    self._rel_coords, reflected_coords,
                    self._atom_labels, self._atom_labels, tolerance
                ):
                    symmetry_planes.append({
                        'normal': plane_normal,
                        'point': self._com,
                        'type': 'molecular_plane',
                        'defining_atoms': 'all'
                    })
            except ValueError:
                pass

        # 2. Test candidate planes from atom pairs
        tested_normals = []
        
        for i, j in combinations(range(n_atoms), 2):
            normal_candidate = self._rel_coords[i] - self._rel_coords[j]
            norm = np.linalg.norm(normal_candidate)

            if norm < tolerance:
                continue

            normal_candidate /= norm

            # Check for duplicate directions
            is_duplicate = any(
                np.abs(np.dot(normal_candidate, tested)) > (1 - tolerance)
                for tested in tested_normals
            )

            if is_duplicate:
                continue

            tested_normals.append(normal_candidate)

            # Test reflection
            reflected_coords = np.array([
                self.reflect_point(self._rel_coords[k], normal_candidate)
                for k in range(n_atoms)
            ])

            if self.check_structure_equivalence(
                self._rel_coords, reflected_coords,
                self._atom_labels, self._atom_labels, tolerance
            ):
                symmetry_planes.append({
                    'normal': normal_candidate,
                    'point': self._com,
                    'type': 'reflection_plane',
                    'defining_atoms': (i, j)
                })

        return symmetry_planes

    def check_rotation_axis(
        self,
        axis: np.ndarray,
        order: int,
        tolerance: float = 1e-5
    ) -> bool:
        """
        Checks if an axis is a C_n rotation axis.
        
        Args:
            axis: Direction vector
            order: Rotation order
            tolerance: Distance tolerance
            
        Returns:
            True if axis is a C_n symmetry axis
        """
        self._initialize_coordinates()
        n_atoms = len(self._atom_labels)
        
        angle = 2 * np.pi / order

        rotated_coords = np.array([
            self.rodrigues_rotation(self._coords[i], axis, angle, self._com)
            for i in range(n_atoms)
        ])

        return self.check_structure_equivalence(
            self._coords, rotated_coords,
            self._atom_labels, self._atom_labels, tolerance
        )

    def find_c2_axes(self, tolerance: float = 1e-5) -> List[Dict[str, Any]]:
        """
        Finds all C2 rotation axes.
        
        Direction vector: d = (P_i + P_j) / ||P_i + P_j||
        
        Returns:
            List of C2 axes
        """
        
        self._initialize_coordinates()
        n_atoms = len(self._atom_labels)
        
        c2_axes = []
        tested_axes = []

        for i, j in combinations(range(n_atoms), 2):
            d = self._coords[i] + self._coords[j]
            norm = np.linalg.norm(d)

            if norm < tolerance:
                continue

            d /= norm

            # Check duplicate
            is_duplicate = any(
                np.abs(np.dot(d, tested)) > (1 - tolerance)
                for tested in tested_axes
            )

            if is_duplicate:
                continue

            tested_axes.append(d)

            if self.check_rotation_axis(d, 2, tolerance):
                c2_axes.append({
                    'axis': d,
                    'order': 2,
                    'point': self._com,
                    'type': 'C2',
                    'defining_atoms': (i, j)
                })

        return c2_axes

    def find_cn_axes(
        self,
        max_order: int = 8,
        tolerance: float = 1e-5
    ) -> List[Dict[str, Any]]:
        """
        Finds all Cn rotation axes (n > 2).
        
        Direction vector: d = (P_i - P_j) × (P_j - P_k) / ||...||
        
        Args:
            max_order: Maximum order to test
            tolerance: Distance tolerance
            
        Returns:
            List of Cn axes
        """
        
        self._initialize_coordinates()
        n_atoms = len(self._atom_labels)
        
        cn_axes = []
        tested_axes = []

        # Test linear molecule axis
        if self.linear:
            try:
                axis = self.get_molecular_axis(tolerance)
                for order in range(3, max_order + 1):
                    if self.check_rotation_axis(axis, order, tolerance):
                        cn_axes.append({
                            'axis': axis,
                            'order': order,
                            'point': self._com,
                            'type': f'C{order}',
                            'defining_atoms': 'molecular_axis'
                        })
                return cn_axes
            except ValueError:
                pass

        # Test planar molecule perpendicular axis
        if self.planar:
            try:
                axis = self.get_plane_normal(tolerance)
                for order in range(3, max_order + 1):
                    if self.check_rotation_axis(axis, order, tolerance):
                        cn_axes.append({
                            'axis': axis,
                            'order': order,
                            'point': self._com,
                            'type': f'C{order}',
                            'defining_atoms': 'plane_normal'
                        })
            except ValueError:
                pass

        # Test axes from atom triplets
        for i, j, k in combinations(range(n_atoms), 3):
            v1 = self._coords[i] - self._coords[j]
            v2 = self._coords[j] - self._coords[k]

            cross_prod = np.cross(v1, v2)
            norm = np.linalg.norm(cross_prod)

            if norm < tolerance:
                continue

            d = cross_prod / norm

            # Check duplicate
            is_duplicate = any(
                np.abs(np.dot(d, tested)) > (1 - tolerance)
                for tested in tested_axes
            )

            if is_duplicate:
                continue

            tested_axes.append(d)

            # Test orders from high to low
            for order in range(max_order, 2, -1):
                if self.check_rotation_axis(d, order, tolerance):
                    cn_axes.append({
                        'axis': d,
                        'order': order,
                        'point': self._com,
                        'type': f'C{order}',
                        'defining_atoms': (i, j, k)
                    })
                    break  # Only record highest order

        return cn_axes

    def find_all_rotation_axes(
        self,
        max_order: int = 8,
        tolerance: float = 1e-5
    ) -> Dict[str, List[Dict[str, Any]]]:
        """
        Finds all proper rotation axes.
        
        Returns:
            Dict with 'C2' and 'Cn' keys containing lists of axes
        """

        c2_axes = self.find_c2_axes(tolerance)
        cn_axes = self.find_cn_axes(max_order, tolerance)

        return {'C2': c2_axes, 'Cn': cn_axes}

    def check_improper_rotation_axis(
        self,
        axis: np.ndarray,
        order: int,
        tolerance: float = 1e-5
    ) -> bool:
        """
        Checks if an axis is an S_n improper rotation axis.
        
        Args:
            axis: Direction vector
            order: Rotation order
            tolerance: Distance tolerance
            
        Returns:
            True if axis is an S_n axis
        """
        self._initialize_coordinates()
        n_atoms = len(self._atom_labels)

        transformed_coords = np.array([
            self.improper_rotation(self._coords[i], axis, order, self._com)
            for i in range(n_atoms)
        ])

        return self.check_structure_equivalence(
            self._coords, transformed_coords,
            self._atom_labels, self._atom_labels, tolerance
        )

    def find_improper_rotation_axes(
        self,
        max_order: int = 8,
        tolerance: float = 1e-5
    ) -> List[Dict[str, Any]]:
        """
        Finds all improper rotation axes (Sn).
        
        Returns:
            List of Sn axes
        """
        
        self._initialize_coordinates()
        n_atoms = len(self._coords)
        
        sn_axes = []
        tested_axes = []

        # Test molecular axis for linear molecules
        if self.linear:
            try:
                axis = self.get_molecular_axis(tolerance)
                for order in range(2, max_order + 1):
                    if self.check_improper_rotation_axis(axis, order, tolerance):
                        sn_axes.append({
                            'axis': axis,
                            'order': order,
                            'point': self._com,
                            'type': f'S{order}',
                            'defining_atoms': 'molecular_axis'
                        })
                return sn_axes
            except ValueError:
                pass

        # Test plane normal for planar molecules
        if self.planar:
            try:
                axis = self.get_plane_normal(tolerance)
                for order in range(2, max_order + 1):
                    if self.check_improper_rotation_axis(axis, order, tolerance):
                        sn_axes.append({
                            'axis': axis,
                            'order': order,
                            'point': self._com,
                            'type': f'S{order}',
                            'defining_atoms': 'plane_normal'
                        })
            except ValueError:
                pass

        # Test axes from atom pairs
        for i, j in combinations(range(n_atoms), 2):
            d = self._coords[i] + self._coords[j]
            norm = np.linalg.norm(d)

            if norm < tolerance:
                continue

            d /= norm

            is_duplicate = any(
                np.abs(np.dot(d, tested)) > (1 - tolerance)
                for tested in tested_axes
            )

            if is_duplicate:
                continue

            tested_axes.append(d)

            for order in range(max_order, 1, -1):
                if self.check_improper_rotation_axis(d, order, tolerance):
                    sn_axes.append({
                        'axis': d,
                        'order': order,
                        'point': self._com,
                        'type': f'S{order}',
                        'defining_atoms': (i, j)
                    })
                    break

        # Test axes from atom triplets
        for i, j, k in combinations(range(n_atoms), 3):
            v1 = self._coords[i] - self._coords[j]
            v2 = self._coords[j] - self._coords[k]

            cross_prod = np.cross(v1, v2)
            norm = np.linalg.norm(cross_prod)

            if norm < tolerance:
                continue

            d = cross_prod / norm

            is_duplicate = any(
                np.abs(np.dot(d, tested)) > (1 - tolerance)
                for tested in tested_axes
            )

            if is_duplicate:
                continue

            tested_axes.append(d)

            for order in range(max_order, 1, -1):
                if self.check_improper_rotation_axis(d, order, tolerance):
                    sn_axes.append({
                        'axis': d,
                        'order': order,
                        'point': self._com,
                        'type': f'S{order}',
                        'defining_atoms': (i, j, k)
                    })
                    print(f"S{order} Axis: Atoms {i}-{j}-{k}")
                    break

        return sn_axes

    # ==================== Point Group Determination ====================

    def get_symmetry_elements(self, tolerance: float = 1e-5, max_order: int = 8) -> Dict[str, Any]:
        """   
        Function that determines all the possible symmetry elements of the molecule.
        """
        self._initialize_coordinates()
        self.check_linearity(tolerance)
        self.check_planarity(tolerance)
        self.check_inversion_center(tolerance)

        symmetry_elements = {
            'linear': self.linear,
            'planar': self.planar,
            'inversion_center': self.has_inversion,
            'symmetry_planes': self.find_symmetry_planes(tolerance),
            'rotation_axes': self.find_all_rotation_axes(max_order, tolerance),
            'improper_rotation_axes': self.find_improper_rotation_axes(max_order, tolerance)
        }

        return symmetry_elements


    def test_analysis_speed(self, tolerance: float = 1e-5, max_order: int = 8) -> None:
        """
        Tests and prints the time taken for full symmetry analysis.
        
        Args:
            tolerance: Numerical tolerance
            max_order: Maximum order for rotation axes
        """
        # Initialize Random Molecules with N_atoms = 3, 5, 10, 20, 30, 50, 100
        atom_counts = list(range(1, 25))
        times = []
        for n_atoms in atom_counts:
            # Generate random coordinates
            coords = np.random.rand(n_atoms, 3) * 10.0  # Random positions in 10x10x10 cube
            atom_labels = random.choices(['H', 'C', 'N', 'O', 'F'], k=n_atoms)

            # Create Molecule object
            test_molecule = Molecule.from_labels_and_coords(atom_labels, coords) 

            # Create SymmetryAnalyzer
            analyzer = SymmetryAnalyzer(molecule=test_molecule)

            # Time the analysis
            start_time = time.time()
            analyzer.analyze_full_symmetry(tolerance=tolerance, max_order=max_order)
            end_time = time.time()
            elapsed_time = end_time - start_time
            times.append(elapsed_time)
        
        # Apply a polyfit to the data to determine the O(n^k) scaling
        log_atoms = np.log(atom_counts)
        log_times = np.log(times)
        coeffs = np.polyfit(log_atoms, log_times, 1)
        scaling_exponent = coeffs[0]
        print(f"Estimated time complexity: O(n^{scaling_exponent:.2f})")

        # Plot the results
        plt.figure(figsize=(8, 6))  
        plt.plot(atom_counts, times, marker='o')
        plt.xlabel('Number of Atoms')
        plt.ylabel('Analysis Time (s)')
        plt.title('Symmetry Analysis Time vs Number of Atoms, Estimated O(n^{:.2f})'.format(scaling_exponent))
        plt.grid(True) 
        plt.savefig("figures/symmetry_analysis_time.png")
        plt.close()
