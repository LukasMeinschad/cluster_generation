import numpy as np
from molecule_class import Molecule
from modules.logger import Logger
from typing import List, Dict, Any
from itertools import combinations



""" 
The symmetry detection module presented here is based on the paper 
https://doi.org/10.1016/j.cpc.2017.01.019 - SYVA: A program to analyze symmetry of molecules based on vector algebra
Authors: László Gyevi-Nagy, Gyula Tasi 
"""

class SymmetryAnalyzer:
    def __init__(self, molecule: Molecule):
        """   
        A class to analyze the symmetry properties of a molecule
        """
        self.molecule = molecule 
        self.point_group = None 
        self.linear = False 

    @staticmethod
    def remove_digits(symbols: List[str]) -> List[str]:
        """   
        Removes digits from a list of atom symbols
        E.g. ['C1', 'H2', 'O'] -> ['C', 'H', 'O']
        """
        return [''.join(filter(str.isalpha, symbol)) for symbol in symbols]

    def check_linearity(self, tolerance: float = 1e-5):
        """   
        Checks if a given molecule is linear

        For this:

        + Calculate the position vectors of all atoms relative to the center of mass P_i - P_com
        + The molecule is linear when the difference vectors betwenn all vectors are parallel

        |(P_i - P_com) * (P_j - P_com) / ||P_i - P_com|| * ||P_j - P_com|| | = 1 for all i,j
        """
        coords = self.molecule.get_coordinates_array() # Get coordinates as numpy array
        com = self.molecule.center_of_mass()
        rel_coords = coords - com  # Position vectors relative to center of mass
        n_atoms = rel_coords.shape[0]
        for i,j in combinations(range(n_atoms), 2):
            vec_i = rel_coords[i] 
            vec_j = rel_coords[j]
            norm_i = np.linalg.norm(vec_i)
            norm_j = np.linalg.norm(vec_j)
            if norm_i < tolerance or norm_j < tolerance:
                continue  # Skip zero-length vectors
            cos_angle = np.dot(vec_i, vec_j) / (norm_i * norm_j)
            if not (np.isclose(cos_angle, 1.0, atol=tolerance) or np.isclose(cos_angle, -1.0, atol=tolerance)):
                self.linear = False
                print("Non-Linear Molecule Detected.")
                return self.linear
        self.linear = True
        print("Linear Molecule Detected.")
        return self.linear
    
    def check_inversion_center(self,tolerance: float = 1e-5) -> bool:
        """
        Checks if a given molecule has a inversion center at the com 

        We have an inversion center if every atom at position P_i, there exists an
        corresponding atom of the same element at position -P_i (relative to com)
        """
        coords = self.molecule.get_coordinates_array() 
        com = self.molecule.center_of_mass()
        rel_coords = coords - com

        # Get atom symbols and remove digits
        atom_labels = self.remove_digits(self.molecule.get_atom_labels())
        
        n_atoms = len(atom_labels)
        matched = [False] * n_atoms  # Track matched atoms

        for i in range(n_atoms):
            if matched[i]:
                continue 
            inv_pos = -rel_coords[i] 
            found_partner = False 

            # Look for matching atom at inverted position
            for j in range(n_atoms):
                if i == j or matched[j]:
                    continue 

                # Check if same element type 
                if atom_labels[i] != atom_labels[j]:
                    continue

                # Check if positions match within tolerance
                distance = np.linalg.norm(rel_coords[j] - inv_pos) 
                if distance < tolerance:
                    matched[i] = True
                    matched[j] = True
                    found_partner = True
                    break
            
            # If atom is at center of mass it matches itself
            if np.linalg.norm(rel_coords[i]) < tolerance:
                matched[i] = True
                found_partner = True
            if not found_partner:
                print("No Inversion Center Detected.")
                return False
            
        print("Inversion Center Detected.")
        return all(matched)
    
    def check_planarity(self,tolerance: float = 1e-5) -> bool:
        """  
        Checks if all atoms of the molecule lie in the same plane 

        Method:
        1. Calculate position vectors relative to the center of mass (P_i - P_com)
        2. Generate normal vector n as cross product of two differnce vectors 
              n = (P_2 - P_1) x (P_3 - P_1)
        3. Check if all difference vectors are orthogonal to n:
              |(P_i - P_1) . n| / ||P_i - P_1|| * ||n|| = 0 for all i
        """
        coords = self.molecule.get_coordinates_array() 
        com = self.molecule.center_of_mass() 
        rel_coords = coords - com 
        n_atoms = rel_coords.shape[0] 

        if n_atoms < 3:
            print("Molecule has less than 3 atoms, considered planar.")
            return True  # Less than 3 atoms are always planar
        
        # Find two non-parallel vectors to define the plane
        normal_vector = None
        for i in range(n_atoms):
            if np.linalg.norm(rel_coords[i]) < tolerance:
                continue # Skip zero-length vectors
            for j in range(i+1, n_atoms):
                if np.linalg.norm(rel_coords[j]) < tolerance:
                    continue 

                # Check if vectors are not parallel
                vec_i = rel_coords[i]
                vec_j = rel_coords[j]
                cross_prod = np.cross(vec_i, vec_j)
                cross_norm = np.linalg.norm(cross_prod)

                if cross_norm >= tolerance:
                    # Found two non-parallel vectors, use their cross product as normal vector 
                    normal_vector = cross_prod / cross_norm 
                    break

            if normal_vector is not None:
                break

        if normal_vector is None:
            print("All atoms are collinear, considered planar.")
            return True  # All atoms are collinear, considered planar
        
        # Check if all atoms are orthogonal to the normal vector
        for i in range(n_atoms):
            vec_i = rel_coords[i] 

            # Skip zero-length vectors (atoms at center of mass)
            if np.linalg.norm(vec_i) < tolerance: 
                continue    

            # Calculate the dot product
            dot_product = np.dot(vec_i, normal_vector)

            if not np.isclose(dot_product, 0.0, atol=tolerance):
                print("Non-Planar Molecule Detected.")
                return False
        print("Planar Molecule Detected.")
        return True

    def get_plane_normal(self, tolerance: float = 1e-5) -> np.ndarray:
        """  
        Calculates the normal vector of the molecular plane for a planar molecule
        """
        coords = self.molecule.get_coordinates_array()
        com = self.molecule.center_of_mass()
        rel_coords = coords - com
        n_atoms = len(rel_coords)

        # Find two non-parallel vectors to define the plane
        for i in range(n_atoms):
            if np.linalg.norm(rel_coords[i]) < tolerance:
                continue 

            for j in range(i+1, n_atoms):
                if np.linalg.norm(rel_coords[j]) < tolerance:
                    continue 

                # Check if vectors are not parallel
                vec_i = rel_coords[i]
                vec_j = rel_coords[j]
                cross_prod = np.cross(vec_i, vec_j)
                cross_norm = np.linalg.norm(cross_prod)

                if cross_norm >= tolerance:
                    # Found two non-parallel vectors, use their cross product as normal vector 
                    normal_vector = cross_prod / cross_norm 
                    return normal_vector
        raise ValueError("Cannot determine plane normal: all atoms may be collinear.")
    
    def reflect_point(self, point: np.ndarray, normal: np.ndarray, plane_point: np.ndarray = None) -> np.ndarray:
        """  
        Reflects a point across a plane defined by a normal vector and a point on a plane 

        P' = P - 2 * proj_n(P - P_0)
        where proj_n(v) = (v cdot n) /(n cdot n) * n
        """
        if plane_point is None:
            plane_point = np.zeros(3)

        # Normalize the normal vector
        n = normal / np.linalg.norm(normal)
        # vector from plane point to the point to be reflected
        v = point - plane_point

        # Projection of v onto n 
        proj = np.dot(v, n) * n

        # Reflected point
        reflected_point = point - 2 * proj

        return reflected_point

    def check_structure_equivalence(self, 
                                    coords1: np.ndarray,
                                    coords2: np.ndarray,
                                    labels1: List[str],
                                    labels2: List[str],
                                    tolerance: float = 1e-5) -> bool:
        """   
        Checks if two molecular structures are equivalent (same atoms at same positions)
        
        Args:
            coords1: Coordinates of first structure
            coords2: Coordinates of second structure
            labels1: Atom labels of first structure
            labels2: Atom labels of second structure
            tolerance: Distance tolerance for matching atoms
        
        Returns:
            True if structures are equivalent, False otherwise
        """
        n_atoms = len(labels1)
        if n_atoms != len(labels2):
            return False  # Different number of atoms
        
        matched = [False] * n_atoms  # Track matched atoms

        for i in range(n_atoms):
            found_match = False
            for j in range(n_atoms):
                if matched[j]:
                    continue 

                # Check if same element type
                if labels1[i] != labels2[j]:
                    continue

                # check if positions match 
                distance = np.linalg.norm(coords1[i] - coords2[j])
                if distance < tolerance:
                    matched[j] = True
                    found_match = True
                    break
            if not found_match:
                return False  # No matching atom found
            
        return all(matched) 
    
    def find_symmetry_planes(self, tolerance: float = 1e-5) -> List[Dict[str, Any]]:
        """
        Finds all symmetry planes of the molecule

        A symmetry plane either:
        1. Contains all atoms (molecular plane for planar molecules)
        2. Contains some atoms and reflects pairs of atoms into each other

        Method:
        - For planaer molecules: Check if the molecule plane is a symmetry plane
        - For all molecules: Generate candidate planes from pairs of atoms with same element using normal vector n = P_i - P_j
        - Test each plane by reflecting all atoms and checking structure equivalence
        """
        coords = self.molecule.get_coordinates_array()  
        com = self.molecule.center_of_mass() 
        rel_coords = coords - com 
        atom_labels = self.remove_digits(self.molecule.get_atom_labels())
        n_atoms = len(atom_labels)

        symmetry_planes = []
        
        #1. Check if molecular plane is symmetry plane for planar molecules
        if self.check_planarity(tolerance=tolerance):
            try:
                plane_normal = self.get_plane_normal(tolerance=tolerance)
                # Reflect all atoms across the molecular plane
                reflected_coords = np.array([
                    self.reflect_point(rel_coords[i], plane_normal,)
                    for i in range(n_atoms)
                ])
                # Check structure equivalence
                if self.check_structure_equivalence(rel_coords, reflected_coords, atom_labels, atom_labels, tolerance=tolerance):
                    symmetry_planes.append({
                        "normal": plane_normal,
                        "point": com,
                        "type": "molecular_plane"
                    })
                    print(f"Symmetry Plane Detected: Molecular Plane")
            except ValueError:
                pass  # Cannot determine plane normal

        #2. Generate candiate planes from pairs of atoms
        tested_normals = []
        for i,j in combinations(range(n_atoms),2):
            #Generate normal vector as difference between position vectors
            normal_candidate = rel_coords[i] - rel_coords[j]
            norm = np.linalg.norm(normal_candidate)
            if norm < tolerance:
                continue  # Skip zero-length normals

            normal_candidate /= norm  # Normalize
            # Check if this direction of opposite direction has already been tested 
            is_duplicate = False 
            for tested_normal in tested_normals: 
                dot_product = abs(np.dot(normal_candidate, tested_normal))
                if np.isclose(dot_product, 1.0, atol=tolerance):
                    is_duplicate = True
                    break
            
            if is_duplicate:
                continue

            tested_normals.append(normal_candidate)

            # Test reflection through plane
            reflected_coords = np.array([
                self.reflect_point(rel_coords[k], normal_candidate) 
                for k in range(n_atoms)
            ])

            # Check if reflected structure matches original
            if self.check_structure_equivalence(rel_coords, reflected_coords, atom_labels, atom_labels, tolerance=tolerance):
                symmetry_planes.append({
                    "normal": normal_candidate,
                    "point": com, 
                    "type": "reflection_plane",
                    "defining_atoms": (i, j)
                })
                print(f"Symmetry Plane Detected: Atoms {i} and {j}")

        print(f"Total Symmetry Planes Found: {len(symmetry_planes)}")
        return symmetry_planes
    
    def rodrigues_rotation(self, point: np.ndarray, axis: np.ndarray, angle: float, axis_point: np.ndarray = None) -> np.ndarray:
        """    
        Rotates a point around an axis using Rodrigues' rotation formula 

        P' = P + sin(theta) * d x P + (1- cos(theda)) * d x (d x P)
        where d is the unit vector along the rotation axis

        Args:
            point: Point to be rotated
            axis: Rotation axis (unit vector)
            angle: Rotation angle in radians
            axis_point: A point on the rotation axis (default is origin)

        Returns:
            Rotated point as numpy array
        """
        if axis_point is None:
            # Default to origin
            axis_point = np.zeros(3)
        
        # Normalize axis direction 
        d = axis / np.linalg.norm(axis) 

        # Shit to origin if needed 
        p = point - axis_point 

        # Apply Rodgrigues Formula 
        sin_phi = np.sin(angle)
        cos_phi = np.cos(angle) 
        cross1 = np.cross(d,p)
        cross2 = np.cross(d, cross1) 

        p_rotated = p + sin_phi * cross1 + (1- cos_phi) * cross2 

        # Shift back
        return p_rotated + axis_point

    def check_rotation_axis(self, axis: np.ndarray, order:int, 
                            axis_point: np.ndarray = None,
                            tolerance: float = 1e-5) -> bool:
        """  
        Checks if a given axis is a symmetry axis of given order 
        
        Args:
            axis: Direction vector of axis 
            order: Order of rotation (n for C_n axis)
            axis_point: Point on the axis (default is center of mass)
            tolerance: Distance tolerance for positon matching
        """
        if axis_point is None:
            axis_point = self.molecule.center_of_mass()
        
        coords = self.molecule.get_coordinates_array()
        atom_labels = self.remove_digits(self.molecule.get_atom_labels())
        n_atoms = len(atom_labels) 

        # Rotation angle for C_n axis
        angle = 2 * np.pi / order

        # Rotate all atoms
        rotated_coords = np.array([ 
            self.rodrigues_rotation(coords[i], axis, angle, axis_point) 
            for i in range(n_atoms)
        ])

        # Check if rotated structure matches original
        return self.check_structure_equivalence(coords, rotated_coords, atom_labels, atom_labels, tolerance=tolerance)

    def get_molecular_axis(self,tolerance: float = 1e-5) -> np.ndarray:
        """  
        Returns the direction vector of the molecular axis for linear molecules 

        Returns:
            Direction vector as numpy array
        """    
        if not self.linear:
            raise ValueError("Molecule is not linear, molecular axis undefined.")
        coords = self.molecule.get_coordinates_array()
        com = self.molecule.center_of_mass()
        rel_coords = coords - com

        # Find the first non zero vector as reference 
        for vec in rel_coords:
            if np.linalg.norm(vec) >= tolerance:
                axis = vec / np.linalg.norm(vec)
                return axis
        raise ValueError("Cannot determine molecular axis: all atoms may be at center of mass.")



    def find_c2_axes(self, tolerance: float = 1e-5) -> List[Dict[str,Any]]:
        """
        Finds all possible C2 axes

        The direction vectors of C2 axes from atoms i and j:
        d = (P_i + P_j) / ||P_i + P_j||

        Returns:
            List of detected C2 axes with their direction and a point on the axis
        """
        coords = self.molecule.get_coordinates_array()
        atom_labels = self.remove_digits(self.molecule.get_atom_labels()) 
        n_atoms = len(atom_labels)

        c2_axes = []
        tested_axes = [] 
        for i,j in combinations(range(n_atoms),2):
            # Direction vector for C2 axis
            d = coords[i] + coords[j]
            norm = np.linalg.norm(d)
            if norm < tolerance:
                continue  # Skip zero-length vectors
            
            d /= norm  # Normalize

            # Check if this direction is already tested 
            is_duplicate = False
            for tested_axis in tested_axes:
                dot_product = abs(np.dot(d, tested_axis))
                if np.isclose(dot_product, 1.0, atol=tolerance):
                    is_duplicate = True
                    break

            if is_duplicate:
                continue

            tested_axes.append(d)

            # Check if this is a C2 axis
            if self.check_rotation_axis(d, order=2, tolerance=tolerance):
                c2_axes.append({
                    "axis": d,
                    "order": 2,
                    "point": self.molecule.center_of_mass(),
                    "defining_atoms": (i,j),
                    "type": "C2"
                })
                print(f"C2 Axis Detected: Atoms {i} and {j}")

        print(f"Total C2 Axes Found: {len(c2_axes)}")
        return c2_axes
    
    def find_cn_axes(self, max_order: int=8, tolerance: float = 1e-5) -> List[Dict[str,Any]]:
        """
        Finds all possible Cn axes up to a given maximum order

        Direction vectors for three atoms (i,j,k):
        d = (P_i - P_j) x (P_j - P_k) / ||(P_i - P_j) x (P_j - P_k)||   

        Consider cases where atoms i,j,k correspond to adjacent vertices of a regular polygon 

        Args:
            max_order: Maximum order of rotation axes to search for
            tolerance: Distance tolerance for position matching
        Returns:
            List of detected Cn axes with their direction, order, and a point on the axis
        """
        coords = self.molecule.get_coordinates_array()
        atom_labels = self.remove_digits(self.molecule.get_atom_labels()) 
        n_atoms = len(atom_labels)

        cn_axes = []
        tested_axes = []
        # Test linear molecule axis:
        if self.linear:
            try:
                axis = self.get_molecular_axis(tolerance=tolerance)
                # Linear molecules have C_infinity symmetry, we test higher orders 
                for order in range(3, max_order + 1):
                    if self.check_rotation_axis(axis, order, tolerance=tolerance):
                        cn_axes.append({
                            "axis": axis,
                            "order": order,
                            "point": self.molecule.center_of_mass(),
                            "type": f"C{order}",
                            "defining_atoms": "molecular_axis"
                        })
                        print(f"C{order} Axis Detected: Linear Molecule Axis")
                return cn_axes  # No need to search further
            except ValueError:
                pass  # Cannot determine molecular axis
        # Test planar molecule perpendicular axes 
        if self.check_planarity(tolerance=tolerance):
            try:
                axis = self.get_plane_normal(tolerance=tolerance)
                for order in range(3, max_order + 1):
                    if self.check_rotation_axis(axis, order, tolerance=tolerance):
                        cn_axes.append({
                            "axis": axis,
                            "order": order,
                            "point": self.molecule.center_of_mass(),
                            "type": f"C{order}",
                            "defining_atoms": "plane_normal"
                        })
                        print(f"C{order} Axis Detected: Planar Molecule Normal")
            except ValueError:
                pass  # Cannot determine plane normal

        # Test axes from triplets of atoms
        for i,j,k in combinations(range(n_atoms),3):
            # Calculate the direction vector 
            v1 = coords[i] - coords[j]
            v2 = coords[j] - coords[k]

            cross_prod = np.cross(v1, v2)
            norm = np.linalg.norm(cross_prod)
            if norm < tolerance:
                continue  # Skip zero-length normals
            d = cross_prod / norm  # Normalize

            # Check if duplicate
            is_duplicate = False
            for tested_axis in tested_axes:
                dot_product = abs(np.dot(d, tested_axis))
                if np.isclose(dot_product, 1.0, atol=tolerance):
                    is_duplicate = True
                    break

            if is_duplicate:
                continue

            tested_axes.append(d)

            # Test different orders
            for order in range(3, max_order + 1):
                if self.check_rotation_axis(d, order, tolerance=tolerance):
                    cn_axes.append({
                        "axis": d,
                        "order": order,
                        "point": self.molecule.center_of_mass(),
                        "type": f"C{order}",
                        "defining_atoms": (i,j,k)
                    })
                    print(f"C{order} Axis Detected: Atoms {i}, {j}, {k}")
        print(f"Total Cn Axes Found: {len(cn_axes)}")
        return cn_axes
    
    def find_all_rotation_axes(self, max_order: int=8, tolerance: float = 1e-5) -> Dict[str, List[Dict[str,Any]]]:
        """  
        Finds all rotation axes (Cn) up to a given maximum order

        Args:
            max_order: Maximum order of rotation axes to search for
            tolerance: Distance tolerance for position matching
        Returns:
            Dictionary with lists of detected rotation axes by order
        """
        print("\n ======== Searching for Rotation Axes ======== \n")

        # find C2 axes
        print("Searching for C2 Axes...")
        c2_axes = self.find_c2_axes(tolerance=tolerance)

        # Find higher order Cn axes 
        cn_axes = self.find_cn_axes(max_order=max_order, tolerance=tolerance)

        all_axes = {
            "C2": c2_axes,
            "Cn": cn_axes
        }
        return all_axes
    
    # TODO Implement Improper Rotation Axes shit and complete this lines of pure evil 

    def determine_point_group(self,tolerance: float=1e-5) -> str:
        """  
        Determines the point group of the molecule based on the symmetry elements 
 
        Classification scheme:
        1. Linear molecules: C∞v or D∞h
        2. Special groups: T, O, I (tetrahedral, octahedral, icosahedral)
        3. High symmetry groups with principal axis
        4. Low symmetry groups: Cs, Ci, C1        
        """
        print("\n ============== Determining Point group ============== \n")
