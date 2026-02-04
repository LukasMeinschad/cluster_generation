import numpy as np
from typing import List, Tuple, Optional, Union, Protocol
import matplotlib.pyplot as plt
from dataclasses import dataclass
import time
import matplotlib.pyplot as plt

@dataclass 
class ReferenceFrame:
    """  
    Dataclass for Reference Frame
    """
    axes: np.ndarray # 3x3 array where each column is an axis
    origin: np.ndarray # 3 element array for origin point
    eigenvalues: Optional[np.ndarray] = None # Optional eigenvalues associated with the reference frame axes

    @property 
    def x_axis(self) -> np.ndarray:
        """ Returns the x-axis of the reference frame """
        return self.axes[:,0]
    @property 
    def y_axis(self) -> np.ndarray:
        """ Returns the y-axis of the reference frame """
        return self.axes[:,1]
    @property 
    def z_axis(self) -> np.ndarray:
        """ Returns the z-axis of the reference frame """
        return self.axes[:,2]

@dataclass
class Quaternion:
    """
    Quaternion class to represent:
    q = w + xi + yj + zk
    """
    w: float # scalar part
    x: float # i component 
    y: float # j component 
    z: float # k component

    @staticmethod 
    def from_axis_angle(axis: np.ndarray, angle_rad: float) -> "Quaternion":
        """  
        Creates a Quaternion from an axis and angle representation 
        q = cos(θ/2) + sin(θ/2)(xi + yj + zk)
        """
        axis_normalizes = GeometryOps.normalize_vector(axis)
        half_angle = angle_rad / 2.0
        sin_half = np.sin(half_angle)
        cos_half = np.cos(half_angle)
        return Quaternion(
            w=cos_half,
            x=axis_normalizes[0] * sin_half,
            y=axis_normalizes[1] * sin_half,
            z=axis_normalizes[2] * sin_half
        )
    
    @staticmethod 
    def from_two_vectors(v1: np.ndarray, v2: np.ndarray) -> "Quaternion":
        """
        Creates a quaternion that rotates vector v1 to vector v2 

        + First normalize both vectors 
        + The ncompute dot and cross product d = v1 * v2 and c = v1 x v2
        + If d =-1 the vectors are opposite and we have to choose rotation axis manually
        + The quaternion that carries out this rotation by:
                + w = sqrt((1 + d) / 2)
                + (x,y,z) = c / (2w)
        """
        v1_norm = GeometryOps.normalize_vector(v1)
        v2_norm = GeometryOps.normalize_vector(v2)
        
        # Parallel case 
        if np.allclose(v1_norm, v2_norm):
            return Quaternion(1.0, 0.0, 0.0, 0.0) # Idenity Quaternion
        
        # Antiparallel case 
        if np.allclose(v1_norm, -v2_norm):
            # Find the orthogonal vector
            if abs(v1_norm[0]) > 1e-8 or abs(v1_norm[1]) > 1e-8:
                ortho = np.array([-v1_norm[1], v1_norm[0], 0.0])
            else:
                ortho = np.array([0.0, -v1_norm[2], v1_norm[1]])
            ortho = GeometryOps.normalize_vector(ortho)
            return Quaternion(0.0, ortho[0], ortho[1], ortho[2])
        
        # General Case
        d = np.dot(v1_norm,v2_norm)
        c = np.cross(v1_norm,v2_norm)
        w = np.sqrt((1 + d) / 2.0)
        xyz = c / (2.0 * w)

        return Quaternion(w, xyz[0], xyz[1], xyz[2])

    def normalize(self) -> "Quaternion":
        """  
        Normalizes the quaternion to unit length
        """
        norm = np.sqrt(self.w**2 + self.x**2 + self.y**2 + self.z**2)
        if norm < 1e-8:
            raise ValueError("Cannot normalize zero-length quaternion")
        return Quaternion(
            w=self.w / norm,
            x=self.x / norm,
            y=self.y / norm,
            z=self.z / norm
        )
    
    def conjugate(self) -> "Quaternion":
        """  
        Returns the conjugate quaternion 
        q = w - xi - yj - zk
        """
        return Quaternion(
            w=self.w,
            x=-self.x,
            y=-self.y,
            z=-self.z
        )
    
    def multiply(self, other: "Quaternion") -> "Quaternion":
        """
        Multiplies two quaternions
        if t = rs then we have components
        + t_0 = (r_0s_0 - r_1s_1 - r_2s_2 - r_3s_3)
        + t_1 = (r_0s_1 + r_1s_0 + r_2s_3 - r_3s_2)
        + t_2 = (r_0s_2 + r_1s_3 + r_2s_0 - r_3s_1)
        + t_3 = (r_0s_3 - r_1s_2 + r_2s_1 + r_3s_0)
        """
        w = (self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z)
        x = (self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y)
        y = (self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x)
        z = (self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w)
        return Quaternion(w, x, y, z)
    
    def to_rotation_matrix(self) -> np.ndarray:
        """  
        Converts a quaternion to a rotation matrix
        """
        w,x,y,z = self.w, self.x, self.y, self.z

        return np.array([
            [1- 2*y**2 - 2*z**2, 2*x*y - 2*w*z, 2*x*z + 2*w*y],
            [2*x*y + 2*w*z, 1 - 2*x**2 - 2*z**2, 2*y*z - 2*w*x],
            [2*x*z - 2*w*y, 2*y*z + 2*w*x, 1- 2*x**2 - 2*y**2]
        ])
    
    def rotate_point(self, point: np.ndarray) -> np.ndarray:
        """  
        Rotates a point using quaternion rotation:
        1. Convert the point to be rotated into a pure quaternion p = 0 + xi + yj + zk
        2. Perform the rutation using qpq^-1
        """ 
    
        # Convert point to pure quaternion
        p = Quaternion(0.0, point[0], point[1], point[2])
        # Rotate using qpq^-1
        q_conj = self.conjugate()
        p_rotated = self.multiply(p).multiply(q_conj)
        return np.array([p_rotated.x, p_rotated.y, p_rotated.z])

    def rotate_points(self, points: np.ndarray) -> np.ndarray:
        """   
        Rotates multiple points using vectorized quaternion rotation
        We use the fast Method given by John D. Cook 
        https://www.johndcook.com/blog/2021/06/16/faster-quaternion-rotations/

        Let t = 2(q_1, q_2,q_3) x (v_1, v_2, v_3) = (t_1, t_2, t_3)
        Then the rotated vector v' is given by:
        v = (v_1, v_2, v_3) + q_0 t + q x t
        """ 
        q0, q1, q2, q3 = self.w, self.x, self.y, self.z
        q_vec = np.array([q1, q2, q3])
        t = 2.0 * np.cross(q_vec, points)

        v_rotated = points + q0 * t + np.cross(q_vec, t)
        return v_rotated

class GeometryOps:
    """   
    Class for static geometric operations
    """
    @staticmethod
    def centroid(coords: np.ndarray) -> np.ndarray: 
        """   
        Calculates the geometric centroid of a set of coordinates 
        """
        return np.mean(coords, axis=0)
    
    @staticmethod 
    def center_of_mass(coords: np.ndarray, masses: np.ndarray) -> np.ndarray:
        """   
        Calculates the Center of Mass of a set of coordinates with given masses
        """
        if masses.ndim == 1: 
            masses = masses[:, np.newaxis]
        return np.sum(masses * coords, axis=0) / np.sum(masses)
    
    @staticmethod 
    def inertia_tenstor(coords: np.ndarray, masses: np.ndarray) -> np.ndarray:
        """    
        Calculates the moment of inertia tensor for a set of coordinates with given masses
        
        I_ij = sum_k m_k (r_k^2 delta_ij - r_ki r_kj)
        """ 
        x,y,z = coords[:,0], coords[:,1], coords[:,2]
        m = masses

        # Vectorized calculation
        Ixx = np.sum(m * (y**2 + z**2))
        Iyy = np.sum(m * (x**2 + z**2))
        Izz = np.sum(m * (x**2 + y**2))
        Ixy = -np.sum(m * x * y)
        Ixz = -np.sum(m * x * z)
        Iyz = -np.sum(m * y * z)

        return np.array([[Ixx, Ixy, Ixz],
                         [Ixy, Iyy, Iyz],
                         [Ixz, Iyz, Izz]])
    
    @staticmethod 
    def normalize_vector(v: np.ndarray) -> np.ndarray: 
        """   
        Normalizes a vector to unit length
        """
        norm = np.linalg.norm(v)
        if norm < 1e-8:
            raise ValueError("Cannot normalize zero-length vector")
        return v / norm
    
    @staticmethod 
    def rotation_matrix_rodrigues(axis: np.ndarray, angle_rad: float) -> np.ndarray: 
        """   
        Construct the Rotation Matrix using Rodrigues' rotation formula

        R = I + sin(θ)K + (1 - cos(θ))K^2
        where K is the cross-product matrix of the unit axis vector
        """
        k = GeometryOps.normalize_vector(axis)
        kx, ky, kz = k

        # Cross-Product Matrix K
        K = np.array([[0, -kz, ky],
                      [kz, 0, -kx],
                      [-ky, kx, 0]])
        
        sin_theta = np.sin(angle_rad)
        cos_theta = np.cos(angle_rad)

        K_squared = K @ K 
        return np.eye(3) + sin_theta * K + (1 - cos_theta) * K_squared
    
    @staticmethod  
    def align_vectors(a: np.ndarray, b: np.ndarray) -> np.ndarray: 
        """  
        Constructs a rotation matrix that aligns vector a to vector b
        """ 
        a_norm = GeometryOps.normalize_vector(a)
        b_norm = GeometryOps.normalize_vector(b)

        # Handle the parallel case
        if np.allclose(a_norm, b_norm):
            return np.eye(3)
        
        # Antiparallel case: 180° rotation orthogonal axis 
        if np.allclose(a_norm, -b_norm): 
            # Find orthogonal vector 
            if abs(a_norm[0]) > 1e-8 or abs(a_norm[1]) > 1e-8:
                v = np.array([-a_norm[1], a_norm[0], 0.0])
            else:
                v = np.array([0.0, -a_norm[2], a_norm[1]])
            v = GeometryOps.normalize_vector(v)
            # Make Householder reflection matrix R = 2vv^T - I
            return 2 * np.outer(v, v) - np.eye(3)
        
        # General case with Rodrigues' rotation formula
        v = np.cross(a_norm, b_norm)
        s = np.linalg.norm(v)
        c = np.dot(a_norm, b_norm)

        # Cross-Product Matrix of v
        v_x = np.array([[0, -v[2], v[1]],
                        [v[2], 0, -v[0]],
                        [-v[1], v[0], 0]])
        
        # Rodrigues Formula
        return np.eye(3) + v_x + (v_x @ v_x) * ((1 - c) / (s**2))
    
    @staticmethod 
    def rotation_matrix_quaternion(axis: np.ndarray, angle_rad: float) -> np.ndarray:
        """  
        Constructs a rotation matrix from an axis and angle using quaternions
        """
        q = Quaternion.from_axis_angle(axis, angle_rad)
        return q.to_rotation_matrix()
    
    @staticmethod 
    def align_vectors_quaternion(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """  
        Constructs a rotation matrix that aligns vector a to vector b using quaternions
        """
        q = Quaternion.from_two_vectors(a, b)
        return q.to_rotation_matrix()


    
class Transformation:
    """ 
    Optimizes Class for 3D-Transformation of Molecules
    """ 
    def __init__(self, name: str = "Transformation"):
        """  
        Initializes the Transformation class
        """
        self.name = name 

    # =========================== Center Point Calculations ===========================

    def get_center(
            self,
            molecule: Union["Molecule", "Submolecule"],
            method: str = "centroid"
        ) -> np.ndarray:
        """ 
        Determines the center of a molecules 

        Args:
            molecule: Molecule or Submolecule object
            method: Method to calculate center point ("centroid" or "com")
        Returns:
            np.ndarray: Center point coordinates
        """
        if method == "centroid":
            return GeometryOps.centroid(molecule.coordinates)
        if method == "com":
            return GeometryOps.center_of_mass(molecule.coordinates, molecule.masses)
        else: 
            raise ValueError("Unknown method '{method}'. Method must be 'centroid' or 'com'")
        
    def get_center_coords( 
            self, 
            coords: np.ndarray,
            masses: Optional[np.ndarray] = None,
            method: str = "centroid"
        ) -> np.ndarray:
        """   
        Determines the center point for a given set of coordinates
        """
        if method == "centroid":
            return GeometryOps.centroid(coords)
        if method == "com":
            if masses is None:
                raise ValueError("Masses must be provided to calculate center of mass")
            return GeometryOps.center_of_mass(coords, masses)
        else:
            raise ValueError("Unknown method '{method}'. Method must be 'centroid' or 'com'")
        
    # ========================= Basic Transformations ==========================

    def translate(
            self, 
            molecule: Union["Molecule", "Submolecule"],
            vector: Union[List[float], np.ndarray],
            in_place: bool = True
        ) -> np.ndarray:
        """ 
        Translates a molecule with a vector T_v(p) = p + v

        Args:
            molecule: Molecule or Submolecule object
            vector: Translation vector
            in_place: Whether to modify the molecule in place or return new coordinates 
        """
        vec = np.asarray(vector, dtype=np.float64)
        if vec.shape != (3,):
            raise ValueError("Translation vector must be of shape (3,)")
        
        new_coords = molecule.coordinates + vec

        if in_place:
            # Modify molecule coordinates in place
            molecule.coordinates = new_coords
            return None
        else:
            return new_coords
        
    def rotate(
            self, 
            molecule: Union["Molecule", "Submolecule"],
            axis: Union[List[float], np.ndarray],
            angle_degree: float,
            in_place: bool = True 
        ) -> Optional[np.ndarray]:
        """   
        Rotates the molecule around a given axis using Rodrigues rotation formula
        
        Args:
            molecule: Molecule or Submolecule object
            axis: Rotation axis
            angle_degree: Rotation angle in degrees
            in_place: Whether to modify the molecule in place or return new coordinates 
        """
        axis_array = np.asarray(axis, dtype=np.float64)
        angle_rad = np.radians(angle_degree)

        R = GeometryOps.rotation_matrix_rodrigues(axis_array, angle_rad)
        new_coords = molecule.coordinates @ R.T

        if in_place:
            molecule.coordinates = new_coords
            return None
        else:
            return new_coords
        
    def rotate_quaternion(
            self, 
            molecue: Union["Molecule", "Submolecule"],
            axis: Union[List[float], np.ndarray],
            angle_degree: float, 
            in_place: bool = True
        ) -> Optional[np.ndarray]:
        """   
        Rotates the molecule using quaternion rotation
        
        Args:
            molecue: Molecule or Submolecule object
            axis: Rotation axis
            angle_degree: Rotation angle in degrees
            in_place: Whether to modify the molecule in place or return new coordinates 
        """
        axis_array = np.asarray(axis, dtype=np.float64)
        angle_rad = np.radians(angle_degree)

        q = Quaternion.from_axis_angle(axis_array, angle_rad)
        q = q.normalize()

        
        # Still use rotation matrix because of linear algebra speed
        R = q.to_rotation_matrix()
        new_coords = molecue.coordinates @ R.T 

        if in_place:
            molecue.coordinates = new_coords
            return None
        else:
            return new_coords
        
    def center_molecule( 
            self, 
            molecule: Union["Molecule", "Submolecule"],
            method: str = "centroid",
            in_place: bool = True
        ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Centers a molecule to a center point

        Args:
            molecule: Molecule or Submolecule object
            method: Method to calculate center point ("centroid" or "com")
            in_place: Whether to modify the molecule in place or return new coordinates

        Returns:
            Tuple of (center, new_coords) if in_place is False
            Otherwise returns center_point only
        """
        center = self.get_center(molecule, method=method)
        new_coords = molecule.coordinates - center
        if in_place:
            molecule.coordinates = new_coords 
            return center, None
        else:
            return center, new_coords
        
    # ======================= Reference Frame Definitions =======================

    def compute_reference_frame( 
            self, 
            coords: np.ndarray,
            masses: np.ndarray,
            center: Optional[np.ndarray] = None
        ) -> ReferenceFrame:
        """
        Computes the Reference Frame based on the inertia tensor of a set of coordinates

        Args:
            coords: Coordinates of the molecule or submolecule
            masses: Masses of the atoms
            center: Optional center point to translate coordinates before computing inertia tensor
        Returns:
            ReferenceFrame object
        """
        if center is None:
            center = GeometryOps.center_of_mass(coords, masses)
        

        coords_new = coords.copy()
        coords_new -= center
        # Compute inertia tensor
        I = GeometryOps.inertia_tenstor(coords_new, masses)
        eigenvalues, eigenvectors = np.linalg.eigh(I)

        # Eigenvectors are columns of the Reference Frame axes
        return ReferenceFrame(
            axes=eigenvectors,
            origin=center,
            eigenvalues=eigenvalues
        )
    
    def set_reference_frame(
            self,
            molecule: Union["Molecule", "Submolecule"],
            method: str = "centroid",
            print_info: bool = False,
        ) -> ReferenceFrame:
        """   
        Sets the reference frame to a center of a molecule and its principal axes

        Args: 
            molecule: Molecule or Submolecule object
            method: Method to calculate center point ("centroid" or "com")
            print_info: Whether to print information about the reference frame
        """
        if len(molecule.coordinates) < 3:
            raise ValueError("Molecule must have at least 3 atoms to define a reference frame")
        
        # Center Molecule
        center = self.get_center(molecule, method=method)
        molecule.coordinates -= center

        ref_frame = self.compute_reference_frame(
            molecule.coordinates,
            molecule.masses,
            center=np.zeros(3) # Already centered
        )
        if print_info:
            print(f"Molecule center: {center}")
            print(f"Inertia tensor eigenvalues: {ref_frame.eigenvalues}")
            print(f"Reference frame axes:")
            print(f"X: {ref_frame.x_axis}")
            print(f"Y: {ref_frame.y_axis}")
            print(f"Z: {ref_frame.z_axis}")
        
        return ref_frame
    
    def set_reference_frame_submolecule(
            self,
            submolecule: "Submolecule",
            method: str = "centroid",
            print_info: bool = False,
        ) -> ReferenceFrame:
        """   
        Computes the Reference Frame based on a Submolecule

        Further aligns all atoms in the parent Molecule accordingly
        """
        if len(submolecule.coordinates) < 3:
            raise ValueError("Submolecule must have at least 3 atoms to define a reference frame")
        # Center Submolecule
        center = self.get_center(submolecule, method=method)
        submolecule.coordinates -= center
        ref_frame = self.compute_reference_frame(
            submolecule.coordinates,
            submolecule.masses,
            center=np.zeros(3) # Already centered
        )
        # Translate Parent Molecule
        if hasattr(submolecule, 'parent') and submolecule.parent is not None: 
            parent_mol = submolecule.parent
            parent_mol.coordinates -= center
        if print_info:
            print(f"Submolecule center: {center}")
            print(f"Inertia tensor eigenvalues: {ref_frame.eigenvalues}")
            print(f"Reference frame axes:")
            print(f"X: {ref_frame.x_axis}")
            print(f"Y: {ref_frame.y_axis}")
            print(f"Z: {ref_frame.z_axis}")
        return ref_frame
    
    def set_reference_frame_from_atoms(
            self,
            molecule: Union["Molecule", "Submolecule"],
            atom_labels: List[int],
            method: str = "centroid",
            print_info: bool = False,
        ) -> ReferenceFrame:
        """   
        Computes the Reference Frame based on a set of atom labels

        Args:
            molecule: Molecule or Submolecule object
            atom_labels: List of atom labels to define the reference frame
            method: Method to calculate center point ("centroid" or "com")
            print_info: Whether to print information about the reference frame
        """
        if len(atom_labels) < 3:
            raise ValueError("At least 3 atom labels must be provided to define a reference frame")

        # Extract coordinates and masses of specified atoms
        selected_coords = []
        selected_masses = []
        for label in atom_labels:
            coord = molecule.get_coords_by_label(label)
            if coord.ndim > 1:
                raise ValueError(f"Atom label {label} corresponds to multiple atoms")
            selected_coords.append(coord)
            selected_masses.append(molecule.get_mass_by_label(label))

        selected_coords = np.array(selected_coords)
        selected_masses = np.array(selected_masses)

        center = self.get_center_coords(selected_coords, masses=selected_masses, method=method)

        # Center the total Molecule
        molecule.coordinates -= center
        selected_coords -= center

        ref_frame = self.compute_reference_frame(
            selected_coords,
            selected_masses,
            center=np.zeros(3) # Already centered
        )
        if print_info:
            print(f"Molecule center (based on labels {atom_labels}): {center}")
            print(f"Inertia tensor eigenvalues: {ref_frame.eigenvalues}")
            print(f"Reference frame axes:")
            print(f"X: {ref_frame.x_axis}")
            print(f"Y: {ref_frame.y_axis}")
            print(f"Z: {ref_frame.z_axis}")

        return ref_frame



    # ======================== Functions for Stochastic Surface Walking ========================

    @staticmethod
    def mb_vector(mass: float, temperature: float = 300):
        """   
        Generates a normalized random vector following the Maxwell-Boltzmann velocity distribution

        Args:
            mass: Atomic mass in atomic mass units (amu)
            temperature: Temperature in Kelvin

        Returns:
            np.ndarray
        """
        # Boltzmann Konstant eV/K
        k_B = 8.617333262e-5 
        sigma = np.sqrt(k_B * temperature / mass)
        vector = np.random.normal(0,sigma,size=3) # 3 normally distributed components

        # Normalize
        norm = np.linalg.norm(vector)
        if norm < 1e-8:
            # Random unit vector if near zero
            vector = np.random.rand(3)
            norm = np.linalg.norm(vector)

        return vector / norm 
    
    def find_non_neighboring_atoms(
            self,
            molecule: Union["Molecule", "Submolecule"],
            distance_threshold: float = 3.0
    ) -> List[Tuple[int,int]]:
        """    
        Finds all pairs of atoms that are not in close contact (distance > threshold)

        Args:
            molecule: Molecule or Submolecule object
            distance_threshold: Minimum distance in Angstrom for atoms to be non-neighbouring
        """
        n_atoms = molecule.get_number_of_atoms()
        non_neighboring_pairs = []

        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                distance = np.linalg.norm(
                    molecule.coordinates[i] - molecule.coordinates[j]
                )
                if distance > distance_threshold:
                    non_neighboring_pairs.append((i,j))
        
        return non_neighboring_pairs
    
    def select_random_non_neighboring_pair(
            self,
            molecule: Union["Molecule", "Submolecule"],
            distance_threshold: float = 3.0
    ):
        """   
        Randomnly selects a pair of non_neighboruing atoms

        Returns:
            Tuple of (atom_A_idx,atom_B_idx) or None if no pairs exists
        """
        non_neighboring_pairs = self.find_non_neighboring_atoms(
            molecule,
            distance_threshold
        )
        if not non_neighboring_pairs:
            return None 
        
        r_int = np.random.randint(len(non_neighboring_pairs))
        return non_neighboring_pairs[r_int]

    def compute_bond_formation_vector(
            self,
            molecule: Union["Molecule", "Submolecule"],
            atom_A_idx: int,
            atom_B_idx: int,
            normalize: bool = True
    ) -> np.ndarray:
        """   
        Computes the bond formation mode N_i^l between two atoms

        N_i^l = (q_B - q_A) / |q_b - q_A|
        """
        q_A = molecule.coordinates[atom_A_idx]
        q_B = molecule.coordinates[atom_B_idx]

        v = q_B - q_A 

        if normalize:
            norm = np.linalg.norm(v)
            if norm < 1e-8:
                raise ValueError(f"Atoms {atom_A_idx} and {atom_B_idx} are in the same position")
            return v / norm 
        
        return v

    def generate_ssw_direction(
            self,
            molecule: Union["Molecule", "Submolecule"],
            lambda_mix: Optional[float] = None,
            temperature: float = 300.0,
            use_total_mass: bool = True,
            distance_threshold: float = 3.0
    ) -> np.ndarray:
        """   
        Generates the initial random direction N_i^0 for SSW Climing

        N_i^0 = (N_i^g + λ*N_i^l) / |N_i^g + λ*N_i^l|

        + N_i^g: Random vector from MW-Distribution at Temp T
        - N_i^l: Soft Local mode, currently --> Bond formation of non-neighboring atoms
        - lambda = Random Mixing Parameter
        """
        # Generate random lambda if not provides
        if lambda_mix is None:
            lambda_mix = np.random.uniform(0.1,1.5)

        # Get total mass
        if use_total_mass:
            effective_mass = molecule.get_total_mass()
        else:
            effective_mass = np.mean(molecule.masses)

        # Generate N_i^g 
        N_g = self.mb_vector(effective_mass, temperature)

        # Generate Random-Pair
        atom_pair = self.select_random_non_neighboring_pair(molecule,distance_threshold)
        atom_A_idx, atom_B_idx = atom_pair 
        N_l = self.compute_bond_formation_vector(
            molecule,
            atom_A_idx,
            atom_B_idx,
            normalize=True
        )
        N_0 = N_g + lambda_mix * N_l 
        norm_0 = np.linalg.norm(N_0)

        return N_0 / norm_0
    

    # ================== Dimer Method ================

    def setup_dimer_images(
            self,
            R0: np.ndarray,
            N: np.ndarray,
            delta_R: float = 0.005
    ) -> Tuple [np.ndarray, np.ndarray]:
        """   
        Sets up the two dimer images R1 and R2 seperated by 2* delta R
        uses the local mode N created beforehand

        R1 = R0 + N * ∆R
        R2 = R= - N * ∆R
        """
        N_norm = GeometryOps.normalize_vector(N) # Ensure normalization

        R1 = R0 + N_norm * delta_R
        R2 = R0 - N_norm * delta_R 

        return R1,R2 


    # ================

    
        



    
    # ======================= Alignment Functions =========================

    def align_molecule_to_target(
            self,
            molecule: Union["Molecule", "Submolecule"],
            ref_frame: ReferenceFrame,
            target_point: np.ndarray,
        ) -> np.ndarray:
        """
        Rotates the molecule so that the RefFrame z-axis aligns with the vector from origin to target_point

        Args:
            molecule: Molecule or Submolecule object
            ref_frame: ReferenceFrame object of the molecule
            target_point: Target point to align the z-axis towards
        Returns:
            np.ndarray: New coordinates of the aligned molecule
        """ 
        # Directions to target 
        direction_to_target = target_point - ref_frame.origin
        target_direction = GeometryOps.normalize_vector(direction_to_target)

        # Rotation matrix: z-axis --> target_direction
        R = GeometryOps.align_vectors(ref_frame.z_axis, target_direction)

        new_coords = molecule.coordinates @ R.T
        return new_coords
    
    def align_to_vector(
            self,
            molecule: Union["Molecule", "Submolecule"],
            current_vector: np.ndarray,
            target_vector: np.ndarray,
            in_place: bool = True
        ) -> Optional[np.ndarray]:
        """  
        Rotates a molecule such that current_vector aligns with target_vector
        
        Args:
            molecule: Molecule or Submolecule object
            current_vector: Current vector in the molecule
            target_vector: Target vector to align to
            in_place: Whether to modify the molecule in place or return new coordinates
        """
        R = GeometryOps.align_vectors(current_vector, target_vector)
        new_coords = molecule.coordinates @ R.T
        if in_place:
            molecule.coordinates = new_coords
            return None
        else:
            return new_coords
        
    # ======================= Selective Transformations =======================

    def translate_submolecule(
            self, 
            submolecule: "Submolecule",
            target_point: Union[List[float], np.ndarray],
            parent_molecule: Optional["Molecule"] = None,
            method: str = "centroid",
        ) -> None:
        """
        Translates the Submolecule to a target point and updates the parent 
        
        Args:
            submolecule: Submolecule object to translate
            target_point: Target point to translate the submolecule to
            parent_molecule: Optional parent Molecule object to update coordinates
            method: Method to calculate center point ("centroid" or "com")
        """
        # Translation vector
        center = self.get_center(submolecule, method=method)
        translation = np.asarray(target_point) - center

        # Update Submolecule
        submolecule.coordinates += translation

        # Update Parent Molecule if provided
        if parent_molecule is not None and hasattr(submolecule, 'get_atom_lables'):
            submol_labels = submolecule.get_atom_labels()

            # Mask
            mask = np.isin(parent_molecule.atom_labels, list(submol_labels))
            parent_molecule.coordinates[mask] += translation

    # ===================== Utility Functions ======================

    def bond_vector(
            self,
            molecule: Union["Molecule", "Submolecule"],
            label1: str,
            label2: str,
            normalize: bool = True
        ) -> np.ndarray:
        """
        Computes the bond vector from atom with label1 to atom with label2
        """
        coord1 = molecule.get_coords_by_label(label1)
        coord2 = molecule.get_coords_by_label(label2)

        if coord1.ndim > 1 or coord2.ndim > 1:
            raise ValueError("Atom labels must correspond to single atoms")
        vector = coord2 - coord1

        if normalize:
            vector = GeometryOps.normalize_vector(vector)
        return vector
 
    # ======================== Speed Testing Functions ======================== --IGNORE---
    
    def test_translate_speed(
            self,
            molecule: Union["Molecule", "Submolecule"]
        ) -> None:
        """   
        Tests the speed of the translate function for a random vector and N=[1000,5000,10000,50000,100000] iterations

        Returns a plot of time vs iterations
        """
        iterations = [1000, 5000, 10000, 50000, 100000]
        times = []
        random_vector = np.random.rand(3)
        for N in iterations:
            start_time = time.time()
            for _ in range(N):
                self.translate(molecule, random_vector, in_place=False)
            end_time = time.time()
            times.append(end_time - start_time)
        # Plotting
        plt.figure()
        plt.plot(iterations, times, marker='o')
        plt.xlabel("Number of Iterations")
        plt.ylabel("Time (s)")
        plt.title("Translation Speed Test")
        plt.grid()
        # Save in the figues path
        plt.savefig("figures/translation_speed_test.png")
        plt.close()

    def test_rotate_speed(
            self,
            molecule: Union["Molecule", "Submolecule"]
        ) -> None:
        """
        Test the speed of the rotate functions for a random axis and angle over N=[1000,5000,10000,50000,100000] iterations
        """
        iterations = [1000, 5000, 10000, 50000, 100000]
        times = []
        random_axis = np.random.rand(3)
        random_angle = np.random.uniform(0, 360)
        for N in iterations:
            start_time = time.time()
            for _ in range(N):
                self.rotate(molecule, random_axis, random_angle, in_place=False)
            end_time = time.time()
            times.append(end_time - start_time)
        # Plotting
        plt.figure()
        plt.plot(iterations, times, marker='o')
        plt.xlabel("Number of Iterations")
        plt.ylabel("Time (s)")
        plt.title("Rotation Speed Test")
        plt.grid()
        # Save in the figues path
        plt.savefig("figures/rotation_speed_test.png")
        plt.close()

    def compare_rotate_methods_speed(
            self,
            molecule: Union["Molecule", "Submolecule"]
        ) -> None:
        """
        Compares the speed of the rotate and rotate_quaternion methods over N=[1000,5000,10000,50000,100000] iterations
 
        """
        iterations = [1000, 5000, 10000, 50000, 100000]
        times_rodrigues = []
        times_quaternion = []
        random_axis = np.random.rand(3)
        random_angle = np.random.uniform(0, 360)

        for N in iterations:
            # Rodrigues method
            start_time = time.time()
            for _ in range(N):
                self.rotate(molecule, random_axis, random_angle, in_place=False)
            end_time = time.time()
            times_rodrigues.append(end_time - start_time)

            # Quaternion method
            start_time = time.time()
            for _ in range(N):
                self.rotate_quaternion(molecule, random_axis, random_angle, in_place=False)
            end_time = time.time()
            times_quaternion.append(end_time - start_time)
        


        

        # Plotting
        plt.figure()
        plt.plot(iterations, times_rodrigues, marker='o', label='Rodrigues')
        plt.plot(iterations, times_quaternion, marker='o', label='Quaternion')
        plt.xlabel("Number of Iterations")
        plt.ylabel("Time (s)")
        plt.title("Rotation Method Speed Comparison")
        plt.legend()
        plt.grid()
        # Save in the figues path
        plt.savefig("figures/rotation_method_speed_comparison.png")
        plt.close()