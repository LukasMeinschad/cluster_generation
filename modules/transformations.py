import numpy as np
from typing import List, Tuple, Optional, Union, Protocol, Dict
import matplotlib.pyplot as plt
from dataclasses import dataclass
import time
import random
import copy


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
    def random_so3() -> np.ndarray:
        """  
        Samples uniformly over the space of rotations.
        https://lavalle.pl/planning/node198.html
        """
        u1 = random.random()
        u2 = random.random()
        u3 = random.random()

        w = np.sqrt(1 - u1) * np.sin(2 * np.pi * u2)
        x = np.sqrt(1- u1) * np.cos(2 * np.pi * u2)
        y = np.sqrt(u1) * np.sin(2 * np.pi * u3)
        z = np.sqrt(u1) * np.cos(2 * np.pi * u3)

        q = Quaternion(w,x,y,z)
        return q.to_rotation_matrix()




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


    def from_rotation_matrix(R: np.ndarray) -> "Quaternion":
        """  
        Converts a rotation matrix to quaternion
        uses the method of the homepage  https://danceswithcode.net/engineeringnotes/quaternions/quaternions.html

        Is the numerically stable Shepperd Method

        Step 1: Compute the magnitude of each quaternion component. This leaves the sign of each component undetermined
        Step 2: Find the largest component, assume its sign is positive,
                then compute t he ramaining components from the of diagonal elements
                Taking the largest magnitude avoids division by small numbers
        """
        # Step 1 compute the magnitude from diagonal elements
        q0_sq = (1 + R[0,0] + R[1,1] + R[2,2]) / 4.0
        q1_sq = (1 + R[0,0] - R[1,1] - R[2,2]) / 4.0
        q2_sq = (1 - R[0,0] + R[1,1] - R[2,2]) / 4.0 
        s3_sq = (1 - R[0,0] - R[1,1] + R[2,2]) / 4.0

        # Pick largest to avoid division by small numbers
        idx = np.argmax([q0_sq, q1_sq, q2_sq, s3_sq])
        if idx == 0:
            q0 = np.sqrt(max(q0_sq, 0))
            q1 = (R[2,1] - R[1,2]) / (4 * q0)
            q2 = (R[0,2] - R[2,0]) / (4 * q0)
            q3 = (R[1,0] - R[0,1]) / (4 * q0)
        elif idx == 1:
            q1 = np.sqrt(max(q1_sq, 0))
            q0 = (R[2,1] - R[1,2]) / (4 * q1)
            q2 = (R[0,1] + R[1,0]) / (4 * q1)
            q3 = (R[1,2] + R[2,1]) / (4 * q1)
        elif idx == 2:
            q2 = np.sqrt(max(q2_sq, 0))
            q0 = (R[0,2] - R[2,0]) / (4 * q2)
            q1 = (R[0,1] + R[1,0]) / (4 * q2)
            q3 = (R[1,2] + R[2,1]) / (4 * q2)
        else:
            q3 = np.sqrt(max(s3_sq, 0))
            q0 = (R[1,0] - R[0,1]) / (4 * q3)
            q1 = (R[0,2] + R[2,0]) / (4 * q3)
            q2 = (R[1,2] + R[2,1]) / (4 * q3)
        return Quaternion(q0, q1, q2, q3)
    

    
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
    def inertia_tensor(coords: np.ndarray, masses: np.ndarray) -> np.ndarray:
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
    def kabsch_rotation(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
        """  
        Computes the optimal rotation matrix R that minimizes the RMSD between two centered
        coordinate sets using the Kabsch algorithm

        Solves min_R || P - R @ Q || via SVD of Cross-Covariand H = P Q
        
        Args:
            P: Target coordinates (Nx3)
            Q: Source coordinates (Nx3)
        
        Returns:
            R: 3x3 rotation matrix that best aligns Q to P
        """
        H = P.T @ Q
        U, S, Vt = np.linalg.svd(H)

        # Ensure proper rotation (det = +1) not reflection
        d = np.linalg.det(U @ Vt)
        sign_matrix = np.diag([1, 1, np.sign(d)])
        R = U @ sign_matrix @ Vt
        return R
    



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
            raise ValueError(f"Unknown method '{method}'. Method must be 'centroid' or 'com'")
        
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
            raise ValueError(f"Unknown method '{method}'. Method must be 'centroid' or 'com'")
        
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
            molecule: Union["Molecule", "Submolecule"],
            axis: Union[List[float], np.ndarray],
            angle_degree: float, 
            in_place: bool = True
        ) -> Optional[np.ndarray]:
        """   
        Rotates the molecule using quaternion rotation
        
        Args:
            molecule: Molecule or Submolecule object
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
        new_coords = molecule.coordinates @ R.T 

        if in_place:
            molecule.coordinates = new_coords
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
        I = GeometryOps.inertia_tensor(coords_new, masses)
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
    
    def set_reference_frame_from_indices(
            self,
            molecule: Union["Molecule", "Submolecule"],
            atom_indices: List[int],
            method: str = "centroid",
            print_info: bool = False,
        ) -> ReferenceFrame:
        """   
        Computes the Reference Frame based on a set of atom indices
        """
        if len(atom_indices) < 3:
            raise ValueError("At least 3 atom indices must be provided to define a reference frame")
        # Extract coordinates and masses of specified atoms
        selected_coords = molecule.coordinates[atom_indices]
        selected_masses = molecule.masses[atom_indices]
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
            print(f"Molecule center (based on indices {atom_indices}): {center}")
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
        if parent_molecule is not None and hasattr(submolecule, 'get_atom_labels'):
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
    

    def kabsch_interpolate(
            self,
            mol_a: "Molecule",
            mol_b: "Molecule",
            lambdas: Optional[List[float]] = None,
            n_points: int = 3
        ) -> List["Molecule"]:
        """
        Generate interpolated structures between two molecules using Kabsch alignment

        1. Center both structures at the origin
        2. Compute optimal rotation R (Kabsch) to align B onto A
        3. For each lambda in [0,1]: r_new = (1-lambda) * r_A + lambda * (R @ r_B)

        Args:
            mol_a: First Molecule object
            mol_b: Second Molecule object (must have same number of atoms)
            lambdas: Optional list of interpolation parameters between 0 and 1
            n_points: Number of interpolated points to generate if lambdas not provided
        """
        if len(mol_a.coordinates) != len(mol_b.coordinates):
            raise ValueError("Molecules must have the same number of atoms for interpolation")
        if lambdas is None:
            lambdas = np.linspace(0, 1, n_points + 2)[1:-1]

        # Center both structures 
        coords_a = mol_a.coordinates.copy()
        coords_b = mol_b.coordinates.copy()
        centroid_a = np.mean(coords_a, axis=0)
        centroid_b = np.mean(coords_b, axis=0)
        P = coords_a - centroid_a
        Q = coords_b - centroid_b
        # Kabsch Rotation
        R = GeometryOps.kabsch_rotation(P,Q)
        Q_aligned = (R @ Q.T).T

        interpolated = []
        for lam in lambdas:
            coords_new = (1- lam) * P + lam * Q_aligned
            # Shift back to centroid of A
            coords_new += centroid_a
            new_mol = mol_a.copy()
            new_mol.coordinates = coords_new
            interpolated.append(new_mol)

        return interpolated


class MolecularOperators:
    """   
    Base class for molecular geometry operators.
    Uses Transformation and GeometryOps for efficient operations
    """
    def __init__(self, simulation_box: Optional["SimulationBox"] = None):
        """
        Initialize Molecular Operators

        Args:
            simulation_box: Optional SimulationBox object to constrain operations
        """
        self.simulation_box = simulation_box
        self.transformer = Transformation()

    def _apply_box_constraints(self, molecule: "Molecule", original: "Molecule") -> "Molecule":
        """
        Apply simulation box constraints to the molecule after transformation

        Args:
            molecule: Transformed Molecule object
            original: Original Molecule object before transformation (returned if constraints violated)
        """ 
        if self.simulation_box is not None:
            new_coords, is_valid = self.simulation_box.apply_boundary_conditions(
                molecule.coordinates, method="reflect"
            )
            if not is_valid:
                return original.copy()  
            molecule.coordinates = new_coords
        return molecule

    def _compute_scaling_factor(
            self,
            molecule: "Molecule",
            operator_type: str = "local"
    ) -> float:
        """  
        Computes a normalized scaling factor
        """ 


        if self.simulation_box is None:
            return 1.0
        coords = molecule.coordinates
        from scipy.spatial.distance import pdist

        if len(coords) > 1:
            pairwise_distances = pdist(coords)
            mol_diameter = np.max(pairwise_distances)
        else:
            mol_diameter = 2.0 * molecule.covalent_radii.max()
        # Add average covalent radius to account for size of atoms
        avg_radius = np.mean(molecule.covalent_radii)
        effective_mol_size = mol_diameter + 2 * avg_radius


        # Get characteristic length
        if self.simulation_box.box_type == "cubic":
            box_length = self.simulation_box.box_dimensions
        elif self.simulation_box.box_type == "sphere":
            box_length = self.simulation_box.radius * 2
        else:
            raise ValueError(f"Unknown box type '{self.simulation_box.box_type}'")
        
        free_space = max(box_length - effective_mol_size, 0.1)
        size_ratio = effective_mol_size / free_space



        if operator_type == "local":
            # Here one can simply change scalings around
            base_scale = 0.05
            max_scale = 0.4
            scale = base_scale + (max_scale - base_scale) * (1 - np.exp(-size_ratio))
        if operator_type == "nonlocal":
            base_scale = 0.2
            max_scale = 0.8
            scale = base_scale + (max_scale - base_scale) * (1 - np.exp(-size_ratio * 0.5))
        return np.clip(scale, base_scale, max_scale)

                                         
    @staticmethod
    def _random_unit_vetor() -> np.ndarray:
        """   
        Generates a random unit vector on the unit sphere
        Uses spherical coordinates for uniform distribution

        Returns:
            Random unit vector (3,)
        """
        phi = np.random.uniform(0, 2 * np.pi)
        costheta = np.random.uniform(-1, 1)
        theta = np.arccos(costheta)

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        return np.array([x, y, z])

    @staticmethod
    def _select_random_submolecule(submolecule_indices: List[List[int]]) -> List[int]:
        """   
        selects a random submolecule from the list
        """
        return random.choice(submolecule_indices)
    
    @staticmethod 
    def _select_two_submolecules(submolecule_indices: List[List[int]]) -> Tuple[List[int], List[int]]:
        """   
        Selects two random submolecules from the list
        """
        if len(submolecule_indices) < 2:
            raise ValueError("At least two submolecules must be provided to select two")
        return random.sample(submolecule_indices, 2)

class LocalOperators(MolecularOperators):
    """   
    Local Operators for small perturbations of the molecule
    """ 
    def random_displacement_submolecule(self,
                                        molecule: "Molecule",
                                        submolecule_indices: List[List[int]],
                                        delta_range: Tuple[float, float] = (-0.5, 0.5),
                                        adaptive: bool = True) -> "Molecule":
        """  
        Randomly displaces a submolecule by a random vector

        Args:
            molecule: Molecule object to modify
            submolecule_indices: List of lists of atom indices defining submolecules
            delta_range: Range for random displacement in Angstrom (default is (-0.3, 0.3))
            adaptive: If True, scale range by box/molecule size ratio
        """
        mol_copy = molecule.copy()
        n_atoms = len(mol_copy.coordinates)

        # select random submolecule
        submol_idx = self._select_random_submolecule(submolecule_indices)
        
        # Validate
        if any(idx >= n_atoms for idx in submol_idx):
            raise ValueError("Submolecule indices must be within the range of molecule atoms")
        
        # Apply adaptive scaling
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="local")
            effective_range = (delta_range[0] * scale, delta_range[1] * scale)
        else:
            effective_range = delta_range
        
        # Generate random displacement vector
        delta_value = random.uniform(*effective_range)
        direction = self._random_unit_vetor()
        displacement = direction * delta_value
        # Apply displacement to selected submolecule atoms
        mol_copy.coordinates[submol_idx] += displacement
        return self._apply_box_constraints(mol_copy, molecule)

    def random_rotation_submolecule(
            self,
            molecule: "Molecule",
            submolecule_indices: List[List[int]],
            angle_range: Tuple[float, float] = (-40, 40),
            adaptive: bool = True) -> "Molecule":
        """  
        Randomly rotate a submolecule around its center of mass using a random axis and a small angle.
        Analogous to the nonlocal twist but with local-scale angle

        Args:
            molecule: Molecule object to modify
            submolecule_indices: List of lists of atom indices defining submolecules
            angle_range: Range for random rotation angle in degrees (default is (-30, 30))
            adaptive: If True, scale angle range by box/molecule size ratio
        """
        mol_copy = molecule.copy()

        submol_idx = self._select_random_submolecule(submolecule_indices)
        submol_coords = mol_copy.coordinates[submol_idx]
        submol_masses = mol_copy.masses[submol_idx] 
        center = GeometryOps.center_of_mass(submol_coords, submol_masses)

        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="local")
            effective_range = (angle_range[0] * scale, angle_range[1] * scale)
        else:
            effective_range = angle_range

        angle = np.radians(random.uniform(*effective_range))
        axis = self._random_unit_vetor()

        R = GeometryOps.rotation_matrix_rodrigues(axis, angle)
        coords_centered = mol_copy.coordinates[submol_idx] - center
        coords_rotated = coords_centered @ R.T
        mol_copy.coordinates[submol_idx] = coords_rotated + center
        return self._apply_box_constraints(mol_copy, molecule)
    
    def local_twist_submolecule(
            self,
            molecule: "Molecule",
            submolecule_indices: List[List[int]],
            angle_range: Tuple[float, float] = (-15, 15),
            adaptive: bool = True) -> "Molecule":
        """  
        Local Twist: Rotate one submolecule around the axis connecting two submolecule
        centers of mass by a small angle

        This is physically motivated - the twist axis is the intermolecular COM-to-COM direction, so it
        presevers the intermolecular distance while changing the relative orientation

        Algorithm:
            1. Select two submolecules (reference + rotating)
            2. Compute the COM-to-COM axis
            3. Rotate the second submolecule around that axis by a small angle

        Args:
            molecule: Molecule object to modify
            submolecule_indices: List of lists of atom indices defining submolecules
            angle_range: Range for random rotation angle in degrees (default is (-15, 15))
            adaptive: If True, scale angle range by box/molecule size ratio
        """
        mol_copy = molecule.copy()
        ref_submol, rot_submol = self._select_two_submolecules(submolecule_indices)
        ref_center = GeometryOps.center_of_mass(mol_copy.coordinates[ref_submol], mol_copy.masses[ref_submol])
        rot_center = GeometryOps.center_of_mass(mol_copy.coordinates[rot_submol], mol_copy.masses[rot_submol])
        
        # Twist axis: COM to COM
        twist_axis = rot_center - ref_center
        if np.linalg.norm(twist_axis) < 1e-8:
            return mol_copy # Cannot define twist because submolecules overlap
        axis = GeometryOps.normalize_vector(twist_axis)

        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="local")
            effective_range = (angle_range[0] * scale, angle_range[1] * scale)
        else:
            effective_range = angle_range
        angle = np.radians(random.uniform(*effective_range))
        R = GeometryOps.rotation_matrix_rodrigues(axis, angle)
        coords_centered = mol_copy.coordinates[rot_submol] - rot_center
        coords_rotated = coords_centered @ R.T
        mol_copy.coordinates[rot_submol] = coords_rotated + rot_center
        return self._apply_box_constraints(mol_copy, molecule)

    def correlated_displacement(
        self,
        molecule: "Molecule",
        submolecule_indices: List[List[int]],
        delta_range: Tuple[float, float] = (-0.3, 0.3),
        adaptive: bool = True) -> "Molecule":
        """ 
        Displace two or more submolecules simultaneously in a correlated manner

        Args:
            molecule: Molecule object to modify
            submolecule_indices: List of lists of atom indices defining submolecules
            delta_range: Range for random displacement in Angstrom (default is (-0.3, 0.3))
            adaptive: If True, scale range by box/molecule size ratio
        """
        mol_copy = molecule.copy()

        n_select = random.randint(2, len(submolecule_indices))
        selected = random.sample(submolecule_indices, n_select)
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="local")
            effective_range = (delta_range[0] * scale, delta_range[1] * scale)
        else:
            effective_range = delta_range
        magnitude = random.uniform(*effective_range)

        direction = self._random_unit_vetor()
        displacement = direction * magnitude
        for submol_idx in selected:
            mol_copy.coordinates[submol_idx] += displacement
        return self._apply_box_constraints(mol_copy, molecule)

    def small_principal_axis_rotation(
            self,
            molecule: "Molecule",
            submolecule_indices: List[List[int]],
            angle_range: Tuple[float, float] = (-15,15),
            adaptive: bool = True) -> "Molecule":
        """  
        Rotate a single submolecule around one of its own principal axes
        computed from the inertia tensor by a small angle

        Algorithm:
            1. Select random submolecule
            2. Compute its reference frame and principal axes
            3. Select random principal axis
            4. Rotate submolecule around that axis by a small angle

        Args:            
            molecule: Molecule object to modify
            submolecule_indices: List of lists of atom indices defining submolecules
            angle_range: Range for random rotation angle in degrees (default is (-15, 15))
            adaptive: If True, scale angle range by box/molecule size ratio
        """
        mol_copy = molecule.copy()
        submol_idx = self._select_random_submolecule(submolecule_indices)
        coords = mol_copy.coordinates[submol_idx]
        masses = mol_copy.masses[submol_idx]
        center = GeometryOps.center_of_mass(coords, masses)

        ref_frame = self.transformer.compute_reference_frame(coords, masses, center=center)
        principal_axes = [ref_frame.x_axis, ref_frame.y_axis, ref_frame.z_axis]
        axis = random.choice(principal_axes)

        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="local")
            effective_range = (angle_range[0] * scale, angle_range[1] * scale)
        else:
            effective_range = angle_range
        angle = np.radians(random.uniform(*effective_range))
        R = GeometryOps.rotation_matrix_rodrigues(axis, angle)
        coords_centered = mol_copy.coordinates[submol_idx] - center
        coords_rotated = coords_centered @ R.T
        mol_copy.coordinates[submol_idx] = coords_rotated + center
        return self._apply_box_constraints(mol_copy, molecule)

    def diagnose_perturbations(
        self,
        molecule: "Molecule",
        submolecule_indices: List[List[int]],
        n_samples: int = 200,
        adaptive: bool = True
    ) -> Dict:
        """
        Diagnose the actual perturbation magnitudes produced by each local operator.
        
        Applies each operator n_samples times and reports statistics on:
          - RMSD between original and perturbed structures
          - Max single-atom displacement
          - COM displacement of the perturbed submolecule
          - Effective scaling factor and parameter ranges
          
        Args:
            molecule: Molecule to test operators on
            submolecule_indices: Submolecule index lists
            n_samples: Number of random applications per operator
            adaptive: Whether to use adaptive scaling
            
        Returns:
            Dict with operator names as keys, each containing displacement statistics
        """
        from scipy.spatial.distance import pdist
            
        # First report the raw scaling info
        scale = self._compute_scaling_factor(molecule, operator_type="local")
        coords = molecule.coordinates
        if len(coords) > 1:
            mol_diameter = np.max(pdist(coords))
        else:
            mol_diameter = 2.0 * molecule.covalent_radii.max()
        avg_radius = np.mean(molecule.covalent_radii)
        effective_mol_size = mol_diameter + 2 * avg_radius
        
        if self.simulation_box is not None:
            if self.simulation_box.box_type == "sphere":
                box_length = self.simulation_box.radius * 2
            else:
                box_length = self.simulation_box.box_dimensions
            free_space = max(box_length - effective_mol_size, 0.1)
            size_ratio = effective_mol_size / free_space
        else:
            box_length = None
            free_space = None
            size_ratio = None
            
        print("=" * 70)
        print("LOCAL OPERATOR PERTURBATION DIAGNOSTICS")
        print("=" * 70)
        print(f"  Molecule diameter:     {mol_diameter:.3f} Å")
        print(f"  Effective mol size:    {effective_mol_size:.3f} Å")
        if box_length is not None:
            print(f"  Box length (diameter): {box_length:.3f} Å")
            print(f"  Free space:            {free_space:.3f} Å")
            print(f"  Size ratio:            {size_ratio:.3f}")
        print(f"  Adaptive scale factor: {scale:.4f}")
        print(f"  Adaptive enabled:      {adaptive}")
        print()
        
        # Define operators with their default ranges and the effective ranges
        operators_info = {
            "random_displacement": {
                "func": self.random_displacement_submolecule,
                "default_range": (-0.3, 0.3),
                "unit": "Å",
                "range_type": "displacement"
            },
            "random_rotation": {
                "func": self.random_rotation_submolecule,
                "default_range": (-30, 30),
                "unit": "°",
                "range_type": "angle"
            },
            "local_twist": {
                "func": self.local_twist_submolecule,
                "default_range": (-15, 15),
                "unit": "°",
                "range_type": "angle"
            },
            "correlated_displacement": {
                "func": self.correlated_displacement,
                "default_range": (-0.3, 0.3),
                "unit": "Å",
                "range_type": "displacement"
            },
            "small_principal_axis_rotation": {
                "func": self.small_principal_axis_rotation,
                "default_range": (-15, 15),
                "unit": "°",
                "range_type": "angle"
            },
        }
        
        results = {}
        
        for op_name, info in operators_info.items():
            default_range = info["default_range"]
            if adaptive:
                effective_range = (default_range[0] * scale, default_range[1] * scale)
            else:
                effective_range = default_range
                
            rmsds = []
            max_displacements = []
            box_rejections = 0
            
            for _ in range(n_samples):
                original_coords = molecule.coordinates.copy()
                perturbed = info["func"](molecule, submolecule_indices, adaptive=adaptive)
                
                # Check if box rejected (returned original)
                if np.allclose(perturbed.coordinates, original_coords):
                    box_rejections += 1
                    continue
                
                diff = perturbed.coordinates - original_coords
                per_atom_disp = np.linalg.norm(diff, axis=1)
                rmsd = np.sqrt(np.mean(per_atom_disp**2))
                
                rmsds.append(rmsd)
                max_displacements.append(np.max(per_atom_disp))
            
            rmsds = np.array(rmsds) if rmsds else np.array([0.0])
            max_displacements = np.array(max_displacements) if max_displacements else np.array([0.0])
            
            results[op_name] = {
                "effective_range": effective_range,
                "rmsd_mean": np.mean(rmsds),
                "rmsd_std": np.std(rmsds),
                "max_disp_mean": np.mean(max_displacements),
                "max_disp_std": np.std(max_displacements),
                "box_rejection_rate": box_rejections / n_samples,
            }
            
            print(f"  {op_name}")
            print(f"    Default range:    ({default_range[0]}, {default_range[1]}) {info['unit']}")
            print(f"    Effective range:  ({effective_range[0]:.4f}, {effective_range[1]:.4f}) {info['unit']}")
            print(f"    RMSD:             {np.mean(rmsds):.4f} ± {np.std(rmsds):.4f} Å")
            print(f"    Max atom disp:    {np.mean(max_displacements):.4f} ± {np.std(max_displacements):.4f} Å")
            print(f"    Box rejections:   {box_rejections}/{n_samples} ({box_rejections/n_samples*100:.1f}%)")
            print()
        
        print("=" * 70)
        print("REFERENCE: Typical chemical scales")
        print(f"  O-H bond length:        0.96 Å")
        print(f"  H-bond (O···H):         1.8-2.0 Å")
        print(f"  O···O in water dimer:    2.9 Å")
        print(f"  Useful displacement:     0.05-0.5 Å")
        print(f"  Useful rotation:         5-45°")
        print("=" * 70)
        
        return results


class NonLocalOperators(MolecularOperators):
    """   
    Non-local operators for large structural perturbations of the molecule
    These operators make significant changes to explore different basins of the PES
    """

    def twist_operator(self,
                       molecule: "Molecule",
                       submolecule_indices: List[List[int]],
                       rotation_angle_range: Tuple[float, float] = (0, 360),
                       adaptive: bool = True) -> "Molecule":
        """
        Twist operator: Rotate one submolecule relative to others
        
        Args:
            molecule: Input Molecule object
            submolecule_indices: List of lists of atom indices defining submolecules
            rotation_angle_range: Range for random rotation angle in degrees
            adaptive: If True, scale angle range by box/molecule size ratio
        """ 
        mol_copy = molecule.copy()

        ref_submol, rot_submol = self._select_two_submolecules(submolecule_indices)

        # Calculate reference_frame
        ref_coords = mol_copy.coordinates[ref_submol]
        ref_masses = mol_copy.masses[ref_submol]
        ref_center = GeometryOps.center_of_mass(ref_coords, ref_masses)

        # Apply adaptive scaling
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="nonlocal")
            effective_range = (rotation_angle_range[0] * scale, rotation_angle_range[1] * scale)
        else:
            effective_range = rotation_angle_range

        angle = np.radians(random.uniform(*effective_range))
        axis = self._random_unit_vetor()
        R = GeometryOps.rotation_matrix_rodrigues(axis, angle)

        coords_centered = mol_copy.coordinates[rot_submol] - ref_center
        coords_rotated = coords_centered @ R.T
        mol_copy.coordinates[rot_submol] = coords_rotated + ref_center
        return self._apply_box_constraints(mol_copy, molecule)
    
    def large_displacement(self,
                           molecule: "Molecule",
                           submolecule_indices: List[List[int]],
                           displacement_range: Tuple[float, float] = (1.0, 3.0),
                           adaptive: bool = True) -> "Molecule":
        """   
        Large displacement operation

        Args:
            molecule: Input Molecule object
            submolecule_indices: List of lists of atom indices defining submolecules
            displacement_range: Range for random displacement in Angstrom (default is (1.0, 3.0))
            adaptive: If True, scale range by box/molecule size ratio
        """
        mol_copy = molecule.copy()

        # Select random submolecule
        submol_idx = self._select_random_submolecule(submolecule_indices)

        # Apply adaptive scaling
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="nonlocal")
            effective_range = (displacement_range[0] * scale, displacement_range[1] * scale)
        else:
            effective_range = displacement_range

        magnitude = random.uniform(*effective_range)
        direction = self._random_unit_vetor()
        displacement = direction * magnitude

        mol_copy.coordinates[submol_idx] += displacement
        return self._apply_box_constraints(mol_copy, molecule)
    
    def mirror_operator(self,
                        molecule: "Molecule",
                        submolecule_indices: List[List[int]]) -> "Molecule":
        """  
        Mirror operator: Reflect submolecule through a plane
        Uses reference frame construction and reflection matrix

        Note: This operator is not scaled as reflection is a discrete transformation

        Algorithm:
        1. Construct reference frame from the submolecule to mirror
        2. Select coordinate plane (XY, YZ or XZ) 
        3. Build Reflection matrix H = I - 2 * n * n^T 
        4. Apply reflection to the submolecule coordinates

        Args:
            molecule: Input Molecule object
            submolecule_indices: List of lists of atom indices defining submolecules
        """
        mol_copy = molecule.copy()

        ref_submol, mir_submol = self._select_two_submolecules(submolecule_indices)

        # Construct reference frame from the submolecule to mirror
        ref_frame = self.transformer.set_reference_frame_from_indices(
            mol_copy,
            atom_indices=ref_submol,
        )
        plane = np.random.choice(['XY', 'YZ', 'XZ'])
        if plane == 'XY':
            normal_vector = ref_frame.z_axis
        elif plane == 'YZ':
            normal_vector = ref_frame.x_axis 
        else: # XZ
            normal_vector = ref_frame.y_axis
        n = normal_vector / np.linalg.norm(normal_vector)
        H = np.eye(3) - 2 * np.outer(n, n)
        coords_centered = mol_copy.coordinates[mir_submol] - ref_frame.origin
        coords_reflected = coords_centered @ H.T
        mol_copy.coordinates[mir_submol] = coords_reflected + ref_frame.origin
        return self._apply_box_constraints(mol_copy, molecule)

    def random_so3_operator(self,
                            molecule: "Molecule",
                            submolecule_indices: List[List[int]]) -> "Molecule":
        """ 
        Applies a uniform Random SO(3) rotation to a submolecule around its own center of mass
        
        Note: This operator is not scaled as it samples uniformly over all rotations
        """
        mol_copy = molecule.copy()

        # Select submolecule
        submol_idx = self._select_random_submolecule(submolecule_indices)
        coords = mol_copy.coordinates[submol_idx]
        masses = mol_copy.masses[submol_idx]
        # Calc com
        com = GeometryOps.center_of_mass(coords, masses)
        # Sample Rotation Quaternion
        R = Quaternion.random_so3()
        coords_centered = coords - com
        coords_rotated = (R @ coords_centered.T).T

        mol_copy.coordinates[submol_idx] = coords_rotated + com
        return self._apply_box_constraints(mol_copy, molecule)

    def principal_axis_rotation_operator(self,
                                         molecule: "Molecule",
                                         submolecule_indices: List[List[int]],
                                         angle_range: Tuple[float, float] = (0, 360),
                                         adaptive: bool = True) -> "Molecule":
        """  
        Rotate a single submolecule around one of its own principal axes (from inertia tensor)

        Similar to random_so3 this operates on a single submolecule around its COM

        Algorithm:
        1. Select a random submolecule
        2. Compute its reference frame to get the principal axes
        3. Pick one random principal axes (x,y,z)
        4. Rotate the submolecule around that axis by a random angle

        Args:
            molecule: Input Molecule object
            submolecule_indices: List of lists of atom indices defining submolecules
            angle_range: Range for random rotation angle in degrees
            adaptive: If True, scale angle range by box/molecule size ratio
        """
        mol_copy = molecule.copy()
        submol_idx = self._select_random_submolecule(submolecule_indices)
        coords = mol_copy.coordinates[submol_idx]
        masses = mol_copy.masses[submol_idx]
        com = GeometryOps.center_of_mass(coords, masses)

        # Compute the reference frame to get the principal axes
        ref_frame = self.transformer.compute_reference_frame(coords, masses, center=com)
        
        # Pick a random principal axis
        axis = ref_frame.axes[:, np.random.choice(3)]
        
        # Apply adaptive scaling
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="nonlocal")
            effective_range = (angle_range[0] * scale, angle_range[1] * scale)
        else:
            effective_range = angle_range

        angle = np.radians(random.uniform(*effective_range))
        R = GeometryOps.rotation_matrix_rodrigues(axis, angle)

        coords_centered = coords - com 
        coords_rotated = coords_centered @ R.T
        mol_copy.coordinates[submol_idx] = coords_rotated + com
        
        return self._apply_box_constraints(mol_copy, molecule)
    




    def roto_reflection_operator(self,
                                 molecule: "Molecule",
                                 submolecule_indices: List[List[int]],
                                 angle_range: Tuple[float, float] = (0, 360),
                                 adaptive: bool = True) -> "Molecule":
        """  
        Roto-reflection operator: Combined rotation and reflection
        
        Args:
            molecule: Input Molecule object
            submolecule_indices: List of lists of atom indices defining submolecules
            angle_range: Range for random rotation angle in degrees
            adaptive: If True, scale angle range by box/molecule size ratio
        """
        mol_copy = molecule.copy()
        ref_submol, rot_ref_submol = self._select_two_submolecules(submolecule_indices)
        ref_frame = self.transformer.set_reference_frame_from_indices(
            mol_copy,
            atom_indices=ref_submol,
        )
        plane = np.random.choice(['XY', 'YZ', 'XZ'])
        if plane == 'XY':
            normal_vector = ref_frame.z_axis
        elif plane == 'YZ':
            normal_vector = ref_frame.x_axis 
        else: # XZ
            normal_vector = ref_frame.y_axis
        n = normal_vector / np.linalg.norm(normal_vector)
        H = np.eye(3) - 2 * np.outer(n, n)

        # Apply adaptive scaling
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="nonlocal")
            effective_range = (angle_range[0] * scale, angle_range[1] * scale)
        else:
            effective_range = angle_range

        angle = np.radians(random.uniform(*effective_range))
        R = GeometryOps.rotation_matrix_rodrigues(n, angle)
        M = R @ H
        coords_centered = mol_copy.coordinates[rot_ref_submol] - ref_frame.origin
        coords_transformed = coords_centered @ M.T
        mol_copy.coordinates[rot_ref_submol] = coords_transformed + ref_frame.origin
        return self._apply_box_constraints(mol_copy, molecule)
    
    def exchange_operator(self,
                          molecule: "Molecule",
                          submolecule_indices: List[List[int]]) -> "Molecule":
        """  
        Exchange Operator: Swaps the positions of two submolecules
        Uses vectorized center of mass calculation
        
        Note: This operator is not scaled as it's a discrete exchange operation
        """
        mol_copy = molecule.copy()
        submol_1, submol_2 = self._select_two_submolecules(submolecule_indices)
        coords_1 = mol_copy.coordinates[submol_1]
        masses_1 = mol_copy.masses[submol_1]
        center_1 = GeometryOps.center_of_mass(coords_1, masses_1)
        coords_2 = mol_copy.coordinates[submol_2]
        masses_2 = mol_copy.masses[submol_2]
        center_2 = GeometryOps.center_of_mass(coords_2, masses_2)
        translation_1 = center_2 - center_1
        translation_2 = center_1 - center_2
        mol_copy.coordinates[submol_1] += translation_1
        mol_copy.coordinates[submol_2] += translation_2
        return self._apply_box_constraints(mol_copy, molecule)


# ========================= Testing Utilities for Operators =========================


# ========================= Testing Utilities for Operators =========================

if __name__ == "__main__":

    from molecule_class import Molecule
    from box import SimulationBox

    # ===================== Setup: H2O Dimer =====================
    xyz = """
    6
    Coordinates from ORCA-job h2o_2 E -152.102897751726
    O           0.20131422818946     -0.13419863189991     -0.38118207664628
    H           1.10826926774619      0.03054527004232     -0.16898516572928
    H          -0.26386315635905     -0.06924940339370      0.43390806815981
    O           3.06513444516272      0.52224610599392      0.05902763481169
    H           3.25817767296947      1.29105665797100     -0.44996761253291
    H           3.66908154229119     -0.14039999871364     -0.23001884806303
    """
    molecule = Molecule.from_xyz(xyz)
    submol_indices = [[0, 1, 2], [3, 4, 5]]


    def animate_operator(molecule,submolecule_indices, operator_func, operator_name,
                         n_frames=60, n_applications=3, save_path=None, adaptive=None):
        
        """  
        Animate a nonlocal operator by interpolation between original and transformed coordinates

        Args:
            molecule: Input Molecule object
            submolecule_indices: List of lists of atom indices defining submolecules
            operator_func: Non-local operator function to apply
            operator_name: Name of the operator (for title)
            n_frames: Number of frames in the animation
            n_applications: Number of times to apply the operator for animation
            save_path: Optional path to save the animation as a GIF
            adaptive: Whether to use adaptive scaling for the operator
        """
        import matplotlib.animation as animation
        from pathlib import Path

        # Collect keyframes: original + each operator application
        keyframes = [molecule.copy()]
        current = molecule.copy()
        for _ in range(n_applications):
            if adaptive is not None:
                new = operator_func(current, submolecule_indices, adaptive=adaptive)
            else:
                new = operator_func(current, submolecule_indices)
            keyframes.append(new)
            current = new.copy()

        # Build interpolated frames between each pair of keyframes
        all_frames = []
        for k in range(len(keyframes) - 1):
            start_coords = keyframes[k].coordinates
            end_coords = keyframes[k + 1].coordinates
            for t in range(n_frames):
                frac = t / n_frames
                interp = start_coords + frac * (end_coords - start_coords)
                all_frames.append((interp.copy(), k+1))
            
            # Hold on final position
            for _ in range(n_frames // 3):
                all_frames.append((end_coords.copy(), k+1))
        
        # Precompute colors and sizes per submolecule
        sub_edge_colors = plt.cm.Set2(np.linspace(0, 1, max(len(submolecule_indices), 1)))
        atom_colors = molecule.get_atom_colors()
        atom_sizes = molecule.get_atom_sizes()

        # Compute fixed axis limits for all keyframes
        all_kf_coords = np.vstack([kf.coordinates for kf in keyframes])
        center = all_kf_coords.mean(axis=0)
        max_range = np.ptp(all_kf_coords, axis=0).max() / 2 + 1.5

        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111, projection='3d')

        def update(frame_idx):
            coords, app_num = all_frames[frame_idx]
            ax.cla()

            for sub_i, indices in enumerate(submolecule_indices):
                sub_coords = coords[indices]
                colors = [atom_colors[i] for i in indices]
                sizes = [atom_sizes[i] for i in indices]

                ax.scatter(
                    sub_coords[:, 0], sub_coords[:, 1], sub_coords[:, 2],
                    c=colors, s=sizes, edgecolors=sub_edge_colors[sub_i],
                    linewidths=2, alpha=0.9, depthshade=True
                )
                # Draw covalent bonds within submolecule
                for i in range(len(indices)):
                    for j in range(i + 1, len(indices)):
                        dist = np.linalg.norm(sub_coords[i] - sub_coords[j])
                        if dist < 1.3:
                            ax.plot(
                                [sub_coords[i, 0], sub_coords[j, 0]],
                                [sub_coords[i, 1], sub_coords[j, 1]],
                                [sub_coords[i, 2], sub_coords[j, 2]],
                                color='dimgray', linewidth=2, alpha=0.7
                            )

            ax.set_title(f"{operator_name}  (step {app_num}/{n_applications})",
                     fontsize=13, fontweight='bold')
            ax.set_xlabel('X (Å)')
            ax.set_ylabel('Y (Å)')
            ax.set_zlabel('Z (Å)')
            ax.set_xlim(center[0] - max_range, center[0] + max_range)
            ax.set_ylim(center[1] - max_range, center[1] + max_range)
            ax.set_zlim(center[2] - max_range, center[2] + max_range)
            ax.view_init(elev=20, azim=30 + frame_idx * 0.3)
            return []

        ani = animation.FuncAnimation(
            fig, update, frames=len(all_frames), interval=33, blit=False
        )

        if save_path:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            if save_path.endswith('.gif'):
                ani.save(save_path, writer='pillow', fps=30, dpi=150)
            else:
                ani.save(save_path, writer='ffmpeg', fps=30, dpi=150)
            print(f"Saved: {save_path}")
            plt.close()
        else:
            plt.show()





    covalent_radii = [0.66, 0.31, 0.31, 0.66, 0.31, 0.31]
    box = SimulationBox.from_covalent_radii(covalent_radii, 6, scale_factor=5.0)
    nonlocal_ops = NonLocalOperators(simulation_box=box)

    # ===================== Generate Animations =====================
    # (operator_name, function, adaptive_kwarg)
    #   adaptive=None means the operator doesn't take an adaptive parameter
    operators_to_animate = [
        ("Twist",                    nonlocal_ops.twist_operator,                      True),
        ("Large Displacement",       nonlocal_ops.large_displacement,                  True),
        ("Mirror",                   nonlocal_ops.mirror_operator,                     None),
        ("Random SO(3)",             nonlocal_ops.random_so3_operator,                 None),
        ("Principal Axis Rotation",  nonlocal_ops.principal_axis_rotation_operator,    True),
        ("Roto-Reflection",          nonlocal_ops.roto_reflection_operator,            True),
        ("Exchange",                 nonlocal_ops.exchange_operator,                   None),
    ]

    for op_name, op_func, adaptive in operators_to_animate:
        safe_name = op_name.lower().replace(" ", "_").replace("(", "").replace(")", "")
        save_path = f"figures/operator_{safe_name}.gif"
        print(f"Animating: {op_name} -> {save_path}")
        animate_operator(
            molecule, submol_indices, op_func, op_name,
            n_frames=100, n_applications=5,
            save_path=save_path,
            adaptive=adaptive,
        )

    print("All animations done.")

