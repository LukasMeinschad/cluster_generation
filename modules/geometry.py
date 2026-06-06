import numpy as np
from dataclasses import dataclass
from typing import Optional, List, Tuple 
from scipy.optimize import linear_sum_assignment
import random
from numba import njit, prange

@dataclass 
class ReferenceFrame:
    """  
    A class representing a possible reference frame for a given system
    """
    axes: np.ndarray # 3x3 array where each column is an axis
    origin: np.ndarray # 3 element array for origin point
    eigenvalues: Optional[np.ndarray] = None # Optional eigenvalues associated with the reference frame axes

    @property 
    def x_axis(self) -> np.ndarray:
        return self.axes[:,0]
    
    @property 
    def y_axis(self) -> np.ndarray:
        return self.axes[:,1]
    
    @property 
    def z_axis(self) -> np.ndarray:
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

    @staticmethod
    def from_rotation_matrix(R: np.ndarray) -> "Quaternion":
        """  
        Converts a rotation matrix to quaternion
        uses the method of the homepage  https://danceswithcode.net/engineeringnotes/quaternions/quaternions.html

        Is the numerically stable Shepperd Method

        Step 1: Compute the magnitude of each quaternion component. This leaves the sign of each component undetermined
        Step 2: Find the largest component, assume its sign is positive,
                then compute the remaining components from the off diagonal elements
                Taking the largest magnitude avoids division by small numbers
        """
        # Step 1 compute the magnitude from diagonal elements
        q0_sq = (1 + R[0,0] + R[1,1] + R[2,2]) / 4.0
        q1_sq = (1 + R[0,0] - R[1,1] - R[2,2]) / 4.0
        q2_sq = (1 - R[0,0] + R[1,1] - R[2,2]) / 4.0 
        q3_sq = (1 - R[0,0] - R[1,1] + R[2,2]) / 4.0

        # Pick largest to avoid division by small numbers
        idx = np.argmax([q0_sq, q1_sq, q2_sq, q3_sq])
        if idx == 0:
            q0 = np.sqrt(max(q0_sq, 0))
            q1 = (R[2,1] - R[1,2]) / (4 * q0)
            q2 = (R[0,2] - R[2,0]) / (4 * q0)
            q3 = (R[1,0] - R[0,1]) / (4 * q0)
        elif idx == 1:
            q1 = np.sqrt(max(q1_sq, 0))
            q0 = (R[2,1] - R[1,2]) / (4 * q1)
            q2 = (R[0,1] + R[1,0]) / (4 * q1)
            q3 = (R[0,2] + R[2,0]) / (4 * q1)
        elif idx == 2:
            q2 = np.sqrt(max(q2_sq, 0))
            q0 = (R[0,2] - R[2,0]) / (4 * q2)
            q1 = (R[0,1] + R[1,0]) / (4 * q2)
            q3 = (R[1,2] + R[2,1]) / (4 * q2)
        else:
            q3 = np.sqrt(max(q3_sq, 0))
            q0 = (R[1,0] - R[0,1]) / (4 * q3)
            q1 = (R[0,2] + R[2,0]) / (4 * q3)
            q2 = (R[1,2] + R[2,1]) / (4 * q3)
        return Quaternion(q0, q1, q2, q3).normalize()

    @staticmethod
    def is_valid_rotation_matrix(R: np.ndarray) -> bool:
        """  
        Checks if the given matrix is truly a valid rotation
        we check 
        1. Orthonormality: R^T R = I
        2. Proper Rotation: det(R) = +1
        """
        should_be_identity = R.T @ R
        I = np.eye(3)
        return np.allclose(should_be_identity, I) and np.isclose(np.linalg.det(R), 1.0)

    def slerp(q1: "Quaternion", q2: "Quaternion", t: float) -> "Quaternion":
        """ 
        Spherical Linear Interpolation (SLERP) between two quaternions.

        Given two unit quaternions q1 and q2 and a parameter t in [0,1]:
            slerp(q_1,q_2,t) = sin((1-t)*theta) / sin(theta) * q_1 + sin(t*theta) / sin(theta) * q_2

        Properties:
            - slerp(q1, q2, 0) = q1
            - slerp(q1, q2, 1) = q2
            - Constant angular velocity along shortest arc on S3

        Args:
            q1: Starting quaternion (unit)
            q2: Ending quaternion (unit)
            t: Interpolation parameter in [0,1]
        """
        # Dot product of quaternion components
        dot = q1.w * q2.w + q1.x * q2.x + q1.y * q2.y + q1.z * q2.z

        # Ensure shortest path: if dot < 0, negate q2
        if dot < 0.0:
            q2 = Quaternion(-q2.w, -q2.x, -q2.y, -q2.z)
            dot = -dot

        # If quaternions are very close, fall back to normalized linear interpolation
        if dot > 0.9995:
            return Quaternion(
                q1.w + t * (q2.w - q1.w),
                q1.x + t * (q2.x - q1.x),
                q1.y + t * (q2.y - q1.y),
                q1.z + t * (q2.z - q1.z),
            ).normalize()

        theta = np.arccos(np.clip(dot, -1.0, 1.0))
        sin_theta = np.sin(theta)

        a = np.sin((1 - t) * theta) / sin_theta
        b = np.sin(t * theta) / sin_theta

        return Quaternion(
            a * q1.w + b * q2.w,
            a * q1.x + b * q2.x,
            a * q1.y + b * q2.y,
            a * q1.z + b * q2.z,
        ).normalize()
    
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
    @njit(fastmath=True)
    def inertia_tensor(coords: np.ndarray, masses: np.ndarray) -> np.ndarray:
        """    
        Calculates the moment of inertia tensor for a set of coordinates with given masses
        
        I_ij = sum_k m_k (r_k^2 delta_ij - r_ki r_kj)
        """ 
        Ixx = Iyy = Izz = Ixy = Ixz = Iyz = 0.0
        for i in range(coords.shape[0]):
            x, y, z = coords[i, 0], coords[i, 1], coords[i, 2]
            m = masses[i]

            x2,y2,z2 = x*x, y*y, z*z
            Ixx += m * (y2 + z2)
            Iyy += m * (x2 + z2)
            Izz += m * (x2 + y2)
            Ixy -= m * x * y
            Ixz -= m * x * z
            Iyz -= m * y * z

        return np.array([[Ixx, Ixy, Ixz],
                         [Ixy, Iyy, Iyz],
                         [Ixz, Iyz, Izz]])

    @staticmethod
    def compute_optimal_correspondence_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
        """
        Computes the optimal-correspondence RMSD between two structures:
        1. Center both (remove translation)
        2. Hungarian assignment on centered coords (remove labelling ambiguity)
        3. Kabsch rotation (remove rotation)
        4. Compute RMSD on the aligned, reordered coords
        """
        P_c = P - P.mean(axis=0)
        Q_c = Q - Q.mean(axis=0)
        Q_reordered = GeometryOps.find_optimal_correspondence(P_c, Q_c)
        R = GeometryOps.kabsch_rotation(P_c, Q_reordered)
        Q_aligned = Q_reordered @ R.T
        diff = P_c - Q_aligned
        return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

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
    @njit(fastmath=True, parallel=True)
    def _fast_cost_matrix(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
        """  
        Computes the cost matrix using numba
        """
        n = P.shape[0]
        m = Q.shape[0]
        cost = np.zeros((n,m), dtype=np.float64)
        for i in prange(n):
            for j in range(m):
                dx = P[i,0] - Q[j,0]
                dy = P[i,1] - Q[j,1]
                dz = P[i,2] - Q[j,2]
                cost[i,j] = dx*dx + dy*dy + dz*dz
        return cost

    @staticmethod
    def find_optimal_correspondence(P: np.ndarray, Q: np.ndarray):
        """ 
        Finds the optimal correspondence between two sets of coordinates P and Q using the
        Hungarian algorithm to solve the linear assignment problem.

        Args:
            P: Target coordinates (Nx3)
            Q: Source coordinates (Nx3)  
        """
        # Compute cost matrix as squared distances
        cost_matrix = GeometryOps._fast_cost_matrix(P, Q)
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        return Q[col_ind]

    
    @staticmethod 
    def pairwise_distance_matrix(coords: np.ndarray) -> np.ndarray:
        """  
        Computes the pairwise distance matrix for a set of coordinates
        D_ij = || r_i - r_j ||
        """
        n = coords.shape[0]
        dist_matrix = np.zeros((n,n), dtype=np.float64)
        for i in range(n):
            for j in range(i+1, n):
                dx = coords[i,0] - coords[j,0]
                dy = coords[i,1] - coords[j,1]
                dz = coords[i,2] - coords[j,2]
                dist = np.sqrt(dx*dx + dy*dy + dz*dz)
                dist_matrix[i,j] = dist
                dist_matrix[j,i] = dist
        return dist_matrix


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
    @njit(fastmath=True)
    def rotation_matrix_rodrigues(axis: np.ndarray, angle_rad: float) -> np.ndarray: 
        """   
        Construct the Rotation Matrix using Rodrigues' rotation formula

        R = I + sin(θ)K + (1 - cos(θ))K^2
        where K is the cross-product matrix of the unit axis vector
        """
        norm = np.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
        kx, ky, kz = axis[0] / norm, axis[1] / norm, axis[2] / norm
        

        sin_theta = np.sin(angle_rad)
        cos_theta = np.cos(angle_rad)
        one_minus_cos = 1 - cos_theta

        R = np.empty((3,3), dtype=np.float64)
        R[0,0] = cos_theta + kx*kx*one_minus_cos
        R[0,1] = kx*ky*one_minus_cos - kz*sin_theta
        R[0,2] = kx*kz*one_minus_cos + ky*sin_theta

        R[1,0] = ky*kx*one_minus_cos + kz*sin_theta
        R[1,1] = cos_theta + ky*ky*one_minus_cos
        R[1,2] = ky*kz*one_minus_cos - kx*sin_theta

        R[2,0] = kz*kx*one_minus_cos - ky*sin_theta
        R[2,1] = kz*ky*one_minus_cos + kx*sin_theta
        R[2,2] = cos_theta + kz*kz*one_minus_cos

        return R

    
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

    @staticmethod
    def distance(coords1: np.ndarray, coords2: np.ndarray) -> float:
        """  
        Computes the Euclidean distance between two points
        """
        return np.linalg.norm(coords1 - coords2)
    
  
