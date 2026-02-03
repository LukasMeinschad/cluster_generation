import numpy as np 
from typing import List, Tuple, Optional, Union 
from dataclasses import dataclass

@dataclass
class CoordinateSpace:
    """   
    Represents a partitioned coordinate space with internal + external coordinates.
    """
    U_ext: np.ndarray  # External coordinate basis vectors
    U_int: np.ndarray  # Internal coordinate basis vectors
    n_atoms: int   # Number of atoms in the molecule
    n_ext: int    # Number of external DOFs
    n_int: int    # Number of internal DOFs
    linear: bool # Indicates if the molecule is linear
    mass_weighted: bool # Indicates if the coordinates are mass-weighted

    @property
    def P_int(self) -> np.ndarray:
        """   
        Projector onto the internal coordinate space = I - U_ext U_ext^T
        """
        return np.eye(3*self.n_atoms) - self.U_ext @ self.U_ext.T
    
    def P_ext(self) -> np.ndarray:
        """   
        Projector onto external space = U_ext U_ext^T 
        """
        return self.U_ext @ self.U_ext.T
    
class CoordinateProjector:
    """
    Projects Cartesian forces and Hessians to internal coordinate space.  
    """

    def __init__(self, tolerance: float = 1e-10):
        """
        Tolerance is numerical tolerance for orthogonalization
        """
        self.tol = tolerance

    # ================= Mass-Weighting ================

    def to_mass_weighted(
            self,
            coords: np.ndarray,
            masses: np.ndarray
    ) -> np.ndarray:
        """   
        Converts Cartesian coordinates to mass-weighted coordinates.
        """ 
        sqrt_masses = np.sqrt(masses)
        return coords * sqrt_masses[:, np.newaxis]
    
    def gradient_to_mass_weighted(
        self,
        gradient: np.ndarray,
        masses: np.ndarray
    ):
        """   
        Converts Cartesian gradient to mass-weighted gradient

           F_q = M^(-1/2) * F_x
        """ 
        original_shape = gradient.shape
        gradient_flat = gradient.flatten()

        inv_sqrt_masses = 1.0 / np.sqrt(masses)
        mass_matrix = np.repeat(inv_sqrt_masses, 3)
        gradient_mw = gradient_flat * mass_matrix
        return gradient_mw.reshape(original_shape)
    
    def hessian_to_mass_weighted(
            self,
            hessian: np.ndarray,
            masses: np.ndarray,
    ) -> np.ndarray:
        """   
        Convert Cartesian Hessian to Mass-weighted

        H_q = M^(-1/2) * H_x * M^(-1/2)
        """
        inv_sqrt_masses = 1.0 / np.sqrt(masses)
        mass_matrix = np.repeat(inv_sqrt_masses, 3)
        hessian_mw = hessian * mass_matrix[:, np.newaxis] * mass_matrix[np.newaxis, :]
        return hessian_mw
    
    # ================= Translational Basis =================
    def build_translational_basis(
            self,
            n_atoms: int,
            masses: np.ndarray
    ) -> np.ndarray:
        """    
        Constructs 3 orthonormal translation vectors in mass-weighted space
        t_x = (sqrt(m1), 0, 0, sqrt(m2), 0, 0, ..., sqrt(mN), 0, 0) and so on
        """
        n_coords = 3 * n_atoms
        T = np.zeros((n_coords, 3))
        sqrt_masses = np.sqrt(masses)
        for i in range(n_atoms):
            idx = 3 * i
            T[idx:idx+3,:] = np.eye(3) * sqrt_masses[i]
        
        # Normalize Each collumn
        for col in range(3):
            norm = np.linalg.norm(T[:,col])
            T[:,col] /= norm
        return T

    # ================= Build Rotational Basis =================

    def build_rotation_basis(
        self,
        coords: np.ndarray,
        masses: np.ndarray,
        linear: bool = False
    ):
        """   
        Constructs the 2 or 3 orthonormal rotation vectors

        First we compute the center of mass and shift the coordinates
        Then for rotation about axis w, displacement of atom i
        delta r_i = w x r_I 

        Mass weighted delta q_i = sqrt(m_i) * (w x r_i)
        """
        n_atoms = len(coords)
        n_coords = 3 * n_atoms 

        # Rotation axes
        axes = np.eye(3)

        n_rot = 2 if linear else 3

        if linear:
            raise NotImplementedError("Linear molecule rotation basis not implemented yet.")
        
        R = np.zeros((n_coords, n_rot))
        sqrt_masses = np.sqrt(masses)
        for rot_idx in range(n_rot):
            omega = axes[rot_idx]
            for atom_idx in range(n_atoms):
                r_i = coords[atom_idx]
                delta_r = np.cross(omega, r_i)
                delta_q = sqrt_masses[atom_idx] * delta_r
                coord_idx = 3 * atom_idx
                R[coord_idx:coord_idx+3, rot_idx] = delta_q

        return R
    
    # ================= Build the External Basis =================
    def build_external_basis(
            self,
            coords: np.ndarray,
            masses: np.ndarray,
            linear: bool = False 
    ) -> np.ndarray:
        """    
        Constructs complete orthonormal external basis (translations + rotations)
        """
        n_atoms = len(coords)

        # Center at com
        com = np.sum(coords * masses[:, np.newaxis], axis=0) / np.sum(masses)
        shifted_coords = coords - com[np.newaxis, :]

        T = self.build_translational_basis(n_atoms, masses)
        R = self.build_rotation_basis(shifted_coords, masses, linear=linear)

        U_ext = np.hstack((T, R))

        # Gram Schmit
        U_ext, _ = np.linalg.qr(U_ext)
        return U_ext


    # =============== Build Internal Basis ==================
    def build_internal_basis(
            self,
            U_ext: np.ndarray,
            n_coords: int
    ) -> np.ndarray:
        """   
        Constructs the orthonormal internal basis (complemen of external basis)

        Uses SVD to find the null space of U_ext^T
        """
        U,S,Vt = np.linalg.svd(U_ext, full_matrices=True)
        # Columns beyond rank span null space
        rank = np.sum(S > self.tol)
        U_int = U[:, rank:]
        return U_int

    # =============== Set up Coordinate Space ==================

    def setup_coordinate_space(
        self,
        coords: np.ndarray,
        masses: np.ndarray,
        linear: bool = False,
    ) -> CoordinateSpace:
        """   
        Complete Setup of the external and internal coordinate spaces
        """
        n_atoms = len(coords)
        n_coords = 3 * n_atoms
        n_ext = 5 if linear else 6
        n_int = n_coords - n_ext

        U_ext = self.build_external_basis(coords, masses, linear=linear)
        U_int = self.build_internal_basis(U_ext, n_coords)

        coord_space = CoordinateSpace(
            U_ext=U_ext,
            U_int=U_int,
            n_atoms=n_atoms,
            n_ext=n_ext,
            n_int=n_int,
            linear=linear,
            mass_weighted=True
        )
        return coord_space

    # ============== Project Forces and Hessians =================
    def project_gradient(
            self,
            gradient: np.ndarray,
            coord_space: CoordinateSpace,
            full_space: bool = True
    ) -> np.ndarray:
        """   
        Projects gradient to internal coordinates

        Option 1 (full_space=True) --> F_int = P_int * F_q (stays 3N)
        Option 2 (full_space=False) --> F_int = U_int^T * F_q (reduces to 3N-6)
        """
        original_shape = gradient.shape
        grad_flat = gradient.flatten()
        if full_space:
            F_int = coord_space.P_int @ grad_flat
        else:
            F_int = coord_space.U_int.T @ grad_flat
        return F_int.reshape(original_shape) if full_space else F_int
    
    def project_hessian(
            self,
            hessian: np.ndarray,
            coord_space: CoordinateSpace,
            full_space: bool = True
    ) -> np.ndarray:
        """   
        Projects Hessian to internal coordinates

        Option 1 (full_space=True) --> H_int = P_int * H_q * P_int (stays 3N x 3N)
        Option 2 (full_space=False) --> H_int = U_int^T * H_q * U_int (reduces to (3N-6) x (3N-6))
        """
        if full_space:
            H_int = coord_space.P_int @ hessian @ coord_space.P_int
        else:
            H_int = coord_space.U_int.T @ hessian @ coord_space.U_int
        return H_int