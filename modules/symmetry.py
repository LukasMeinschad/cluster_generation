import numpy as np
from molecule_class import Molecule
from modules.logger import Logger



class MoleculeSymmetry:
    """
    Helper Class to determine the symmetry properties of a molecular system
    """
    def __init__(self, molecule: Molecule, logger: Logger):
        self.molecule = molecule
        self.logger = logger

    def rotor_classification(self):
        """
        Function that classifies the molecule as linear, symmetric, spherical or assymetric tops

        The definition of rigid rotor stems from classical mechanics. In classical mechanics the kinetic energy of rotation of a rigid rotor is linear 
        in a quanity called inertia tensot. This is a real symmetric 3x3 matrix which has three real eigenvalues (the principal moments of intertia)

        For linear rotos:

        + I_A << I_B = I_C (Normally I_A = 0)

        For symmetric tops:

        + I_A = I_B < I_C (frisbee or disc shape) or I_A < I_B = I_C (cigar shape)

        For spherical tops:

        Special case of symmetric tops with equal moment of inertia about all three axes I_A = I_B = I_C
        
        Special case 
        """

        # Get inertia tensor
        inertia_tensor = self._get_inertia_tensor()
        principal_moments, principal_axes = self._get_principal_moments(inertia_tensor)

        #print("Principal Moments of Inertia: ", principal_moments)


        if np.isclose(principal_moments[0], 0, atol=1e-3) and np.isclose(principal_moments[1], principal_moments[2], atol=1e-3):
            rotor_type = "Linear Rotor"
        elif np.isclose(principal_moments[0], principal_moments[1], atol=1e-3) and np.isclose(principal_moments[1], principal_moments[2], atol=1e-3):
            rotor_type = "Spherical Top"
        elif np.isclose(principal_moments[0], principal_moments[1], atol=1e-3) or np.isclose(principal_moments[1], principal_moments[2], atol=1e-3):
            rotor_type = "Symmetric Top"
        else:
            rotor_type = "Asymmetric Top"

        if self.logger:
            message = f"Molecule Symmetry Analysis:\n"
            message += f"Principal Moments of Inertia (amu·Å²): {principal_moments[0]:.4f}, {principal_moments[1]:.4f}, {principal_moments[2]:.4f}\n"
            message += f"Rotor Classification: {rotor_type}\n"
            self.logger.write_message_block(message)

    def _get_inertia_tensor(self):
        """
        Function that computes the inertia tensor for the given molecules
        """
        coords = self.molecule.coordinates
        m = self.molecule.masses

        x,y,z = coords[:,0], coords[:,1], coords[:,2]

        I = np.zeros((3,3))
        
        Ixx = np.sum(m * (y**2 + z**2))
        Iyy = np.sum(m * (x**2 + z**2))
        Izz = np.sum(m * (x**2 + y**2))
        Ixy = -np.sum(m * x * y)
        Ixz = -np.sum(m * x * z)
        Iyz = -np.sum(m * y * z)

        I = np.array([[Ixx, Ixy, Ixz],
                      [Ixy, Iyy, Iyz],
                      [Ixz, Iyz, Izz]])
        
        return I

    def _get_principal_moments(self,inertia_tensor: np.ndarray):
        """
        Function that computes the principal moments of inertia via diagonalization of the inertia tensor
        """
        eigenvalues, eigenvectors = np.linalg.eigh(inertia_tensor)  
        return eigenvalues, eigenvectors


    def distance_matrix(self):
        """
        Function that computes the distance matrix of a given molecule

        Example 

        H  O  H 

        Distance Matrix:
        [[0.     0.9572 0.9572]
         [0.9572 0.     1.5145]
         [0.9572 1.5145 0.    ]]
        """
        coords = self.molecule.coordinates
        num_atoms = coords.shape[0]
        dist_matrix = np.zeros((num_atoms, num_atoms))

        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                dist = np.linalg.norm(coords[i] - coords[j])
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist  # Symmetric matrix

        if self.logger:
            message = f"Molecule Distance Matrix:\n{dist_matrix}\n"
            self.logger.write_message_block(message)

        return dist_matrix


