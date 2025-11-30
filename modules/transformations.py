import numpy as np
from mendeleev.fetch import fetch_table
import timeit
from plotting import Plotting
from molecule_class import Molecule
from itertools import combinations
from typing import List, Tuple, Optional, Union


class Transformation:
    """
    Class that performs transformation of molecules in R^3
    """
    def __init__(self, name: str = "Transformation", logger: Optional["Logger"] = None):
        """ 
        Initializes the Transformation class

        Optionally one can provide a Logger object to log information about the transformations
        """
        self.name = None
        self.logger = logger


    def center_point(self,molecule: Union["Molecule", "Submolecule"], method: str = "centroid") -> np.ndarray:
        """ 
        Determine the center point for a molecule or a submolecule
        """

        coords = molecule.coordinates

        if method == "centroid":
            center_point = np.mean(coords,axis=0)
        elif method == "com":
            center_point = self.com(molecule)
        else:
            raise ValueError("Method must be 'centroid' or 'com'")

        return center_point
    
    def center_point_coords(self,coords: np.ndarray, masses: Optional[np.ndarray] = None, method: str = "centroid") -> np.ndarray:
        """ 
        Determine the center point for a set of coordinates 
        """
        
        if method == "centroid":
            center_point = np.mean(coords,axis=0)
        elif method == "com":
            if masses is None:
                raise ValueError("Masses must be provided to calculate center of mass")
            center_point = np.sum(masses[:, np.newaxis] * coords, axis=0) / np.sum(masses)
        else:
            raise ValueError("Method must be 'centroid' or 'com'")

        return center_point
    
    def com(self,molecule: Union["Molecule", "Submolecule"]) -> np.ndarray:
        """ 
         Calculates the center of mass for a given molecule or submolecule
        """
        coords = molecule.coordinates
        masses = molecule.masses

        if masses is None:
            raise ValueError("Molecule must have masses defined to calculate center of mass")
         
        
        com = np.sum(masses[:, np.newaxis] * coords, axis=0) / np.sum(masses)
        return com

    def translate(self,molecule: Union["Molecule", "Submolecule"],vector: Union[List[float], np.ndarray], print_info: bool = False) -> None:
        """
        Translates a molecule with a vector T_v(p) = p + v

        Note that the coordinates of the molecule are modified in place
        """

        vec = np.asarray(vector, dtype=float)
        if print_info:
            print("Performing translation of the molecule")
            print(f"Initial coordinates of the molecule: {molecule.coordinates}")
            print(f"Translation vector: {vec}")
        
        molecule.coordinates = molecule.coordinates  + vec

        if print_info:
            print("New coordinates of the molecule")
            print(molecule.coordinates)



    def test_translation_speed(self,molecule, num_trials=1000):
        """ 
        Function to test the speed of the translation function

        For this num_trials random translations are performed and the average and total time is reported
        """            

        total_time = 0.0
        for _ in range(num_trials):
            random_vector = np.random.rand(3)  # Random translation vector
            start_time = timeit.default_timer()
            self.translate(molecule, random_vector, print_info=False)
            end_time = timeit.default_timer()
            total_time += (end_time - start_time)
        
        average_time = total_time / num_trials
        print(f"Total time for {num_trials} translations: {total_time:.6f} seconds")
        print(f"Average time per translation: {average_time:.6f} seconds")

    def rodrigues_rotation(self,molecule: Union["Molecule", "SubMolecule"],
                           axis: Union[List[float]],
                           angle_degree: float,
                           print_info: bool =False) -> None:
        """ 
        Performs a rodrigues rotation of a molecule around a given axis and a given angle

        If v is a vector in R^3 and k is a unit vector and θ is an angle in degrees, the rotation of v around k by θ is given by:

        v_rot = v cos(θ) + (k × v) sin(θ) + k (k · v)(1 - cos(θ)
        
        As a rotation matrix R, this can be written as:

        K = |  0   -kz   ky |
            |  kz   0   -kx |
            | -ky   kx   0  |

        R = I + sin(θ)K + (1 - cos(θ))K^2
        """
        theta = np.radians(angle_degree)
        k = axis / np.linalg.norm(axis)  # Ensure k is a unit vector
        kx,ky,kz = k
        # Build the K matrix
        K = np.array([[0, -kz, ky],
                      [kz, 0, -kx],
                      [-ky, kx, 0]])
        # Build the Rotation Matrix R
        R = np.eye(3) + np.sin(theta) * K + (1- np.cos(theta)) * np.dot(K,K)

        if print_info:
            print(f"Performing Rodrigues rotation")
            print(f"Rotation axis (unit vector): {k}")
            print(f"Rotation angle (degrees): {angle_degree}")
            print(f"Rotation matrix R:\n{R}")
            print("Initial coordinates of the molecule:")
            print(molecule.coordinates)
        # Apply rotation to all coordinates
        molecule.coordinates = np.dot(molecule.coordinates, R.T)

    def test_rodrigues_speed(self,molecule,num_trials=1000):
        """  
        Function to test the speed of the rodrigues rotation function

        For this num_trials random rotations are performed, aka we draw random axis and angle
        """
        total_time = 0.0
        for _ in range(num_trials):
            random_axis = np.random.rand(3)  # Random rotation axis
            random_angle = np.random.uniform(0, 360)  # Random rotation angle in degrees
            start_time = timeit.default_timer()
            self.rodrigues_rotation(molecule, random_axis, random_angle, print_info=False)
            end_time = timeit.default_timer()
            total_time += (end_time - start_time)
        average_time = total_time / num_trials
        print(f"Total time for {num_trials} rotations: {total_time:.6f} seconds")
        print(f"Average time per rotation: {average_time:.6f} seconds")


    def set_reference_frame(self,
                            molecule: Union["Molecule", "SubMolecule"],
                            method: str = "centroid",
                            print_info: bool = False,
                            plot_frame: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """ 
        Sets the reference frame to a molecules center and principal axes

        Returns:
            Tuple of (reference_frame, center_point)
        """
        if len(molecule.coordinates) < 3:
            raise ValueError("Molecule must have at least 3 atoms to define a reference frame")
        
        # Calculate the center and translate molecule
        center_point = self.center_point(molecule, method=method)
        centered_coords = molecule.coordinates - center_point 
        # Update molecule coordinates
        molecule.coordinates = centered_coords

        # Calculate inertia tensor and principal axes
        I = self.inertia_tensor(centered_coords,molecule.masses)
        eigval_I, eigvec_I = np.linalg.eigh(I)

        ref_frame = np.array([eigvec_I[:,0], eigvec_I[:,1], eigvec_I[:,2]]).T

        if print_info:
            print(f"Molecule center: {center_point}")
            print(f"Inertia tensor eigenvalues: {eigval_I}")
            print(f"Reference frame axes:")
            print(f"X: {ref_frame[:,0]}")
            print(f"Y: {ref_frame[:,1]}") 
            print(f"Z: {ref_frame[:,2]}")
        if plot_frame:
            Plotting().plot_reference_frame(molecule, reference_frame=ref_frame, reference_frame_origin=np.zeros(3))
        
        if self.logger is not None:
            message = f"Set reference frame based on molecule with center at {center_point} using method '{method}'"
            message += f"\nInertia tensor eigenvalues: {eigval_I}"
            message += f"\nReference frame axes:\nX: {ref_frame[:,0]}\nY: {ref_frame[:,1]}\nZ: {ref_frame[:,2]}"
            self.logger.write_message_block(message)

        
        return ref_frame, center_point
    

    def set_reference_frame_labels(self,
                                   molecule: Union["Molecule", "SubMolecule"],
                                   atom_labels=List[int],
                                   method: str = "centroid",
                                   print_info: bool = False,
                                   plot_frame: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """ 
        Sets the reference frame to a molecules center and principal axes based on given atom labels
        """ 
        if len(atom_labels) < 3:
            raise ValueError("At least 3 atom labels must be provided to define a reference frame")
        # Extract coordinates of the specified atom labels
        label_coords = []
        label_masses = []
        for label in atom_labels:
            coord = molecule.get_coords_by_label(label)
            label_coords.append(coord)
            label_masses.append(molecule.get_mass_by_label(label))
        label_coords = np.array(label_coords)
        label_masses = np.array(label_masses)
        if method == "centroid":
            center_point = self.center_point_coords(label_coords, method=method)
        elif method == "com":
            center_point = self.center_point_coords(label_coords, masses=label_masses, method=method)
        
        centered_coords = label_coords - center_point
        # Update molecule coordinates
        molecule.coordinates = molecule.coordinates - center_point
        # Calculate inertia tensor and principal axes
        I = self.inertia_tensor(centered_coords,label_masses)
        eigval_I, eigvec_I = np.linalg.eigh(I)
        ref_frame = np.array([eigvec_I[:,0], eigvec_I[:,1], eigvec_I[:,2]]).T
        if print_info:
            print(f"Molecule center (based on labels {atom_labels}): {center_point}")
            print(f"Inertia tensor eigenvalues: {eigval_I}")
            print(f"Reference frame axes:")
            print(f"X: {ref_frame[:,0]}")
            print(f"Y: {ref_frame[:,1]}") 
            print(f"Z: {ref_frame[:,2]}")
        if plot_frame:
            Plotting().plot_reference_frame(molecule, reference_frame=ref_frame, reference_frame_origin=np.zeros(3))
        return ref_frame, center_point

        

    @classmethod
    def inertia_tensor(cls, centered_coords,masses):
        """ 
        Calculates the Inertia Tensor for a given molecule

        I = sum_i m_i r_i^2 where r_i i s the dinstance of each atom to the center of mass
        """
        I = np.zeros((3,3))
        x,y,z = centered_coords[:,0], centered_coords[:,1], centered_coords[:,2]
        m = masses

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
    
    def set_reference_frame_submolecule(self,
                                        submolecule: "SubMolecule",
                                        parent_molecule: Optional["Molecule"] = None,
                                        method: str = "centroid",
                                        print_info: bool = False,
                                        plot_frame: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        
        """ 
        Defines a reference frame based on a submolecule's geometry

        This centers both the submolecule and the parent molecule to the submolecule center.
        Further it defines axes based on the principal axes of the submolecule

        Args:
            submolecules: List of submolecule objects
            submol_index: Index of the submolecule to define the reference frame
            parent_molecule: The full parent molecule object
            method: Method to calculate center point ("centroid" or "com")
            print_info: Whether to print information about the reference frame
            plot_frame: Whether to plot the reference frame
        """ 
        
        if len(submolecule.coordinates) < 3:
            raise ValueError("Submolecule must have at least 3 atoms to define a reference frame")
        
        # Calculate the center of the submolecule
        center_point = self.center_point(submolecule, method=method)

        if print_info:
            print(f"Submolecule center: {center_point}")
        
        # Center the submolecule coordinates
        submolecule.coordinates = submolecule.coordinates - center_point

        # If parent molecule is provided, center its coordinates as well
        if parent_molecule is not None:
            if print_info:
                print("Centering parent molecule coordinates to submolecule center")
            parent_molecule.coordinates = parent_molecule.coordinates - center_point

        # Calculate inertia tensor and principal axes of the submolecule
        I = self.inertia_tensor(submolecule.coordinates, submolecule.masses)
        eigval_I, eigvec_I = np.linalg.eigh(I)
        ref_frame = np.array([eigvec_I[:,0], eigvec_I[:,1], eigvec_I[:,2]]).T

        if print_info:
            print(f"Submolecule inertia tensor eigenvalues: {eigval_I}")
            print(f"Reference frame axes:")
            print(f"X: {ref_frame[:,0]}")
            print(f"Y: {ref_frame[:,1]}") 
            print(f"Z: {ref_frame[:,2]}")
        if plot_frame:
            Plotting().plot_reference_frame(submolecule, reference_frame=ref_frame, reference_frame_origin=np.zeros(3))

        if self.logger is not None:
            message = f"Set reference frame based on submolecule with center at {center_point} using method '{method}'"
            message += f"\nInertia tensor eigenvalues: {eigval_I}"
            message += f"\nReference frame axes:\nX: {ref_frame[:,0]}\nY: {ref_frame[:,1]}\nZ: {ref_frame[:,2]}"
            self.logger.write_message_block(message)
        return ref_frame, center_point
    

    
    def rotate_molecule_to_point(self,
                                 molecule: Union["Molecule", "SubMolecule"],
                                 center_point: np.ndarray,
                                 ref_frame: np.ndarray,
                                 target_point: np.ndarray) -> np.ndarray:
        """ 
        Rotate molecule so that its reference z-axis points towards the target point
        """

        target_direction = target_point - center_point 
        target_direction = target_direction / np.linalg.norm(target_direction)

        # Align molecules z-axis to target direction
        current_z_axis = ref_frame[:,2]
        R = self.align_vectors(current_z_axis, target_direction)
        rotated_coords = np.dot(molecule.coordinates, R.T)
        
        return rotated_coords

    def hbond_configuration_ref_frame(self,
                                      configuration: "SubMolecule",
                                      parent_molecule: Optional["Molecule"] = None,
                                      method: str = "centroid",
                                      ref_type: str = "Inertia",
                                      print_info: bool = False,
                                      plot_frame: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """ 
        Defines a local reference frame based on a hydrogen bond donor configuration.

        The method can be "centroid", "com" or "donor_atom" to define the center point.

        The ref_type can be "Inertia" to use the principal axes of the configuration or "bond" to align the z-axis along the hydrogen bond donor bond vector.
        """
        if len(configuration.coordinates) < 2:
            raise ValueError("Configuration must have at least 2 atoms to define a reference frame")
        # TODO implement from here


    @staticmethod
    def align_vectors(a: np.ndarray,b: np.ndarray) -> np.ndarray:
        """ 
        Construct rotation matrix that alings vector a to vector b
        """


        a = a / np.linalg.norm(a)
        b = b / np.linalg.norm(b)

        # Handle the parallel case
        if np.allclose(a, b):
            return np.eye(3)
        
        if np.allclose(a, -b):
            # 180 Degree rotationn
            if abs(a[0]) > 1e-8 or abs(a[1]) > 1e-8:
                v = np.array([-a[1], a[0], 0])
            else:
                v = np.array([0, -a[2], a[1]])
            v = v / np.linalg.norm(v)
            return 2 * np.outer(v, v) - np.eye(3)
        

        v = np.cross(a,b)
        s = np.linalg.norm(v) 
        c = np.dot(a,b) 

        
        v_x = np.array([
            [0, -v[2], v[1]],
            [v[2], 0, -v[0]],
            [-v[1], v[0], 0]
        ])
        I = np.eye(3)

        R = I + v_x + np.dot(v_x, v_x) * ((1 - c) / (s**2))
        return R
    
    @staticmethod
    def axis_along_bond(molecule: "Molecule", atom_label1: int, atom_label2: int) -> np.ndarray:
        """ 
        Returns the unit vector along the bond defined by atom_index1 and atom_index2
        """
        coord1 = molecule.get_coords_by_label(atom_label1)
        coord2 = molecule.get_coords_by_label(atom_label2)

        bond_vector = coord2 - coord1
        bond_unit_vector = bond_vector / np.linalg.norm(bond_vector)

        return bond_unit_vector
    
    def translate_only_submolecule(self,
                                   submolecule: "SubMolecule",
                                   parent_molecule: Optional["Molecule"],
                                   target_point: np.ndarray) -> None:
        """ 
        Translates only the submolecule to the target point
        """
        submolecule_center = self.center_point(submolecule, method="centroid")
        translation_vector = target_point - submolecule_center
        submolecule.coordinates += translation_vector
        
        submolecule_labels = submolecule.get_atom_labels()
        if parent_molecule is not None:
            for i, label in enumerate(parent_molecule.atom_labels):
                if label in submolecule_labels:
                    parent_molecule.coordinates[i] += translation_vector
        
