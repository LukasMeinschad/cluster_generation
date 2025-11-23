import numpy as np
from mendeleev.fetch import fetch_table
import timeit
from plotting import Plotting
from molecule_class import Molecule

class Transformation:
    """
    Class that performs transformation of molecules in R^3
    """
    def __init__(self):
        self.name = None
    def center_point(self,molecule,method="centroid"):
        """ 
        Function that determines a center in space 
        """
        coords = molecule.coordinates

        if method == "centroid":
            # Calculates th e mean point in space
            center_point = np.mean(coords,axis=0)

        elif method == "com":
            # If the mass of the molecule is empty run the respective function
            # Calculates the center of mass
            center_point = self.com(molecule)

        else:
            raise ValueError("Method must be 'centroid' or 'com'")

        return center_point
    
    
            

    def com(self,molecule):
        """ 
        Calculates the center of mass for a given molecule
        """
        coords = molecule.coordinates
        masses = molecule.masses 
        
        com = np.sum(masses[:, np.newaxis] * coords, axis=0) / np.sum(masses)
        return com

    def translate(self,molecule,vector, print_info=False):
        """
        Makes a translation by a vector T_v(p) = p + v 
        """
        vec = np.asarray(vector, dtype=float)
        if print_info:
            print("Performing translation operation")
            print("Initial coordinates of the molecule")
            print(molecule.coordinates)
            print(f"Translation vector: {vec}")

        molecule.coordinates = molecule.coordinates  + vec
        
        if print_info:
            print("New coordinates of the molecule")
            print(molecule.coordinates)

    def center_to_submol(self,submolecule, molecule,method="centroid",print_info=False):
        """ 
        Function that centers the molecule with respect to a given molecule
        
        Args:
            submolecule: Submolecule object to center to
            molecule: Molecule object to be centered
            method: Method to determine the center ("centroid" or "com")
            print: If True, prints the transformation details
        """
        center_point = None

        print(f"Centering molecule to submolecule with labels  {submolecule.atom_labels}")
        print(f"Using method: {method}")
        
        if method == "centroid":
            center_point = self.center_point(submolecule, method="centroid")
        elif method == "com":
            center_point = self.center_point(submolecule, method="com")
        else:
            raise ValueError("Method must be 'centroid' or 'com'")

        # Translate Molecule to origin  
        translation_vector = -center_point
        if print_info:
            print(f"Centering molecule to submolecule using {method} method")
            print(f"Center point of submolecule: {center_point}")
            print(f"Translation vector: {translation_vector}")
        self.translate(molecule, translation_vector, print_info=print_info)

    def rodrigues_rotation(self,molecule,axis,angle_degree,print_info=False):
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

    def set_reference_frame_molecule(self,molecule,method="centroid",print_info=False,plot_frame=False):
        """ 
        Sets the reference frame to a molecules center of mass or centroid

        1. Translate the molecule so that it is centered at the origin
        2. Define orthogonal, x,y,z axes based on atoms
        3. Returns transformed coordinates and the new reference frame
        """

        if len(molecule.coordinates) < 3:
            raise ValueError("Molecule must have at least 3 atoms to define a reference frame")
        
        # Calculate the new center and translate molecules
        if method=="centroid":
            center_point = self.center_point(molecule,method="centroid")
        elif method == "com":
            center_point = self.center_point(molecule,method="com")

        centered_coords = molecule.coordinates - center_point
        I = self.inertia_tensor(centered_coords,molecule.masses)
        eigval_I, eigvec_I = np.linalg.eigh(I)
        # Obtain the reference frame from this diagonalization

        x_axis = eigvec_I[:,0]
        y_axis = eigvec_I[:,1]
        z_axis = eigvec_I[:,2]

        ref_frame = np.array([x_axis,y_axis,z_axis]).T
        if plot_frame:
            Plotting().plot_reference_frame(molecule, reference_frame=ref_frame, reference_frame_origin=center_point)


    @classmethod
    def inertia_tensor(cls, centered_coords,masses):
        """ 
        Calculates the Inertia Tensor for a given molecule

        I = sum_i m_i r_i^2 where r_i i s the dinstance of each atom to the center of mass
        """
        I = np.zeros((3,3))
        for i, (x,y,z) in enumerate(centered_coords):
            m = masses[i]
            I[0,0] += m * (y**2 + z**2) # Ixx
            I[1,1] += m * (x**2 + z**2) # Iyy
            I[2,2] += m * (x**2 + y**2) # Ixx

            # Negative off diagonal elements
            I[0,1] -= m * x * y # Ixy
            I[0,2] -= m * x * z # Ixz
            I[1,2] -= m * y * z # Iyz

        #TODO this motherfucker needs forever
        I = np.where(I,I,I.T) 
        return I
        
    
    def set_reference_frame_submolecule(self,submolecules, submol_index,molecule, method = "centroid", print_info=False, plot_frame=False):
        """ 
        For a given submolecule define a reference frame on this submolecule using the inertia tensor

        This centers the reference frame to the submolecule and defines axes based on the submolecules geometry 
        """ 
         
        # Select the submolecule
        submolecule = submolecules[submol_index]
        
        if len(submolecule.coordinates) < 3:
            raise ValueError("Submolecule needs at least 3 atoms to define a reference frame")



        # First calculate the center of the submolecule with one of the methods
        if method == "centroid":
            center_point = self.center_point(submolecule,method="centroid")
        if method == "com":
            center_point = self.center_point(submolecule,method="com")

        # Center the submolecule coordinates
        centered_subcoords = submolecule.coordinates - center_point

        # Apply the transformation to all submolecules
        for submol in submolecules:
            submol.coordinates = submol.coordinates - center_point

        # Calculate the inertia tensor for the submolecule
        I = self.inertia_tensor(centered_subcoords, submolecule.masses)
        eigval_I, eigvec_I = np.linalg.eigh(I)

        # Obtain the reference frame from the diagonalization
        x_axis = eigvec_I[:,0]
        y_axis = eigvec_I[:,1]
        z_axis = eigvec_I[:,2]

        ref_frame = np.array([x_axis,y_axis,z_axis]).T

        # We can now transform to the new reference frame
        centered_full_coords = molecule.coordinates - center_point

        molecule.coordinates = centered_full_coords # Overwrite the molecule coords

        if print_info:
            print(f"Submolecule center: {center_point}")
            print(f"Inertia tensor eigenvalues: {eigval_I}")
            print(f"Reference frame axes:")
            print(f"X: {x_axis}")
            print(f"Y: {y_axis}") 
            print(f"Z: {z_axis}")
        if plot_frame:
            Plotting().plot_reference_frame(molecule,reference_frame=ref_frame, reference_frame_origin=np.zeros(3))

        return ref_frame, center_point
    
class ConfigSampler:
    """
    This class is used to sample configurations of molecules within a given reference frame
    """
    def __init__(self,molecule,reference_frame, center_point):
        self.molecule = molecule
        self.reference_frame = reference_frame
        self.center_point = center_point 
        self.x_axis = reference_frame[:,0]
        self.y_axis = reference_frame[:,1]
        self.z_axis = reference_frame[:,2]

        # sample storage of molecules
        self.sampled_molecules = []
        # This sampling volume will store information about the geometric object used for sampling
        self.sampling_volumes = {}

    def create_cube_volume(self, size=5, plot_cube=False):
        """ 
        Function that takes the reference frame as the middle of a cube and creates 
        a cube volume around it  
        """
        half_size = size / 2
        corners = []
        for dx in [-half_size, half_size]:
            for dy in [-half_size, half_size]:
                for dz in [-half_size, half_size]:
                    corner = (self.center_point + 
                              dx * self.x_axis + 
                              dy * self.y_axis + 
                              dz * self.z_axis)
                    corners.append(corner)
        
        cube_cornes = np.array(corners)
        self.sampling_volumes["cube"] = {
            "type": "cube",
            "size": size,
            "axes": (self.x_axis, self.y_axis, self.z_axis),
            "corners": cube_cornes
        }

        if plot_cube:
            #Plotting().plot_cube(corners)  --- IGNORE ---
            Plotting().plot_cube_molecule(self.molecule, np.array(corners), ref_frame=self.reference_frame) 
        return cube_cornes

    def create_spherical_volume(self, radius=5):
        """
        Function that creates a spherical volume around the center point
        """
        # This function can be implemented in the future
        pass

    def place_submolecule_uniformly(self,submolecules,submol_trans_id, submol_fixed_id,num_samples=10, volume_type="cube", plot_samples=False):
        """ 
        Places a submolecule uniformly within the sampling volume
        
        Args:
            submolecules: List of submolecule objects to be placed
            submol_trans: Submolecule to be translated
            submol_fixed: Submolecule to remain fixed
            num_samples: Number of samples to generate
            volume_type: Type of sampling volume ("cube" or "sphere")
        """
        volume = self.sampling_volumes.get(volume_type)

        submol_trans = submolecules[submol_trans_id]
        submol_fixed = submolecules[submol_fixed_id]

        sampled_mols = []

        # Workflow:
        # For each sample for the pick uniformly a random position in the volume
        # Then translate the submolecule to this position
        # Then make a new molecule object that contains both submolecules
        for i in range(num_samples):
            # TODO this is not optimal maybe there is a better option so work with this sampling
            random_pos = self.random_position_cube(volume)
            # Create a copy of the submolecule to be translated
            submol_copy = Molecule(name=f"{submol_trans.name}_sample_{i+1}")
            submol_copy.atom_labels = submol_trans.atom_labels.copy()
            submol_copy.coordinates = submol_trans.coordinates.copy()
            submol_copy.masses = submol_trans.masses.copy()
            self.place_submolecule_position(submol_copy, random_pos)
            # Create a new molecule object that contains both submolecules
            new_mol = Molecule(name=f"sampled_molecule_{i+1}")
            new_mol.add_atoms_batch(submol_fixed.atom_labels, submol_fixed.coordinates, submol_fixed.masses)
            new_mol.add_atoms_batch(submol_copy.atom_labels, submol_copy.coordinates, submol_copy.masses)
            sampled_mols.append(new_mol)

        print(sampled_mols)

        if plot_samples:
            cube_cornes = volume["corners"]
            Plotting().plot_uniform_sampling(sampled_mols,cube_cornes)

        return sampled_mols



    def random_position_cube(self, cube_volume):
        """ 
        Generates a random position within the cube volume
        
        For this we obtain the min and max corners of the cube and take for each axis a random value between these two
        """
        corners = cube_volume["corners"]
        min_corner = np.min(corners, axis=0)
        max_corner = np.max(corners, axis=0)

        random_x = np.random.uniform(min_corner[0], max_corner[0])
        random_y = np.random.uniform(min_corner[1], max_corner[1])
        random_z = np.random.uniform(min_corner[2], max_corner[2])

        return np.array([random_x, random_y, random_z])
    
    def place_submolecule_position(self,submolecule,position):
        """
        Places a submolecule at a given position in space
        """
        # TODO maybe use com here
        translation_vector = position - np.mean(submolecule.coordinates, axis=0)
        submolecule.coordinates += translation_vector





    
        

        
        

    
        



        

        

        
        
