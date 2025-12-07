import numpy as np
from typing import List, Optional, Tuple, Union
import matplotlib.pyplot as plt
from plotting import Plotting
from molecule_class import Molecule
import timeit
from transformations import Transformation


class ConfigSampler:
    """
    Simplified Class for sampling molecular configurations
    """
    def __init__(self,reference_frame: np.ndarray, center_point: np.ndarray, molecule, logger: Optional["Logger"] = None):
        """ 
        Initializes the ConfigSampler with a reference frame and center point

        Optionally one can provide a Logger object to log sampling information
        """

        self.reference_frame = reference_frame  
        self.center_point = center_point
        self.sampled_configurations = []

        # Optional logger
        self.logger = logger

        # Store the sampling region parameters for r
        self.sampling_region = { 
            "form": None, # "sphere", "cone", "rectangle" etc.
            "parameters": {} # Dictionary to hold relevant parameters
        }


    def sample_rectangle(self,
                         center_point: np.ndarray,
                         dir_vector1: np.ndarray,
                         dir_vector2: np.ndarray,
                         length1: float,
                         length2, float,
                         num_points: int = 100,
                         plot_sampling: bool = False) -> np.ndarray:
        """
        Sample a rectangular region defined by two directions vectors and the respective lenghts  

        length1: Half-length along dir_vector1, full length is 2*length1
        length2: Half-length along dir_vector2, full length is 2*length2
        """
        dir1 = dir_vector1 / np.linalg.norm(dir_vector1)
        dir2 = dir_vector2 / np.linalg.norm(dir_vector2)

        # Ensure that dir1 and dir2 are orthogonal
        dir2 = dir2 - np.dot(dir2, dir1) * dir1
        dir2 /= np.linalg.norm(dir2)

        # Generate random parameters in [-1,1] range
        u = np.random.uniform(-1, 1, num_points)
        v = np.random.uniform(-1, 1, num_points)

        # Scale parameters by the half lengths
        u_scaled = u * length1 
        v_scaled = v * length2

        # Calculate sampled points
        sampled_points = center_point + np.outer(u_scaled, dir1) + np.outer(v_scaled, dir2)

        if plot_sampling:
            raise NotImplementedError("Rectangle sampling plot not implemented yet")

        return sampled_points

    def sample_submol_rectangle(self,
                                submolecule: "SubMolecule",
                                center_point: np.ndarray,
                                dir_vector1: np.ndarray,
                                dir_vector2: np.ndarray,
                                length1: float,
                                length2: float,
                                num_points: int = 100,
                                rotation: bool = False,
                                rotational_grid_dim: int = 4,
                                plot_sampling: bool = False) -> List:
        """  
        Samples a submolecule within a rectangular volume, the other atoms in the molecule remain fixed
        """
        
        # Save parameters of sampling region
        self.sampling_region["form"] = "rectangle"
        self.sampling_region["parameters"] = {
            "center_point": center_point,
            "dir_vector1": dir_vector1,
            "dir_vector2": dir_vector2,
            "length1": length1,
            "length2": length2,
            "num_points": num_points
        }
        
        
        mol = submolecule.parent
        submol_indices = submolecule.get_index_in_parent()
        submol_center = np.mean(submolecule.coordinates, axis=0)
        # Sample points already at self.center_point
        sampled_points = self.sample_rectangle(center_point, dir_vector1, dir_vector2, length1, length2, num_points)
        # Create template molecule without coordinates (shared properties)
        template_mol = self._generate_template_molecule(mol)
        # Calculate translation vectors for submolecule
        translation_vectors = sampled_points - submol_center
        # Use optimized batch generation
        if rotation:
            # Generate rotational grid
            rotation_unit_vectors = self.generate_spherical_grid(grid_dim=rotational_grid_dim)
            # Compute reference frame once
            trans_ref_frame, submol_center = Transformation().set_reference_frame_submolecule(
                submolecule, mol, method="com", print_info=False, plot_frame=True
            )

            # Allocation of total list, much faster than appending
            total_mols = len(rotation_unit_vectors) * len(translation_vectors)
            sampled_mols = [None] * total_mols
            mol_counter = 0
            # Precompute base coordinates and mask
            base_coords = mol.coordinates.copy()
            mask = np.zeros(base_coords.shape[0], dtype=bool)
            mask[submol_indices] = True

            mol_counter = 0
            for rot_idx, rot_vec in enumerate(rotation_unit_vectors):
                # Rotate submolecule once per rotation vector
                rotated_coords = Transformation().rotate_molecule_to_point(
                    submolecule,
                    center_point=submol_center,
                    ref_frame=trans_ref_frame,
                    target_point=rot_vec + submol_center
                )

                # Apply all translations for this rotation
                for trans_idx, translation_vector in enumerate(translation_vectors):
                    new_mol = Molecule(name=f"rectangle_sample_rot{rot_idx}_trans{trans_idx}")
                    new_mol.atom_labels = template_mol.atom_labels
                    new_mol.masses = template_mol.masses
                    # Efficient coordinate assignment
                    new_coords = base_coords.copy()
                    new_coords[mask] = rotated_coords + translation_vector
                    new_mol.coordinates = new_coords
                    sampled_mols[mol_counter] = new_mol
                    mol_counter += 1
        else:
            sampled_mols = self._generate_molecules_batch_optimized(
                template_mol, mol.coordinates, submol_indices, translation_vectors
            )

        if plot_sampling:
            raise NotImplementedError("Rectangle submolecule sampling plot not implemented yet")
        
        return sampled_mols

    def sample_cone(self, 
                    apex: np.ndarray,
                    axis: np.ndarray,
                    height: float,
                    base_radius: float,
                    num_points: int = 100,
                    plot_sampling: bool = False) -> np.ndarray:
        """  
        Generates a uniform distribution of points within a cone volume.

        Point inside the cone can be written with cylindrical coordinates:

        + (x,y,z) = (r*cos(theta), h, rsin(theta))
        """
        apex = np.asarray(apex, dtype=float)
        axis = np.asarray(axis, dtype=float)
        axis = axis / np.linalg.norm(axis)  # Normalize the axis vector

        # Build orthonormal frame {axis_unit, u,v}
        # Pick vector not parallel to axis
        guess = np.array([0,0,1], dtype=float)
        if abs(np.dot(axis, guess)) > 0.9:
            guess = np.array([0,1,0], dtype=float)
        u = np.cross(axis, guess)
        u /= np.linalg.norm(u)
        v = np.cross(axis, u)

        # Cube-root distribution for height
        h = height * np.random.random(num_points) ** (1/3)

        # Radius at that height
        r_max = (base_radius / height) * h
        
        # Uniform distribution inside disk
        r = r_max * np.sqrt(np.random.random(num_points))
        theta = 2 * np.pi * np.random.random(num_points)
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)

        # Local coordinates in cone form
        pts_local = (apex + np.outer(h, axis) +
                     np.outer(r * cos_t, u) +
                     np.outer(r * sin_t, v))
        
        if plot_sampling:
            Plotting().plot_sampling_cone(pts_local, apex, axis, height, base_radius, u, v)
        return pts_local

    def sample_cylinder(self,
                        base_point: np.ndarray,
                        axis: np.ndarray,
                        height: float,
                        radius: float,
                        num_points: int = 100,
                        plot_sampling: bool = False) -> np.ndarray:
        """ 
        Uniform samples of points within a cylinder volume

        We take the base_point as the lower base center of the cylinder and sample points along the axis up to height
        """
        base_point = np.asarray(base_point, dtype=float)
        axis = np.asarray(axis, dtype=float)
        axis = axis / np.linalg.norm(axis)  # Normalize the axis vector

        # Build orthonormal frame {axis_unit, u,v}
        # Pick vector not parallel to axis
        guess = np.array([0,0,1], dtype=float)
        if abs(np.dot(axis, guess)) > 0.9:
            guess = np.array([0,1,0], dtype=float)
        u = np.cross(axis, guess)
        u /= np.linalg.norm(u)
        v = np.cross(axis, u)

        # Uniform height distribution
        h = height * np.random.random(num_points)
        # Uniform distribution inside disk
        r = radius * np.sqrt(np.random.random(num_points))
        theta = 2 * np.pi * np.random.random(num_points)
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)
        # Local coordinates in cylinder form
        pts_local = (base_point + np.outer(h, axis) +
                     np.outer(r * cos_t, u) +
                     np.outer(r * sin_t, v))
        
        if plot_sampling:
            Plotting().plot_sampling_cylinder(pts_local, base_point, axis, height, radius)
        return pts_local
    

    def sample_sphere(self,
                      radius: float,
                      num_points: int = 100,
                      plot_sampling: bool = False) -> np.ndarray:
        """ 
        Uniformly samples points inside a sphere volume

        For this we sample in spherical coordinates (with distribution correction for uniformity) and then convert to cartesian
        """

        # Precompute random points
        random_data = np.random.random((num_points, 3))

        phi = 2 * np.pi * random_data[:,0]  # azimuthal angle
        u = 2 * random_data[:,1] - 1  # cos(polar angle)
        r = radius *  np.cbrt(random_data[:,2])  # cube-root for uniform distribution in volume

        # Precompute sqrt(1-u^2) for sin(theta)
        sin_theta = np.sqrt(1 - u**2)

        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)

        x = r * sin_theta * cos_phi
        y = r * sin_theta * sin_phi
        z = r * u

        pts_local = np.column_stack((x, y, z)) + self.center_point


        if plot_sampling:
            raise NotImplementedError("3D sphere sampling plot not implemented yet")

        return pts_local 
    
    def test_sampling_speed(self, sampling_type: str):
        """  
        Helper function to test the speed of the different sampling functions implemented
        """
        if sampling_type == "sphere":
            pts = [1000, 10000, 100000,]
            for n in pts:
                start_time = timeit.default_timer()
                self.sample_sphere(radius=5.0, num_points=n)
                end_time = timeit.default_timer()
                print(f"Sphere sampling of {n} points took {end_time - start_time:.6f} seconds")

    def sample_plane_along_vector(self,
                                  center_point: np.ndarray,
                                  normal_vector: np.ndarray,
                                  direction_vector: np.ndarray,
                                  s_length: int,
                                  t_params: np.array = np.array([0,5]),
                                  num_points: int = 100,
                                  plot_sampling: bool = False
                                  ):
        """
        Function aims to sample points on a plane given by 

        + a center point
        + a normal vector
        + a direction vector inside the plane
        + a second direction vector orthogonal to the first one
        + length parameters are given in the form [s,t] which say how long we go along the direction vector

        Equation for the plane is given by:
        p(s,t) = center_point + s * dir_vector1 + t * dir_vector2
        """
        d = direction_vector / np.linalg.norm(direction_vector)
        n = normal_vector / np.linalg.norm(normal_vector)
        # Compute second in-plane direction vector
        d2 = np.cross(n, d)
        norm_d2 = np.linalg.norm(d2)
        if norm_d2 < 1e-6:
            raise ValueError("Direction vector cannot be parallel to normal vector")
        d2 /= norm_d2

        # Starting point is center_point and we sample along t and s
        s_values = np.random.uniform(0, s_length, num_points) 
        t_values = np.random.uniform(t_params[0], t_params[1], num_points)     


        sampled_points = center_point + np.outer(s_values, d) + np.outer(t_values, d2)

        if plot_sampling:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(sampled_points[:, 0], sampled_points[:, 1], sampled_points[:, 2], c='b', marker='o')
            ax.set_title('Plane Sampling')

            # Draw center point and vectors
            ax.quiver(center_point[0], center_point[1], center_point[2],
                      normal_vector[0], normal_vector[1], normal_vector[2],
                      length=1.0, color='r', label='Normal Vector')
            ax.quiver(center_point[0], center_point[1], center_point[2],
                      direction_vector[0], direction_vector[1], direction_vector[2],
                        length=1.0, color='g', label='Direction Vector')
            ax.quiver(center_point[0], center_point[1], center_point[2],
                      d2[0], d2[1], d2[2],
                      length=1.0, color='m', label='Second Direction Vector')
            
            
            
            plt.legend()
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.savefig('plane_sampling.png')
            plt.close()
        
        return sampled_points



    def generate_spherical_grid(self,grid_dim=4):
        """ 
        Generates a spherical grid of points which are later used for sampling the rotational space of a submolecule//molecule

        returns:
        A list of unit vectors representing points on the sphere
        """  
        phi = np.linspace(0, np.pi, grid_dim)
        theta = np.linspace(0, 2 * np.pi, grid_dim)
        phi_grid, theta_grid = np.meshgrid(phi, theta)
        x = np.sin(phi_grid) * np.cos(theta_grid)
        y = np.sin(phi_grid) * np.sin(theta_grid)
        z = np.cos(phi_grid)
        unit_vectors = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
        return unit_vectors
        
 

    def _generate_template_molecule(self, mol):
        """ 
        Creates a template molecule with shared properties for sampling
        """
        template_mol = Molecule(name="template_molecule")
        template_mol.atom_labels = mol.atom_labels.copy()
        template_mol.masses = mol.masses.copy()
        return template_mol


    def sample_submol_sphere(self, 
                         submolecule: "SubMolecule",
                         radius: float,
                         num_points: int = 100,
                         rotation: bool = True,
                         rotational_grid_dim: int = 4,
                         plot_sampling: bool = False) -> List:
        """ 
        Samples a submolecule within a spherical volume, the other atoms in the molecule remain fixed

        If rotation is True, for each translation a set of rotations will be generated
        """ 

        # TODO: Check this function again here

        # Save parameters of sampling region
        self.sampling_region["form"] = "sphere"
        self.sampling_region["parameters"] = {
            "radius": radius,
            "num_points": num_points
        }

    
        mol = submolecule.parent
        submol_indices = submolecule.get_index_in_parent()
        submol_center = np.mean(submolecule.coordinates, axis=0)

        # Sample points already at self.center_point 
        sampled_points = self.sample_sphere(radius, num_points)
    
        # Create template molecule without coordinates (shared properties)
        template_mol = self._generate_template_molecule(mol)

        # Calculate translation vectors for submolecule
        translation_vectors = sampled_points - submol_center

        if rotation:
            # Generate rotational grid
            rotation_unit_vectors = self.generate_spherical_grid(grid_dim=rotational_grid_dim)

            # Compute reference frame once
            trans_ref_frame, submol_center = Transformation().set_reference_frame_submolecule(
                submolecule, mol, method="com", print_info=False, plot_frame=True
            )

            # Allocation of total list, much faster than appending
            total_mols = len(rotation_unit_vectors) * len(translation_vectors)
            sampled_mols = [None] * total_mols

            # Precompute base coordinates and mask
            base_coords = mol.coordinates.copy()
            mask = np.zeros(base_coords.shape[0], dtype=bool)
            mask[submol_indices] = True

            mol_counter = 0
            for rot_idx, rot_vec in enumerate(rotation_unit_vectors):
                # Rotate submolecule once per rotation vector
                rotated_coords = Transformation().rotate_molecule_to_point(
                    submolecule,
                    center_point=submol_center,
                    ref_frame=trans_ref_frame,
                    target_point=rot_vec + submol_center
                )

                # Apply all translations for this rotation
                for trans_idx, translation_vector in enumerate(translation_vectors):
                    new_mol = Molecule(name=f"sphere_sample_rot{rot_idx}_trans{trans_idx}")
                    new_mol.atom_labels = template_mol.atom_labels
                    new_mol.masses = template_mol.masses

                    # Efficient coordinate assignment
                    new_coords = base_coords.copy()
                    new_coords[mask] = rotated_coords + translation_vector
                    new_mol.coordinates = new_coords

                    sampled_mols[mol_counter] = new_mol
                    mol_counter += 1
        else:
            # Use optimized batch generation
            sampled_mols = self._generate_molecules_batch_optimized(
                template_mol, mol.coordinates, submol_indices, translation_vectors
            )

        if plot_sampling:
            Plotting().plot_sampled_molecules(
                sampled_mols, sampling_type="Sphere", center_point=self.center_point, radius=radius
            )

        if self.logger is not None:
            message = f"Sampled {len(sampled_mols)} configurations of submolecule '{submolecule.name}' within sphere of radius {radius} Å, centered at {self.center_point}."
            if rotation:
                message += f" Each translation included {len(rotation_unit_vectors)} rotations."
            self.logger.write_message_block(message)
        
        return sampled_mols
    
    def sample_submol_cone(self,
                           submolecule: "SubMolecule",
                           apex: np.ndarray,
                            axis: np.ndarray,
                            height: float,
                            base_radius: float,
                            num_points: int = 100,
                            rotation: bool = True,
                            rotational_grid_dim: int = 4,
                            plot_sampling: bool = False) -> List:
        """ 
        Implementation of the submolecule sampling within a cone volume

        For this one needs to provide the apex point of the cone and the axis vector. If rotation is true,
        """
        mol = submolecule.parent
        submol_indices = submolecule.get_index_in_parent()
        submol_center = np.mean(submolecule.coordinates, axis=0)
        
        # Sample points already at self.center_point 
        sampled_points = self.sample_cone(apex, axis, height, base_radius, num_points)
        
        # Create template molecule without coordinates (shared properties)
        template_mol = self._generate_template_molecule(mol)
        
        # Calculate translation vectors for submolecule
        translation_vectors = sampled_points - submol_center

        if rotation:
            # Generate rotational grid
            rotation_unit_vectors = self.generate_spherical_grid(grid_dim=rotational_grid_dim)

            # Compute reference frame once
            trans_ref_frame, submol_center = Transformation().set_reference_frame_submolecule(
                submolecule, mol, method="com", print_info=False, plot_frame=True
            )
            # Allocation of total list, much faster than appending
            total_mols = len(rotation_unit_vectors) * len(translation_vectors)
            sampled_mols = [None] * total_mols
            mol_counter = 0
            # Precompute base coordinates and mask
            base_coords = mol.coordinates.copy()
            mask = np.zeros(base_coords.shape[0], dtype=bool)
            mask[submol_indices] = True
            for rot_idx, rot_vec in enumerate(rotation_unit_vectors):
                # Rotate submolecule once per rotation vector
                rotated_coords = Transformation().rotate_molecule_to_point(
                    submolecule,
                    center_point=submol_center,
                    ref_frame=trans_ref_frame,
                    target_point=rot_vec + submol_center
                )
                # Apply all translations for this rotation
                for trans_idx, translation_vector in enumerate(translation_vectors):
                    new_mol = Molecule(name=f"cone_sample_rot{rot_idx}_trans{trans_idx}")
                    new_mol.atom_labels = template_mol.atom_labels
                    new_mol.masses = template_mol.masses
                    
                    # Efficient coordinate assignment
                    new_coords = base_coords.copy()
                    new_coords[mask] = rotated_coords + translation_vector
                    new_mol.coordinates = new_coords
                    sampled_mols[mol_counter] = new_mol
                    mol_counter += 1
        else:
            # Use optimized batch generation
            sampled_mols = self._generate_molecules_batch_optimized(
                template_mol, mol.coordinates, submol_indices, translation_vectors
            )

        
        return sampled_mols 
            



    def _generate_molecules_batch_optimized(self, template, base_coords, submol_indices, translation_vectors):
        """ 
        Optimized helper function for batch molecule generation
        """
        num_points = len(translation_vectors)
    
        # Create mask for efficient coordinate updates
        mask = np.zeros(base_coords.shape[0], dtype=bool)
        mask[submol_indices] = True
    
        # Extract submolecule coordinates once
        submolecule_coords = base_coords[submol_indices]
    
        # Pre-allocate list
        sampled_mols = [None] * num_points
    
        for i in range(num_points):
            new_mol = Molecule(name=f"sphere_sample_{i}")

            # Copy shared properties (these are references, not deep copies)
            new_mol.atom_labels = template.atom_labels
            new_mol.masses = template.masses

            # Efficient coordinate assignment using the mask
            new_coords = base_coords.copy()
            new_coords[mask] = submolecule_coords + translation_vectors[i]
            new_mol.coordinates = new_coords

            sampled_mols[i] = new_mol
    
        return sampled_mols


    # TODO Delete From here

    def generate_spherical_grid_centers(self, center_point,ref_frame,radius,grid_dim=4, plot_sampling=False):
        """
        Generates points on a sphere surface using 3 x 3 grid in spherical coordinates.

        We then sample each mid point of the grid cell 
        """
        # Generate 3x3 grid
        phi_edges = np.linspace(0,np.pi,grid_dim)  # 4 edges = 3 cells
        theta_edges = np.linspace(0, 2*np.pi, grid_dim) # Edges = 3 cells

        # Calculate mid points of each cell
        phi_centers = (phi_edges[:-1] + phi_edges[1:]) / 2
        theta_centers = (theta_edges[:-1] + theta_edges[1:]) / 2

        # Create meshgrid of mid points
        phi_grid, theta_grid = np.meshgrid(phi_centers, theta_centers)

        # Convert to cartesian coordinates in local frame
        x_local = radius * np.sin(phi_grid) * np.cos(theta_grid)
        y_local = radius * np.sin(phi_grid) * np.sin(theta_grid)
        z_local = radius * np.cos(phi_grid)

        # Reshape into list of points and transform to global coordinates
        local_points = np.column_stack([x_local.ravel(), y_local.ravel(), z_local.ravel()])
        global_points = center_point + local_points @ ref_frame.T

        if plot_sampling:
            Plotting().plot_sampling_sphere(phi_grid.ravel(), theta_grid.ravel())
            Plotting().plot_sampling_sphere_3d(global_points, ref_frame=ref_frame)

        return global_points
    

    
    def sample_spherical_points(self,molecule,target_points,ref_frame,center_point,plot_sampling=False):
        """
        For a given list of target points on the sphere rotate the molecule and save all the configurations
        """
        sampled_mols = []

     #   # Calculate angle in original molecule
     #   for atom1,atom2,atom3 in combinations(molecule.coordinates,3):
     #       vec1 = atom2 - atom1
     #       vec2 = atom3 - atom2
     #       cos_angle = np.dot(vec1,vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
     #       angle = np.arccos(cos_angle) * 180 / np.pi
     #       print(f"Original Molecule: Angle between atoms: {angle:.4f} degrees")

        for i, point in enumerate(target_points):
            rotated_coords = Transformation().rotate_molecule_to_point(molecule.coordinates, center_point, ref_frame, point) 
            new_mol = Molecule(name=f"sphere_sample_{i}")
            new_mol.atom_labels = molecule.atom_labels.copy()
            new_mol.coordinates = rotated_coords
            new_mol.masses = molecule.masses.copy()
            sampled_mols.append(new_mol)
        self.sampled_molecules = sampled_mols

             

        if plot_sampling:
            Plotting().plot_spherical_mol_sampling(sampled_mols,ref_frame=ref_frame, target_points=target_points, center_point=center_point)
        return sampled_mols
        


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
            "corners": cube_cornes,
            "min_corner": np.min(cube_cornes, axis=0),
            "max_corner": np.max(cube_cornes, axis=0)
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

    def generate_rotation_grid(self,grid_dim=4):
        """
        Generates a grid of unit vectors on on the sphere using spherical coordinates, this are then used for rotations as xyz grid points
        """
        phi = np.linspace(0, np.pi, grid_dim)  # polar angle [0. pi]
        theta = np.linspace(0, 2 * np.pi, grid_dim)  # azimuthal angle [0, 2pi]
        phi_grid, theta_grid = np.meshgrid(phi, theta)
        x = np.sin(phi_grid) * np.cos(theta_grid)
        y = np.sin(phi_grid) * np.sin(theta_grid)
        z = np.cos(phi_grid)
        unit_vectors = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
        return unit_vectors


    def place_submolecule_uniformly_rotation(self,submolecules,submol_trans_id, submol_fixed_id, num_samples=10, volume_type="cube", max_attempts = None, min_distance=0.8, rot_grid=4, plot_samples=False):
        """  
        Places a submolecule uniformly within the sampling volume 

        For each position within the sampling volume, a set of rotations will be generated
        """
        if volume_type not in self.sampling_volumes:
            raise ValueError(f"Volume type '{volume_type}' not found in sampling volumes")

        if max_attempts is None:
            max_attempts = num_samples * 10

        volume = self.sampling_volumes[volume_type]

        submol_trans = submolecules[submol_trans_id] # Submolecule to be translated
        submol_fixed = submolecules[submol_fixed_id] # Submolecule to remain fixed
        sampled_mols = []
        sucessful_samples = 0
        attempts = 0
        while sucessful_samples < num_samples and attempts < max_attempts:
            attempts += 1

            # Generate a random position
            if volume_type == "cube":
                position = self._get_uniform_position_cube(volume)

            else:
                raise NotImplementedError("Only 'cube' volume type is implemented")
            
            # Check minimum distance
            if not self._check_minimum_distance(position, submol_trans, submol_fixed, min_distance=min_distance):
                print(f"Sample {attempts}: Minimum distance not respected, trying again.")
                continue

            # Generate rotation grid
            rotation_vectors = self.generate_rotation_grid(grid_dim=rot_grid)
            for rot_vec in rotation_vectors:
                # Get Reference Frame of submolecule
                trans_ref_frame, trans_center = Transformation().set_reference_frame_submolecule(submolecules, submol_trans_id, self.molecule, method="com", print_info=False, plot_frame=False)
                rotated_submol = Molecule(name=f"{submol_trans.name}_rotated_sample_{sucessful_samples}")
                rotated_submol.atom_labels = submol_trans.atom_labels.copy()
                rotated_submol.masses = submol_trans.masses.copy()
                rotated_submol.coordinates = submol_trans.coordinates.copy()
                # Rotate submolecule to rotation vector
                rotated_coords = Transformation().rotate_molecule_to_point(rotated_submol.coordinates,center_point= trans_center, ref_frame=trans_ref_frame, target_point= rot_vec +  trans_center)
                rotated_submol.coordinates = rotated_coords
                # Create new molecule with placed submolecule
                new_mol = self._create_placed_molecule(rotated_submol, submol_fixed, position, sucessful_samples)

                sampled_mols.append(new_mol)
            sucessful_samples += 1

            print(f"Sample {attempts}: Successfully placed submolecule at {position} with {len(rotation_vectors)} rotations")   

        if sucessful_samples < num_samples:
            print(f"Warning: Only {sucessful_samples} samples were generated out of requested {num_samples} after {attempts} attempts.")

        if plot_samples:
            cube_corners = volume.get("corners", None)
            Plotting().plot_uniform_sampling(sampled_mols, cube_corners=cube_corners)
        return sampled_mols 




    def place_submolecule_uniformly(self,submolecules,submol_trans_id, submol_fixed_id,num_samples=10, volume_type="cube",max_attempts=None, min_distance=0.8, rotation=True,plot_samples=False):
        """ 
        Places a submolecule uniformly within the sampling volume
        
        Args:
            submolecules: List of submolecule objects to be placed
            submol_trans: Submolecule to be translated
            submol_fixed: Submolecule to remain fixed
            num_samples: Number of samples to generate
            volume_type: Type of sampling volume ("cube" or "sphere")

        If rotation is true, for each position a set of rotations will be generated
        based on a spherical grid with dimension rot_grid x rot_grid
        """
        if volume_type not in self.sampling_volumes:
            raise ValueError(f"Volume type '{volume_type}' not found in sampling volumes")

        if max_attempts is None:
            max_attempts = num_samples * 10

        volume = self.sampling_volumes[volume_type]
        submol_trans = submolecules[submol_trans_id] # Submolecule to be translated
        submol_fixed = submolecules[submol_fixed_id] # Submolecule to remain fixed

        sampled_mols = []
        sucessful_samples = 0
        attempts = 0



        while sucessful_samples < num_samples and attempts < max_attempts:
            attempts += 1

            # Generate a random position
            if volume_type == "cube":
                position = self._get_uniform_position_cube(volume)
            else:
                raise NotImplementedError("Only 'cube' volume type is implemented")
            
            # Check minimum distance
            if not self._check_minimum_distance(position, submol_trans, submol_fixed, min_distance=min_distance):
                print(f"Sample {attempts}: Minimum distance not respected, trying again.")
                continue

            # Create new molecule with placed submolecule
            new_mol = self._create_placed_molecule(submol_trans, submol_fixed, position, sucessful_samples)
            sampled_mols.append(new_mol)
            sucessful_samples += 1
            print(f"Sample {attempts}: Successfully placed submolecule at {position}")
        
        if sucessful_samples < num_samples:
            print(f"Warning: Only {sucessful_samples} samples were generated out of requested {num_samples} after {attempts} attempts.")
        
        if plot_samples:
            cube_corners = volume.get("corners", None)
            Plotting().plot_uniform_sampling(sampled_mols, cube_corners=cube_corners)

        return sampled_mols
    


    def _get_uniform_position_cube(self,volume):
        """ 
        Uniformly samples a position within a cube volume
        """
        if volume["type"] != "cube":
            raise ValueError("Volume type must be 'cube' for this method")
        
        min_corner = volume["min_corner"]
        max_corner = volume["max_corner"]


        in_cube_volume=False
        while not in_cube_volume:
            rand_x = np.random.uniform(min_corner[0], max_corner[0])
            rand_y = np.random.uniform(min_corner[1], max_corner[1])
            rand_z = np.random.uniform(min_corner[2], max_corner[2])
            position = np.array([rand_x, rand_y, rand_z])
            in_cube_volume = self._check_position_in_cube(position, volume)
        return position

    def _check_position_in_cube(self,position,volume):
        """ 
        Checks if a given position is inside a cube volume
        """
        if volume["type"] != "cube":
            raise ValueError("Volume type must be 'cube' for this method")
        min_corner = volume["min_corner"]
        max_corner = volume["max_corner"]
        return np.all(position >= min_corner) and np.all(position <= max_corner)

    def _check_minimum_distance(self,position,submol_trans,submol_fixed,min_distance=0.8):
        """ 
        Checks if the minimum distance between two submolecules is respected
        """ 
        # TODO Implement also Center of Mass here
        trans_center = np.mean(submol_trans.coordinates, axis=0)
        translation_vector =  position - trans_center
        translated_coords = submol_trans.coordinates + translation_vector

        for fixed_atom in submol_fixed.coordinates:
            for trans_atom in translated_coords:
                distance = np.linalg.norm(fixed_atom - trans_atom)
                if distance < min_distance:
                    return False
                
        return True
    
    def _create_placed_molecule(self,submol_trans,submol_fixed,position,sample_id):
        """ 
        Creates a new Molecule object with the placed submolecule
        """ 
        submol_copy = Molecule(name=f"{submol_trans.name}_sample_{sample_id}")
        submol_copy.atom_labels = submol_trans.atom_labels.copy()
        submol_copy.coordinates = submol_trans.coordinates.copy()
        submol_copy.masses = submol_trans.masses.copy()

        # Place the copy at the target position
        self._place_submolecule_at_position(submol_copy, position)

        # Create new combined molecule
        new_mol = Molecule(name=f"sampled_molecule_{sample_id}")
        new_mol.add_atoms_batch(submol_fixed.atom_labels, submol_fixed.coordinates)
        new_mol.add_atoms_batch(submol_copy.atom_labels, submol_copy.coordinates)
        return new_mol
    
    def _place_submolecule_at_position(self,submolecule, position):
        """
        Places the submolecule at the given position
        """
        # TODO Implement center of mass here as well
        submol_center = np.mean(submolecule.coordinates, axis=0)
        translation_vector = position - submol_center
        submolecule.coordinates += translation_vector