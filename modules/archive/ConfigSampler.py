import numpy as np 
from typing import List, Optional, Tuple, Union, Protocol, Dict 
import matplotlib.pyplot as plt 
from dataclasses import dataclass
import time
import random

from cluster_generation.modules.archive.plotting import Plotting
from molecule_class import Molecule 
from symmetry import SymmetryAnalyzer   
from transformations import Transformation, ReferenceFrame, Quaternion, GeometryOps
from scipy.spatial.distance import cdist




    


class GeometricSampler: 
    """  
    Class that is used for geometric Sampling operations
    """
    @staticmethod 
    def sample_sphere( 
        center: np.ndarray, 
        radius: float,
        num_points: int
        ) -> np.ndarray:
        """
        Uniformly samples points inside a sphere volume
        """  
        # Precompute random points
        random_data = np.random.random((num_points, 3)) 

        phi = 2 * np.pi * random_data[:,0] # azimuthal angle
        u = 2 * random_data[:,1] - 1 # cos(polar angle)
        r = radius * np.cbrt(random_data[:,2]) # cube-root for uniform distribution in volume

        sin_theta = np.sqrt(1 - u**2)

        x = r * sin_theta * np.cos(phi)
        y = r * sin_theta * np.sin(phi)
        z = r * u

        

        return np.column_stack((x,y,z)) + center
    
    @staticmethod
    def sample_cone( 
        apex: np.ndarray,
        axis: np.ndarray, 
        height: float,
        base_radius: float,
        num_points: int
        ) -> np.ndarray: 
        """  
        Uniform sampling inside of a cone volume
        
        Args:
            apex (np.ndarray): Apex point of the cone
            axis (np.ndarray): Axis vector of the cone
            height (float): Height of the cone
            base_radius (float): Base radius of the cone
            num_points (int): Number of points to sample
        """
        apex = np.asarray(apex, dtype=float)
        axis = np.asarray(axis, dtype=float)
        axis = axis / np.linalg.norm(axis)  # Normalize the axis vector

        # Construct Orthonormal Basis
        u, v = GeometricSampler._build_orthonormal_basis(axis)

        # Cube-root distribution for height
        h = height * np.cbrt(np.random.random(num_points))

        # Radius at that height
        r_max = (base_radius / height) * h

        # Uniform distribution inside Disk
        r = r_max * np.sqrt(np.random.random(num_points))
        theta = 2 * np.pi * np.random.random(num_points)

        cos_t = np.cos(theta) 
        sin_t = np.sin(theta)

        # Local coordinates
        points = (
            apex + np.outer(h, axis) + np.outer(r * cos_t, u) + np.outer(r * sin_t, v)
        )

        return points 
    
    @staticmethod 
    def sample_cylinder(
        base_point: np.ndarray,
        axis: np.ndarray,
        height: float,
        radius: float,
        num_points: int 
        ) -> np.ndarray:
        """  
        Uniform sampling inside of a cylinder volume

        Args:
            base_point (np.ndarray): Base center point of the cylinder
            axis (np.ndarray): Axis vector of the cylinder
            height (float): Height of the cylinder
            radius (float): Radius of the cylinder
            num_points (int): Number of points to sample
        """
        base_point = np.asarray(base_point, dtype=float)
        axis = np.asarray(axis, dtype=float)
        axis = axis / np.linalg.norm(axis)  # Normalize the axis vector

        # Construct Orthonormal Basis
        u, v = GeometricSampler._build_orthonormal_basis(axis)

        # Uniform sampling along height
        h = height * np.random.random(num_points)

        # Uniform distribution inside disk 
        r = radius * np.sqrt(np.random.random(num_points))
        theta = 2 * np.pi * np.random.random(num_points)

        cos_t = np.cos(theta)
        sin_t = np.sin(theta)

        # Local coordinates
        points = ( 
            base_point + np.outer(h, axis) + np.outer(r * cos_t, u) + np.outer(r * sin_t, v)
        )
        return points
    
    @staticmethod
    def sample_rectangle(
        center: np.ndarray, 
        dir1: np.ndarray,
        dir2: np.ndarray,
        length1: float,
        length2: float,
        num_points: int 
        ) -> np.ndarray:
        """  
        Uniform sampling inside of a rectangular volume
        
        Args: 
            center (np.ndarray): Center point of the rectangle
            dir1 (np.ndarray): Direction vector 1
            dir2 (np.ndarray): Direction vector 2
            length1 (float): Half-length along dir1
            length2 (float): Half-length along dir2
            num_points (int): Number of points to sample
        """

        # Normalize direction vectors
        d1 = dir1 / np.linalg.norm(dir1)
        # Ensure orthogonality
        d2 = dir2 - np.dot(dir2, d1) * d1
        d2 /= np.linalg.norm(d2)

        # Uniform sampling in [-1, 1] range 
        u = np.random.uniform(-1, 1, num_points)
        v = np.random.uniform(-1,1, num_points)

        # Scale 
        points = center + np.outer(u * length1, d1) + np.outer(v * length2, d2)

        return points
    


    @staticmethod
    def _build_orthonormal_basis(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """ 
        Helper function to build an orthonormal basis given a single vector
        """
        # Choose non parallel vector
        if abs(np.dot(axis, [0,0,1])) > 0.9:
            guess = np.array([0,1,0], dtype=float)
        else:
            guess = np.array([0,0,1], dtype=float)

        u = np.cross(axis, guess)
        u = u / np.linalg.norm(u)

        v = np.cross(axis, u)

        return u, v
    
class RotationSampler:
    """   
    Class that is used for the rotational sampling operations
    """
    @staticmethod 
    def spherical_grid(grid_dim: int = 4) -> np.ndarray: 
        """  
        Generates spherical Grid of unit vectors for rotational sampling
        
        Args:
            grid_dim (int): Grid dimension along each spherical

        Returns:
            np.ndarray: Array of unit vectors on the sphere
        """
        phi = np.linspace(0, np.pi, grid_dim)
        theta = np.linspace(0, 2 * np.pi, grid_dim)
        phi_grid, theta_grid = np.meshgrid(phi, theta)

        # Convert to cartesian coordinates
        x = np.sin(phi_grid) * np.cos(theta_grid)
        y = np.sin(phi_grid) * np.sin(theta_grid)
        z = np.cos(phi_grid)

        return np.column_stack([x.ravel(), y.ravel(), z.ravel()])

    def circle_grid(grid_dim: int = 4) -> np.ndarray:
        """   
        Samples points in a unit circle grid in 2D
        Args:
            grid_dim (int): Grid dimension along each axis
        Returns:
            np.ndarray: Array of points in the circle
        """
        theta = np.linspace(0, 2 * np.pi, grid_dim)
        r = np.linspace(0, 1, grid_dim)
        r_grid, theta_grid = np.meshgrid(r, theta)
        x = r_grid * np.cos(theta_grid)
        y = r_grid * np.sin(theta_grid)
        return np.column_stack([x.ravel(), y.ravel()]) 


    def fibonacci_sphere(n_points: int) -> np.ndarray: 
        """   
        Generates uniformly distributed points on the surface of a sphere using the Fibonacci lattice method
        
        Args:
            n_points (int): Number of points to generate

        Returns: 
            np.ndarray: Array of points on the sphere surface
        """
        golden_ratio = (1 + np.sqrt(5)) / 2
        indices = np.arange(n_points)

        theta = 2 * np.pi * indices / golden_ratio
        phi = np.arccos(1 - 2*(indices + 0.5) / n_points)

        # Convert to Cartesian Coordinates
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)

        return np.column_stack((x, y, z))
    
    def fibonacci_circle(n_points: int) -> np.ndarray:
        """  
        Generates points inside a circle using the Fibonacci lattice method

        t_i = (x_i, y_i) with (floor(i / pihi), i/n) for i = 0,..., n-1

        where phi = (1 + sqrt(5)) / 2 is the golden ratio
        """
        golden_ratio = (1 + np.sqrt(5)) / 2
        indices = np.arange(n_points)

        # Obtain Sampled Polar Coordinates
        theta = 2 * np.pi * indices / golden_ratio
        r = np.sqrt(indices) / np.sqrt(n_points - 1)

        # Convert to Cartesian Coordinates
        x = r * np.cos(theta)
        y = r * np.sin(theta)

        return np.column_stack((x, y))
    
    def canonical_circle(n_points: int) -> np.ndarray:
        """   
        Generates points inside a unit circle using the canonical method

        t_i(theta, phi) = (i / varphi,(i + 0.5)/n) for i = 0,..., n-1

        where varphi is the golden ratio

        This method divide the unit interval into n equal parts and places a point at the
        midpoint of each interval
        """
        golden_ratio = (1 + np.sqrt(5)) / 2
        indices = np.arange(n_points)
        theta = 2 * np.pi * indices / golden_ratio
        r = np.sqrt((indices + 0.5) / n_points)
        # Convert to Cartesian Coordinates
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        return np.column_stack((x, y))
    
    def canonical_sphere(n_points: int) -> np.ndarray:
        """  
        Generates points on the surface of a unit sphere using the canonical method

        t_i(theta, phi) = (i / varphi, (i + 0.5)/n) for i = 0,..., n-1

        where varphi is the golden ratio
        """
        golden_ratio = (1 + np.sqrt(5)) / 2
        indices = np.arange(n_points)
        theta = 2 * np.pi * indices / golden_ratio
        phi = np.arccos(1 - 2*(indices + 0.5) / n_points)
        # Convert to Cartesian Coordinates
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)
        return np.column_stack((x, y, z))

    
class ConfigSampler:
    """   
    Optimized Class for Configuration Sampling
    """
    def __init__(
            self,
            reference_frame: Optional[ReferenceFrame] = None,
            logger: Optional[object] = None
        ):
        """  
        Initializes the ConfigSampler with an optional reference frame and logger
        """
        self.reference_frame = reference_frame
        self.logger = logger
        self.sampled_configurations: List[Molecule] = []
        self.sampling_region: Optional[SamplingRegion] = None

        # Transformation Helper
        self._transformer = Transformation()

        self.sampling_statistics = {
            "num_samples": 0,
            "avg_distance": 0.0,
            "sampling_time": 0.0,
            "total_sampling_volume": 0.0,
            "avg_distance_to_center": 0.0,
            "total_molecular_volume": 0.0,
            "volume_coverage_ratio": 0.0
        }

    # ==================== Sampling Statistics Methods ===================== #
    def calculate_sampling_statistics(self, sampled_mols: List[Molecule], submolecule: "Submolecule") -> None: 
        """  
        Calculate basic sampling statistics and store them in self.sampling_statistics


        The following statistics of the sampling are calculated:

        + Number of Samples 
        + Total Sampling Volume 
        + Average Distance to Center
        + Nearest Neighbour Distances
        + Minimum Pairwise Distances
        
        Args:
            sampled_mols (List[Molecule]): List of sampled molecules
            submolecule (SubMolecule): The submolecule that was sampled
        
        """
        if not sampled_mols:
            raise ValueError("No sampled molecules provided for statistics calculation")
        if self.sampling_region is None:
            raise ValueError("Sampling region not defined")
        num_samples = len(sampled_mols)
        self.sampling_statistics["num_samples"] = num_samples

        # Calculate Total Sampling Volume
        total_volume = self._get_sampling_region_volume()
        self.sampling_statistics["total_sampling_volume"] = total_volume

        # Calculate Average Distance to Center
        center = self.sampling_region.center
        distances_to_center = []
        for mol in sampled_mols:
            submol_coords = mol.coordinates[submolecule.get_index_in_parent()]
            submol_center = np.mean(submol_coords, axis=0)
            dist = np.linalg.norm(submol_center - center)
            distances_to_center.append(dist)
        avg_distance_to_center = np.mean(distances_to_center)
        self.sampling_statistics["avg_distance_to_center"] = avg_distance_to_center

        # Calculate Pairwise Distances between sampled configurations
        submol_centers = []
        for mol in sampled_mols:
            submol_coords = mol.coordinates[submolecule.get_index_in_parent()]
            submol_center = np.mean(submol_coords, axis=0)
            submol_centers.append(submol_center)
        submol_centers = np.array(submol_centers)
        pairwise_dists = cdist(submol_centers, submol_centers)
        # Set diagonal to infinity to ignore self-distances
        np.fill_diagonal(pairwise_dists, np.inf)
        min_pairwise_dist = np.min(pairwise_dists)
        self.sampling_statistics["min_nn_distance"] = min_pairwise_dist
        # Calculate average nearest neighbour distance
        nearest_neighbour_dists = np.min(pairwise_dists, axis=1)
        avg_nearest_neighbour_dist = np.mean(nearest_neighbour_dists)
        self.sampling_statistics["avg_nearest_neighbour_distance"] = avg_nearest_neighbour_dist 
        

    def _get_sampling_region_volume(self) -> float:
        """   
        Helper Function to calculate the volume of the current sampling region
        """
        if self.sampling_region is None: 
           raise ValueError("No sampling region defined")

        region = self.sampling_region

        if region.shape == "sphere":
            radius = region.parameters["radius"]
            return (4/3) * np.pi * radius**3

    # =================== Sphere Sampling Methods =================== #
    def sample_sphere(
            self,
            center: np.ndarray,
            radius: float,
            num_points: int = 100,
            plot: bool = False
        ) -> np.ndarray:
        """  
        Samples points inside a sphere volume

        Args:
            center (np.ndarray): Center point of the sphere
            radius (float): Radius of the sphere
            num_points (int): Number of points to sample
            plot (bool): Whether to plot the sampled points
        
        Retunrs:
            np.ndarray: Sampled points inside the sphere
        """
        points = GeometricSampler.sample_sphere(center, radius, num_points)
        if plot:
            raise NotImplementedError("Sphere sampling plot not implemented yet")
        
        return points 
    
    def sample_submol_sphere(
            self,
            submolecule: "SubMolecule",
            center: np.ndarray,
            radius: float, 
            num_points: int = 100,
            rotation: bool = True,
            rotation_grid_dim: int = 4,
            rotation_method: str = "spherical", # "spherical" or "fibonacci"
            plot: bool = False 
        ) -> List[Molecule]:
        """  
        Sample a submolecule inside a spherical volume with optional rotations

        Args:
            submolecule (SubMolecule): Submolecule to sample
            center (np.ndarray): Center point of the sphere
            radius (float): Radius of the sphere
            num_points (int): Number of translation points to sample
            rotation (bool): Whether to include rotational sampling
            rotation_grid_dim (int): Grid dimension for rotational sampling
            rotation_method (str): Method for rotational sampling ("spherical" or "fibonacci")
            plot (bool): Whether to plot the sampled configurations
        """
        start_time = time.time()

        self.sampling_region = SamplingRegion(
            shape = "sphere",
            center = center,
            parameters={"radius": radius, "num_points": num_points}
        )


        # Sample the Translation Points
        sampled_points = self.sample_sphere(center, radius, num_points)


        # Generate Molecules 
        sampled_mols = self._generate_translated_rotated_molecules(
            submolecule=submolecule,
            translation_points = sampled_points,
            rotation = rotation,
            rotation_grid_dim = rotation_grid_dim,
            rotation_method = rotation_method,
            name_prefix = "sphere"
        )

        if plot:
            raise NotImplementedError("Sphere submolecule sampling plot not implemented yet")
        
        end_time = time.time()
        self.sampling_statistics["sampling_time"] = end_time - start_time
        
        return sampled_mols
        

    # =================== Cone Sampling Methods =================== #

    def sample_cone(
        self,
        apex: np.ndarray,
        axis: np.ndarray,
        height: float,
        base_radis: float,
        num_points: int = 100,
        plot: bool = False
        ) -> np.ndarray:
        """
        Samples points inside of a cone volume
        """ 
        points = GeometricSampler.sample_cone(apex, axis, height, base_radis, num_points)

        if plot:
            raise NotImplementedError("Cone sampling plot not implemented yet")
        return points
    
    def sample_submol_cone(
        self, 
        submolecule: "SubMolecule",
        apex: np.ndarray,
        axis: np.ndarray,
        height: float,
        base_radius: float,
        num_points: int = 100,
        rotation: bool = True,
        rotation_grid_dim: int = 4,
        rotation_method: str = "spherical", # "spherical" or "fibonacci"
        plot: bool = False
        ) -> List[Molecule]:
        """
        Samples a submolecule inside of a cone volume with optional rotations
        
        
        Args:
            submolecule (SubMolecule): Submolecule to sample
            apex (np.ndarray): Apex point of the cone
            axis (np.ndarray): Axis vector of the cone
            height (float): Height of the cone
            base_radius (float): Base radius of the cone
            num_points (int): Number of translation points to sample
            rotation (bool): Whether to include rotational sampling
            rotation_grid_dim (int): Grid dimension for rotational sampling
            rotation_method (str): Method for rotational sampling ("spherical" or "fibonacci")
            plot (bool): Whether to plot the sampled configurations
        """
        start_time = time.time()
        self.sampling_region = SamplingRegion(
            shape = "cone",
            center=apex,
            parameters={
                "axis": axis,
                "height": height,
                "base_radius": base_radius,
                "num_points": num_points
            }
        )

        sampled_points = self.sample_cone(apex, axis, height, base_radius, num_points)
        sampled_mols = self._generate_translated_rotated_molecules(
            submolecule = submolecule,
            translation_points = sampled_points,
            rotation = rotation,
            rotation_grid_dim = rotation_grid_dim,
            rotation_method = rotation_method,
            name_prefix = "cone"
        )

        if plot:
            raise NotImplementedError("Cone submolecule sampling plot not implemented yet")
        
        end_time = time.time()
        self.sampling_statistics["sampling_time"] = end_time - start_time
        return sampled_mols

    # ================== Cylinder Sampling Methods ================== #

    def sample_cylinder(
        self,
        base_point: np.ndarray,
        axis: np.ndarray,
        height: float,
        radius: float,
        num_points: int = 100,
        plot: bool = False
        ) -> np.ndarray:
        """
        Samples points inside of a cylinder volume
        """
        points = GeometricSampler.sample_cylinder(base_point, axis, height, radius, num_points)
        if plot:
            raise NotImplementedError("Cylinder sampling plot not implemented yet")
        return points

    def sample_submol_cylinder(
        self,
        submolecule: "SubMolecule",
        base_point: np.ndarray,
        axis: np.ndarray,
        height: float,
        radius: float,
        num_points: int = 100,
        rotation: bool = True,
        rotation_grid_dim: int = 4,
        rotation_method: str = "spherical", # "spherical" or "fibonacci"
        plot: bool = False
        ) -> List[Molecule]:
        """
        Samples a submolecule inside of a cylinder volume with optional rotations

        Args:
            submolecule (SubMolecule): Submolecule to sample
            base_point (np.ndarray): Base center point of the cylinder
            axis (np.ndarray): Axis vector of the cylinder
            height (float): Height of the cylinder
            radius (float): Radius of the cylinder
            num_points (int): Number of translation points to sample
            rotation (bool): Whether to include rotational sampling
            rotation_grid_dim (int): Grid dimension for rotational sampling
            rotation_method (str): Method for rotational sampling ("spherical" or "fibonacci")
            plot (bool): Whether to plot the sampled configurations
        """
        start_time = time.time()
        self.sampling_region = SamplingRegion(
            shape= "cylinder",
            center = base_point,
            parameters= {
                "axis": axis,
                "height": height,
                "radius": radius,
                "num_points": num_points
            }
        )

        sampled_points = self.sample_cylinder(base_point, axis, height, radius, num_points)
        
        sampled_mols = self._generate_translated_rotated_molecules( 
            submolecule = submolecule,
            translation_points = sampled_points,
            rotation = rotation,
            rotation_grid_dim = rotation_grid_dim,
            rotation_method = rotation_method,
            name_prefix = "cylinder"
        )

        if plot:
            raise NotImplementedError("Cylinder submolecule sampling plot not implemented yet")
        
        end_time = time.time()
        self.sampling_statistics["sampling_time"] = end_time - start_time
        return sampled_mols

    # ==================== Rectangle Sampling Methods ==================== #

    def sample_rectangle(
        self,
        center: np.ndarray,
        dir_vector1: np.ndarray,
        dir_vector2: np.ndarray,
        length1: float,
        length2: float,
        num_points: int = 100,
        plot: bool = False
        ) -> np.ndarray:
        """
        Samples points inside of a rectangular volume
        """
        points = GeometricSampler.sample_rectangle(
            center, dir_vector1, dir_vector2, length1, length2, num_points
        )

        if plot:
            raise NotImplementedError("Rectangle sampling plot not implemented yet")
        return points

    # TODO finish this method
    def sample_along_bond_coordinate(
            self,
            submolecule: "SubMolecule",
            bond_vector: np.ndarray,
            start_distance: float,
            end_distance: float,
            step_size: float,
            name_prefix: str = "bond_coord",
        ) -> List[Molecule]:
        """   
        Samples a submolecule along a bond coordinate and determines the Potential Energy Surface (PES)
        """
        start_time = time.time()

        parent = submolecule.parent
        submol_indices = submolecule.get_index_in_parent()
        template = self._create_template_molecule(parent)

        # Normalize the bond vector 
        bond_vector = np.asarray(bond_vector, dtype=float)
        bond_vector /= np.linalg.norm(bond_vector)

        # Generate distance points
        distances = np.arange(start_distance, end_distance + step_size, step_size)
        num_points = len(distances)

        # Store Sampling Region info
        submol_center = np.mean(submolecule.coordinates, axis=0)
        self.sampling_region = SamplingRegion(
           shape = "line",
           center = submol_center,
           parameters={
               "bond_vector": bond_vector,
               "start_distance": start_distance,
               "end_distance": end_distance,
               "step_size": step_size,
               "num_points": num_points
           }
        )

        # Pre-allocate molecule list 
        sampled_mols = [None] * num_points 

        # Get original submolecule center
        original_submol_center = np.mean(submolecule.coordinates, axis=0)     

        # Pre compute submolecule coordinates relative to center
        submol_coords_centered = submolecule.coordinates - original_submol_center 
        
        for i, dist in enumerate(distances):
            new_mol = Molecule(name=f"{name_prefix}_d{dist:.3f}_sample_{i}")
            new_mol.atom_labels = template.atom_labels.copy()
            new_mol.masses = template.masses.copy()

            # Calculate new position along bond vector 
            displacement = dist * bond_vector
            new_submol_center = original_submol_center + displacement

            # Build new coordinates
            new_coords = parent.coordinates.copy()
            new_coords[submol_indices] = submol_coords_centered + new_submol_center
            new_mol.coordinates = new_coords
            sampled_mols[i] = new_mol
        end_time = time.time()
        self.sampling_statistics["sampling_time"] = end_time - start_time

        return sampled_mols    
    
    def sample_submol_rectangle(
        self,
        submolecule: "SubMolecule",
        center: np.ndarray,
        dir_vector1: np.ndarray,
        dir_vector2: np.ndarray,
        length1: float,
        length2: float,
        num_points: int = 100,
        rotation: bool = True,
        rotation_grid_dim: int = 4,
        rotation_method: str = "spherical", # "spherical" or "fibonacci"
        plot: bool = False
        ) -> List[Molecule]:
        """
        Samples a submolecule inside of a rectangular volume with optional rotations
        """ 
        time_start = time.time()
        self.sampling_region = SamplingRegion(
            shape = "rectangle",
            center = center,
            parameters={
                "dir_vector1": dir_vector1,
                "dir_vector2": dir_vector2,
                "length1": length1,
                "length2": length2,
                "num_points": num_points
            }
        )
        sampled_points = self.sample_rectangle(
            center, dir_vector1, dir_vector2, length1, length2, num_points
        )

        sampled_mols = self._generate_translated_rotated_molecules(
            submolecule = submolecule,
            translation_points = sampled_points,
            rotation = rotation,
            rotation_grid_dim = rotation_grid_dim,
            rotation_method = rotation_method,
            name_prefix = "rectangle"
        )

        if plot:
            raise NotImplementedError("Rectangle submolecule sampling plot not implemented yet")
        
        end_time = time.time()
        self.sampling_statistics["sampling_time"] = end_time - time_start
        return sampled_mols
    
    # ==================== Core Molecule Generation Methods ==================== #
    def _generate_translated_rotated_molecules(
        self,
        submolecule: "SubMolecule",
        translation_points: np.ndarray,  # Absolute Positionen!
        rotation: bool,
        rotation_grid_dim: int,
        rotation_method: str,
        name_prefix: str
        ) -> List[Molecule]:
        """
        Main Method to generate translated and rotated molecules based on sampled points
    
        Args:
            submolecule (SubMolecule): Submolecule to be moved
            translation_points (np.ndarray): ABSOLUTE target positions (not vectors!)
            rotation (bool): Whether to include rotational sampling
            rotation_grid_dim (int): Grid dimension for rotational sampling
            rotation_method (str): Method for rotational sampling ("spherical" or "fibonacci")
            name_prefix (str): Prefix for naming the generated molecules
        """
        parent = submolecule.parent
        submol_indices = submolecule.get_index_in_parent()

        # Template-Molecule 
        template = self._create_template_molecule(parent)

        if not rotation:
            # Just translations - simpler path
            return self._generate_translated_molecules_simple(
                template=template,
                parent_coords=parent.coordinates,
                submolecule=submolecule,
                submol_indices=submol_indices,
                target_positions=translation_points,
                name_prefix=name_prefix
            )

        # With rotations
        return self._generate_translated_rotated_molecules_batch(
            submolecule=submolecule,
            template=template,
            parent_coords=parent.coordinates,
            submol_indices=submol_indices,
            target_positions=translation_points,  # Direct positions!
            rotation_grid_dim=rotation_grid_dim,
            rotation_method=rotation_method,
            name_prefix=name_prefix
        )
    
    def _generate_translated_molecules_simple(
        self,
        template: Molecule,
        parent_coords: np.ndarray,
        submolecule: "SubMolecule",
        submol_indices: List[int],
        target_positions: np.ndarray,
        name_prefix: str
        ) -> List[Molecule]:
        """
        Generate molecules with only translations (no rotation)
        """
        num_points = len(target_positions)
        sampled_mols = [None] * num_points

        # Current submolecule center
        submol_center = np.mean(submolecule.coordinates, axis=0)

        # Mask for efficient updates
        mask = np.zeros(parent_coords.shape[0], dtype=bool)
        mask[submol_indices] = True

        # Submolecule coordinates relative to center
        submol_coords_centered = submolecule.coordinates - submol_center

        for i, target_pos in enumerate(target_positions):
            new_mol = Molecule(name=f"{name_prefix}_sample_{i}")
            new_mol.atom_labels = template.atom_labels.copy()
            new_mol.masses = template.masses.copy()

            # Efficient coordinate assignment
            new_coords = parent_coords.copy()
            new_coords[mask] = submol_coords_centered + target_pos
            new_mol.coordinates = new_coords

            sampled_mols[i] = new_mol

        return sampled_mols

    def _generate_translated_rotated_molecules_batch(
        self,
        submolecule: "SubMolecule",
        template: Molecule,
        parent_coords: np.ndarray,
        submol_indices: List[int],
        target_positions: np.ndarray,  # Changed from translation_vectors!
        rotation_grid_dim: int,
        rotation_method: str,
        name_prefix: str
        ) -> List[Molecule]:
        """  
        Generates molecules with Translations and Rotations

        Args:
            target_positions: Absolute 3D positions where submolecule center should be placed
        """
        # Get rotation vectors
        if rotation_method == "spherical":
            rotation_vectors = RotationSampler.spherical_grid(rotation_grid_dim)
        elif rotation_method == "fibonacci":
            n_rotations = rotation_grid_dim ** 2
            rotation_vectors = RotationSampler.fibonacci_sphere(n_rotations)
        else:
            raise ValueError(f"Unknown rotation method: {rotation_method}")

        # Calculate original submolecule center and center it
        submol_center = np.mean(submolecule.coordinates, axis=0)
        centered_submol_coords = submolecule.coordinates - submol_center

        # Compute reference frame for centered submolecule
        ref_frame = self._transformer.compute_reference_frame(
            centered_submol_coords,
            submolecule.masses,
            center=np.zeros(3)
        )

        # Pre-allocation
        total_mols = len(rotation_vectors) * len(target_positions)
        sampled_mols = [None] * total_mols

        # Mask for efficient coordinate assignment
        mask = np.zeros(parent_coords.shape[0], dtype=bool)
        mask[submol_indices] = True

        mol_counter = 0

        # For each rotation
        for rot_idx, rot_vec in enumerate(rotation_vectors):
            # Create temporary submolecule object for rotation
            from types import SimpleNamespace
            temp_submol = SimpleNamespace(
                coordinates=centered_submol_coords.copy(),
                masses=submolecule.masses
            )

            # Rotate submolecule (centered at origin)
            rotated_coords = self._transformer.align_molecule_to_target(
                temp_submol,
                ref_frame,
                rot_vec  # Direction vector from origin
            )

            # For each target position
            for trans_idx, target_pos in enumerate(target_positions):
                new_mol = Molecule(name=f"{name_prefix}_rot{rot_idx}_trans{trans_idx}")
                new_mol.atom_labels = template.atom_labels.copy()
                new_mol.masses = template.masses.copy()

                # Build coordinates: parent + rotated submolecule at target position
                new_coords = parent_coords.copy()
                new_coords[mask] = rotated_coords + target_pos
                new_mol.coordinates = new_coords

                sampled_mols[mol_counter] = new_mol
                mol_counter += 1

        return sampled_mols

    
    # ==================== Helper Methods ==================== #

    def _create_template_molecule(self, molecule: "Molecule") -> "Molecule":
        """  
        Creates a Template Molecules with Labels and Masses only
        """
        template = Molecule(name="template")
        template.atom_labels = molecule.atom_labels.copy()
        template.masses = molecule.masses.copy()
        return template

    # ================== Testing and Debugging Methods ================== #
    def plot_and_test_fibonacci_circle(self) -> None:
        """  
        Plots the Fibonacci Circle Sampling for n = [10, 20, 50, 100], determines the speed 
        and calculates average distance between points. Further plots the points
        """
        grid_dim = [8,10,20,40]
        n_points_list = [dim**2 for dim in grid_dim]
        fig, axs = plt.subplots(2,2, figsize=(10,10))
        axs = axs.ravel()
        for i, n_points in enumerate(n_points_list):
            start_time = time.time()
            points = RotationSampler.fibonacci_circle(n_points)
            end_time = time.time()
            elapsed_time = end_time - start_time
            # Calculate average distance between points
            dists = []
            for j in range(n_points):
                for k in range(j+1, n_points):
                    dist = np.linalg.norm(points[j] - points[k])
                    dists.append(dist)
            avg_dist = np.mean(dists)
            axs[i].scatter(points[:,0], points[:,1])
            axs[i].set_title(f"Fibonacci Circle n={n_points}\nTime: {elapsed_time:.4f}s, Avg Dist: {avg_dist:.4f}")
            axs[i].set_xlabel("X")
            axs[i].set_ylabel("Y")
            axs[i].axis("equal")
        plt.tight_layout()
        plt.savefig("figures/fibonacci_circle_test.png")
        plt.close()

    def plot_and_test_circle_grid(self) -> None:
        """   
        Tests and plots the circle grid sampling method for N = [10, 20, 50, 100]
        """
        n_points_list = [8, 10, 20, 40]
        fig, axs = plt.subplots(2,2, figsize=(10,10))
        axs = axs.ravel()
        for i, n_points in enumerate(n_points_list):
            start_time = time.time()
            points = RotationSampler.circle_grid(n_points)
            end_time = time.time()
            elapsed_time = end_time - start_time
            # Calculate average distance between points
            dists = []
            for j in range(len(points)):
                for k in range(j+1, len(points)):
                    dist = np.linalg.norm(points[j] - points[k])
                    dists.append(dist)
            avg_dist = np.mean(dists)
            axs[i].scatter(points[:,0], points[:,1])
            axs[i].set_title(f"Circle Grid n={n_points}\nTime: {elapsed_time:.4f}s, Avg Dist: {avg_dist:.4f}")
            axs[i].set_xlabel("X")
            axs[i].set_ylabel("Y")
            axs[i].axis("equal")
        plt.tight_layout()
        plt.savefig("figures/circle_grid_test.png")
        plt.close()

    def plot_and_test_canonical_circle(self) -> None:
        """   
        Tests and plots the canonical circle sampling method for grid_dim = [8, 10, 20, 40]

        Note that we need to square the grid_dim to get the number of points to get equal results
        """
        grid_dims = [8, 10, 20, 40]
        n_points_list = [dim**2 for dim in grid_dims]
        fig, axs = plt.subplots(2,2, figsize=(10,10))
        axs = axs.ravel()
        for i, n_points in enumerate(n_points_list):
            start_time = time.time()
            points = RotationSampler.canonical_circle(n_points)
            end_time = time.time()
            elapsed_time = end_time - start_time
            # Calculate average distance between points
            dists = []
            for j in range(len(points)):
                for k in range(j+1, len(points)):
                    dist = np.linalg.norm(points[j] - points[k])
                    dists.append(dist)


            avg_dist = np.mean(dists)
            axs[i].scatter(points[:,0], points[:,1])
            axs[i].set_title(f"Canonical Circle n={n_points}\nTime: {elapsed_time:.4f}s, Avg Dist: {avg_dist:.4f}")
            axs[i].set_xlabel("X")
            axs[i].set_ylabel("Y")
            axs[i].axis("equal")
        plt.tight_layout()
        plt.savefig("figures/canonical_circle_test.png")
        plt.close()

    def plot_and_test_spherical_grid(self) -> None:
        """  
        Tests and plots the spherical grid sampling method for grid_dim = [8, 10, 20, 40]
        """
        grid_dims = [8, 10, 20, 40]
        fig = plt.figure(figsize=(16,16))
        axs = [fig.add_subplot(2,2,i+1, projection='3d') for i in range(4)]
        for i, grid_dim in enumerate(grid_dims):
            start_time = time.time()
            points = RotationSampler.spherical_grid(grid_dim)
            end_time = time.time()
            elapsed_time = end_time - start_time
            # Calculate average distance between points
            dists = []
            for j in range(len(points)):
                for k in range(j+1, len(points)):
                    dist = np.linalg.norm(points[j] - points[k])
                    dists.append(dist)
            avg_dist = np.mean(dists)
            # Draw unit sphere as wireframe
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            x = np.outer(np.cos(u), np.sin(v))
            y = np.outer(np.sin(u), np.sin(v))
            z = np.outer(np.ones(np.size(u)), np.cos(v))
            axs[i].plot_wireframe(x, y, z, color='lightgray', alpha=0.5)

            axs[i].scatter(points[:,0], points[:,1], points[:,2], color="r", s=15, )
            axs[i].set_title(f"Spherical Grid dim={grid_dim}\nTime: {elapsed_time:.4f}s, Avg Dist: {avg_dist:.4f}")
            axs[i].set_xlabel("X")
            axs[i].set_ylabel("Y")
            axs[i].set_zlabel("Z")
            axs[i].set_box_aspect([1,1,1])
        plt.tight_layout()
        plt.savefig("figures/spherical_grid_test.png")
        plt.close()

    def plot_and_test_fibonacci_sphere(self) -> None:
        """
        Plots the Fibonacci Sphere Sampling for n = [10, 20, 50, 100], determines the speed
        """
        grid_dim = [8, 10, 20, 40]
        n_points = [dim**2 for dim in grid_dim]
        fig = plt.figure(figsize=(16,16))
        axs = [fig.add_subplot(2,2,i+1, projection='3d') for i in range(4)]
        for i, n in enumerate(n_points):
            start_time = time.time()
            points = RotationSampler.fibonacci_sphere(n)
            end_time = time.time()
            elapsed_time = end_time - start_time
            # Calculate average distance between points
            dists = []
            for j in range(len(points)):
                for k in range(j+1, len(points)):
                    dist = np.linalg.norm(points[j] - points[k])
                    dists.append(dist)
            avg_dist = np.mean(dists)
            # Draw unit sphere as wireframe
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            x = np.outer(np.cos(u), np.sin(v))
            y = np.outer(np.sin(u), np.sin(v))
            z = np.outer(np.ones(np.size(u)), np.cos(v))
            axs[i].plot_wireframe(x, y, z, color='lightgray', alpha=0.5)

            axs[i].scatter(points[:,0], points[:,1], points[:,2], color="b", s=15)
            axs[i].set_title(f"Fibonacci Sphere n={n}\nTime: {elapsed_time:.4f}s, Avg Dist: {avg_dist:.4f}")
            axs[i].set_xlabel("X")
            axs[i].set_ylabel("Y")
            axs[i].set_zlabel("Z")
            axs[i].set_box_aspect([1,1,1])
        plt.tight_layout()
        plt.savefig("figures/fibonacci_sphere_test.png")
        plt.close()

    def plot_and_test_canonical_sphere(self) -> None:
        """   
        Plots the Canonical Sphere Sampling for grid_dim = [8, 10, 20, 40], determines the speed
        and calculates average distance between points. Further plots the points
        """
        grid_dims = [8,10,20,40]
        fig = plt.figure(figsize=(16,16))
        axs = [fig.add_subplot(2,2,i+1, projection='3d') for i in range(4)]
        for i, grid_dim in enumerate(grid_dims):
            n_points = grid_dim **2
            start_time = time.time()
            points = RotationSampler.canonical_sphere(n_points)
            end_time = time.time()
            elapsed_time = end_time - start_time
            # Calculate average distance between points
            dists = []
            for j in range(len(points)):
                for k in range(j+1, len(points)):
                    dist = np.linalg.norm(points[j] - points[k])
                    dists.append(dist)
            avg_dist = np.mean(dists)
            # Draw unit sphere as wireframe
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            x = np.outer(np.cos(u), np.sin(v))
            y = np.outer(np.sin(u), np.sin(v))
            z = np.outer(np.ones(np.size(u)), np.cos(v))
            axs[i].plot_wireframe(x, y, z, color='lightgray', alpha=0.5)

            axs[i].scatter(points[:,0], points[:,1], points[:,2], color="g", s=15)
            axs[i].set_title(f"Canonical Sphere dim={grid_dim}\nTime: {elapsed_time:.4f}s, Avg Dist: {avg_dist:.4f}")
            axs[i].set_xlabel("X")
            axs[i].set_ylabel("Y")
            axs[i].set_zlabel("Z")
            axs[i].set_box_aspect([1,1,1])
        plt.tight_layout()
        plt.savefig("figures/canonical_sphere_test.png")
        plt.close()

    def test_sampling_speed(self,
                            submolecule: "SubMolecule",
                            method="sphere", 
                            rotation: bool = False, 
                            rotation_grid_dim=4) -> None:
        """  
        Tests the sampling speed for different methods and configurations

        Args:
            method (str): Sampling method to test ("sphere", "cone", "cylinder", "rectangle")
            rotation (bool): Whether to include rotational sampling
            rotation_grid_dim (int): Grid dimension for rotational sampling
        """
        if method == "sphere":
            center = np.array([0.0, 0.0, 0.0])  # Set center to origin
            radius = 5.0
            points_list = [100, 500, 1000, 5000, 10000, 20000, 50000]
            times = []
            for num_points in points_list:
                start_time = time.time()
                self.sample_submol_sphere(
                    submolecule=submolecule,
                    center=center,
                    radius=radius,
                    num_points=num_points,
                    rotation=rotation,
                    rotation_grid_dim=rotation_grid_dim,
                    rotation_method="spherical",
                    plot=False
                )
                end_time = time.time()
                elapsed_time = end_time - start_time
                times.append(elapsed_time)
            fig = plt.figure(figsize=(8,6))
            plt.plot(points_list, times, marker='o')
            plt.xlabel("Number of Sampled Points")
            plt.ylabel("Sampling Time (s)")
            plt.title(f"Sampling Speed Test - Sphere Method\nRotation: {rotation}, Grid Dim: {rotation_grid_dim}")
            plt.xlabel("Number of Sampled Points")
            plt.ylabel("Sampling Time (s)")
            plt.tight_layout()
            plt.grid()
            plt.savefig("figures/sampling_speed_sphere.png")
            plt.close()
    
    def test_convergence_metrics(
        self,
        submolecule: "SubMolecule",
        center: np.ndarray = np.array([0.0, 0.0, 0.0]),
        radius: float = 5.0,
        method="sphere",
        rotation=False,
        rotation_grid_dim=4,
        max_points=1e6
        ) -> None:
        """
        Test the convergence of the average pairwise distance 
        for different sampling regions and methods

        Metrics:
            - Average Nearest Neighbor Distance
            - Average Distance to Sampling Center
            - Minimum pairwise distance 
            - Standard Deviation of Distances to Center

        """ 
        start_time = time.time()
        print("Starting Convergence Test for Pairwise Distance...")
        if method == "sphere":

            # Storage for different metrics
            avg_nn_distances = [] # Nearest Neighbor
            avg_center_distances = [] # Distance to Center
            min_pairwise_distances = [] # Minimum Pairwise Distance
            std_center_distances = [] # Std Dev of Distances to Center

            points_list = [100, 500, 1000, 5000, 10000, 20000, 30000, 50000]
            avg_distances = []
            
            for num_points in points_list:
                print(f"Sampling {num_points} points...")

                sampled_mols = self.sample_submol_sphere(
                    submolecule=submolecule,
                    center=center,
                    radius=radius,
                    num_points=num_points,
                    rotation=rotation,
                    rotation_grid_dim=rotation_grid_dim,
                    rotation_method="spherical",
                    plot=False
                )

                # Extract Submolecule Centers
                submol_indices = submolecule.get_index_in_parent()
                centers = np.array([
                    np.mean(mol.coordinates[submol_indices], axis=0) for mol in sampled_mols
                ])
                # Nearest Neighbour Distances
                # We use cdist this computes the distance between all pairs efficiently
                dist_matrix = cdist(centers, centers)
                # Set Diagonal to infinity to ignore self-distances
                np.fill_diagonal(dist_matrix, np.inf)

                # Nearest neighbor for each point
                nn_distances = np.min(dist_matrix, axis=1)
                avg_nn_dist = np.mean(nn_distances)
                avg_nn_distances.append(avg_nn_dist)

                # Distances to Sampling Center
                center_dists = np.linalg.norm(centers - center, axis=1)
                avg_center_dist = np.mean(center_dists)
                std_center_dist = np.std(center_dists)
                avg_center_distances.append(avg_center_dist)
                std_center_distances.append(std_center_dist)

                # Minimum pairwise distance
                min_dist = np.min(dist_matrix)
                min_pairwise_distances.append(min_dist)

            # Plotting the metrics
            fig, axs = plt.subplots(2,2, figsize=(16,12))
            axs = axs.ravel()
            # Set bigger font size
            plt.rcParams.update({'font.size': 16})


            axs[0].plot(points_list, avg_nn_distances, marker='o')
            axs[0].set_title("Average Nearest Neighbor Distance")
            axs[0].set_xlabel("Number of Sampled Points", fontsize=14)
            axs[0].set_ylabel("Distance (Å)", fontsize=14)
            axs[0].grid()
            axs[1].plot(points_list, avg_center_distances, marker='o', color='orange')
            axs[1].set_title("Average Distance to Sampling Center")
            axs[1].set_xlabel("Number of Sampled Points", fontsize=14)
            axs[1].set_ylabel("Distance (Å)", fontsize=14)
            axs[1].grid()
            axs[2].plot(points_list, min_pairwise_distances, marker='o', color='green')
            axs[2].set_title("Minimum Pairwise Distance")
            axs[2].set_xlabel("Number of Sampled Points", fontsize=14)
            axs[2].set_ylabel("Distance (Å)", fontsize=14)
            axs[2].grid()
            axs[3].plot(points_list, std_center_distances, marker='o', color='red')
            axs[3].set_title("Std Dev of Distances to Center")
            axs[3].set_xlabel("Number of Sampled Points", fontsize=14)  
            axs[3].set_ylabel("Distance (Å)", fontsize=14)
            axs[3].grid()

            plt.tight_layout()
            
            plt.savefig("figures/convergence_pairwise_distance.png")
            plt.close()

            