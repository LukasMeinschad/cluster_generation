"""  
Module for the Generation of the Simulation Box for the BHMC Algorithm
"""

import numpy as np
from typing import Optional, List, Tuple, Union
import matplotlib.pyplot as plt

class SimulationBox:
    """   
    Defines a simulation box to constrain molecular structures during BHMC

    Supports both spherical and cubic boundary conditions
    """
    def __init__(self,
                 box_type: str = "sphere",
                 radius: Optional[float] = None,
                 box_dimensions: Optional[np.ndarray] = None,
                 center: Optional[np.ndarray] = None):
        """  
        Initializes the Simulation Box

        Args:
            box_type: Type of box ("sphere" or "cube")
            radius: Radius for spherical box
            box_dimensions: Dimensions for cubic box (lengths in x, y, z)
            center: Center of the box (default is origin)
        """
        self.box_type = box_type.lower()
        self.radius = radius
        self.box_dimensions = np.array(box_dimensions) if box_dimensions is not None else None
        self.center = np.array(center) if center is not None else np.zeros(3)
        
    def _validate_parameters(self):
        """   
        Validates the box parameters
        """
        if self.box_type not in ["sphere", "cube"]:
            raise ValueError("box_type must be 'sphere' or 'cube'")
        if self.box_type == "sphere" and self.radius is None:
            raise ValueError("Radius must be provided for spherical box")
        if self.box_type == "cube" and self.box_dimensions is None:
            raise ValueError("Box dimensions must be provided for cubic box")
        if self.box_type == "cube" and np.any(self.box_dimensions <= 0):
            raise ValueError("Box dimensions must be positive")
        
    @staticmethod 
    def from_covalent_radii(covalent_radii: Union[List[float], np.ndarray],
                            n_atoms: int,
                            box_type: str = "sphere",
                            scale_factor: float = 2.0) -> 'SimulationBox':
        """   
        Creates a SimulationBox using the formula
        R = [(3/(4π)) * (2*Rc + 1)^2 * N]^(1/3)

        Args:
            covalent_radii: List or array of covalent radii for each atom
            n_atoms: Number of atoms in the molecule
            box_type: Type of box ("sphere" or "cube")
            scale_factor: Scaling factor to ensure sufficient space (default is 2.0)
        """
        # Average Covalent Radius
        R_c = np.mean(covalent_radii)
        R_sphere = ((3.0 / (4.0 * np.pi)) * (2 * R_c + 1)**2 * n_atoms)**(1/3) * scale_factor

        if box_type.lower() == "sphere":
            return SimulationBox(box_type="sphere", radius=R_sphere)
        elif box_type.lower() == "cube":
            # Cubic box with the Volume V = 4pi R^3 / 3
            side_length = (4.0 * np.pi * R_sphere**3 / 3.0)**(1/3)
            dimensions = np.array([side_length, side_length, side_length])
            return SimulationBox(box_type="cube", box_dimensions=dimensions)
        else:
            raise ValueError("box_type must be 'sphere' or 'cube'")
        
    def is_inside(self, coordinates: np.ndarray) -> bool:
        """   
        Checks if all coordinates are inside the box
        
        Args:
            coordinates: Array of shape (N, 3) representing atomic coordinates
        """
        coords = np.atleast_3d(coordinates)
        if self.box_type == "sphere":
            dist = np.linalg.norm(coords - self.center, axis=-1)
            return np.all(dist <= self.radius)
        
        elif self.box_type == "cube":
            half_dims = self.box_dimensions / 2.0
            rel_coords = coords - self.center
            inside = np.all(np.abs(rel_coords) <= half_dims, axis=-1)
            return np.all(inside)
        
        return False
    
    def check_atoms_inside(self, coordinates: np.ndarray) -> np.ndarray:
        """   
        Checks which atoms are inside the box

        Args:
            coordinates: Array of shape (N, 3) representing atomic coordinates
        """
        coords = np.atleast_3d(coordinates)

        if self.box_type == "sphere":
            dist = np.linalg.norm(coords - self.center, axis=-1)
            return dist <= self.radius
        elif self.box_type == "cube":
            half_dims = self.box_dimensions / 2.0
            rel_coords = coords - self.center
            inside = np.all(np.abs(rel_coords) <= half_dims, axis=-1)
            return inside
                
        return np.zeros(len(coords), dtype=bool)

    def apply_boundary_conditions(self,
                                  coordinates: np.ndarray,
                                  method: str = "reject") -> Tuple[np.ndarray, bool]:
        """   
        Applys the boundary condition to coordinates

        Args:
            coordinates: Atomic coordinates (N,3)
            method: Method to apply ("reject", "reflect", "wrap")
                - "reject": Rejects coordinates outside the box
                - "reflect": Reflects coordinates back into the box
                - "wrap": Wraps coordinates around the box (periodic)
        """
        if self.is_inside(coordinates):
            return coordinates.copy(), True
        if method == "reject":
            return coordinates.copy(), False
        elif method == "reflect":
            return self._reflect_coordinates(coordinates)
        elif method == "wrap":
            #TODO Implement this someday
            raise NotImplementedError("Wrap method is not implemented yet")

    def _reflect_coordinates(self, coordinates: np.ndarray) -> Tuple[np.ndarray, bool]:
        """   
        Reflects the coordinates back into the box
        """
        coords = coordinates.copy()
        if self.box_type == "sphere":
            for i, coord in enumerate(coords):
                # Vector from center to coordinate
                rel_pos = coord - self.center
                dist = np.linalg.norm(rel_pos)
                if dist > self.radius:
                    # Reflect back into sphere
                    direction = rel_pos / dist
                    excess = dist - self.radius
                    coords[i] = coord - 2 * excess * direction

        elif self.box_type == "cube":
            half_dims = self.box_dimensions / 2.0 
            rel_coords = coords - self.center

            # Reflect at each boundary
            for dim in range(3):
                mask_high = rel_coords[:, dim] > half_dims[dim]
                mask_low = rel_coords[:, dim] < -half_dims[dim]
                
                # Reflect atoms beyond upper boundary
                excess_high = rel_coords[mask_high, dim] - half_dims[dim]
                rel_coords[mask_high, dim] = half_dims[dim] - excess_high

                # Reflect atoms beyond lower boundary
                excess_low = -half_dims[dim] - rel_coords[mask_low, dim]
                rel_coords[mask_low, dim] = -half_dims[dim] + excess_low

            coords = rel_coords + self.center

        return coords, self.is_inside(coords)
    
    def get_random_position(self, n_points: int=1) -> np.ndarray:
        """    
        Generate Random positions inside the simulation box

        Args:
            n_points: Number of random points to generate
        """
        positions = np.zeros((n_points, 3))

        if self.box_type == "sphere":
            # Uniform distribution in a sphere
            for i in range(n_points):
                while True:
                    pos = np.random.uniform(-1,1,3) * self.radius
                    if np.linalg.norm(pos) <= self.radius:
                        positions[i] = self.center + pos
                        break
        elif self.box_type == "cube":
            # Uniform distribution in the cube
            half_dims = self.box_dimensions / 2.0
            positions = np.random.uniform(-half_dims, half_dims, (n_points, 3)) 
            positions += self.center

        return positions if n_points > 1 else positions[0]
    
    def get_volume(self) -> float:
        """  
        Calculate volume of the simulation box
        """
        if self.box_type == "sphere":
            return (4.0 / 3.0) * np.pi * self.radius**3
        elif self.box_type == "cube":
            return np.prod(self.box_dimensions)
        return 0.0
    
    def get_surface_area(self) -> float:
        """  
        Calculate surface area of the simulation box
        """
        if self.box_type == "sphere":
            return 4.0 * np.pi * self.radius**2 
        elif self.box_type == "cubic":
            return 2.0 * np.sum(self.box_dimensions[[0,0,1]] * self.box_dimensions[[1,2,2]])
        return 0.0
    
    