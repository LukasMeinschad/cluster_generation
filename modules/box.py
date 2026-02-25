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
        self._validate_parameters()
        
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
        Reflects the coordinates back into the box.
        Uses iterative reflection to handle overshoots.
        """
        coords = coordinates.copy()
        max_iterations = 10  # Prevent infinite loops
        
        if self.box_type == "sphere":
            for iteration in range(max_iterations):
                all_inside = True
                for i, coord in enumerate(coords):
                    rel_pos = coord - self.center
                    dist = np.linalg.norm(rel_pos)
                    if dist > self.radius:
                        all_inside = False
                        # Reflect: fold the excess distance back inside
                        direction = rel_pos / dist
                        # Use modular reflection to handle any overshoot
                        excess = dist - self.radius
                        # Fold excess back: if excess > 2*radius, mod it
                        excess = excess % (2 * self.radius)
                        if excess <= self.radius:
                            coords[i] = self.center + (self.radius - excess) * direction
                        else:
                            coords[i] = self.center + (excess - self.radius) * direction
                if all_inside:
                    break

        elif self.box_type == "cube":
            half_dims = self.box_dimensions / 2.0 
            rel_coords = coords - self.center

            # Use modular reflection (like a zigzag function) for each dimension
            for dim in range(3):
                L = 2.0 * half_dims[dim]  # Full box length
                # Shift so that box goes from 0 to L
                shifted = rel_coords[:, dim] + half_dims[dim]
                # Apply zigzag (triangle wave) reflection
                # Period is 2L, fold into [0, 2L], then reflect
                shifted = shifted % (2 * L)
                mask = shifted > L
                shifted[mask] = 2 * L - shifted[mask]
                rel_coords[:, dim] = shifted - half_dims[dim]

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
            # Proper spherical sampling: r = R * cbrt(U), theta/phi uniform
            for i in range(n_points):
                u = np.random.uniform(0, 1)
                r = self.radius * np.cbrt(u)
                phi = np.random.uniform(0, 2 * np.pi)
                costheta = np.random.uniform(-1, 1)
                sintheta = np.sqrt(1 - costheta**2)
                positions[i] = self.center + r * np.array([
                    sintheta * np.cos(phi),
                    sintheta * np.sin(phi),
                    costheta
                ])
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
        elif self.box_type == "cube":
            return 2.0 * np.sum(self.box_dimensions[[0,0,1]] * self.box_dimensions[[1,2,2]])
        return 0.0
    

    def to_dict(self) -> dict:
        """  
        Serialize the SimulationBox to a dictionary
        """
        return {
            "box_type": self.box_type,
            "radius": self.radius,
            "box_dimensions": self.box_dimensions.tolist() if self.box_dimensions is not None else None,
            "center": self.center.tolist() if self.center is not None else None
        }
    @staticmethod
    def from_dict(data: dict) -> 'SimulationBox':
        """  
        Deserialize a SimulationBox from a dictionary
        """
        return SimulationBox(
            box_type=data.get("box_type", "sphere"),
            radius=data.get("radius"),
            box_dimensions=np.array(data["box_dimensions"]) if data.get("box_dimensions") is not None else None,
            center=np.array(data["center"]) if data.get("center") is not None else None
        )
