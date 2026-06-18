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
        
    @staticmethod
    def from_vdw_radii(vdw_radii: Union[List[float], np.ndarray],
                       n_atoms: int,
                       box_type: str = "sphere",
                       eta_factor: float = 0.4,
                       padding: float = 0.0) -> 'SimulationBox':
        """
        Generates a Simulation Box based on vdW Sphere Packing, plus a fixed
        padding distance added on top to leave room for rigid-body moves.

        V_vdw = sum_i^N (4/3 * pi * R_vdw_i^3)
        R_eff = (V_vdw / eta)^(1/3) + padding
        """
        R_vdw = np.array(vdw_radii)
        V_vdw = np.sum((4.0 / 3.0) * np.pi * R_vdw**3)
        if box_type.lower() == "sphere":
            R_eff = (V_vdw / eta_factor)**(1/3) + padding
            return SimulationBox(box_type="sphere", radius=R_eff)
        elif box_type.lower() == "cube":
            side_length = (V_vdw / eta_factor)**(1/3) + 2.0 * padding
            dimensions = np.array([side_length, side_length, side_length])
            return SimulationBox(box_type="cube", box_dimensions=dimensions)
        else:
            raise ValueError("box_type must be 'sphere' or 'cube'")
        
        
    def is_inside(self, coordinates: np.ndarray) -> bool:
        """   
        Checks if the centroid of the coordinates is inside the box
        
        Args:
            coordinates: Array of shape (N, 3) representing atomic coordinates
        """
        coords = np.atleast_2d(coordinates)

        inside_mask = self.check_atoms_inside(coords)
        
        return bool(np.all(inside_mask))
    
    def check_atoms_inside(self, coordinates: np.ndarray) -> np.ndarray:
        """   
        Checks which atoms are inside the box

        Args:
            coordinates: Array of shape (N, 3) representing atomic coordinates
        """
        coords = np.atleast_2d(coordinates)
        

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
                                  method: str = "reflect") -> Tuple[np.ndarray, bool]:
        """   
        Applys the boundary condition to coordinates.

        Further uses the center of mass for the reflection method
        """
        coords = np.atleast_2d(coordinates.copy())
        if self.is_inside(coords):
            return coords, True
        
        if method == "reject":
            return coordinates.copy(), False
        elif method == "reflect":
            return self._reflect_coordinates(coords)
        elif method == "wrap":
            #TODO Implement this someday
            raise NotImplementedError("Wrap method is not implemented yet")
        elif method == "rescale":
            # Rescales the coordinates along its COM vector to fit inside the box
            return self._rescale_coordinates(coords)
        
    
    def _rescale_coordinates(self, 
                             coordinates: np.ndarray,
                             scale_factor: float = 0.5) -> Tuple[np.ndarray, bool]:
        """ 
        Rescales the coordinates along its COM vector to fit inside the box

        Algorithm:
        1. Compute the geometric center of the molecule
        2. Rescale the COM so that r_new = r_old * (R - epsilon) / r_old

        Args:
            coordinates: Array of shape (N, 3) representing atomic coordinates
            scale_factor: Percentage factor of box radius to rescale to (default is 0.5, meaning 50% inside the box)
        """
        coords = coordinates.copy()
        centroid = np.mean(coords, axis=0)
        com_vec = centroid - self.center
        dist = np.linalg.norm(com_vec)
        if dist == 0:
            # If the COM is exactly at the center, we can just place it at a random point inside
            return self.get_random_position(n_points=len(coords)), True
        direction = com_vec / dist
        if self.box_type == "sphere":
            r_epsilon = self.radius * scale_factor
            new_dist = min(dist, self.radius - r_epsilon)
            new_centroid = self.center + direction * new_dist
            translation = new_centroid - centroid
            coords += translation
            return coords, self.is_inside(coords)
        elif self.box_type == "cube":
            # For cube we need to rescale each coordinate separately
            half_dims = self.box_dimensions / 2.0
            rel_coords = coords - self.center
            for i in range(3):
                if abs(rel_coords[:, i]).max() > half_dims[i] * scale_factor:
                    scale = (half_dims[i] * scale_factor) / abs(rel_coords[:, i]).max()
                    rel_coords[:, i] *= scale
            new_coords = self.center + rel_coords
            return new_coords, self.is_inside(new_coords)
        else:
            raise ValueError("box_type must be 'sphere' or 'cube'")
        


    def _reflect_coordinates(self, coordinates: np.ndarray) -> Tuple[np.ndarray, bool]:
        """   
        Reflects the molecule into the box as a rigid body

        Algorithm:
        1. Compute the geometric center of the molecule
        2. If centroid is outside relfect back
        3. Translate entire molecule by the same displacement
        4. Repeat until inside or max iterations reached
        """
        coords = coordinates.copy()
        centroid = np.mean(coords, axis=0)

        max_iterations = 10
        for _ in range(max_iterations):
            if self.is_inside(coords):
                return coords, True 
            
            # recompute centroid and displacement
            centroid = np.mean(coords, axis=0)
            new_centroid = self._reflect_point(centroid)
            translation = new_centroid - centroid
            coords += translation

        return coords, self.is_inside(coords)
        
        
    def _reflect_point(self, point: np.ndarray) -> np.ndarray:
        """  
        Reflects a single point back into the box
        """
        if self.box_type == "sphere":
            rel_pos = point - self.center
            dist = np.linalg.norm(rel_pos)
            if dist <= self.radius:
                return point.copy()
            
            # Relfect: new_dist = radius - (dist - radius) = 2*radius - dist
            direction = rel_pos / dist
            excess = dist - self.radius
            excess_mod = excess % (2 * self.radius)
            if excess_mod <= self.radius:
                new_dist = self.radius - excess_mod
            else:
                new_dist = excess_mod - self.radius
            return self.center + direction * new_dist

        elif self.box_type == "cube":
            raise NotImplementedError("Reflection for cube is not implemented yet")
    
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


def _run_box_self_tests() -> None:
    print("\n" + "=" * 70)
    print("Simulation Box Self Tests")
    print("=" * 70)

    sphere = SimulationBox(box_type="sphere", radius=5.0, center=np.array([0.0, 0.0, 0.0]))

    # Generate coords outside the sphere
    coords_outside = np.array([
        [7.0, 0.0, 0.0],
        [7.5, 0.2, -0.1],
        [6.8, -0.3, 0.1],
    ])
    print("\n [Sphere] Outside check:", sphere.is_inside(coords_outside))  # Should be False
    
    # Test reject
    rejected_coords, ok_reject = sphere.apply_boundary_conditions(coords_outside, method="reject")
    print(" [Sphere] Reject method:", ok_reject)  # Should be False
    print(" [Sphere]Rejected coords unchanged:", np.allclose(rejected_coords, coords_outside))  # Should be True

    # Test Reflect should move rigid body inside
    reflected_coords, ok_reflect = sphere.apply_boundary_conditions(coords_outside, method="reflect")
    print(" [Sphere] Reflect method:", ok_reflect)  # Should be True
    print(" [Sphere] Reflected coords inside:", sphere.is_inside(reflected_coords))

    # Rigid-body check: pairwise distances should be preserved
    def pairwise_distances(c):
        n = len(c)
        d = []
        for i in range(n):
            for j in range(i + 1, n):
                d.append(np.linalg.norm(c[i] - c[j]))
        return np.array(d)
    d_before = pairwise_distances(coords_outside)
    d_after = pairwise_distances(reflected_coords)
    print(" [Sphere] Rigid-body preserved:", np.allclose(d_before, d_after, atol=1e-8))

    # Test rescaling should move the COM inside
    rescaled_coords, ok_rescale = sphere.apply_boundary_conditions(coords_outside, method="rescale")
    print(" [Sphere] Rescale method:", ok_rescale)  # Should be True
    print(" [Sphere] Rescaled coords inside:", sphere.is_inside(rescaled_coords))

    # Test if rescale preserves rigid_body
    d_rescaled = pairwise_distances(rescaled_coords)
    print(" [Sphere] Rescale rigid-body preserved:", np.allclose(d_before, d_rescaled, atol=1e-8))

    # Random points inside the sphere should all be inside
    rnd = sphere.get_random_position(n_points=1000)
    inside_random = sphere.is_inside(rnd)
    print(" [Sphere] Random points inside check:", np.all(inside_random))  #

    # Make a figure of how the reflection works as an X-Y projection
    
    # First make coords that are outside in x but inside in y and z
    coords_outside_x_y = np.array([
        [7.0, 0.0, 0.0],
        [7.5, 0.2, -0.1],
        [6.8, -0.3, 0.1],
    ])
    reflected_coords_xy, _ = sphere.apply_boundary_conditions(coords_outside_x_y, method="reflect")
    # Add rescaled coords for comparison
    rescaled_coords_xy, _ = sphere.apply_boundary_conditions(coords_outside_x_y, method="rescale")


    # Plot the original and reflected points
    plt.figure(figsize=(6, 6))
    plt.scatter(coords_outside_x_y[:, 0], coords_outside_x_y[:, 1], color='red', label='Original (Outside)', s=100)
    plt.scatter(reflected_coords_xy[:, 0], reflected_coords_xy[:, 1], color='blue', label='Reflected (Inside)', s=100)
    plt.scatter(rescaled_coords_xy[:, 0], rescaled_coords_xy[:, 1], color='green', label='Rescaled (Inside)', s=100, marker='x')
    # Plot the circle representing the sphere boundary
    theta = np.linspace(0, 2 * np.pi, 100)
    x_circle = sphere.radius * np.cos(theta) + sphere.center[0]
    y_circle = sphere.radius * np.sin(theta) + sphere.center[1]
    plt.plot(x_circle, y_circle, color='gray', linestyle='--', label='Sphere Boundary')
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Reflection of Points into Sphere")
    plt.legend()
    plt.grid(True)
    plt.axis('equal')
    plt.savefig("figures/simulation_box_reflection_test.png", dpi=300)

if __name__ == "__main__":


    # Compute the Box Size for several molecules and print the results
    from modules.molecule_class import Molecule
    from modules.calculator import EnergyEvaluator

    h2o_dimer_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/water_clusters/h2o_dimer.xyz"
    h2o_trimer_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/water_clusters/h2o_timer.xyz"
    h2o_tetramer_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/water_clusters/h2o_4.xyz"
    h2o_pentamer_file = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/water_clusters/tip4p_literature/TIP4P-5.xyz"
    co2_dimer_file  = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/co2_co2.xyz"
    formic_acid_dimer = "/media/storage_6/lme/master_thesis/cluster_generation/test_molecules/fad.xyz"

    # Load molecules and compute box sizes
    for file in [h2o_dimer_file, h2o_trimer_file, h2o_tetramer_file, h2o_pentamer_file, co2_dimer_file, formic_acid_dimer]:
        with open(file, "r") as f:
            content = f.read()
        mol = Molecule.from_xyz(content)
        box = SimulationBox.from_covalent_radii(covalent_radii=mol.covalent_radii, n_atoms=len(mol.atom_labels), box_type="sphere", scale_factor=1.0)
        print(f"File: {file.split('/')[-1]}, Atoms: {mol.atom_labels}, Box Radius: {box.radius:.2f} Å, Volume: {box.get_volume():.2f} Å³")
    print("\n" + "=" * 70)
    # Load molecules and compute box sizes using vdW radii
    for file in [h2o_dimer_file, h2o_trimer_file, h2o_tetramer_file, h2o_pentamer_file, co2_dimer_file, formic_acid_dimer]:
        with open(file, "r") as f:
            content = f.read()
        mol = Molecule.from_xyz(content)
        box = SimulationBox.from_vdw_radii(vdw_radii=mol.vdw_radii, n_atoms=len(mol.atom_labels), box_type="sphere", eta_factor=0.4)
        print(f"File: {file.split('/')[-1]}, Atoms: {mol.atom_labels}, Box Radius: {box.radius:.2f} Å, Volume: {box.get_volume():.2f} Å³")


    _run_box_self_tests()
