import numpy as np 
from typing import List, Tuple, Optional, Dict, Any
import random

from modules.transformations import Transformation
from modules.geometry import GeometryOps, Quaternion, ReferenceFrame
from scipy.spatial.distance import pdist, squareform


class MolecularOperators:
    """   
    Base class for molecular geometry operators.
    Uses Transformation and GeometryOps for efficient operations
    """
    def __init__(self, simulation_box: Optional["SimulationBox"] = None):
        """
        Initialize Molecular Operators

        Args:
            simulation_box: Optional SimulationBox object to constrain operations
        """
        self.simulation_box = simulation_box
        self.transformer = Transformation()

    def _apply_box_constraints(self, molecule: "Molecule", original: "Molecule") -> "Molecule":
        """
        Apply simulation box constraints to the molecule after transformation

        Args:
            molecule: Transformed Molecule object
            original: Original Molecule object before transformation (returned if constraints violated)
        """ 
        if self.simulation_box is not None:
            new_coords, is_valid = self.simulation_box.apply_boundary_conditions(
                molecule.coordinates, method="reflect"
            )
            if not is_valid:
                return original.copy()  
            molecule.coordinates = new_coords
        return molecule
    
    def _has_inter_submolecule_clashes(self,
                                       molecule: "Molecule",
                                       submolecule_indices: List[List[int]],
                                       scale=0.5) -> bool:
        """  
        Compute wether any inter-submolecule atom pairs are closer than the clash threshold

        Algorithm:
        1. Get the VDW radii of all atoms and apply scaling factor to find clash threshold
        """ 
        coords = molecule.coordinates
        radii = molecule.vdw_radii
        n_atoms = len(coords)
        if n_atoms < 2:
            return False  # No clashes possible with 1 or 0 atoms
        radii = getattr(molecule, "vdw_radii", None)
        if radii is None or len(radii) != n_atoms:
            # Fallback
            radii = molecule.covalent_radii * 1.5
        
        # Map atom -> fragment id
        frag_id = np.full(n_atoms, -1, dtype=int)
        for frag, idxs in enumerate(submolecule_indices):
            frag_id[np.asarray(idxs, dtype=int)] = frag
        
        dists = squareform(pdist(coords))
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                if frag_id[i] != frag_id[j]:  # Only check inter-fragment pairs
                    threshold = scale * (radii[i] + radii[j])
                    if dists[i, j] < threshold:
                        return True
        return False
        



    def _compute_scaling_factor(
            self,
            molecule: "Molecule",
            operator_type: str = "local"
    ) -> float:
        """  
        Computes a normalized scaling factor
        """ 


        if self.simulation_box is None:
            return 1.0
        coords = molecule.coordinates
        from scipy.spatial.distance import pdist

        if len(coords) > 1:
            pairwise_distances = pdist(coords)
            mol_diameter = np.max(pairwise_distances)
        else:
            mol_diameter = 2.0 * molecule.covalent_radii.max()
        # Add average covalent radius to account for size of atoms
        avg_radius = np.mean(molecule.covalent_radii)
        effective_mol_size = mol_diameter + 2 * avg_radius


        # Get characteristic length
        if self.simulation_box.box_type == "cube":
            box_length = float(np.max(self.simulation_box.box_dimensions))
        elif self.simulation_box.box_type == "sphere":
            box_length = self.simulation_box.radius * 2
        else:
            raise ValueError(f"Unknown box type '{self.simulation_box.box_type}'")

        free_space = max(box_length - effective_mol_size, 0.1)
        # Room available per molecule-width: more room -> bigger allowed moves,
        # a cramped box -> smaller ones (previously inverted, which let tight
        # boxes take the largest moves and clash constantly).
        room_ratio = free_space / effective_mol_size

        if operator_type == "local":
            base_scale = 0.2
            max_scale = 0.8
            scale = base_scale + (max_scale - base_scale) * (1 - np.exp(-room_ratio))
        elif operator_type == "nonlocal":
            base_scale = 0.3
            max_scale = 1.0
            scale = base_scale + (max_scale - base_scale) * (1 - np.exp(-room_ratio * 0.5))
        else:
            scale = 1.0
        return np.clip(scale, base_scale, max_scale)

                                         
    @staticmethod
    def _random_unit_vetor() -> np.ndarray:
        """   
        Generates a random unit vector on the unit sphere
        Uses spherical coordinates for uniform distribution

        Returns:
            Random unit vector (3,)
        """
        phi = np.random.uniform(0, 2 * np.pi)
        costheta = np.random.uniform(-1, 1)
        theta = np.arccos(costheta)

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        return np.array([x, y, z])

    @staticmethod
    def _select_random_submolecule(submolecule_indices: List[List[int]]) -> List[int]:
        """   
        selects a random submolecule from the list
        """
        return random.choice(submolecule_indices)
    
    @staticmethod 
    def _select_two_submolecules(submolecule_indices: List[List[int]]) -> Tuple[List[int], List[int]]:
        """   
        Selects two random submolecules from the list
        """
        if len(submolecule_indices) < 2:
            raise ValueError("At least two submolecules must be provided to select two")
        return random.sample(submolecule_indices, 2)

class LocalOperators(MolecularOperators):
    """   
    Local Operators for small perturbations of the molecule
    """ 
    def random_displacement_submolecule(self,
                                        molecule: "Molecule",
                                        submolecule_indices: List[List[int]],
                                        delta_range: Tuple[float, float] = (-0.5, 0.5),
                                        adaptive: bool = True) -> "Molecule":
        """  
        Randomly displaces a submolecule by a random vector

        Args:
            molecule: Molecule object to modify
            submolecule_indices: List of lists of atom indices defining submolecules
            delta_range: Range for random displacement in Angstrom (default is (-0.3, 0.3))
            adaptive: If True, scale range by box/molecule size ratio
        """
        mol_copy = molecule.copy()
        n_atoms = len(mol_copy.coordinates)

        # select random submolecule
        submol_idx = self._select_random_submolecule(submolecule_indices)
        
        # Validate
        if any(idx >= n_atoms for idx in submol_idx):
            raise ValueError("Submolecule indices must be within the range of molecule atoms")
        
        # Apply adaptive scaling
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="local")
            effective_range = (delta_range[0] * scale, delta_range[1] * scale)
        else:
            effective_range = delta_range
        
        # Generate random displacement vector
        delta_value = random.uniform(*effective_range)
        direction = self._random_unit_vetor()
        displacement = direction * delta_value
        # Apply displacement to selected submolecule atoms
        mol_copy.coordinates[submol_idx] += displacement
        return self._apply_box_constraints(mol_copy, molecule)

    def random_rotation_submolecule(
            self,
            molecule: "Molecule",
            submolecule_indices: List[List[int]],
            angle_range: Tuple[float, float] = (-10, 10),
            adaptive: bool = True) -> "Molecule":
        """  
        Randomly rotate a submolecule around its center of mass using a random axis and a small angle.
        Analogous to the nonlocal twist but with local-scale angle

        Args:
            molecule: Molecule object to modify
            submolecule_indices: List of lists of atom indices defining submolecules
            angle_range: Range for random rotation angle in degrees (default is (-30, 30))
            adaptive: If True, scale angle range by box/molecule size ratio
        """
        mol_copy = molecule.copy()

        submol_idx = self._select_random_submolecule(submolecule_indices)
        submol_coords = mol_copy.coordinates[submol_idx]
        submol_masses = mol_copy.masses[submol_idx] 
        center = GeometryOps.center_of_mass(submol_coords, submol_masses)

        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="local")
            effective_range = (angle_range[0] * scale, angle_range[1] * scale)
        else:
            effective_range = angle_range

        angle = np.radians(random.uniform(*effective_range))
        axis = self._random_unit_vetor()

        R = GeometryOps.rotation_matrix_rodrigues(axis, angle)
        coords_centered = mol_copy.coordinates[submol_idx] - center
        coords_rotated = coords_centered @ R.T
        mol_copy.coordinates[submol_idx] = coords_rotated + center
        return self._apply_box_constraints(mol_copy, molecule)
    
    def local_twist_submolecule(
            self,
            molecule: "Molecule",
            submolecule_indices: List[List[int]],
            angle_range: Tuple[float, float] = (-15, 15),
            adaptive: bool = True) -> "Molecule":
        """  
        Local Twist: Rotate one submolecule around the axis connecting two submolecule
        centers of mass by a small angle

        This is physically motivated - the twist axis is the intermolecular COM-to-COM direction, so it
        presevers the intermolecular distance while changing the relative orientation

        Algorithm:
            1. Select two submolecules (reference + rotating)
            2. Compute the COM-to-COM axis
            3. Rotate the second submolecule around that axis by a small angle

        Args:
            molecule: Molecule object to modify
            submolecule_indices: List of lists of atom indices defining submolecules
            angle_range: Range for random rotation angle in degrees (default is (-15, 15))
            adaptive: If True, scale angle range by box/molecule size ratio
        """
        mol_copy = molecule.copy()
        ref_submol, rot_submol = self._select_two_submolecules(submolecule_indices)
        ref_center = GeometryOps.center_of_mass(mol_copy.coordinates[ref_submol], mol_copy.masses[ref_submol])
        rot_center = GeometryOps.center_of_mass(mol_copy.coordinates[rot_submol], mol_copy.masses[rot_submol])
        
        # Twist axis: COM to COM
        twist_axis = rot_center - ref_center
        if np.linalg.norm(twist_axis) < 1e-8:
            return mol_copy # Cannot define twist because submolecules overlap
        axis = GeometryOps.normalize_vector(twist_axis)

        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="local")
            effective_range = (angle_range[0] * scale, angle_range[1] * scale)
        else:
            effective_range = angle_range
        angle = np.radians(random.uniform(*effective_range))
        R = GeometryOps.rotation_matrix_rodrigues(axis, angle)
        coords_centered = mol_copy.coordinates[rot_submol] - rot_center
        coords_rotated = coords_centered @ R.T
        mol_copy.coordinates[rot_submol] = coords_rotated + rot_center
        return self._apply_box_constraints(mol_copy, molecule)

    def correlated_displacement(
        self,
        molecule: "Molecule",
        submolecule_indices: List[List[int]],
        delta_range: Tuple[float, float] = (-0.3, 0.3),
        adaptive: bool = True) -> "Molecule":
        """ 
        Displace two or more submolecules simultaneously in a correlated manner

        Args:
            molecule: Molecule object to modify
            submolecule_indices: List of lists of atom indices defining submolecules
            delta_range: Range for random displacement in Angstrom (default is (-0.3, 0.3))
            adaptive: If True, scale range by box/molecule size ratio
        """
        mol_copy = molecule.copy()

        n_select = random.randint(2, len(submolecule_indices))
        selected = random.sample(submolecule_indices, n_select)
        #print(f"Selected {n_select} submolecules for correlated displacement: {selected}")  
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="local")
            effective_range = (delta_range[0] * scale, delta_range[1] * scale)
        else:
            effective_range = delta_range
        # for each selected submolecule sample a random displacement and apply to the submolecule
        random_displacements = []
        for submol_idx in selected:
            direction = self._random_unit_vetor()
            magnitude = random.uniform(*effective_range)
            displacement = direction * magnitude
            random_displacements.append(displacement)
            mol_copy.coordinates[submol_idx] += displacement
        
        # Check for inter-submolecule clashes and reject if any are found
        if self._has_inter_submolecule_clashes(mol_copy, submolecule_indices):
            return molecule.copy()  # Reject move by returning original
        return self._apply_box_constraints(mol_copy, molecule)

    def small_principal_axis_rotation(
            self,
            molecule: "Molecule",
            submolecule_indices: List[List[int]],
            angle_range: Tuple[float, float] = (-15,15),
            adaptive: bool = True) -> "Molecule":
        """  
        Rotate a single submolecule around one of its own principal axes
        computed from the inertia tensor by a small angle

        Algorithm:
            1. Select random submolecule
            2. Compute its reference frame and principal axes
            3. Select random principal axis
            4. Rotate submolecule around that axis by a small angle

        Args:            
            molecule: Molecule object to modify
            submolecule_indices: List of lists of atom indices defining submolecules
            angle_range: Range for random rotation angle in degrees (default is (-15, 15))
            adaptive: If True, scale angle range by box/molecule size ratio
        """
        mol_copy = molecule.copy()
        submol_idx = self._select_random_submolecule(submolecule_indices)
        coords = mol_copy.coordinates[submol_idx]
        masses = mol_copy.masses[submol_idx]
        center = GeometryOps.center_of_mass(coords, masses)

        ref_frame = self.transformer.compute_reference_frame(coords, masses, center=center)
        principal_axes = [ref_frame.x_axis, ref_frame.y_axis, ref_frame.z_axis]
        axis = random.choice(principal_axes)

        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="local")
            effective_range = (angle_range[0] * scale, angle_range[1] * scale)
        else:
            effective_range = angle_range
        angle = np.radians(random.uniform(*effective_range))
        R = GeometryOps.rotation_matrix_rodrigues(axis, angle)
        coords_centered = mol_copy.coordinates[submol_idx] - center
        coords_rotated = coords_centered @ R.T
        mol_copy.coordinates[submol_idx] = coords_rotated + center
        return self._apply_box_constraints(mol_copy, molecule)

    def diagnose_perturbations(
        self,
        molecule: "Molecule",
        submolecule_indices: List[List[int]],
        n_samples: int = 200,
        adaptive: bool = True
    ) -> Dict:
        """
        Diagnose the actual perturbation magnitudes produced by each local operator.
        
        Applies each operator n_samples times and reports statistics on:
          - RMSD between original and perturbed structures
          - Max single-atom displacement
          - COM displacement of the perturbed submolecule
          - Effective scaling factor and parameter ranges
          
        Args:
            molecule: Molecule to test operators on
            submolecule_indices: Submolecule index lists
            n_samples: Number of random applications per operator
            adaptive: Whether to use adaptive scaling
            
        Returns:
            Dict with operator names as keys, each containing displacement statistics
        """
        from scipy.spatial.distance import pdist
            
        # First report the raw scaling info
        scale = self._compute_scaling_factor(molecule, operator_type="local")
        coords = molecule.coordinates
        if len(coords) > 1:
            mol_diameter = np.max(pdist(coords))
        else:
            mol_diameter = 2.0 * molecule.covalent_radii.max()
        avg_radius = np.mean(molecule.covalent_radii)
        effective_mol_size = mol_diameter + 2 * avg_radius
        
        if self.simulation_box is not None:
            if self.simulation_box.box_type == "sphere":
                box_length = self.simulation_box.radius * 2
            else:
                box_length = self.simulation_box.box_dimensions
            free_space = max(box_length - effective_mol_size, 0.1)
            size_ratio = effective_mol_size / free_space
        else:
            box_length = None
            free_space = None
            size_ratio = None
            
        print("=" * 70)
        print("LOCAL OPERATOR PERTURBATION DIAGNOSTICS")
        print("=" * 70)
        print(f"  Molecule diameter:     {mol_diameter:.3f} Å")
        print(f"  Effective mol size:    {effective_mol_size:.3f} Å")
        if box_length is not None:
            print(f"  Box length (diameter): {box_length:.3f} Å")
            print(f"  Free space:            {free_space:.3f} Å")
            print(f"  Size ratio:            {size_ratio:.3f}")
        print(f"  Adaptive scale factor: {scale:.4f}")
        print(f"  Adaptive enabled:      {adaptive}")
        print()
        
        # Define operators with their default ranges and the effective ranges
        operators_info = {
            "random_displacement": {
                "func": self.random_displacement_submolecule,
                "default_range": (-0.3, 0.3),
                "unit": "Å",
                "range_type": "displacement"
            },
            "random_rotation": {
                "func": self.random_rotation_submolecule,
                "default_range": (-30, 30),
                "unit": "°",
                "range_type": "angle"
            },
            "local_twist": {
                "func": self.local_twist_submolecule,
                "default_range": (-15, 15),
                "unit": "°",
                "range_type": "angle"
            },
            "correlated_displacement": {
                "func": self.correlated_displacement,
                "default_range": (-0.3, 0.3),
                "unit": "Å",
                "range_type": "displacement"
            },
            "small_principal_axis_rotation": {
                "func": self.small_principal_axis_rotation,
                "default_range": (-15, 15),
                "unit": "°",
                "range_type": "angle"
            },
        }
        
        results = {}
        
        for op_name, info in operators_info.items():
            default_range = info["default_range"]
            if adaptive:
                effective_range = (default_range[0] * scale, default_range[1] * scale)
            else:
                effective_range = default_range
                
            rmsds = []
            max_displacements = []
            box_rejections = 0
            
            for _ in range(n_samples):
                original_coords = molecule.coordinates.copy()
                perturbed = info["func"](molecule, submolecule_indices, adaptive=adaptive)
                
                # Check if box rejected (returned original)
                if np.allclose(perturbed.coordinates, original_coords):
                    box_rejections += 1
                    continue
                
                diff = perturbed.coordinates - original_coords
                per_atom_disp = np.linalg.norm(diff, axis=1)
                rmsd = np.sqrt(np.mean(per_atom_disp**2))
                
                rmsds.append(rmsd)
                max_displacements.append(np.max(per_atom_disp))
            
            rmsds = np.array(rmsds) if rmsds else np.array([0.0])
            max_displacements = np.array(max_displacements) if max_displacements else np.array([0.0])
            
            results[op_name] = {
                "effective_range": effective_range,
                "rmsd_mean": np.mean(rmsds),
                "rmsd_std": np.std(rmsds),
                "max_disp_mean": np.mean(max_displacements),
                "max_disp_std": np.std(max_displacements),
                "box_rejection_rate": box_rejections / n_samples,
            }
            
            print(f"  {op_name}")
            print(f"    Default range:    ({default_range[0]}, {default_range[1]}) {info['unit']}")
            print(f"    Effective range:  ({effective_range[0]:.4f}, {effective_range[1]:.4f}) {info['unit']}")
            print(f"    RMSD:             {np.mean(rmsds):.4f} ± {np.std(rmsds):.4f} Å")
            print(f"    Max atom disp:    {np.mean(max_displacements):.4f} ± {np.std(max_displacements):.4f} Å")
            print(f"    Box rejections:   {box_rejections}/{n_samples} ({box_rejections/n_samples*100:.1f}%)")
            print()
        
        print("=" * 70)
        print("REFERENCE: Typical chemical scales")
        print(f"  O-H bond length:        0.96 Å")
        print(f"  H-bond (O···H):         1.8-2.0 Å")
        print(f"  O···O in water dimer:    2.9 Å")
        print(f"  Useful displacement:     0.05-0.5 Å")
        print(f"  Useful rotation:         5-45°")
        print("=" * 70)
        
        return results


class NonLocalOperators(MolecularOperators):
    """   
    Non-local operators for large structural perturbations of the molecule
    These operators make significant changes to explore different basins of the PES
    """

    def twist_operator(self,
                       molecule: "Molecule",
                       submolecule_indices: List[List[int]],
                       rotation_angle_range: Tuple[float, float] = (-20, 20),
                       adaptive: bool = True) -> "Molecule":
        """
        Twist operator: Rotate one submolecule relative to others
        
        Args:
            molecule: Input Molecule object
            submolecule_indices: List of lists of atom indices defining submolecules
            rotation_angle_range: Range for random rotation angle in degrees
            adaptive: If True, scale angle range by box/molecule size ratio
        """ 
        mol_copy = molecule.copy()

        ref_submol, rot_submol = self._select_two_submolecules(submolecule_indices)

        # Calculate reference_frame
        ref_coords = mol_copy.coordinates[ref_submol]
        ref_masses = mol_copy.masses[ref_submol]
        ref_center = GeometryOps.center_of_mass(ref_coords, ref_masses)

        # Apply adaptive scaling
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="nonlocal")
            effective_range = (rotation_angle_range[0] * scale, rotation_angle_range[1] * scale)
        else:
            effective_range = rotation_angle_range

        angle = np.radians(random.uniform(*effective_range))
        axis = self._random_unit_vetor()
        R = GeometryOps.rotation_matrix_rodrigues(axis, angle)

        coords_centered = mol_copy.coordinates[rot_submol] - ref_center
        coords_rotated = coords_centered @ R.T
        mol_copy.coordinates[rot_submol] = coords_rotated + ref_center


        return self._apply_box_constraints(mol_copy, molecule)
    
    def large_displacement(self,
                           molecule: "Molecule",
                           submolecule_indices: List[List[int]],
                           displacement_range: Tuple[float, float] = (1.0, 3.0),
                           adaptive: bool = True) -> "Molecule":
        """   
        Large displacement operation

        Args:
            molecule: Input Molecule object
            submolecule_indices: List of lists of atom indices defining submolecules
            displacement_range: Range for random displacement in Angstrom (default is (1.0, 3.0))
            adaptive: If True, scale range by box/molecule size ratio
        """
        mol_copy = molecule.copy()

        # Select random submolecule
        submol_idx = self._select_random_submolecule(submolecule_indices)

        # Apply adaptive scaling
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="nonlocal")
            effective_range = (displacement_range[0] * scale, displacement_range[1] * scale)
        else:
            effective_range = displacement_range

        magnitude = random.uniform(*effective_range)
        direction = self._random_unit_vetor()
        displacement = direction * magnitude

        mol_copy.coordinates[submol_idx] += displacement
        return self._apply_box_constraints(mol_copy, molecule)

    def com_com_approach_operator(self,
                                  molecule: "Molecule",
                                  submolecule_indices: List[List[int]],
                                  approach_distance_range: Tuple[float, float] = (0.05, 0.3),
                                  adaptive: bool = True) -> "Molecule":
        """   
        Implement a com_com_approach operator:

        This moves the com_com distance between two submolecules closer by a random amount within the specified range, while keeping the direction fixed.

        Algorithm:
        1. Select two submolecules
        2. Compute their centers of mass and the vector between them
        3. Compute the current distance and the desired approach distance
        4. Move the second submolecule towards the first along the com_com vector by the difference between current and desired distance

        Args:
            molecule: Input Molecule object
            submolecule_indices: List of lists of atom indices defining submolecules
            approach_distance_range: Range for random approach distance in Angstrom (default is (0.1, 0.8))
            adaptive: If True, scale range by box/molecule size ratio
        """
        mol_copy = molecule.copy()
        submol_1, submol_2 = self._select_two_submolecules(submolecule_indices)

        # Compute com of submol1
        coords_1 = mol_copy.coordinates[submol_1]
        masses_1 = mol_copy.masses[submol_1]
        center_1 = GeometryOps.center_of_mass(coords_1, masses_1)

        # compute com of submol2
        coords_2 = mol_copy.coordinates[submol_2]
        masses_2 = mol_copy.masses[submol_2]
        center_2 = GeometryOps.center_of_mass(coords_2, masses_2)

        com_vector = center_2 - center_1
        current_distance = np.linalg.norm(com_vector)
        if current_distance < 1e-8:
            return mol_copy  # Cannot define approach because submolecules overlap
        
        # Find direction vector
        direction = com_vector / current_distance
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="nonlocal")
            effective_range = (approach_distance_range[0] * scale, approach_distance_range[1] * scale)
        else:
            effective_range = approach_distance_range
        
        delta = random.uniform(*effective_range) 
        new_distance = max(current_distance - delta, 0.1)  # Avoid collapse
        displacement_magnitude = current_distance - new_distance

        displacement_vector = direction * displacement_magnitude

        # Move submol2 towards submol1
        mol_copy.coordinates[submol_2] -= displacement_vector
        return self._apply_box_constraints(mol_copy, molecule)
    
    def com_com_separation_operator(self,
                                      molecule: "Molecule",
                                      submolecule_indices: List[List[int]],
                                      separation_distance_range: Tuple[float, float] = (0.05, 0.1),
                                      adaptive: bool = True) -> "Molecule":
          """   
          Implement a com_com_separation operator:
    
          This moves the com_com distance between two submolecules farther by a random amount within the specified range, while keeping the direction fixed.
    
          Algorithm:
          1. Select two submolecules
          2. Compute their centers of mass and the vector between them
          3. Compute the current distance and the desired separation distance
          4. Move the second submolecule away from the first along the com_com vector by the difference between desired and current distance
    
          Args:
                molecule: Input Molecule object
                submolecule_indices: List of lists of atom indices defining submolecules
                separation_distance_range: Range for random separation distance in Angstrom (default is (0.1, 0.8))
                adaptive: If True, scale range by box/molecule size ratio
          """
          mol_copy = molecule.copy()
          submol_1, submol_2 = self._select_two_submolecules(submolecule_indices)
    
          # Compute com of submol1
          coords_1 = mol_copy.coordinates[submol_1]
          masses_1 = mol_copy.masses[submol_1]
          center_1 = GeometryOps.center_of_mass(coords_1, masses_1)
    
          # compute com of submol2
          coords_2 = mol_copy.coordinates[submol_2]
          masses_2 = mol_copy.masses[submol_2]
          center_2 = GeometryOps.center_of_mass(coords_2, masses_2)
    
          com_vector = center_2 - center_1
          current_distance = np.linalg.norm(com_vector)
          if current_distance < 1e-8:
                return mol_copy  # Cannot define separation because submolecules overlap
          
          # Find direction vector
          direction = com_vector / current_distance
          if adaptive:
                scale = self._compute_scaling_factor(molecule, operator_type="nonlocal")
                effective_range = (separation_distance_range[0] * scale, separation_distance_range[1] * scale)
          else:
                effective_range = separation_distance_range
          
          delta = random.uniform(*effective_range)
          new_distance = current_distance + delta
          displacement_magnitude = new_distance - current_distance

          displacement_vector = direction * displacement_magnitude

          # Move submol2 away from submol1
          mol_copy.coordinates[submol_2] += displacement_vector
          return self._apply_box_constraints(mol_copy, molecule)
          


    def mirror_operator(self,
                        molecule: "Molecule",
                        submolecule_indices: List[List[int]]) -> "Molecule":
        """  
        Mirror operator: Reflect submolecule through a plane
        Uses reference frame construction and reflection matrix

        Note: This operator is not scaled as reflection is a discrete transformation

        Algorithm:
        1. Construct reference frame from the submolecule to mirror
        2. Select coordinate plane (XY, YZ or XZ) 
        3. Build Reflection matrix H = I - 2 * n * n^T 
        4. Apply reflection to the submolecule coordinates

        Args:
            molecule: Input Molecule object
            submolecule_indices: List of lists of atom indices defining submolecules
        """
        mol_copy = molecule.copy()

        ref_submol, mir_submol = self._select_two_submolecules(submolecule_indices)

        
        ref_coords = mol_copy.coordinates[ref_submol].copy()
        ref_masses = mol_copy.masses[ref_submol].copy()
        ref_center = GeometryOps.center_of_mass(ref_coords, ref_masses)
        ref_frame = self.transformer.compute_reference_frame(ref_coords, 
                                                             ref_masses, 
                                                             center=ref_center)
        




        plane = np.random.choice(['XY', 'YZ', 'XZ'])
        if plane == 'XY':
            normal_vector = ref_frame.z_axis
        elif plane == 'YZ':
            normal_vector = ref_frame.x_axis 
        else: # XZ
            normal_vector = ref_frame.y_axis
        n = normal_vector / np.linalg.norm(normal_vector)
        H = np.eye(3) - 2 * np.outer(n, n)
        coords_centered = mol_copy.coordinates[mir_submol] - ref_frame.origin
        coords_reflected = coords_centered @ H.T
        mol_copy.coordinates[mir_submol] = coords_reflected + ref_frame.origin
        return self._apply_box_constraints(mol_copy, molecule)

    def random_so3_operator(self,
                            molecule: "Molecule",
                            submolecule_indices: List[List[int]]) -> "Molecule":
        """ 
        Applies a uniform Random SO(3) rotation to a submolecule around its own center of mass
        
        Note: This operator is not scaled as it samples uniformly over all rotations
        """
        mol_copy = molecule.copy()

        # Select submolecule
        submol_idx = self._select_random_submolecule(submolecule_indices)
        coords = mol_copy.coordinates[submol_idx]
        masses = mol_copy.masses[submol_idx]
        # Calc com
        com = GeometryOps.center_of_mass(coords, masses)
        # Sample Rotation Quaternion
        R = Quaternion.random_so3()
        coords_centered = coords - com
        coords_rotated = (R @ coords_centered.T).T

        mol_copy.coordinates[submol_idx] = coords_rotated + com
        return self._apply_box_constraints(mol_copy, molecule)

    def principal_axis_rotation_operator(self,
                                         molecule: "Molecule",
                                         submolecule_indices: List[List[int]],
                                         angle_range: Tuple[float, float] = (0, 360),
                                         adaptive: bool = True) -> "Molecule":
        """  
        Rotate a single submolecule around one of its own principal axes (from inertia tensor)

        Similar to random_so3 this operates on a single submolecule around its COM

        Algorithm:
        1. Select a random submolecule
        2. Compute its reference frame to get the principal axes
        3. Pick one random principal axes (x,y,z)
        4. Rotate the submolecule around that axis by a random angle

        Args:
            molecule: Input Molecule object
            submolecule_indices: List of lists of atom indices defining submolecules
            angle_range: Range for random rotation angle in degrees
            adaptive: If True, scale angle range by box/molecule size ratio
        """
        mol_copy = molecule.copy()
        submol_idx = self._select_random_submolecule(submolecule_indices)
        coords = mol_copy.coordinates[submol_idx]
        masses = mol_copy.masses[submol_idx]
        com = GeometryOps.center_of_mass(coords, masses)

        # Compute the reference frame to get the principal axes
        ref_frame = self.transformer.compute_reference_frame(coords, masses, center=com)
        
        # Pick a random principal axis
        axis = ref_frame.axes[:, np.random.choice(3)]
        
        # Apply adaptive scaling
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="nonlocal")
            effective_range = (angle_range[0] * scale, angle_range[1] * scale)
        else:
            effective_range = angle_range

        angle = np.radians(random.uniform(*effective_range))
        R = GeometryOps.rotation_matrix_rodrigues(axis, angle)

        coords_centered = coords - com 
        coords_rotated = coords_centered @ R.T
        mol_copy.coordinates[submol_idx] = coords_rotated + com
        
        return self._apply_box_constraints(mol_copy, molecule)
    




    def roto_reflection_operator(self,
                                 molecule: "Molecule",
                                 submolecule_indices: List[List[int]],
                                 angle_range: Tuple[float, float] = (0, 360),
                                 adaptive: bool = True) -> "Molecule":
        """  
        Roto-reflection operator: Combined rotation and reflection
        
        Args:
            molecule: Input Molecule object
            submolecule_indices: List of lists of atom indices defining submolecules
            angle_range: Range for random rotation angle in degrees
            adaptive: If True, scale angle range by box/molecule size ratio
        """
        mol_copy = molecule.copy()
        ref_submol, rot_ref_submol = self._select_two_submolecules(submolecule_indices)
        
        ref_coords = mol_copy.coordinates[ref_submol].copy()
        ref_masses = mol_copy.masses[ref_submol].copy()
        ref_center = GeometryOps.center_of_mass(ref_coords, ref_masses)
        ref_frame = self.transformer.compute_reference_frame(ref_coords,
                                                                ref_masses, 
                                                                center=ref_center)
    



        plane = np.random.choice(['XY', 'YZ', 'XZ'])
        if plane == 'XY':
            normal_vector = ref_frame.z_axis
        elif plane == 'YZ':
            normal_vector = ref_frame.x_axis 
        else: # XZ
            normal_vector = ref_frame.y_axis
        n = normal_vector / np.linalg.norm(normal_vector)
        H = np.eye(3) - 2 * np.outer(n, n)

        # Apply adaptive scaling
        if adaptive:
            scale = self._compute_scaling_factor(molecule, operator_type="nonlocal")
            effective_range = (angle_range[0] * scale, angle_range[1] * scale)
        else:
            effective_range = angle_range

        angle = np.radians(random.uniform(*effective_range))
        R = GeometryOps.rotation_matrix_rodrigues(n, angle)
        M = R @ H
        coords_centered = mol_copy.coordinates[rot_ref_submol] - ref_frame.origin
        coords_transformed = coords_centered @ M.T
        mol_copy.coordinates[rot_ref_submol] = coords_transformed + ref_frame.origin
        return self._apply_box_constraints(mol_copy, molecule)
    
    def exchange_operator(self,
                          molecule: "Molecule",
                          submolecule_indices: List[List[int]]) -> "Molecule":
        """  
        Exchange Operator: Swaps the positions of two submolecules
        Uses vectorized center of mass calculation
        
        Note: This operator is not scaled as it's a discrete exchange operation
        """
        mol_copy = molecule.copy()
        submol_1, submol_2 = self._select_two_submolecules(submolecule_indices)
        coords_1 = mol_copy.coordinates[submol_1]
        masses_1 = mol_copy.masses[submol_1]
        center_1 = GeometryOps.center_of_mass(coords_1, masses_1)
        coords_2 = mol_copy.coordinates[submol_2]
        masses_2 = mol_copy.masses[submol_2]
        center_2 = GeometryOps.center_of_mass(coords_2, masses_2)
        translation_1 = center_2 - center_1
        translation_2 = center_1 - center_2
        mol_copy.coordinates[submol_1] += translation_1
        mol_copy.coordinates[submol_2] += translation_2
        return self._apply_box_constraints(mol_copy, molecule)

if __name__ == "__main__":

    from molecule_class import Molecule

    h2o_2 = """ 
    6
    Coordinates from ORCA-job h2o_2 E -152.102897751726
    O           0.20131422818946     -0.13419863189991     -0.38118207664628
    H           1.10826926774619      0.03054527004232     -0.16898516572928
    H          -0.26386315635905     -0.06924940339370      0.43390806815981
    O           3.06513444516272      0.52224610599392      0.05902763481169
    H           3.25817767296947      1.29105665797100     -0.44996761253291
    H           3.66908154229119     -0.14039999871364     -0.23001884806303
    """
    mol = Molecule.from_xyz(h2o_2)
    submolecule_indices = [[0,1,2], [3,4,5]]
    # Check the clash function
    operators = NonLocalOperators()
    print(operators._has_inter_submolecule_clashes(mol, submolecule_indices))
    # Apply com_com_approach operator
    perturbed = operators.com_com_approach_operator(mol, submolecule_indices, approach_distance_range=(0.1, 0.5), adaptive=False)
    print("Original COM-COM distance:", np.linalg.norm(GeometryOps.center_of_mass(mol.coordinates[0:3], mol.masses[0:3]) - GeometryOps.center_of_mass(mol.coordinates[3:6], mol.masses[3:6])))
    print("Perturbed COM-COM distance:", np.linalg.norm(GeometryOps.center_of_mass(perturbed.coordinates[0:3], perturbed.masses[0:3]) - GeometryOps.center_of_mass(perturbed.coordinates[3:6], perturbed.masses[3:6])))
    # Apply com_com_separation operator
    perturbed_sep = operators.com_com_separation_operator(mol, submolecule_indices, separation_distance_range=(0.1, 0.5), adaptive=False)
    print("Separated COM-COM distance:", np.linalg.norm(GeometryOps.center_of_mass(perturbed_sep.coordinates[0:3], perturbed_sep.masses[0:3]) - GeometryOps.center_of_mass(perturbed_sep.coordinates[3:6], perturbed_sep.masses[3:6])))


