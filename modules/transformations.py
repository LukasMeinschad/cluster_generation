import numpy as np
from typing import List, Tuple, Optional, Union, Protocol, Dict
import matplotlib.pyplot as plt
from dataclasses import dataclass
import time
import random
import copy
from scipy.optimize import linear_sum_assignment

# Module import

from modules.geometry import ReferenceFrame, GeometryOps, Quaternion




class Transformation:
    """ 
    Optimizes Class for 3D-Transformation of Molecules
    """ 
    def __init__(self, name: str = "Transformation"):
        """  
        Initializes the Transformation class
        """
        self.name = name 

    # =========================== Center Point Calculations ===========================

    def get_center(
            self,
            molecule: Union["Molecule", "Submolecule"],
            method: str = "centroid"
        ) -> np.ndarray:
        """ 
        Determines the center of a molecules 

        Args:
            molecule: Molecule or Submolecule object
            method: Method to calculate center point ("centroid" or "com")
        Returns:
            np.ndarray: Center point coordinates
        """
        if method == "centroid":
            return GeometryOps.centroid(molecule.coordinates)
        if method == "com":
            return GeometryOps.center_of_mass(molecule.coordinates, molecule.masses)
        else: 
            raise ValueError(f"Unknown method '{method}'. Method must be 'centroid' or 'com'")
        
    def get_center_coords( 
            self, 
            coords: np.ndarray,
            masses: Optional[np.ndarray] = None,
            method: str = "centroid"
        ) -> np.ndarray:
        """   
        Determines the center point for a given set of coordinates
        """
        if method == "centroid":
            return GeometryOps.centroid(coords)
        if method == "com":
            if masses is None:
                raise ValueError("Masses must be provided to calculate center of mass")
            return GeometryOps.center_of_mass(coords, masses)
        else:
            raise ValueError(f"Unknown method '{method}'. Method must be 'centroid' or 'com'")
        
    # ========================= Basic Transformations ==========================

    def translate(
            self, 
            molecule: Union["Molecule", "Submolecule"],
            vector: Union[List[float], np.ndarray],
            in_place: bool = True
        ) -> np.ndarray:
        """ 
        Translates a molecule with a vector T_v(p) = p + v

        Args:
            molecule: Molecule or Submolecule object
            vector: Translation vector
            in_place: Whether to modify the molecule in place or return new coordinates 
        """
        vec = np.asarray(vector, dtype=np.float64)
        if vec.shape != (3,):
            raise ValueError("Translation vector must be of shape (3,)")
        
        new_coords = molecule.coordinates + vec

        if in_place:
            # Modify molecule coordinates in place
            molecule.coordinates = new_coords
            return None
        else:
            return new_coords
        
    def rotate(
            self, 
            molecule: Union["Molecule", "Submolecule"],
            axis: Union[List[float], np.ndarray],
            angle_degree: float,
            in_place: bool = True 
        ) -> Optional[np.ndarray]:
        """   
        Rotates the molecule around a given axis using Rodrigues rotation formula
        
        Args:
            molecule: Molecule or Submolecule object
            axis: Rotation axis
            angle_degree: Rotation angle in degrees
            in_place: Whether to modify the molecule in place or return new coordinates 
        """
        axis_array = np.asarray(axis, dtype=np.float64)
        angle_rad = np.radians(angle_degree)

        R = GeometryOps.rotation_matrix_rodrigues(axis_array, angle_rad)
        new_coords = molecule.coordinates @ R.T

        if in_place:
            molecule.coordinates = new_coords
            return None
        else:
            return new_coords
        
    def rotate_quaternion(
            self, 
            molecule: Union["Molecule", "Submolecule"],
            axis: Union[List[float], np.ndarray],
            angle_degree: float, 
            in_place: bool = True
        ) -> Optional[np.ndarray]:
        """   
        Rotates the molecule using quaternion rotation
        
        Args:
            molecule: Molecule or Submolecule object
            axis: Rotation axis
            angle_degree: Rotation angle in degrees
            in_place: Whether to modify the molecule in place or return new coordinates 
        """
        axis_array = np.asarray(axis, dtype=np.float64)
        angle_rad = np.radians(angle_degree)

        q = Quaternion.from_axis_angle(axis_array, angle_rad)
        q = q.normalize()

        
        # Still use rotation matrix because of linear algebra speed
        R = q.to_rotation_matrix()
        new_coords = molecule.coordinates @ R.T 

        if in_place:
            molecule.coordinates = new_coords
            return None
        else:
            return new_coords
        
    def center_molecule( 
            self, 
            molecule: Union["Molecule", "Submolecule"],
            method: str = "centroid",
            in_place: bool = True
        ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Centers a molecule to a center point

        Args:
            molecule: Molecule or Submolecule object
            method: Method to calculate center point ("centroid" or "com")
            in_place: Whether to modify the molecule in place or return new coordinates

        Returns:
            Tuple of (center, new_coords) if in_place is False
            Otherwise returns center_point only
        """
        center = self.get_center(molecule, method=method)
        new_coords = molecule.coordinates - center
        if in_place:
            molecule.coordinates = new_coords 
            return center, None
        else:
            return center, new_coords
        
    # ======================= Reference Frame Definitions =======================

    def compute_reference_frame( 
            self, 
            coords: np.ndarray,
            masses: np.ndarray,
            center: Optional[np.ndarray] = None
        ) -> ReferenceFrame:
        """
        Computes the Reference Frame based on the inertia tensor of a set of coordinates

        Args:
            coords: Coordinates of the molecule or submolecule
            masses: Masses of the atoms
            center: Optional center point to translate coordinates before computing inertia tensor
        Returns:
            ReferenceFrame object
        """
        if center is None:
            center = GeometryOps.center_of_mass(coords, masses)
        

        coords_new = coords.copy()
        coords_new -= center
        # Compute inertia tensor
        I = GeometryOps.inertia_tensor(coords_new, masses)
        eigenvalues, eigenvectors = np.linalg.eigh(I)

        # Eigenvectors are columns of the Reference Frame axes
        return ReferenceFrame(
            axes=eigenvectors,
            origin=center,
            eigenvalues=eigenvalues
        )
    
    def set_reference_frame(
            self,
            molecule: Union["Molecule", "Submolecule"],
            method: str = "centroid",
            print_info: bool = False,
        ) -> ReferenceFrame:
        """   
        Sets the reference frame to a center of a molecule and its principal axes

        Args: 
            molecule: Molecule or Submolecule object
            method: Method to calculate center point ("centroid" or "com")
            print_info: Whether to print information about the reference frame
        """
        if len(molecule.coordinates) < 3:
            raise ValueError("Molecule must have at least 3 atoms to define a reference frame")
        
        # Center Molecule
        center = self.get_center(molecule, method=method)
        molecule.coordinates -= center

        ref_frame = self.compute_reference_frame(
            molecule.coordinates,
            molecule.masses,
            center=np.zeros(3) # Already centered
        )
        if print_info:
            print(f"Molecule center: {center}")
            print(f"Inertia tensor eigenvalues: {ref_frame.eigenvalues}")
            print(f"Reference frame axes:")
            print(f"X: {ref_frame.x_axis}")
            print(f"Y: {ref_frame.y_axis}")
            print(f"Z: {ref_frame.z_axis}")
        
        return ref_frame
    
    def set_reference_frame_submolecule(
            self,
            submolecule: "Submolecule",
            method: str = "centroid",
            print_info: bool = False,
        ) -> ReferenceFrame:
        """   
        Computes the Reference Frame based on a Submolecule

        Further aligns all atoms in the parent Molecule accordingly
        """
        if len(submolecule.coordinates) < 3:
            raise ValueError("Submolecule must have at least 3 atoms to define a reference frame")
        # Center Submolecule
        center = self.get_center(submolecule, method=method)
        submolecule.coordinates -= center
        ref_frame = self.compute_reference_frame(
            submolecule.coordinates,
            submolecule.masses,
            center=np.zeros(3) # Already centered
        )
        # Translate Parent Molecule
        if hasattr(submolecule, 'parent') and submolecule.parent is not None: 
            parent_mol = submolecule.parent
            parent_mol.coordinates -= center
        if print_info:
            print(f"Submolecule center: {center}")
            print(f"Inertia tensor eigenvalues: {ref_frame.eigenvalues}")
            print(f"Reference frame axes:")
            print(f"X: {ref_frame.x_axis}")
            print(f"Y: {ref_frame.y_axis}")
            print(f"Z: {ref_frame.z_axis}")
        return ref_frame
    
    def set_reference_frame_from_indices(
            self,
            molecule: Union["Molecule", "Submolecule"],
            atom_indices: List[int],
            method: str = "centroid",
            print_info: bool = False,
        ) -> ReferenceFrame:
        """   
        Computes the Reference Frame based on a set of atom indices
        """
        if len(atom_indices) < 3:
            raise ValueError("At least 3 atom indices must be provided to define a reference frame")
        # Extract coordinates and masses of specified atoms
        selected_coords = molecule.coordinates[atom_indices]
        selected_masses = molecule.masses[atom_indices]
        center = self.get_center_coords(selected_coords, masses=selected_masses, method=method)
        # Center the total Molecule
        molecule.coordinates -= center
        selected_coords -= center
        ref_frame = self.compute_reference_frame(
            selected_coords,
            selected_masses,
            center=np.zeros(3) # Already centered
        )
        if print_info:
            print(f"Molecule center (based on indices {atom_indices}): {center}")
            print(f"Inertia tensor eigenvalues: {ref_frame.eigenvalues}")
            print(f"Reference frame axes:")
            print(f"X: {ref_frame.x_axis}")
            print(f"Y: {ref_frame.y_axis}")
            print(f"Z: {ref_frame.z_axis}")
        return ref_frame
    
    def set_reference_frame_from_atoms(
            self,
            molecule: Union["Molecule", "Submolecule"],
            atom_labels: List[int],
            method: str = "centroid",
            print_info: bool = False,
        ) -> ReferenceFrame:
        """   
        Computes the Reference Frame based on a set of atom labels

        Args:
            molecule: Molecule or Submolecule object
            atom_labels: List of atom labels to define the reference frame
            method: Method to calculate center point ("centroid" or "com")
            print_info: Whether to print information about the reference frame
        """
        if len(atom_labels) < 3:
            raise ValueError("At least 3 atom labels must be provided to define a reference frame")

        # Extract coordinates and masses of specified atoms
        selected_coords = []
        selected_masses = []
        for label in atom_labels:
            coord = molecule.get_coords_by_label(label)
            if coord.ndim > 1:
                raise ValueError(f"Atom label {label} corresponds to multiple atoms")
            selected_coords.append(coord)
            selected_masses.append(molecule.get_mass_by_label(label))

        selected_coords = np.array(selected_coords)
        selected_masses = np.array(selected_masses)

        center = self.get_center_coords(selected_coords, masses=selected_masses, method=method)

        # Center the total Molecule
        molecule.coordinates -= center
        selected_coords -= center

        ref_frame = self.compute_reference_frame(
            selected_coords,
            selected_masses,
            center=np.zeros(3) # Already centered
        )
        if print_info:
            print(f"Molecule center (based on labels {atom_labels}): {center}")
            print(f"Inertia tensor eigenvalues: {ref_frame.eigenvalues}")
            print(f"Reference frame axes:")
            print(f"X: {ref_frame.x_axis}")
            print(f"Y: {ref_frame.y_axis}")
            print(f"Z: {ref_frame.z_axis}")

        return ref_frame
    
    # ======================= Alignment Functions =========================

    def align_molecule_to_target(
            self,
            molecule: Union["Molecule", "Submolecule"],
            ref_frame: ReferenceFrame,
            target_point: np.ndarray,
        ) -> np.ndarray:
        """
        Rotates the molecule so that the RefFrame z-axis aligns with the vector from origin to target_point

        Args:
            molecule: Molecule or Submolecule object
            ref_frame: ReferenceFrame object of the molecule
            target_point: Target point to align the z-axis towards
        Returns:
            np.ndarray: New coordinates of the aligned molecule
        """ 
        # Directions to target 
        direction_to_target = target_point - ref_frame.origin
        target_direction = GeometryOps.normalize_vector(direction_to_target)

        # Rotation matrix: z-axis --> target_direction
        R = GeometryOps.align_vectors(ref_frame.z_axis, target_direction)

        new_coords = molecule.coordinates @ R.T
        return new_coords
    
    def align_to_vector(
            self,
            molecule: Union["Molecule", "Submolecule"],
            current_vector: np.ndarray,
            target_vector: np.ndarray,
            in_place: bool = True
        ) -> Optional[np.ndarray]:
        """  
        Rotates a molecule such that current_vector aligns with target_vector
        
        Args:
            molecule: Molecule or Submolecule object
            current_vector: Current vector in the molecule
            target_vector: Target vector to align to
            in_place: Whether to modify the molecule in place or return new coordinates
        """
        R = GeometryOps.align_vectors(current_vector, target_vector)
        new_coords = molecule.coordinates @ R.T
        if in_place:
            molecule.coordinates = new_coords
            return None
        else:
            return new_coords
        
    # ======================= Selective Transformations =======================

    def translate_submolecule(
            self, 
            submolecule: "Submolecule",
            target_point: Union[List[float], np.ndarray],
            parent_molecule: Optional["Molecule"] = None,
            method: str = "centroid",
        ) -> None:
        """
        Translates the Submolecule to a target point and updates the parent 
        
        Args:
            submolecule: Submolecule object to translate
            target_point: Target point to translate the submolecule to
            parent_molecule: Optional parent Molecule object to update coordinates
            method: Method to calculate center point ("centroid" or "com")
        """
        # Translation vector
        center = self.get_center(submolecule, method=method)
        translation = np.asarray(target_point) - center

        # Update Submolecule
        submolecule.coordinates += translation

        # Update Parent Molecule if provided
        if parent_molecule is not None and hasattr(submolecule, 'get_atom_labels'):
            submol_labels = submolecule.get_atom_labels()

            # Mask
            mask = np.isin(parent_molecule.atom_labels, list(submol_labels))
            parent_molecule.coordinates[mask] += translation

    # ===================== Utility Functions ======================

    def bond_vector(
            self,
            molecule: Union["Molecule", "Submolecule"],
            label1: str,
            label2: str,
            normalize: bool = True
        ) -> np.ndarray:
            """
            Computes the bond vector from atom with label1 to atom with label2
            """
            coord1 = molecule.get_coords_by_label(label1)
            coord2 = molecule.get_coords_by_label(label2)

            if coord1.ndim > 1 or coord2.ndim > 1:
                raise ValueError("Atom labels must correspond to single atoms")
            vector = coord2 - coord1

            if normalize:
                vector = GeometryOps.normalize_vector(vector)
            return vector
    



# ========================= Testing Utilities for Operators =========================

if __name__ == "__main__":

    from molecule_class import Molecule
    from box import SimulationBox
    from scipy.optimize import linear_sum_assignment
    import matplotlib.pyplot as plt
    import time
    



#    # ===================== Setup: H2O Dimer =====================
#    xyz_1 = """
#    6
#    Coordinates from ORCA-job h2o_2 E -152.102897751726
#    O           0.20131422818946     -0.13419863189991     -0.38118207664628
#    H           1.10826926774619      0.03054527004232     -0.16898516572928
#    H          -0.26386315635905     -0.06924940339370      0.43390806815981
#    O           3.06513444516272      0.52224610599392      0.05902763481169
#    H           3.25817767296947      1.29105665797100     -0.44996761253291
#    H           3.66908154229119     -0.14039999871364     -0.23001884806303
#    """
#    # Second molecule with one atom permutated to test the symmetry effect
#    # of the rmsd calculation and the operators
#    xyz_2 = """
#6
#    Coordinates from ORCA-job h2o_2 E -152.102897751726
#    O           0.20131422818946     -0.13419863189991     -0.38118207664628
#    H          -0.26386315635905     -0.06924940339370      0.43390806815981
#    H           1.10826926774619      0.03054527004232     -0.16898516572928
#    O           3.06513444516272      0.52224610599392      0.05902763481169
#    H           3.25817767296947      1.29105665797100     -0.44996761253291
#    H           3.66908154229119     -0.14039999871364     -0.23001884806303
#    """
#
#    # Make n test points to calculate the speed of the operators and the RMSD calculation
#
#
#
#
#    molecule_1 = Molecule.from_xyz(xyz_1)
#    molecule_2 = Molecule.from_xyz(xyz_2)
#    submol_indices = [[0, 1, 2], [3, 4, 5]]
#
#    # Compute Rotation matrix for alignment
#    R = GeometryOps.kabsch_rotation(molecule_1.coordinates, molecule_2.coordinates)
#    print("Rotation matrix to align molecule 2 to molecule 1:")
#    # Compute RMSD after alignment
#    aligned_coords_2 = molecule_2.coordinates @ R.T
#
#    def rmsd(coords1, coords2):
#        return np.sqrt(np.mean(np.sum((coords1 - coords2)**2, axis=1)))
#    print(f"RMSD after alignment: {rmsd(molecule_1.coordinates, aligned_coords_2):.4f} Å") 
#
#    def cost_matrix(coords1, coords2):
#        """
#        Computes a cost matrix where each element (i,j)
#        is the squared norm of the positional vectors a_i, b_j
#        """
#        cost_matrix = np.zeros((len(coords1), len(coords2)))
#        for i in range(len(coords1)):
#            for j in range(len(coords2)):
#                cost_matrix[i, j] = np.linalg.norm(coords1[i] - coords2[j])**2
#        return cost_matrix
#    
#    def find_optimal_correspondence(coords1,coords2):
#        """
#        Find the optimal one-to-one correspondence between atoms
#        using the Hungarian algorithm on the cost matrix
#        """
#        cost_mat= cost_matrix(coords1, coords2)
#        row_ind, col_ind = linear_sum_assignment(cost_mat)
#        # Reorder coords2 to match the optimal correspondence with coords1
#        reordered_coords2 = coords2[col_ind]
#        return reordered_coords2
#
#    # Test the optimal correspondence function
#    reordered_coords_2 = find_optimal_correspondence(molecule_1.coordinates, molecule_2.coordinates)
#    print(f"RMSD after optimal correspondence: {rmsd(molecule_1.coordinates, reordered_coords_2):.4f} Å")
#
#    times_kabsch = []
#    times_opt = [] 
#    dims = [5,10,20,50,100,200,500,1000, 2000,5000]
#
#    for dim in dims:
#        # Speed test the RMSD with Kabsch and with optimal correspondence
#        coords1 = np.random.rand(dim, 3)
#        coords2 = np.random.rand(dim, 3)
#        start = time.time()
#        R = GeometryOps.kabsch_rotation(coords1, coords2)
#        aligned_coords2 = coords2 @ R.T
#        rmsd_kabsch = rmsd(coords1, aligned_coords2)
#        end = time.time()
#        times_kabsch.append(end - start)
#        start = time.time()
#        reordered_coords2 = find_optimal_correspondence(coords1, coords2)
#        rmsd_opt = rmsd(coords1, reordered_coords2)
#        end = time.time()
#        times_opt.append(end - start)
#    
#    plt.figure(figsize=(8,5))
#    plt.plot(dims, times_kabsch, label="Kabsch RMSD", marker='o')
#    plt.plot(dims, times_opt, label="Optimal Correspondence RMSD", marker='o')
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.xlabel("Number of atoms (log scale)")
#    plt.ylabel("Time (seconds, log scale)")
#    plt.title("RMSD Calculation Time: Kabsch vs Optimal Correspondence")
#    plt.legend()
#    plt.grid(True, which="both", ls="--")
#    plt.tight_layout()
#    plt.savefig("rmsd_timing.png", dpi=150)
#    plt.close()
#                                                        
#                                                    
#






#    def animate_operator(molecule,submolecule_indices, operator_func, operator_name,
#                         n_frames=60, n_applications=3, save_path=None, adaptive=None):
#        
#        """  
#        Animate a nonlocal operator by interpolation between original and transformed coordinates
#
#        Args:
#            molecule: Input Molecule object
#            submolecule_indices: List of lists of atom indices defining submolecules
#            operator_func: Non-local operator function to apply
#            operator_name: Name of the operator (for title)
#            n_frames: Number of frames in the animation
#            n_applications: Number of times to apply the operator for animation
#            save_path: Optional path to save the animation as a GIF
#            adaptive: Whether to use adaptive scaling for the operator
#        """
#        import matplotlib.animation as animation
#        from pathlib import Path
#
#        # Collect keyframes: original + each operator application
#        keyframes = [molecule.copy()]
#        current = molecule.copy()
#        for _ in range(n_applications):
#            if adaptive is not None:
#                new = operator_func(current, submolecule_indices, adaptive=adaptive)
#            else:
#                new = operator_func(current, submolecule_indices)
#            keyframes.append(new)
#            current = new.copy()
#
#        # Build interpolated frames between each pair of keyframes
#        all_frames = []
#        for k in range(len(keyframes) - 1):
#            start_coords = keyframes[k].coordinates
#            end_coords = keyframes[k + 1].coordinates
#            for t in range(n_frames):
#                frac = t / n_frames
#                interp = start_coords + frac * (end_coords - start_coords)
#                all_frames.append((interp.copy(), k+1))
#            
#            # Hold on final position
#            for _ in range(n_frames // 3):
#                all_frames.append((end_coords.copy(), k+1))
#        
#        # Precompute colors and sizes per submolecule
#        sub_edge_colors = plt.cm.Set2(np.linspace(0, 1, max(len(submolecule_indices), 1)))
#        atom_colors = molecule.get_atom_colors()
#        atom_sizes = molecule.get_atom_sizes()
#
#        # Compute fixed axis limits for all keyframes
#        all_kf_coords = np.vstack([kf.coordinates for kf in keyframes])
#        center = all_kf_coords.mean(axis=0)
#        max_range = np.ptp(all_kf_coords, axis=0).max() / 2 + 1.5
#
#        fig = plt.figure(figsize=(9, 7))
#        ax = fig.add_subplot(111, projection='3d')
#
#        def update(frame_idx):
#            coords, app_num = all_frames[frame_idx]
#            ax.cla()
#
#            for sub_i, indices in enumerate(submolecule_indices):
#                sub_coords = coords[indices]
#                colors = [atom_colors[i] for i in indices]
#                sizes = [atom_sizes[i] for i in indices]
#
#                ax.scatter(
#                    sub_coords[:, 0], sub_coords[:, 1], sub_coords[:, 2],
#                    c=colors, s=sizes, edgecolors=sub_edge_colors[sub_i],
#                    linewidths=2, alpha=0.9, depthshade=True
#                )
#                # Draw covalent bonds within submolecule
#                for i in range(len(indices)):
#                    for j in range(i + 1, len(indices)):
#                        dist = np.linalg.norm(sub_coords[i] - sub_coords[j])
#                        if dist < 1.3:
#                            ax.plot(
#                                [sub_coords[i, 0], sub_coords[j, 0]],
#                                [sub_coords[i, 1], sub_coords[j, 1]],
#                                [sub_coords[i, 2], sub_coords[j, 2]],
#                                color='dimgray', linewidth=2, alpha=0.7
#                            )
#
#            ax.set_title(f"{operator_name}  (step {app_num}/{n_applications})",
#                     fontsize=13, fontweight='bold')
#            ax.set_xlabel('X (Å)')
#            ax.set_ylabel('Y (Å)')
#            ax.set_zlabel('Z (Å)')
#            ax.set_xlim(center[0] - max_range, center[0] + max_range)
#            ax.set_ylim(center[1] - max_range, center[1] + max_range)
#            ax.set_zlim(center[2] - max_range, center[2] + max_range)
#            ax.view_init(elev=20, azim=30 + frame_idx * 0.3)
#            return []
#
#        ani = animation.FuncAnimation(
#            fig, update, frames=len(all_frames), interval=33, blit=False
#        )
#
#        if save_path:
#            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
#            if save_path.endswith('.gif'):
#                ani.save(save_path, writer='pillow', fps=30, dpi=150)
#            else:
#                ani.save(save_path, writer='ffmpeg', fps=30, dpi=150)
#            print(f"Saved: {save_path}")
#            plt.close()
#        else:
#            plt.show()
#
#
#
#
#
#    covalent_radii = [0.66, 0.31, 0.31, 0.66, 0.31, 0.31]
#    box = SimulationBox.from_covalent_radii(covalent_radii, 6, scale_factor=5.0)
#    nonlocal_ops = NonLocalOperators(simulation_box=box)
#
#    # ===================== Generate Animations =====================
#    # (operator_name, function, adaptive_kwarg)
#    #   adaptive=None means the operator doesn't take an adaptive parameter
#    operators_to_animate = [
#        ("Twist",                    nonlocal_ops.twist_operator,                      True),
#        ("Large Displacement",       nonlocal_ops.large_displacement,                  True),
#        ("Mirror",                   nonlocal_ops.mirror_operator,                     None),
#        ("Random SO(3)",             nonlocal_ops.random_so3_operator,                 None),
#        ("Principal Axis Rotation",  nonlocal_ops.principal_axis_rotation_operator,    True),
#        ("Roto-Reflection",          nonlocal_ops.roto_reflection_operator,            True),
#        ("Exchange",                 nonlocal_ops.exchange_operator,                   None),
#    ]
#
#    for op_name, op_func, adaptive in operators_to_animate:
#        safe_name = op_name.lower().replace(" ", "_").replace("(", "").replace(")", "")
#        save_path = f"figures/operator_{safe_name}.gif"
#        print(f"Animating: {op_name} -> {save_path}")
#        animate_operator(
#            molecule, submol_indices, op_func, op_name,
#            n_frames=100, n_applications=5,
#            save_path=save_path,
#            adaptive=adaptive,
#        )
#
#    print("All animations done.")

