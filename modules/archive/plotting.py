"""
Plotting module for visualization
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Import for 3D plotting  
from typing import List, Optional
from pathlib import Path 




    


class Plotting:
    """
    Different Functionalities for Plotting of a Molecule
    """
    def __init__(self):
        pass

    color_dict = {
            "H": "grey",
            "O": "red",
            "C": "black",
            "N": "blue",
            "S": "yellow"
    }
    size_dict = {
            "H": 50,
            "O": 100,
            "C": 80,
            "N": 90,
            "S": 120
    }

    def plot_molecule_3d(self,molecule,cov_bonds=None,hydrogen_bonds=None):
        """ 
        Plots a 3D representation of the molecule

        Args:
            molecule: Molecule object to be plotted
            cov_bonds: List of covalent bonds as tuples (atom_label1, atom_label2)
            hydrogen_bonds: List of hydrogen bonds as tuples (atom_label1, atom_label2)
        """
        

        coords = molecule.coordinates
        elements = molecule.atom_labels
        # Remove digits from element symbols
        elements = [''.join(filter(str.isalpha, el)) for el in elements]
        colors = [self.color_dict.get(el, "green") for el in elements]
        sizes = [self.size_dict.get(el, 70) for el in elements]

        # Add bonds if provided

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for i, (x, y, z) in enumerate(coords):
            ax.scatter(x, y, z, color=colors[i], s=sizes[i], label=elements[i] if i == 0 or elements[i] != elements[i-1] else "") 
        
        if cov_bonds:
            for bond in cov_bonds:
                idx1 = np.where(molecule.atom_labels == bond[0])[0][0]
                idx2 = np.where(molecule.atom_labels == bond[1])[0][0]
                x_vals = [coords[idx1][0], coords[idx2][0]]
                y_vals = [coords[idx1][1], coords[idx2][1]]
                z_vals = [coords[idx1][2], coords[idx2][2]]
                ax.plot(x_vals, y_vals, z_vals, color='black', linewidth=2)
        if hydrogen_bonds:
            for bond in hydrogen_bonds:
                idx1 = np.where(molecule.atom_labels == bond[0])[0][0]
                idx2 = np.where(molecule.atom_labels == bond[1])[0][0]
                x_vals = [coords[idx1][0], coords[idx2][0]]
                y_vals = [coords[idx1][1], coords[idx2][1]]
                z_vals = [coords[idx1][2], coords[idx2][2]]
                ax.plot(x_vals, y_vals, z_vals, color='lightgrey', linestyle='dashed', linewidth=1.5)

        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        # Make Legend in outside box
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
        plt.title("3D Molecule Plot")
        plt.savefig("molecule_3d_plot.png")
        #plt.show()
        plt.close()


    def plot_reference_frame(self, molecule,reference_frame, reference_frame_origin):
        """ 
        Plots the reference frame of a molecule in 3D
        
        Args:
            molecule: Molecule object to be plotted
            reference_frame: 3x3 numpy array where columns are the x, y, z axes
        """
        
        coords = molecule.coordinates
        elements = molecule.atom_labels
        # Remove digits from element symbols
        elements = [''.join(filter(str.isalpha, el)) for el in elements]
        colors = [self.color_dict.get(el, "green") for el in elements]
        sizes = [self.size_dict.get(el, 70) for el in elements]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for i, (x, y, z) in enumerate(coords):
            ax.scatter(x, y, z, color=colors[i], s=sizes[i], label=elements[i] if i == 0 or elements[i] != elements[i-1] else "")
        
        # Plot reference frame as arrows
        origin = reference_frame_origin
        x_axis = reference_frame[:,0]
        y_axis = reference_frame[:,1]
        z_axis = reference_frame[:,2]
        ax.quiver(origin[0], origin[1], origin[2], x_axis[0], x_axis[1], x_axis[2], color='r', length=1.0, normalize=True, label='X axis')
        ax.quiver(origin[0], origin[1], origin[2], y_axis[0], y_axis[1], y_axis[2], color='g', length=1.0, normalize=True, label='Y axis')
        ax.quiver(origin[0], origin[1], origin[2], z_axis[0], z_axis[1], z_axis[2], color='b', length=1.0, normalize=True, label='Z axis')
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        # Make Legend in outside box
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
        plt.title("3D Molecule with Reference Frame")
        #plt.savefig("molecule_reference_frame_plot.png")

        #plt.show()
        plt.close()


    def plot_cube(self, cube_corners):
        """ 
        Plots a cube define by its corner points in 3D
        
        Args:
            cube_corners: 8x3 numpy array with the coordinates of the cube corners
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Define the 12 edges of the cube
        edges = [
            (0, 1), (1, 3), (3, 2), (2, 0),  # Bottom face
            (4, 5), (5, 7), (7, 6), (6, 4),  # Top face
            (0, 4), (1, 5), (2, 6), (3, 7)   # Vertical edges
        ]

        # Plot each edge
        for edge in edges:
            p1 = cube_corners[edge[0]]
            p2 = cube_corners[edge[1]]
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color='orange')

        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        plt.title("3D Cube Plot")
        plt.savefig("cube_3d_plot.png")
        plt.close()

        
        #plt.show()

    def plot_cube_molecule(self,molecule,cube_corners, ref_frame=None):
        """ 
        Plots a cube and a molecule in 3D
        
        Args:
            molecule: Molecule object to be plotted
            cube_corners: 8x3 numpy array with the coordinates of the cube corners
        """
        coords = molecule.coordinates
        elements = molecule.atom_labels
        # Remove digits from element symbols
        elements = [''.join(filter(str.isalpha, el)) for el in elements]
        colors = [self.color_dict.get(el, "green") for el in elements]
        sizes = [self.size_dict.get(el, 70) for el in elements]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for i, (x, y, z) in enumerate(coords):
            ax.scatter(x, y, z, color=colors[i], s=sizes[i], label=elements[i] if i == 0 or elements[i] != elements[i-1] else "") 
        
        # Define the 12 edges of the cube
        edges = [
            (0, 1), (1, 3), (3, 2), (2, 0),  # Bottom face
            (4, 5), (5, 7), (7, 6), (6, 4),  # Top face
            (0, 4), (1, 5), (2, 6), (3, 7)   # Vertical edges
        ]

        # Plot each edge
        for edge in edges:
            p1 = cube_corners[edge[0]]
            p2 = cube_corners[edge[1]]
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color='orange')

        if ref_frame is not None:
            # Plot reference frame as arrows
            origin = np.mean(cube_corners, axis=0)
            x_axis = ref_frame[:,0]
            y_axis = ref_frame[:,1]
            z_axis = ref_frame[:,2]
            ax.quiver(origin[0], origin[1], origin[2], x_axis[0], x_axis[1], x_axis[2], color='r', length=1.0, normalize=True, label='X axis')
            ax.quiver(origin[0], origin[1], origin[2], y_axis[0], y_axis[1], y_axis[2], color='g', length=1.0, normalize=True, label='Y axis')
            ax.quiver(origin[0], origin[1], origin[2], z_axis[0], z_axis[1], z_axis[2], color='b', length=1.0, normalize=True, label='Z axis')

        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        # Make Legend in outside box
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
        plt.title("3D Molecule and Cube Plot")
        plt.savefig("molecule_cube_3d_plot.png")
        plt.close()
        #plt.show()

    def plot_uniform_sampling(self,sampled_molecules, cube_corners):
        """
        Plots uniformly samples molecules in 3D
        """ 
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plot cube edges
        edges = [
            (0, 1), (1, 3), (3, 2), (2, 0),  # Bottom face
            (4, 5), (5, 7), (7, 6), (6, 4),  # Top face
            (0, 4), (1, 5), (2, 6), (3, 7)   # Vertical edges
        ]

        for edge in edges:
            p1 = cube_corners[edge[0]]
            p2 = cube_corners[edge[1]]
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color='orange')
        


        for mol in sampled_molecules:
            coords = mol.coordinates
            elements = mol.atom_labels
            # Remove digits from element symbols
            elements = [''.join(filter(str.isalpha, el)) for el in elements]
            colors = [self.color_dict.get(el, "green") for el in elements]
            sizes = [self.size_dict.get(el, 70) for el in elements]

            for i, (x, y, z) in enumerate(coords):
                ax.scatter(x, y, z, color=colors[i], s=sizes[i], alpha=0.5)

        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        plt.title("3D Uniformly Sampled Molecules")
        plt.savefig("uniformly_sampled_molecules_3d_plot.png")

        #plt.show()
        plt.close()

    def plot_sampling_sphere(self,phi,theta):
        """ 
        Plots the sampling that was carried out on the sphere.

        For this we use spherical coordinates, which means we sample the phi [0, pi] and theta [0,2pi] space
        """
        phi_deg = np.rad2deg(phi)
        theta_deg = np.rad2deg(theta)

        plt.figure(figsize=(12,6))
        plt.scatter(phi_deg, theta_deg, label="Sampled Points")
        plt.xlabel("$\phi$ (Deg)")
        plt.ylabel("θ (Deg)")
        plt.title("Sampled Space on the Sphere Surface")
        plt.savefig("sampling_points_theta_phi_plot.png")
        plt.close()

    def plot_sampling_sphere_3d(self,global_points, ref_frame=None):
        """ 
        Function to plot the sampling points on a sphere in 3D
        
        args:
            global_points: Nx3 numpy array with the sampled points in Cartesian coordinates
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(global_points[:,0], global_points[:,1], global_points[:,2], color='purple', s=20)

        if ref_frame is not None:
            # Plot reference frame as arrows
            origin = np.array([0.0, 0.0, 0.0])
            x_axis = ref_frame[:,0]
            y_axis = ref_frame[:,1]
            z_axis = ref_frame[:,2]
            ax.quiver(origin[0], origin[1], origin[2], x_axis[0], x_axis[1], x_axis[2], color='r', length=1.0, normalize=True, label='X axis')
            ax.quiver(origin[0], origin[1], origin[2], y_axis[0], y_axis[1], y_axis[2], color='g', length=1.0, normalize=True, label='Y axis')
            ax.quiver(origin[0], origin[1], origin[2], z_axis[0], z_axis[1], z_axis[2], color='b', length=1.0, normalize=True, label='Z axis')

        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        plt.title("3D Sampling Points on Sphere")
        plt.savefig("sampling_points_sphere_3d_plot.png")
        #plt.show()
        plt.close()

    def plot_spherical_mol_sampling(self,sampled_molecules,ref_frame=None, target_points=None, center_point=None):
        """
        Plots the sampled molecules on the sphere in 3D
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plot molecule with bonds 
        for mol in sampled_molecules:
            coords = mol.coordinates
            elements = mol.atom_labels
            # Remove digits from element symbols
            elements = [''.join(filter(str.isalpha, el)) for el in elements]
            colors = [self.color_dict.get(el, "green") for el in elements]
            sizes = [self.size_dict.get(el, 70) for el in elements]

            for i, (x, y, z) in enumerate(coords):
                ax.scatter(x, y, z, color=colors[i], s=sizes[i], alpha=0.5)
        # Plot target points if provided
        if target_points is not None:
            ax.scatter(target_points[:,0], target_points[:,1], target_points[:,2], color='purple', s=50, label="Target Points")
        if ref_frame is not None:
            # Plot reference frame as arrows
            origin = np.array([0.0, 0.0, 0.0])
            x_axis = ref_frame[:,0]
            y_axis = ref_frame[:,1]
            z_axis = ref_frame[:,2]
            ax.quiver(origin[0], origin[1], origin[2], x_axis[0], x_axis[1], x_axis[2], color='r', length=1.0, normalize=True, label='X axis')
            ax.quiver(origin[0], origin[1], origin[2], y_axis[0], y_axis[1], y_axis[2], color='g', length=1.0, normalize=True, label='Y axis')
            ax.quiver(origin[0], origin[1], origin[2], z_axis[0], z_axis[1], z_axis[2], color='b', length=1.0, normalize=True, label='Z axis')

        if all([v is not None for v in [target_points, center_point, ref_frame]]):
            # Plot vector from origin to target points
            for point in target_points:
                ax.plot([center_point[0], point[0]], [center_point[1], point[1]], [center_point[2], point[2]], color='grey', linestyle='dashed', linewidth=1)
            

            
        
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        plt.title("3D Sampled Molecules on Sphere")
        plt.legend()
        plt.savefig("sampled_molecules_sphere_3d_plot.png")
        #plt.show()

    def plot_sampling_cone(self,pts_local, apex,axis,height,base_radius, u ,v):
        """ 
        Plots the sampling inside a cone in 3D
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(pts_local[:,0], pts_local[:,1], pts_local[:,2], color='purple', s=20)

        # Plot cone surface
        num_circle_points = 50
        theta = np.linspace(0, 2 * np.pi, num_circle_points)
        circle_x = base_radius * np.cos(theta)
        circle_y = base_radius * np.sin(theta)
        circle_z = np.zeros(num_circle_points)

        # Create the circle at the base
        circle_points = np.vstack((circle_x, circle_y, circle_z)).T

        # Rotate and translate the circle to the cone base
        cone_base_center = apex + axis * height
        rotation_matrix = np.column_stack((u, v, axis))
        rotated_circle = circle_points @ rotation_matrix.T + cone_base_center

        # Plot lines from apex to base circle
        for i in range(num_circle_points):
            ax.plot([apex[0], rotated_circle[i,0]], [apex[1], rotated_circle[i,1]], [apex[2], rotated_circle[i,2]], color='orange', alpha=0.3)

        # Plot base circle
        ax.plot(rotated_circle[:,0], rotated_circle[:,1], rotated_circle[:,2], color='orange')

        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        plt.title("3D Cone Sampling Plot")
        plt.savefig("cone_sampling_3d_plot.png")
        #plt.show()
        plt.close()
            
    def plot_sampling_cylinder(self,pts_local, base_point, axis, height, radius):
        """ 
        Plots the sampling inside a cylinder in 3D
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(pts_local[:,0], pts_local[:,1], pts_local[:,2], color='purple', s=20)

        # Plot cylinder surface
        num_circle_points = 50
        theta = np.linspace(0, 2 * np.pi, num_circle_points)
        circle_x = radius * np.cos(theta)
        circle_y = radius * np.sin(theta)
        circle_z = np.zeros(num_circle_points)

        # Create the circle at the base
        circle_points = np.vstack((circle_x, circle_y, circle_z)).T
        # Rotate and translate the circle to the cylinder base
        cylinder_top_center = base_point + axis * height
        # Find two orthogonal vectors to the axis
        if np.allclose(axis, [0, 0, 1]):
            u = np.array([1, 0, 0])
        else:
            u = np.cross(axis, [0, 0, 1])
            u /= np.linalg.norm(u)
        v = np.cross(axis, u)
        rotation_matrix = np.column_stack((u, v, axis))
        rotated_circle_base = circle_points @ rotation_matrix.T + base_point
        rotated_circle_top = circle_points @ rotation_matrix.T + cylinder_top_center
        # Plot lines from base to top circle
        for i in range(num_circle_points):
            ax.plot([rotated_circle_base[i,0], rotated_circle_top[i,0]], 
                    [rotated_circle_base[i,1], rotated_circle_top[i,1]], 
                    [rotated_circle_base[i,2], rotated_circle_top[i,2]], color='orange', alpha=0.3)
            
        # Plot base and top circles
        ax.plot(rotated_circle_base[:,0], rotated_circle_base[:,1], rotated_circle_base[:,2], color='orange')
        ax.plot(rotated_circle_top[:,0], rotated_circle_top[:,1], rotated_circle_top[:,2], color='orange')
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        plt.title("3D Cylinder Sampling Plot")
        plt.savefig("cylinder_sampling_3d_plot.png")
        #plt.show()
        plt.close()

    # Write a main plott sampled molecules function which checks the given sampling volume
    def plot_sampled_molecules(self,sampled_molecules,sampling_type="Sphere", radius=None, center_point=None):
        """  
        Plots the sampled molecules and the sampling volume in 3D

        Distinguishes between different sampling types
        """

        if sampling_type == "Sphere":
            # Plot sphere with radius around center point
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            # Create sphere surface
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            x = radius * np.outer(np.cos(u), np.sin(v)) + center_point[0]
            y = radius * np.outer(np.sin(u), np.sin(v)) + center_point[1]
            z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center_point[2]
            ax.plot_surface(x, y, z, color='cyan', alpha=0.3)
            # Plot sampled molecules
            for mol in sampled_molecules:
                coords = mol.coordinates
                elements = mol.atom_labels
                # Remove digits from element symbols
                elements = [''.join(filter(str.isalpha, el)) for el in elements]
                colors = [self.color_dict.get(el, "green") for el in elements]
                sizes = [self.size_dict.get(el, 70) for el in elements]

                for i, (x, y, z) in enumerate(coords):
                    ax.scatter(x, y, z, color=colors[i], s=sizes[i], alpha=0.5)
            ax.set_xlabel('X (Å)')
            ax.set_ylabel('Y (Å)')
            ax.set_zlabel('Z (Å)')
            plt.title("3D Sampled Molecules in Sphere")
            plt.savefig("sampled_molecules_sphere_3d_plot.png")
            #plt.show()
            plt.close()

        if sampling_type == "Cone":
            # Plot cone with given parameters
            # TODO implement this 
            pass


    def plot_ref_frame_hbond_config(self, molecule, ref_frame, center_point):
        """
        Helper function that plots the reference frame of a local hbond donor configuration

        For this we use the center point of this configuration and the corresponding reference frame obtained 
        via the diagonalization of the inertia tensor
        """ 
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        coords = molecule.coordinates
        elements = molecule.atom_labels
        # Remove digits from element symbols
        elements = [''.join(filter(str.isalpha, el)) for el in elements]
        colors = [self.color_dict.get(el, "green") for el in elements]
        sizes = [self.size_dict.get(el, 70) for el in elements]

        for i, (x, y, z) in enumerate(coords):
            ax.scatter(x, y, z, color=colors[i], s=sizes[i], label=elements[i] if i == 0 or elements[i] != elements[i-1] else "") 
        
        # Plot reference frame as arrows
        origin = center_point
        x_axis = ref_frame[:,0]
        y_axis = ref_frame[:,1]
        z_axis = ref_frame[:,2]
        ax.quiver(origin[0], origin[1], origin[2], x_axis[0], x_axis[1], x_axis[2], color='r', length=1.0, normalize=True, label='X axis')
        ax.quiver(origin[0], origin[1], origin[2], y_axis[0], y_axis[1], y_axis[2], color='g', length=1.0, normalize=True, label='Y axis')
        ax.quiver(origin[0], origin[1], origin[2], z_axis[0], z_axis[1], z_axis[2], color='b', length=1.0, normalize=True, label='Z axis')
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        # Make Legend in outside box
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
        plt.title("3D Molecule with Hbond Config Reference Frame")

        plt.show()
        plt.close()