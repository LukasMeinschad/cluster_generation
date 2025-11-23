import numpy as np
import matplotlib.pyplot as plt


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
        plt.show()

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
        plt.savefig("molecule_reference_frame_plot.png")
        plt.show()

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
        plt.show()

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
        plt.show()

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
        plt.show()

