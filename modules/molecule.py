import numpy as np


class Molecule:
    def __init__(self, name: str = "Unnamed Molecule"):
        """ 
        Initializes a Molecule object
        """
        self.name = name
        self.atom_labels = np.array([]) # empty array for element labels
        self.coordinates = np.empty((0,3), dtype=np.float64) # Empty 2D array
