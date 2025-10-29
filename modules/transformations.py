import numpy as np
from mendeleev.fetch import fetch_table


class Transformation:
    """
    Class that performs transformation of molecules in R^3
    """
    def __init__(self):
        self.name= None
    def centre(self,molecule,method="centroid"):
        """ 
        Function that determines a center in space 
        """
        coords = molecule.coordinates

        if method == "centroid":
            # Calculates the mean point in space
            center_point = np.mean(coords,axis=0)

        if method == "com":
            # Calculates the center of mass
            center_point = self.com(molecule)

        else:
            raise ValueError("Method must be 'centroid' or 'com'")

        return centre_point
            

    def com(self,molecule):
        """ 
        Calculates the center of mass for a given molecule
        """
        coords = molecule.coordinates
        masses = molecule.masses

        com = np.sum(masses[:, np.newaxis] * coords, axis=0) / np.sum(masses)
        return com
        
