""" 
Helper Module to construct redundant internal coordinates for a given molecule.

This follows closly the algorithm proposed by 
The efficient optimization of molecular geometries usingredundant internal coordinates Vebjørn Bakken; Trygve Helgaker
https://pubs.aip.org/aip/jcp/article/117/20/9160/464558/The-efficient-optimization-of-molecular-geometries?getftrtoken=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAALIwga8GCSqGSIb3DQEHBqCBoTCBngIBADCBmAYJKoZIhvcNAQcBMB4GCWCGSAFlAwQBLjARBAzfLOP2QTYjtpcdiQMCARCAa1j66AHFHg_edQns_SEkjvcWosKtK_C6eCMxIbCXAgIfLyutXgoQlViHuFdfZKEYp8UPmrbzlyx4L5jxVuGDX3lloAcTAfsqqKm65tge8e-BGFxV-WF4n_JAI7x88IhHlvEj-0rR7i8OH9bf
"""
import numpy as np 
import copy
from typing import List, Tuple, Optional
from itertools import combinations

# Module imports
from geometry import GeometryOps
from molecule_class import Molecule

class RedundantInternalCoordinates:
    """ 
    Class that constructs and stores the redundant internal coordinates for a given molecules. The internal coordinates are stored
    as a list of tuples where each tuple contains the indices of the atoms involved in the internal coordinate
    - regular bonds: (i, j)
    - angles: (i, j, k)
    - dihedrals: (i, j, k ,l)
    
    Further these are grouped either in intramolecular or intermolecular coordinates.
    Intermolecular coordinates are defined between atoms of submolecules within the cluster structure.
    As an example one can consider the water dimer with each water molecule being a submolecule. Then the intramolecular coordinates are defined between atoms of the same water molecule and the intermolecular coordinates are defined between atoms of different water molecules.
    """

    def __init__(self, mol: Molecule, submolecule_indices: Optional[List[List[int]]] = None):
        self.mol = mol
        self.submolecule_indices = submolecule_indices
        
        # Intramolecular coordinates
        self.bonds = []
        self.angles = []
        self.dihedrals = []

        # Intermolecular coordinates
        self.inter_bonds = []
        self.inter_angles = []
        self.inter_dihedrals = []


    def covalent_bonds(self):
        """ 
        Uses the respective method of the molecule class to find covalent bonds
        """
        # Compute bonds 
        self.mol.compute_bonds()
        bond_ls = [(b.atom1, b.atom2) for b in self.mol._covalent_bonds]
        # Convert into list of indices
        bond_idx = [(self.mol.get_index_by_label(b[0]), self.mol.get_index_by_label(b[1])) for b in bond_ls]
        # Append to bonds list
        self.bonds.extend(bond_idx)
    
    def interfragment_distance(self):
        """  
        Constructs bond coordinates along the shortest interfragmental distance between two submolecules
        """
        if self.submolecule_indices is None:
            raise ValueError("Submolecule indices must be provided for interfragment distance calculation")
        
        coords = self.mol.coordinates
        # Loop over all pairs of submolecules
        for submol1, submol2 in combinations(self.submolecule_indices, 2):
            min_distance = np.inf
            closest_pair = None
            # Loop over all pairs of atoms between the two submolecules
            for idx1 in submol1:
                for idx2 in submol2:
                    dist = np.linalg.norm(coords[idx1] - coords[idx2])
                    if dist < min_distance:
                        min_distance = dist
                        closest_pair = (idx1, idx2)
            self.inter_bonds.append(closest_pair)

    def is_inter(self, atoms):
        """ 
        Helper function that checks if atoms span over multiple fragments
        """
        for frag in self.submolecule_indices:
            if all(a in frag for a in atoms):
                return False
        return True
    

    
    def build_angles(self):
        for (i,j), (k,l) in combinations(self.bonds + self.inter_bonds, 2):
            shared = set([i,j]) & set([k,l])
            if len(shared) == 1:
                center = shared.pop()

                a = (set([i,j]) - {center}).pop()
                b = center
                c = (set([k,l]) - {center}).pop()

                angle = (a, b, c)
                if self.is_inter(angle):
                    self.inter_angles.append(angle)
                else:
                    self.angles.append(angle)
                


    def build_dihedrals(self):
        all_angles = self.angles + self.inter_angles
        for a1, a2 in combinations(all_angles, 2):
            shared = set(a1) & set(a2)
            if len(shared) == 2:
                b,c = shared

                a = (set(a1) - {b,c}).pop()
                d = (set(a2) - {b,c}).pop()

                dih = (a, b, c, d)
                if self.is_inter(dih):
                    self.inter_dihedrals.append(dih)
                else:
                    self.dihedrals.append(dih)


    def get_intermolecular_ics(self):
        """ 
        Returns a list of all intermolecular internal coordinates
        """
        q_inter = (
            self.inter_bonds +
            self.inter_angles +
            self.inter_dihedrals
        )
        return q_inter
    
    def build_q_inter(self, q_inter):
        """  
        Builds the ntermolecular coordinate vector q(x)
        """
        coords = self.mol.coordinates
        q = []
        
        # Distances
        for i,j in self.inter_bonds:
            dist = GeometryOps.distance(coords[i], coords[j])
            q.append(dist)
        print(q)

        

    
             


if __name__ == "__main__":
    """ 
    The following code is used for testing of this algorithm
    """

    h2o_dimer = """
    6
    Coordinates from ORCA-job h2o_2 E -152.102897751726
    O           0.20131422818946     -0.13419863189991     -0.38118207664628
    H           1.10826926774619      0.03054527004232     -0.16898516572928
    H          -0.26386315635905     -0.06924940339370      0.43390806815981
    O           3.06513444516272      0.52224610599392      0.05902763481169
    H           3.25817767296947      1.29105665797100     -0.44996761253291
    H           3.66908154229119     -0.14039999871364     -0.23001884806303  
    """
    mol = Molecule.from_xyz(h2o_dimer)
    print("Molecule initialized", mol)
    submolecules = [[0, 1, 2], [3, 4, 5]]

    RIC = RedundantInternalCoordinates(mol, submolecule_indices=submolecules)
    RIC.covalent_bonds()
    RIC.interfragment_distance()
    RIC.build_angles()
    RIC.build_dihedrals()
    q_inter = RIC.get_intermolecular_ics()
    RIC.build_q_inter(q_inter)
