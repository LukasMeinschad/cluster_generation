import matplotlib.pyplot as plt
import numpy as np
from mendeleev.fetch import fetch_table
import itertools
import networkx as nx


def parse_molecule(xyz_file):
    """ 
    Parses in a molecue from a xyz file
    """
    molecule = []
    lines = xyz_file.strip().split("\n")

    for line in lines [2:]: # skip first line
        parts = line.split()
        if len(parts) >= 4:
            atom = {
                'element': parts[0],
                'coords': np.array([float(parts[1]), float(parts[2]), float(parts[3])])
            }
            molecule.append(atom)

    # Count elements then enumerate
    element_count = {}
    for atom in molecule:
        element = atom['element']
        if element not in element_count:
            element_count[element] = 0
        element_count[element] += 1
    
    # now enumerate starting from O1 to O2, H1 to H4 and so on
    for element in element_count:
        count = 1
        for atom in molecule:
            if atom['element'] == element:
                atom['element'] = f"{element}{count}"
                count += 1

        
    return molecule


def get_covalent_radius(element):
    """ 
    Helper function to obtain the covalent radius of an element
    """
    # remove the string number from the element
    ptable = fetch_table('elements')
    element = ''.join(filter(str.isalpha, element))
    row = ptable[ptable['symbol'] == element]
    if not row.empty:
        return row.iloc[0]['covalent_radius_pyykko'] / 100  # convert pm to angstrom
    else:
        raise ValueError(f"Element {element} not found in periodic table.")


def degree_of_covalence(molecule):
    """ 
    Determines the degree of covalence of a bond and groups into covalent and non-covalent bonds
    
    This is given by exp(-[(r_ij)/(C_ij) - 1])
    with r_ij the distance between atoms i and j and C_ij the sum of covalent radii of atoms i and j
    """
    cov_bonds = []
    hydrogen_bonds = []
    for atom1, atom2 in itertools.combinations(molecule, 2):
        coords1 = atom1['coords']
        coords2 = atom2['coords']
        euclidean_distance = np.linalg.norm(coords1 - coords2)
        cov_radius1 = get_covalent_radius(atom1['element'])
        cov_radius2 = get_covalent_radius(atom2['element'])
        C_ij = cov_radius1 + cov_radius2
        degree = np.exp(-((euclidean_distance / C_ij) - 1))
        if degree >= 0.2 and degree < 0.7:
            hydrogen_bonds.append((atom1['element'], atom2['element'], degree))
        elif degree >= 0.7:
            cov_bonds.append((atom1['element'], atom2['element'], degree))

    # Loop through hydrogen bonds to only consider O-H pairs:
    hydrogen_bonds = [bond for bond in hydrogen_bonds if ('H' in bond[0] and 'O' in bond[1]) or ('H' in bond[1] and 'O' in bond[0])]
    return cov_bonds, hydrogen_bonds


def find_submolecules(molecule,bonds):
    """ 
    Finds connected structures using networkx
    """
    G = nx.Graph()

    # Atoms = nodes
    for atom in molecule:
        G.add_node(atom['element'])

    for bond in bonds:
        G.add_edge(bond[0],bond[1])

    submolecules = list(nx.connected_components(G))

    molecules_list = []

    for submol in submolecules:
        submol_atoms = [atom for atom in molecule if atom['element'] in submol]
        molecules_list.append(submol_atoms)

    print(molecules_list)
    return molecules_list
        


if __name__ == "__main__":
    with open('/media/storage_6/lme/master_thesis/initial_tests/h2o_2/h2o_2.xyz', 'r') as file:
        xyz_content = file.read()
   

    molecule = parse_molecule(xyz_content)
    cov_bonds, hydrogen_bonds = degree_of_covalence(molecule)

    find_submolecules(molecule,cov_bonds)


