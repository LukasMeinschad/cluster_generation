import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import time 
from itertools import permutations

class MolecularGraph():
    """    
    Class that transforms molecular structures into graph representations to compute descriptors invariants
    and infer certain properties of the molecules.
    """

    def __init__(self, molecule, bonds, logger):
        """ 
        Initializes the Molecular Graph Class
        """
        
        self.molecule = molecule
        if self.molecule.covalent_bonds is None:
            self.bonds = bonds
        else:
            self.bonds = self.molecule.covalent_bonds
        self.logger = logger

        # Initialize Graph 
        self.graph = self.build_graph()
        self.euclidean_graph = self.build_euclidean_graph()

    def build_graph(self):
        """
        Builds a molecular graph from the molecule's atoms and bonds
        """
        G = nx.Graph()
        # Add atoms as nodes
        for atom in  self.molecule.atom_labels:
            G.add_node(atom)

        # Add bonds and edges
        for bond in self.bonds:
            G.add_edge(bond[0], bond[1])
            
        
        self.graph = G

        return G

    def get_vertices(self):
        """
        Returns the vertices of the molecular graph
        """
        return self.graph.nodes()

    def get_edges(self):
        """
        Returns the edges of the molecular graph
        """
        return self.graph.edges()

    def draw_graph(self):
        """
        Draws the molecular graph using matplotlib
        """
        plt.figure(figsize=(8, 6))
        pos = nx.spring_layout(self.graph)
        nx.draw(self.graph, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=500, font_size=10)
        plt.title("Molecular Graph Representation")
        plt.savefig("molecular_graph.png") 

    def get_adjacency_matrix(self):
        """
        Returns the adjacency matrix of the molecular graph
        """
        return nx.adjacency_matrix(self.graph).todense() 
    
    def find_automorphism_group(self):
        """ 
        Finds the automorphism group of the molecular graph
        """
        start_time  = time.time()
        automorphisms = []
        n = len(self.graph.nodes())
        nodes = list(self.graph.nodes())

        for perm in permutations(range(n)):
            mapping = {nodes[i]: nodes[perm[i]] for i in range(n)}

            # Check if its and automorphism
            is_isomprhic = nx.is_isomorphic(self.graph, nx.relabel_nodes(self.graph, mapping))

            if is_isomprhic:
                automorphisms.append(mapping)
        end_time = time.time()
        print(f"Automorphism group found in {end_time - start_time} seconds")
        return automorphisms

    def build_euclidean_graph(self):
        """ 
        Builds a euclidean weighted graph based on euclidean distances between atoms
        
        # This graph has the following adjacency matrix:

        D_ij = d_ij if i neq j and d_ij is the euclidean distance between atoms i and j
                0   if i = j
        """
        G = nx.Graph()
        coords = self.molecule.coordinates
        atom_labels = self.molecule.atom_labels

        for i, atom_label in enumerate(atom_labels):
            G.add_node(i, label = atom_label, coord=coords[i])
        n_atoms = len(atom_labels)

        for i in range(n_atoms):
            for j in range(i+1,n_atoms):
                dist_ij = np.linalg.norm(coords[i]-coords[j])
                G.add_edge(i,j,weight=dist_ij)
        self.euclidean_graph = G
        return G
    
    def draw_euclidean_graph(self):
        """ 
        Draws the weighted euclidean graph of the molecule
        """

        if self.euclidean_graph is None:
            self.build_euclidean_graph()

        plt.figure(figsize=(8, 6)) 
        pos = nx.spring_layout(self.euclidean_graph)
        
        # Nodes
        nx.draw_networkx_nodes(self.euclidean_graph, pos, node_color='lightgreen', node_size=500)

        # Edges
        nx.draw_networkx_edges(self.euclidean_graph, pos, edgelist=self.euclidean_graph.edges(), width=3)
        nx.draw_networkx_edges(
            self.euclidean_graph,
            pos,
            edgelist=self.euclidean_graph.edges(),
            width=3,
            alpha=0.5,
            edge_color="b",
            style="dashed",
        )

        # Node labels
        nx.draw_networkx_labels(self.euclidean_graph,pos, font_family="sans-serif")

        # Edges with weights
        edge_labels = nx.get_edge_attributes(self.euclidean_graph, 'weight')
        nx.draw_networkx_edge_labels(self.euclidean_graph, pos, edge_labels=edge_labels)


        ax = plt.gca()
        ax.margins(0.1)
        plt.axis("off")
        plt.tight_layout()
        plt.savefig("euclidean_graph.png")

