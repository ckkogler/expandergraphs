"""
This a package for Cayley graphs, which we are building on the groups package.  
We are using Python 3.10.

It contains the following classes:
CayleyGraph: 
"""

import numpy as np
import networkx as nx 
import math
import itertools
import time
from numpy import matrix
from numpy import linalg

class CayleyGraph:
    
    def __init__(self,group,set):
        """
        Given an object in the class of finite_groups and a set of elements in finite_group, we initalize the associated Cayley graph.
        finite_group:   Underlying group.       
        set:            Subset of finite_group that determines the Cayley graph. 
        """

        graph = nx.Graph()
        graph.add_nodes_from(group.elements)
        
        for g in group.elements:
            for s in set:
                graph.add_edge(g,group.operation(s,g))

        self.graph = graph
        self.group = group
        self.set = set
        self.degree = graph.degree(group.elements[0])

    def adjacency_spectrum(self):
        """
        Returns the spectrum of the adjacency matric of the Cayley graph, orderd in decreasing order. For convienience, we round the eigenvalues to 8 significant digits.
        """
        eigenvalues = nx.adjacency_spectrum(self.graph)
        eigenvalues = [round(eig.real,8) for eig in eigenvalues]
        eigenvalues.sort(reverse = True)
        return eigenvalues
    
    def is_bipartite(self):
        """
        Returns whether or not the Cayley graph is bipartite.
        """
        abs_eigenvalues = [abs(eig) for eig in self.adjacency_spectrum()]
        abs_eigenvalues.sort(reverse = True)
        if abs_eigenvalues[1] == self.degree: return True
        return False 
    
    def adjacency_spectral_gap(self):
        """
        Returns the spectral gap of the adjacency matrix. 
        """
        abs_eigenvalues = [abs(eig) for eig in self.adjacency_spectrum()]
        abs_eigenvalues.sort(reverse = True)
        return round(self.degree - abs_eigenvalues[1],8)
    
    def adjacency_strong_spectral_gap(self):
        """
        Returns the strong spectral gap of the adjacency matrix. 
        """
        abs_eigenvalues = [abs(eig) for eig in self.adjacency_spectrum()]
        abs_eigenvalues.sort(reverse = True)
        if self.is_bipartite(): return round(self.degree - abs_eigenvalues[2],8)
        else: return round(self.degree - abs_eigenvalues[1],8)
    
    def laplacian_spectrum(self):
        """
        Returns the spectrum of the Laplacian of the Cayley graph, orderd in increasing order. For convienience, we round the eigenvalues to 8 significant digits.
        """
        eigenvalues = nx.laplacian_spectrum(Cay.graph)
        abs_of_eigenvalues = [round(eig,8) for eig in eigenvalues]
        abs_of_eigenvalues.sort()
        return abs_of_eigenvalues


    def laplacian_spectral_gap(self):
        """
        Returns the first eigenvalue of the Laplacian.
        """
        return self.laplacian_spectrum()[1]