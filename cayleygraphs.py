"""
This a package for Cayley graphs, which we are building on the groups package.  
We are using Python 3.10.

It contains the following classes:
CayleyGraph: An implementation of Cayley graphs, building on the groups package.

Taking a Cayley graph as input we define the following functions:
is_connected: Returns wether or not the input Cayley graph is connected.
is_bipartite: Returns wether or not the input Cayley graph is bipartite.
adjacency_spectral_gap: Calculates the (strong) spectral gap of the normalized adjacency matrix.
laplacian_spectral_gap: Calculates the spectral gap of the normalized Laplacian.
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

        #We furthermore calculate the eigenvalues of our graph. 
        #For convenience we sort the eigenvalues in descending order, normalize them and round the values to 8 digits.
        eigenvalues = list(nx.adjacency_spectrum(self.graph))
        eigenvalues.sort(reverse = True)
        self.degree = round(eigenvalues[0].real)
        self.eigenvalues = [round(eig.real/self.degree,8) for eig in eigenvalues]
    
def is_connected(Cay):
        """
        Returns wether or not the Cayley graph is connected.
        We use the condition that a graph is connected if and only if the second largest eigenvalue of the adjacency matrix is strictly less than 1.
        """
        if Cay.eigenvalues[1] < 1: return True
        return False

def is_bipartite(Cay):
        """
        Returns wether or not a connected Cayley graph is bipartite.
        We use the spectral condition that a connected graph is bipatite if and only if -1 is an eigenvalue of the graph.
        """
        if is_connected(Cay) == False: return "Graph not connected."
        if Cay.eigenvalues[-1] == -1: return True
        return False

def adjacency_spectral_gap(Cay):
        """
        Returns the (strong) spectral gap of the normalized adjacency matrix. 
        """
        return max(1 - Cay.eigenvalues[1], abs(-1 - Cay.eigenvalues[-1]))

def laplacian_spectral_gap(Cay):
        """
        Returns the first eigenvalue of the normalized Laplacian. 
        We note that the spectrum of the normalized Laplacian is simply the 1 - (spectrum of the normalized adjacency matrix).
        """
        return (1 - Cay.eigenvalues[1])