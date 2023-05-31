"""
This a package for Cayley graphs, which we are building on the groups package.  
We are using Python 3.10.

It contains the following classes:
CayleyGraph: An implementation of Cayley graphs, building on the groups package.

Taking a Cayley graph as input we define the following functions:
is_connected: Returns whether or not the input Cayley graph is connected.
is_bipartite: Returns whether or not the input Cayley graph is bipartite.
diameter: Calculates the diameter of the Cayley graph.
adjacency_spectral_gap: Calculates the strong spectral gap of the normalized adjacency matrix.
adjacency_star: Calculates the spectral gap of the adjacency matrix. 
laplacian_spectral_gap: Calculates the spectral gap of the normalized Laplacian.
girth: Calculates the girth of the Cayley graph (i.e. the length of the shortest cycle).
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
        Returns whether or not the Cayley graph is connected.
        We use the condition that a graph is connected if and only if the second largest eigenvalue of the adjacency matrix is strictly less than 1.
        """
        if Cay.eigenvalues[1] < 1: return True
        return False

def is_bipartite(Cay):
        """
        Returns whether or not a connected Cayley graph is bipartite.
        We use the spectral condition that a connected graph is bipatite if and only if -1 is an eigenvalue of the graph.
        """
        if is_connected(Cay) == False: return "Graph not connected."
        if Cay.eigenvalues[-1] == -1: return True
        return False

def diameter(Cay):
        """
        Returns the diameter of the Cayley graph.
        """
        if is_connected(Cay)==False: return "Not connected."
        return nx.diameter(Cay.graph)

def adjacency_spectral_gap(Cay):
        """
        Returns the strong spectral gap of the normalized adjacency matrix of a connected graph. 
        """
        return 1 - max(abs(Cay.eigenvalues[1]),abs(Cay.eigenvalues[-1]))

def adjacency_star(Cay):
        """
        Returns the spectral gap of the normalized adjacency matrix of a connected graph. 
        This is the distance to one of the largest modulus that is not one.  
        """
        if Cay.is_bipartite: 
                return 1 - max(abs(Cay.eigenvalues[1]),abs(Cay.eigenvalues[-2]))
        else: return adjacency_spectral_gap(Cay)

def laplacian_spectral_gap(Cay):
        """
        Returns the first eigenvalue of the normalized Laplacian. 
        We note that the spectrum of the normalized Laplacian is simply the 1 - (spectrum of the normalized adjacency matrix).
        """
        return (1 - Cay.eigenvalues[1])

def girth(Cay):
        """
        Returns girth (length of shortest cycle) of the Cayley graph.
        """
        return min(len(cycle) for cycle in nx.cycle_basis(Cay.graph))

def size_of_ball(d,r):
    if r in {0,1}: return d**r
    else: return d*((d-1)**(r-1))

def injectivity_radius(Cay):
    A = nx.to_numpy_array(Cay.graph)  
    n = len(Cay.group.elements)
    d = Cay.degree
    r = 1
    B = A
    while list(B[0,:]).count(1) == size_of_ball(d,r):
        B = np.matmul(B,A)
        r += 1
    return r-1
