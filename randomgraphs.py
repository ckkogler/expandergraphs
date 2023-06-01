"""
In this script we construct random graphs, with the hope of establishing Ramanujan graphs.
"""
import numpy as np
from numpy import linalg
import networkx as nx

def random_permutation_matrix(n):
    """
    Returns a random (n,n) permutation matrix
    
    Input: 
    n: Positive integer.

    Output: 
    Random (n,n) permutation matrix
    """
    
    perm = np.random.permutation(n)
    A = np.zeros((n,n), dtype=int)
    for k in range(0,n):
        A[k,perm[k]] = 1
    return A

def random_adjacency_matrix(n,d):
    """
    Returns a random adjacecy matrix corresponding to a graph with n vertices and degree 2*d.
    
    Input: 
    n: positive integer.
    d: positive integer.

    Return: 
    Random adjacency matrix corresponding to a graph with n vertices and degree 2*d.
    """

    A = np.zeros((n,n), dtype=int)
    for k in range(d):
        A += random_permutation_matrix(n)
    return A + A.T


def is_simple(A):
    """
    Returns if a random adjacecy matrix is a simple graph.
    
    Input: 
    A: (n,n) matrix.

    Return: 
    Returns if A corresponds to a matrix of a connected graph.
    """
    return A.max()==1


def random_simple_adjacency_matrix(n,d):
    """
    Returns a random adjacecy matrix corresponding to a graph with n vertices and degree 2*d.
    
    Input: 
    n: positive integer.
    d: positive integer.


    Return: 
    Random adjacency matrix corresponding to a graph with n vertices and degree d.
    """

    A = np.zeros((n,n), dtype=int)
    count = 0
    while count < d:
        R = random_permutation_matrix(n)
        A_new = A + R + R.T
        if is_simple(A_new):
            A = A_new
            count += 1
    return A 

def sorted_eigenvalues(A):
    return np.round(np.sort(linalg.eigvalsh(A))[::-1],6)    

def is_connected(A):
    """
    Returns whether or not the Cayley graph is connected.
    We use the condition that a graph is connected if and only if the second largest eigenvalue of the adjacency matrix is strictly less than 1.
    """
    if sorted_eigenvalues(A)[1] < 1: return True
    return False

def is_bipartite(A):
    """
    Returns whether or not a connected Cayley graph is bipartite.
    We use the spectral condition that a connected graph is bipatite if and only if -1 is an eigenvalue of the graph.
    """
    if is_connected(Cay) == False: return "Graph not connected."
    if sorted_eigenvalues(A) == -1: return True
    return False

def diameter(A):
    diameter = 1
    B = A
    while B.min() ==0:
        B = np.matmul(B,A)
        diameter += 1
    return diameter

def adjacency_spectral_gap(A):
    """
    Returns the strong spectral gap of the normalized adjacency matrix of a connected graph. 
    """
    return 1 - max(abs(sorted_eigenvalues(A)[1]),abs(sorted_eigenvalues(A)[-1]))

def adjacency_star(A):
    """
    Returns the spectral gap of the normalized adjacency matrix of a connected graph. 
    This is the distance to one of the largest modulus that is not one.  
    """
    if is_bipartite(A): 
            return 1 - max(abs(sorted_eigenvalues(A)[1]),abs(sorted_eigenvalues(A)[-2]))
    else: return adjacency_spectral_gap(A)

def laplacian_spectral_gap(A):
    """
    Returns the first eigenvalue of the normalized Laplacian. 
    We note that the spectrum of the normalized Laplacian is simply the 1 - (spectrum of the normalized adjacency matrix).
    """
    return (1 - sorted_eigenvalues(A)[1])

def girth(A):
    """
    Returns girth (length of shortest cycle) of the graph associated to A.
    """
    gph = nx.from_numpy_array(A)
    return min(len(cycle) for cycle in nx.cycle_basis(gph))

def size_of_sphere(d,r):
    """
    Calculates the size of the d-sphere around the identity in the free group.
    """
    if r in {0,1}: return d**r
    else: return d*((d-1)**(r-1))

def injectivity_radius(A):
    n = A.shape[0]
    d = A[:,1].sum()
    output = np.zeros((n,2),int) - 1
    B = A
    r = 1
    while output.min() == -1:
        for i in range(0,n):
            if output[i,1] == -1 and list(B[i,:]).count(1) < size_of_sphere(d,r):
                output[i,0] = r-1
                output[i,1] = 0
        r += 1
        B = np.matmul(B,A)
    return output[:,0]

def mean_injectivity_radius(A):
    return np.mean(injectivity_radius(A))
