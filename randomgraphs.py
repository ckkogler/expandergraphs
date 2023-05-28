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

def spectral_gap(A):
    deg = A[0,:].sum()
    ev = np.round(np.linalg.eigvals(A),5)/deg
    lambda1 = 1 - max(ev[1:])
    lambda2 = 1 - abs(min(ev[1:]))
    return min(lambda1,lambda2)

def diameter(A):
    diameter = 1
    B = A
    while B.min() ==0:
        B = np.matmul(B,A)
        diameter += 1
    return diameter

def girth(A):
    G = nx.from_numpy_array(A)
    girth = min(len(cycle) for cycle in nx.cycle_basis(G))
    return girth
        

def ramanujan_bound(d):
    return 1 - 2*np.sqrt(d-1)/d



def find_ramanujan(n,d,t):
    r = 0
    for k in range(1,t):
        r1 = spectral_gap(random_simple_adjacency_matrix(n,d))
        if r1 > r: 
                r = r1
                print(r)
                if r > ramanujan_bound(2*d):
                    print("Is Ramanujan!")


def count_simple_ramanujan(n,d,t):   
    count = 0
    for k in range(1,t):
        r1 = spectral_gap(random_simple_adjacency_matrix(n,d))
        if r1 > ramanujan_bound(2*d): 
            count += 1
            print("Ramanujan Graph found.")
    print("We sampled ", t, "random simple graphs.")
    print("Out of them", count, "were Ramanujan.")      


def find_graph_with_low_grith(n,d,t):   
    girth = 0
    for k in range(0,t):
        girth_new = spectral_gap(random_simple_adjacency_matrix(n,d))
        if r1 > ramanujan_bound(2*d): 
            count += 1
            print("Ramanujan Graph found.")
    print("We sampled ", t, "random simple graphs.")
    print("Out of them", count, "were Ramanujan.") 