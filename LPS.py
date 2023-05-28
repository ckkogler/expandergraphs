"""
In this module we are implementing the Ramanujan graphs constructed by Lubotzky-Phillips-Sarnak im 1988.
The original paper can be found here: https://link.springer.com/article/10.1007/BF02126799.
The used implmentation is also described on Wikipedia: https://en.wikipedia.org/wiki/Ramanujan_graph.

We are building the LPS graphs on the cayleygraphs package. 
"""

import numpy as np
import itertools
import math
import networkx as nx
import groups as gp
import cayleygraphs as cy

#We first implememt some basic number theoric functions that are necessary to constuct the LPS graphs.

def is_prime(n):
    """
    Checks whether the integer n is prime or not.
    """
    if n  == 1: return False
    if n == 2: return True
    for k in range(2,int(np.sqrt(n))+1):
        if n%k == 0: return False
    return True

def is_quadratic_residue(p,m):
    """
    Calculates wheter the integer p is a quadratic residue modulo the integer m.
    For a definition of being a quadratic residue we refer to: https://en.wikipedia.org/wiki/Quadratic_residue.
    """
    quadratic_residues = {(x**2)%m for x in range(1,m) if math.gcd(x,m) == 1}
    if (p%m) in quadratic_residues: return True
    return False

def sum_of_four_squares(p):
    """
    Given an integer p, we calculate all the ways p can be written as a sum of four squares.
    If p is prime, there are 8*(p + 1) such solutions by Jacobi's four square theorem (https://en.wikipedia.org/wiki/Jacobi%27s_four-square_theorem).

    Input:
    p: Integer.

    Output:
    List of all tuples (x_0, x_1, x_2, x_3) satisfying p == x_0**2 + x_1**2 + x_2**2 + x_3**2 with x_i being integers.
    """
    # We first of all tuples (x_0, x_1, x_2, x_3) with 0 <= x_i**2 =< p for all i. 
    list_of_four_tuples = itertools.product(range(-int(np.sqrt(p))-1,int(np.sqrt(p))+1),repeat = 4)
    output = []
    for tuple in list_of_four_tuples:
        if tuple[0] > 0 and tuple[0]%2 == 1 and tuple[1]%2 == 0 and tuple[2]%2 == 0 and tuple[3]%2 == 0:
            if p == (tuple[0]**2 + tuple[1]**2 + tuple[2]**2 + tuple[3]**2):
                output.append(tuple)

    return output



class LPS(cy.CayleyGraph):
    """
    Given two distict primes p and q both satisfying == 1 mod 4, we constuct the (p,q) LPS graph.  
    The original paper by Lubotzky-Phillips-Sarnak can be found here: https://link.springer.com/article/10.1007/BF02126799.
    The used implmentation can be found on Wikipedia: https://en.wikipedia.org/wiki/Ramanujan_graph.

    Input: 
    p: prime %4 == 1.
    q: prime %4 == 1 different from p.

    Return:
    (p,q) LPS graph as object in cy.CayleyGraphs. The resulting graph has degree (p + 1). 
    """
    def __init__(self,p,q):
        """
        Let p and q be distinct primes = 1 mod 4. 
        """
        if is_quadratic_residue(p,q):
            group = gp.PSL(2,q)
            iota = min([i for i in range(0,q) if i**2%q == q-1])
            p_inv = min([i for i in range(0,q) if (i*p)%q == 1])
            quad_res = min([i for i in range(0,q) if (i**2)%q == p_inv])
            set = [(quad_res*(tuple[0] + iota*tuple[1])%q, quad_res*(tuple[2] + iota*tuple[3])%q, quad_res*(-tuple[2] + iota*tuple[3])%q, quad_res*(tuple[0] - iota*tuple[1])%q) for tuple in sum_of_four_squares(p)]

        if is_quadratic_residue(p,q) == False:
            group = gp.PGL(2,q)
            iota = min([i for i in range(0,q) if i**2%q == q-1])
            set = [((tuple[0] + iota*tuple[1])%q, (tuple[2] + iota*tuple[3])%q, (-tuple[2] + iota*tuple[3])%q, (tuple[0] - iota*tuple[1])%q) for tuple in sum_of_four_squares(p)]

        super().__init__(group,set)
