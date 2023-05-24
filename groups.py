"""
This a package for finite groups.
We are using Python 3.10.

It contains the following classes:
FiniteGroups: The main class for this module. Defines a group an several useful methods.
FiniteGroups has the following subclasses:
CyclicGroup: Class of cyclic groups Z/nZ.
GL: Class of groups GL_n(Z/qZ).
SL: Class of groups SL_n(Z/qZ).
PSL: Class of groups PSL_n(Z/qZ).
PGL: Class of groups PGL_n(Z/qZ).
"""

import numpy as np
import math
import itertools
import time
from numpy import matrix
from numpy import linalg

class FiniteGroup:
    """
    This class implements finite groups.

    Initialization:
    elements: Elements of the group. We usually choose the elements to be a set or list of tuples, as this type is hashable. 
    identity: Identity element of the group.
    operation: The operation ouputs an element for an input of two elements.

    Methods:
    inverse: Calculates the inverse of a given element.
    is_group: Checks if multiplication is well-defined and inverses exist.
    """
    
    def __init__(self, elements, identity, operation):
        """
        Initalizes FiniteGroup.

        Arguments:
        elements: Elements of the group. We usually choose the elements to be a set or list of tuples, as this type is hashable. 
        identity: Identity element of the group.
        operation: The operation ouputs an element for an input of two elements.

        Return: 
        The group determined by the above data.
        """

        self.elements = elements
        self.identity = identity
        self.operation = operation
    
    def inverse(self, g):
        """
        g: Find inverse of the group element g. If no inverse exists, we return False.
        """

        for h in self.elements:
            if self.operation(g, h) == self.identity:
                return h
        return False
    
    def is_group(self):
        """
        Checks if multiplication is well-defined and the inverses exist.  
        """

        for e1 in self.elements:
            for e2 in self.elements:
                if self.operation(e1,e2) not in self.elements: return False

        for e1 in self.elements:
            if not self.inverse(e1): return False
            if self.inverse(e1) not in self.elements: return False
        
        return True
    
    
class CyclicGroup(FiniteGroup):
    """This class implements cyclic groups Z/nZ."""

    def __init__(self,n):
        """
        Class of cyclic groups Z/nZ.

        Argument: 
        n: Poistive integer.

        Return:
        The group Z/nZ as an FiniteGroup object.
        """

        set = range(0,n)
        def op(a,b): return (a+b)%n
        super().__init__(set, 0, op)


class GL(FiniteGroup):
    """
    This class implements the groups GL_n(Z/qZ) as a subclass of FiniteGroup.

    Initialization: 
    n: Poistive integer.
    q: Poistive integer.
    """

    def __init__(self,n,q):
        """
        Initialization of GL(n,q).

        Arguments: 
        n: Poistive integer.
        q: Poistive integer.

        Return:
        The group GL_n(Z/qZ) as a FiniteGroup object.
        """

        # We initalize the set of tuples coresponding to matrices of non-zero determinant.
        # A (n,n)-array corresponds to a tuple of length n**2 by writing the rows of the matrix behind each other. 
        set_of_matrices = set()
        set_of_tuples = set(itertools.product(range(0,q),repeat = n**2))
        for t in set_of_tuples: 
            A = np.reshape(t,(n,n))
            if round(linalg.det(A))%q != 0:
                set_of_matrices.add(t) 
    
        # To define the group operation, we reshape the tuples and apply np.matmul.
        def op(A,B): 
            A_matrix = np.reshape(A,(n,n))
            B_matrix = np.reshape(B,(n,n))
            C = np.matmul(A_matrix, B_matrix)%q
            return tuple(C.reshape(n**2))

        super().__init__(set_of_matrices, tuple(np.identity(n,dtype = int).reshape(n**2)), op)

# We assert that we have implemented GL correctly.
# For cardinality formulas of linear groups we refer to https://groupprops.subwiki.org/wiki/Order_formulas_for_linear_groups
G = GL(2,3)
assert G.is_group(), "GL is not a group."
assert len(G.elements) == 48, "GL has the wrong cardinality."
assert G.inverse((1,1,0,1)) == (1,2,0,1), "GL does not calculate inverses correctly."
G = GL(3,2)
assert G.is_group(), "GL is not a group."
assert len(G.elements) == 168, "GL has the wrong cardinality."

class SL(FiniteGroup):
    """
    This class implements the groups SL_n(Z/qZ) as a subclass of FiniteGroup.

    Initialization: 
    n: Poistive integer.
    q: Poistive integer.
    """

    def __init__(self,n,q):
        """
        Initialization of SL(n,q).

        Arguments: 
        n: Poistive integer.
        q: Poistive integer.

        Return:
        The group SL_n(Z/qZ) as an FiniteGroup object.
        """

        # This class has the same initalization as GL(n,q), with the only difference to require the matrices to have determinant to be 1.
        G = GL(n,q)
        set_of_matrices = set()
        for t in G.elements: 
            A = np.reshape(t,(n,n))
            if round(linalg.det(A))%q == 1:
                set_of_matrices.add(t) 

        super().__init__(set_of_matrices, G.identity, G.operation)

# We assert that we have implemented SL correctly.
# For cardinality formulas of linear groups we refer to https://groupprops.subwiki.org/wiki/Order_formulas_for_linear_groups
G = SL(2,3)
assert G.is_group(), "GL is not a group."
assert len(G.elements) == 24, "GL has the wrong cardinality."
G = SL(3,2)
assert G.is_group(), "SL is not a group."
assert len(G.elements) == 168, "SL has the wrong cardinality."

class PGL(FiniteGroup):
    """
    This class implements the groups PGL_n(Z/qZ) as a subclass of FiniteGroup.
    For a definition of PGL_n(Z/qZ) we refer to https://en.wikipedia.org/wiki/Projective_linear_group.

    Initialization: 
    n: poistive integer.
    q: prime number.
    """

    def __init__(self,n,q):
        """
        Initialization of PGL(n,q).

        Arguments: 
        n: poistive integer.
        q: prime number.

        Return:
        The group PGL_n(Z/qZ) as a FiniteGroup object.
        """
        
        # We start with the initialization of GL(n,q) and endow the latter group with an equivalence relation.
        G = GL(n,q)
        def is_equivalent(A,B):
            "Checks if two matrices in GL(n,q) are equivalent in PGL(n,q)."
            for l in range(1,q):
                if tuple([(l*a)%q for a in A]) == B: return True
            return False
        
        # For each element of the resulting equivalence classes, we choose one element as a representative.
        # We implement the resulting projection mapping from GL(n,q) to PGL(n,q) as a dictionary.
        # To make sure that the identity matrix is chosen as the representative of its equivalence class we initialize the dicitionary as follows.
        projection = {G.identity: G.identity}
        for A in G.elements:
            for B in set(projection.values()):
                if is_equivalent(A,B):
                    projection[A] = B
            projection.setdefault(A,A)
        
        def op(A,B): 
            return projection[G.operation(A,B)]

        super().__init__(set(projection.values()),G.identity,op) 

# We assert that we have implemented PGL correctly.
# For cardinality formulas of linear groups we refer to https://groupprops.subwiki.org/wiki/Order_formulas_for_linear_groups
G = PGL(2,3)
assert G.is_group(), "GL is not a group."
assert len(G.elements) == 24, "GL has the wrong cardinality."
G = PGL(3,2)
assert G.is_group(), "SL is not a group."
assert len(G.elements) == 168, "SL has the wrong cardinality."

# PSL_n(Z/qZ) is also implemented. To do so, we need to find the nth roots of unity of Z/qZ.
def roots_of_unity(n,q):
    """
    Returns the set of the nth roots of unity in Z/qZ.
    
    Arguments: 
    n: Poistive integer.
    q: prime number.

    Retrun: 
    Set of nth roots of unity in Z/qZ viewed as a subset of {0,1,2,3,...,q-1}
    """

    roots = {1}
    for k in range(2,q):
        if (k**n)%q ==1:
            roots.add(k)
    return roots

class PSL(FiniteGroup):
    """
    This class implements the groups PSL_n(Z/qZ) as a subclass of FiniteGroup.
    For a definition of PSL_n(Z/qZ) we refer to https://en.wikipedia.org/wiki/Projective_linear_group.

    Initialization: 
    n: Poistive integer.
    q: prime number.
    """

    def __init__(self,n,q):
        """
        Initialization of PSL(n,q).

        Arguments: 
        n: Poistive integer.
        q: prime number.

        Return:
        The group PSL_n(Z/qZ) as a FiniteGroup object.
        """

        # We use the same implementation as PGL(n,q), with the only difference of starting with SL(n,q) and using a different equivalence relation.
        G = SL(n,q)
        roots = roots_of_unity(n,q)
        def is_equivalent(A,B):
            "Checks if two matrices in SL(n,q) are equivalent in PSL(n,q)."
            for l in roots:
                if tuple([(l*a)%q for a in A]) == B: return True
            return False

        projection = {G.identity: G.identity}
        for A in G.elements:
            for B in set(projection.values()):
                if is_equivalent(A,B):
                    projection[A] = B
            projection.setdefault(A,A)
        
        def op(A,B): 
            return projection[G.operation(A,B)]

        super().__init__(set(projection.values()),G.identity,op)   

# We assert that we have implemented PSL correctly.
# For cardinality formulas of linear groups we refer to https://groupprops.subwiki.org/wiki/Order_formulas_for_linear_groups
G = PSL(2,3)
assert G.is_group(), "GL is not a group."
assert len(G.elements) == 12, "GL has the wrong cardinality."
G = PSL(3,2)
assert G.is_group(), "SL is not a group."
assert len(G.elements) == 168, "SL has the wrong cardinality."