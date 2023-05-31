"""
With this script we study the 6-regular Cayley graph on the group SL(2,k).
We use the generating set [(1,\pm 2, 0, 1), (1, 0, \pm 2, 0), (2,0,0,2^{-1}), (2^{-1},0,0,2)].
The parameter k ranges from 3 to 35 and we compute several graph theoretic properties. 
"""

#We first fix the system path such that we import the modules from the parent directory.
import os,sys,inspect
current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

import math
import numpy as np
import groups as gp
import cayleygraphs as cg
import LPS

print("We show the following quantities:")
print("k & Number of elements & Diameter & gamma & lambda_1 & girth & girth ratio & mean injectivity radius")

#In 
def inverse(n,k):
    """
    Calculates the multiplicative inverse of n in Z/kZ. The multiplicative inverse exists if and only if gcd(n,k) = 1.
    """
    if math.gcd(n,k) != 1: 
        return "Inverse does not exist."

    for m in range(1,k):
        if (m*n)%k == 1: 
            return m 

for k in range(3,35):
    #2 needs to be invertible and therefore we only consider odd k.
    if k%2!=0:

        #Construction of Cayley graph. 
        G = gp.SL(2,k)
        set = [(1,2,0,1), (1,0,2,1), (2,0,0,inverse(2,k))]
        Cay = cg.CayleyGraph(G,set)

        #Calculation of graph theoretic properties.
        noe = len(Cay.group.elements)
        diam = cg.diameter(Cay)
        asg = cg.adjacency_spectral_gap(Cay)
        lsg = cg.laplacian_spectral_gap(Cay)
        gth = cg.girth(Cay)
        deg = Cay.degree
        gthratio = gth*np.log(deg - 1)/np.log(noe)
        minj = cg.injectivity_radius(Cay)

        #We output the graph theoretic quantities such that it can be used for a Latex table. 
        print(f"{k} & {noe} & {diam} &  {round(asg,3)} & {round(lsg,3)} & {gth} & {round(gthratio,3)} & {minj} \\\\")
