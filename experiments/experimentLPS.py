"""
With this script we study the Ramanujan graphs constructed by Lubotzky-Phillips-Sarnak (LPS).
"""

#We first fix the system path such that we import the modules from the parent directory.
import os,sys,inspect
current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

import numpy as np
import groups as gp
import cayleygraphs as cg
import LPS

print("We show the following quantities:")
print("(p,q) & Number of elements & Diameter & gamma & lambda_star & lambda_1 & girth & girth ratio & mean injectivity radius")

for p in range(20,50):
    for q in range(3,20):
        if LPS.is_prime(p) and p%4==1:
            if LPS.is_prime(q) and q%4==1 and q!=p: 
                #Construction of Cayley graph. It suffices to take the set [(1,1,0,1), (1,0,1,1)]
                Cay = LPS.LPS(p,q)

                #Calculation of graph theoretic properties.
                noe = len(Cay.group.elements)
                diam = cg.diameter(Cay)
                asg = cg.adjacency_spectral_gap(Cay)
                lst = cg.adjacency_star(Cay)
                lsg = cg.laplacian_spectral_gap(Cay)
                gth = cg.girth(Cay)
                deg = Cay.degree
                gthratio = gth*np.log(deg-1)/np.log(noe)
                minj = cg.injectivity_radius(Cay)

                #We output the graph theoretic quantities such that it can be used for a Latex table.  
                print(f"{p,q} & {noe} & {diam} &  {round(asg,3)} & {round(lst,3)} & {round(lsg,3)} & {gth} & {round(gthratio,3)} & {minj} \\\\")
