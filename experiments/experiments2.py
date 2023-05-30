import os,sys,inspect
current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)


import numpy as np
import groups as gp
import cayleygraphs as cg
import LPS

print("We show the following quantities:")
print("k & Number of elements & Diameter & gamma & lambda_1 & girth & girth ratio & mean injectivity radius")

for k in range(2,35):
    if k%3!=0:
        G = gp.SL(2,k)
        set = [(1,1,0,1), (1,0,1,1), (1,3,0,1), (1,0,3,1)]
        Cay = cg.CayleyGraph(G,set)
        noe = len(Cay.group.elements)
        diam = cg.diameter(Cay)
        asg = cg.adjacency_spectral_gap(Cay)
        lsg = cg.laplacian_spectral_gap(Cay)
        gth = cg.girth(Cay)
        deg = Cay.degree
        gthratio = gth*np.log(deg - 1)/np.log(noe)
        minj = cg.injectivity_radius(Cay)
        print(f"{k} & {noe} & {diam} &  {round(asg,3)} & {round(lsg,3)} & {gth} & {round(gthratio,3)} & {minj} \\\\")
