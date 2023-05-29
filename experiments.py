import numpy as np
import groups as gp
import cayleygraphs as cg
import LPS

for k in range(3,16):
    G = gp.SL(2,k)
    set = [(1,2,0,1), (1,0,2,1)]
    Cay = cg.CayleyGraph(G,set)
    noe = len(Cay.group.elements)
    asg = cg.adjacency_spectral_gap(Cay)
    lsg = cg.laplacian_spectral_gap(Cay)
    gth = cg.girth(Cay)
    gthratio = gth*np.log(3)/np.log(noe)
    print(f"k = {k}, noe = {noe}, girth = {gth}, asg = {asg}, lsg = {lsg}, gthrat = {gthratio}")