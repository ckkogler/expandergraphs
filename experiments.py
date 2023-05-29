import numpy as np
import groups as gp
import cayleygraphs as cg
import LPS

print("We show the following quantities:")
print("k & Number of elements & Diameter & gamma & lambda_1 & girth & girth ratio & mean injectivity radius")

for k in range(10,17):
    G = gp.SL(2,2*k+1)
    set = [(1,2,0,1), (1,0,2,1)]
    Cay = cg.CayleyGraph(G,set)
    noe = len(Cay.group.elements)
    diam = cg.diameter(Cay)
    asg = cg.adjacency_spectral_gap(Cay)
    lsg = cg.laplacian_spectral_gap(Cay)
    gth = cg.girth(Cay)
    gthratio = gth*np.log(3)/np.log(noe)
    minj = cg.injectivity_radius(Cay)
    print(f"{2*k+1} & {noe} & {diam} &  {round(asg,3)} & {round(lsg,3)} & {gth} & {round(gthratio,3)} & {minj} \\\\")
