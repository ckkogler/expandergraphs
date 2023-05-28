import groups as gp
import cayleygraphs as cg
import LPS

for k in range(20,100):
    G = gp.SL(2,k)
    set = [(1,1,0,1), (1,0,1,1)]
    Cay = cg.CayleyGraph(G,set)
    noe = len(Cay.group.elements)
    asg = cg.adjacency_spectral_gap(Cay)
    lsg = cg.laplacian_spectral_gap(Cay)
    gth = cg.girth(Cay)
    print(f"k = {k}, noe = {noe}, girth = {gth}, asg = {asg}, lsg = {lsg}")