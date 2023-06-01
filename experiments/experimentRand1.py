"""
With this script we study random graphs
"""

#We first fix the system path such that we import the modules from the parent directory.
import os,sys,inspect
current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

import numpy as np
import randongraphs as rg

print("We show the following quantities:")
print("Diameter & gamma & lambda_1 & girth & girth ratio & mean injectivity radius")

n = 64
k = 2
m = 1000000

list_of_diam = np.array([])
list_of_asg = np.array([])
list_of_ast = np.array([])
list_of_lsg = np.array([])
list_of_gth = np.array([])
list_of_gthratio = np.array([])
list_of_minj = np.array([])

for l in range(0,m):

    #Sampling of matrix.
    A = rg.random_simple_adjacency_matrix(n,k)

    #Calculation of graph theoretic properties.
    diam = rg.diameter(A)
    asg = rg.adjacency_spectral_gap(A)
    ast = rg.adjacency_star(A)
    lsg = rg.laplacian_spectral_gap(A)
    gth = rg.girth(A)
    gthratio = gth*np.log(2*k-1)/np.log(noe)
    minj = rg.mean_injectivity_radius(A)

    #Add calculated values to list.
    list_of_diam = np.append(list_of_diam, diam)
    list_of_asg = np.append(list_of_asg, asg)
    list_of_ast = np.append(list_of_ast, ast)
    list_of_lsg = np.append(list_of_lsg, lsg)
    list_of_gth = np.append(list_of_gth, gth)
    list_of_gthratio = np.append(list_of_gthratio, gthratio)
    list_of_minj = np.append(list_of_minj, minj)

    #We output the graph theoretic quantities such that it can be used for a Latex table.  
    print(f"Mean: {diam} & {round(asg,3)} & {round(lsg,3)} & {gth} & {round(gthratio,3)} & {minj} \\\\")