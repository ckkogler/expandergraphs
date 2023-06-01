"""
With this script we study random graphs
"""

#We first fix the system path such that we import the modules from the parent directory.
import os,sys,inspect
current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

import numpy as np
import randomgraphs as rg

print("We show the following quantities:")
print("Diameter & gamma & lambda_{star} & lambda_1 & girth & girth ratio & mean injectivity radius")

n = 256
k = 2
m = 100

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
    gthratio = gth*np.log(2*k-1)/np.log(n)
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
print(f"mean & {list_of_diam.mean()} & {round(list_of_asg.mean(),3)} & {round(list_of_ast.mean(),3)} &  \
{round(list_of_lsg.mean(),3)} & {list_of_gth.mean()} & {round(list_of_gthratio.mean(),3)} & {list_of_minj.mean()} \\\\")
print(f"std & {round(list_of_diam.std(),3)} & {round(list_of_asg.std(),3)} & {round(list_of_ast.std(),3)} &  \
{round(list_of_lsg.std(),3)} & {round(list_of_gth.std(),3)} & {round(list_of_gthratio.std(),3)} & {round(list_of_minj.std(),3)} \\\\")
print(f"max & {list_of_diam.min()} & {round(list_of_asg.max(),3)} & {round(list_of_ast.max(),3)} &  \
{round(list_of_lsg.max(),3)} & {list_of_gth.max()} & {round(list_of_gthratio.max(),3)} & {list_of_minj.max()} \\\\")