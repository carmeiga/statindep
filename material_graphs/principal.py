#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 15:52:35 2021

@author: carlos
"""

#http://www.crystallography.net/cod/search.html
#http://hoffman.physics.harvard.edu/materials/CuprateIntro.php
#https://materialsproject.org/materials/mp-5986/#

# sed -r 's/.cif//' idprop > id_prop.csv

import os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import random
from matplotlib.backends.backend_pdf import PdfPages

os.chdir('/home/carlos/repodir/cgcnn-master')
from cgcnn.data import CIFData


dataset = CIFData()

adjacencymat=[]
nomes_materiais=[]

for k in range(len(dataset)):
    try:
        
        mitmat=dataset[k][0][1]
    except ValueError:
        continue
    nomes_materiais.append(dataset[k][2])
    nver=mitmat.shape[0]
    adj=np.zeros((nver,nver))

    for i in range(nver):
        for j in range(12):
            if mitmat[i,j] != i:
                adj[i,mitmat[i,j]]=adj[i,mitmat[i,j]] or 1 
    adj=np.logical_or(adj, adj.T)
    adj=np.triu(adj)
    adjacencymat.append(adj.astype(int))
    
        
   
    

pp = PdfPages('foo.pdf')

def show_graph_with_labels(adjacency_matrix,title):
    plt.close()
    rows, cols = np.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    plt.title(title)
    gr = nx.Graph()
    gr.add_edges_from(edges)
    graplo=nx.draw(gr, node_size=100, edgecolors='black',with_labels=False,node_color='gray')
    pp.savefig(graplo)
    
    
    

for i in range(len(adjacencymat)):
    show_graph_with_labels(adjacencymat[i],nomes_materiais[i])

pp.close()


exec(open('gromov.py').read())

n=len(adjacencymat)

lb=np.zeros((n,n))
ub=np.zeros((n,n))
for i in range(n):
        for j in range(n):
            lb[i,j], ub[i,j] = gromov_hausdorff(adjacencymat[i], adjacencymat[j])
        
kub=np.zeros((n,n))
klb=np.zeros((n,n)) # inicializamos

optimos=[]
xoptimos=[]

for k in range(n):

    for i in range(n):
        for j in range(n):
            kub[i,j]=lb[i,k]+lb[j,k]-lb[i,j]
            
    
    for i in range(n):
        for j in range(n):
            klb[i,j]=ub[i,k]+ub[j,k]-ub[i,j]
    
    import cvxpy as cp
    
    # Define and solve the CVXPY problem.
    # Create a symmetric matrix variable.
    X = cp.Variable((n,n), symmetric=True)
    
    
    # The operator >> denotes matrix inequality.
    constraints = [X >> 0]
    
    prob = cp.Problem(cp.Minimize(cp.atoms.norm(X-klb, p='fro', axis=None)+cp.atoms.norm(X-kub, p='fro', axis=None)),
                      constraints)
    prob.solve()
    
    optimos.append(prob.value)
    xoptimos.append(X.value)


min_value = min(optimos)
minimos=[i for i, x in enumerate(optimos) if x == min_value]
# Print result.
z0=random.choice(minimos)
kernel1=xoptimos[z0]

rho=np.zeros((n,n))
# corollary 16 annals 2013
for i in range(n):
    for j in range(n):
        rho[i,j]=kernel1[i,i]+kernel1[j,j]-2*kernel1[i,j]
        
# rho resulta que e negative type
# facemola strong negative type tomando raiz cadrada r (primeiro paragrafo 3286 lyons)
# enunciamos resultado de que o kernel q sae de aqui é characteristic 
# escollemos un z0 aleatorio e volvemos ao kernel
np.linalg.eigvals(kernel1)

rhostrong=np.sqrt(rho)
np.savetxt("rhostrong.txt", rhostrong, fmt="%s")

# eliximos un kernel da familia
z1=random.choice(range(n))

final=np.zeros((n,n))

for i in range(n):
        for j in range(n):
            final[i,j]=rhostrong[i,z1]+rhostrong[j,z1]-rhostrong[i,j]

np.savetxt("gromov.txt", final, fmt="%s")
np.savetxt("nomesmateriais.txt", nomes_materiais, fmt="%s")


