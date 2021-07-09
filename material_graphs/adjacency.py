#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 19:07:55 2021

@author: carlos
"""

import argparse
import os
import shutil
import sys
import time
import warnings
from random import sample

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn import metrics
from torch.autograd import Variable
from torch.optim.lr_scheduler import MultiStepLR

from cgcnn.data import CIFData
from cgcnn.data import collate_pool, get_train_val_test_loader
from cgcnn.model import CrystalGraphConvNet

dataset = CIFData()

n=len(dataset)

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
    
        
    #http://www.crystallography.net/cod/search.html
    #http://hoffman.physics.harvard.edu/materials/CuprateIntro.php
    
   # https://materialsproject.org/materials/mp-5986/#
    
import matplotlib.pyplot as plt
import networkx as nx

from matplotlib.backends.backend_pdf import PdfPages

pp = PdfPages('foo.pdf')

def show_graph_with_labels(adjacency_matrix):
    plt.close()
    rows, cols = np.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    plt.title('A title')
    gr = nx.Graph()
    gr.add_edges_from(edges)
    graplo=nx.draw(gr, node_size=100, edgecolors='black',with_labels=False,node_color='gray')
    pp.savefig(graplo)
    
    
    

for i in range(len(adjacencymat)):
    show_graph_with_labels(adjacencymat[i])

pp.close()

lb=np.zeros((n,n))
ub=np.zeros((n,n))
for i in range(n):
        for j in range(n):
            lb[i,j], ub[i,j] = gromov_hausdorff(adjacencymat[i], adjacencymat[j])
        
disgh=(ub+lb)/2

ghsim=(disgh+np.transpose(disgh))/2

np.fill_diagonal(ghsim,0)

kernopt=1-ghsim/np.amax(ghsim)

kub=1-lb/np.amax(lb)

klb=1-ub/np.amax(ub)

import cvxpy as cp
import numpy as np

# Generate a random SDP.
n = 3
p = 3
np.random.seed(1)
C = np.random.randn(n, n)
A = []
b = []
for i in range(p):
    A.append(np.random.randn(n, n))
    b.append(np.random.randn())

# Define and solve the CVXPY problem.
# Create a symmetric matrix variable.
X = cp.Variable((n,n), symmetric=True)

#np.linalg.norm(X-A, 'fro')
# The operator >> denotes matrix inequality.
constraints = [X >> 0]
# constraints += [
#     X <= kub
# ]
# constraints += [
#     klb <= X
# ]
prob = cp.Problem(cp.Minimize(cp.atoms.norm(X-klb, p='fro', axis=None)+cp.atoms.norm(X-kub, p='fro', axis=None)),
                  constraints)
prob.solve()

# Print result.
print("The optimal value is", prob.value)
print("A solution X is")
print(X.value)

final=X.value


def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

np.linalg.eigvals(final)

np.savetxt("gromov.txt", final, fmt="%s")



