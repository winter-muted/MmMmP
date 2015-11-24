# Do some global imports
from __future__ import division
import numpy as np

import matplotlib.pyplot as plt
import os
from multiprocessing import Process
import copy

elems = 4
nodes = elems + 1

alpha = 1.0
beta = 1.0

A = np.full(nodes*nodes,0.0).reshape((nodes,nodes))
b = np.zeros(nodes)

# Fix row 1
A[0][0] = -2*beta
A[0][1] = alpha
# general tridiagonal sys
for i in range(1,nodes-1):
    A[i][i-1] = alpha
    A[i][i] = -2*beta
    A[i][i+1] = alpha
    b[i] = 1

# fix end row
A[nodes-1][nodes-2] = alpha
A[nodes-1][nodes-1] = -2*beta

x = np.matmul(np.linalg.inv(A),b)
# print A
# print b
print x

plt.plot(x)
plt.show()
