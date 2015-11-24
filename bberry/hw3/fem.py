# Do some global imports
from __future__ import division
import numpy as np

import matplotlib.pyplot as plt
import os
from multiprocessing import Process
import copy

elems = 16
nodes = elems + 1
h = 5.0/elems
w = .025
print "h=",h

alpha = 2.0/h
beta = -1.0/h
gamma = w*h/2
A = np.full(nodes*nodes,0.0).reshape((nodes,nodes))
b = np.zeros(nodes)

# Fix row 1
A[0][0] = alpha/2
# A[0][1] = beta
# b[0] = gamma
# general tridiagonal sys
for i in range(1,nodes-1):
    A[i][i-1] = beta
    A[i][i] = alpha
    A[i][i+1] = beta
    b[i] = -2*gamma

# fix end row
# A[nodes-1][nodes-2] = beta
A[nodes-1][nodes-1] = alpha/2
# b[nodes-1] = gamma



print "A=\n",A
print "b=\n",b

x = np.matmul(np.linalg.inv(A),b)

print "x=\n",x
print np.min(x)


plt.plot(x)
plt.show()
