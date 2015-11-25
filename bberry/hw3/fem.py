# Do some library imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


elems = 16
nodes = elems + 1
h = 5.0/elems
w = -.025

alpha = 2.0/h
beta = -1.0/h
gamma = w*h/2
A = np.full(nodes*nodes,0.0).reshape((nodes,nodes))
b = np.zeros(nodes)

# Fix first row
A[0][0] = alpha/2
# general tridiagonal sys
for i in range(1,nodes-1):
    A[i][i-1] = beta
    A[i][i] = alpha
    A[i][i+1] = beta
    b[i] = 2*gamma
# fix end row
A[nodes-1][nodes-1] = alpha/2
# solve
x = np.matmul(np.linalg.inv(A),b)

t1 = np.linspace(0,5,nodes)
t2 = np.linspace(0,5,250)  # change to nodes to get proper error evaluation
y = (w/2)*(-t2*t2 + 5*t2)
# residual error and minima
r = np.zeros(nodes)
# r = y - x
# print "error = ",np.sum(r)
# print np.min(x)
# print np.min(y)


# plot

plt.figure(figsize=(8,4.5))
plt.plot(t1,x,marker='o',linestyle='--',color='red',label='FEM')
plt.plot(t2,y,color='blue',label='Analytic')
plt.title('4 elements')
plt.xlabel('Position (ft)')
plt.ylabel('Displacement (ft)')
plt.legend()

# plt.show()
plt.savefig('16-elems.png')
