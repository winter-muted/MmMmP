# A grain growth routine in two dimensions
# uses spectral solving method
# Do some global imports
from __future__ import division
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
from multiprocessing import Process
import copy
import sys

# Do some global defines and initialization


kappa = 1.0   # controls interface width
dx = 1.0
dy = 1.0
dt = .1
L = 1.0

nx = 256
ny = 256
nxy = nx*ny
nsteps = 1000
plot_interval = 100
P = 10 # the number of order parameters to use
sumop = 0.0 # the sum of the order parameters (at each node point?)
maxop = 0.0 # the maximum order parameter (at each node?)

# for i in range(0,P):
# phi = 0.1 + 0.2*(np.random.rand(P,nx,ny) - 0.5)

phi = np.full(P*nxy,0.5)
for i in range(P*nxy):
    a = np.random.rand()
    if a < .5:
        phi[i] -= .1*a
    else:
        phi[i] += .1*a
phi = phi.reshape((P,nx,ny))

df = np.zeros((P,nx,ny))
df_FFT = np.zeros((P,nx,ny))
phi_FFT = np.zeros((P,nx,ny))

# This code generates fourier-type meshing and indexing
# values go from 0 to pi and -pi to zero
kx = (((.5 + np.arange(nx)/nx) % 1) - .5)*(2*np.pi/dx)
ky = (((.5 + np.arange(ny)/ny) % 1) - .5)*(2*np.pi/dy)
# kx1 = mod( .5 + (0:(nx-1))/nx,1) - .5
# kx = kx1*(2*pi/dx)
# ky1 = mod( .5 + (0:(ny-1))/ny,1) - .5
# ky = ky1*(2*pi/dy)

KX, KY = np.meshgrid(kx,ky,sparse=True)
k2 = KX*KX + KY*KY
# k4 = k2**
# k2(1,:,:) = KX.*KX + KY.*KY

def AC_update():
    # calculate df/dphi
    global phi
    global df
    global df_FFT
    global phi_FFT
    for i in range(1,P):
        sumj = np.zeros((1,nx,ny))
        for j in range(1,P):
            if (j != i):
                sumj += phi[j][:][:]**2
            else:
                pass
        df[i][:][:] = .2*12*(phi[i][:][:]**3 - phi[i][:][:]**2 +phi[i][:][:]*sumj)


    for i in range(1,P):
        df_FFT[i][:][:] = np.fft.fft2(np.squeeze(df[i][:][:]))/nxy
        phi_FFT[i][:][:] = np.fft.fft2(np.squeeze(phi[i][:][:]))/nxy

    for i in range(1,P):
        # phi_FFT.shape
        # df_FFT.shape
        # sys.exit(0)
        phi_FFT[i][:][:] = (phi_FFT[i] - dt*L*df_FFT[i]) / (1.0 + dt*k2*L*kappa)
        # df_FFT[i][:][:] = (df_FFT[i][:][:] - dt*L*df_FFT[i][:][:]) / (1.0 + dt*k2*L*kappa)
        phi[i][:][:] = np.fft.ifft2(np.squeeze(phi_FFT[i][:][:]))*nxy





def AC_run():
    for step in range(nsteps):
        AC_update()
        print step
        # plot on the given interval
        if (step == 0) or (step % plot_interval == 0):

            average_grain_size(step)
            visualize_op(step)



# Plotting routines
def average_grain_size():
    # find the index of the maximum order parameter for each node point

    #



def visualize_op(step):

    sumop = np.zeros((nx,ny))

    for i in range(nx):
        for j in range(ny):
            op_sum = 0.0
            for k in range(P):
                op_sum += phi[k][i][j]**2
            sumop[i][j] = op_sum

    filename = 'op-step' + str(step) + '.png'
    plt.imshow(sumop)
    plt.colorbar()
    plt.clim(0,1)
    plt.savefig(filename)
    plt.clf()

# Boilerplate
def main():
    AC_run()

if __name__ == "__main__":
    main()
