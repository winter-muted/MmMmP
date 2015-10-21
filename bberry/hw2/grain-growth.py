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

# Do some global defines and initialization


kappa = 1.0   # controls interface width
dx = 1.0
dy = 1.0
dt = .1
L = 1.0

nx = 128
ny = 128
nxy = nx*ny
nsteps = 1000
plot_invterval = 100
P = 10 # the number of order parameters to use
sumop = 0.0 # the sum of the order parameters (at each node point?)
maxop = 0.0 # the maximum order parameter (at each node?)

# for i in range(0,P):
phi = 0.1 + 0.2*(np.random.rand(P,nx,ny) - 0.5)

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
    for i in range(1,P):
        sumj = np.zeros(1,nx,ny)
        for j in range(1,P):
            if (j != i):
                sumj += phi[j][:][:]**2
            else:
                pass
            df[i][:][:] = .2*12*(phi[i][:][:]**3 - phi[i][:][:]**2 +phi[i][:][:]*sumj)


    for i in range(1,P):
        df_FFT[i][:][:] = fft2(squeeze(df[i][:][:]))/nxy
        phi_FFT = fft2(squeeze(phi[i][:][:]))/nxy

    for i in range(1,P):
        phi_FFT[i][:][:] = (phi_FFT[i][:][:] - dt*L*df_FFT[i][:][:] / (1.0 + dt*k2*L*kappa)
        # df_FFT[i][:][:] = (df_FFT[i][:][:] - dt*L*df_FFT[i][:][:]) / (1.0 + dt*k2*L*kappa)
        phi[i][:][:] = ifft2(squeeze(phi_FFT[i][:][:]))*nxy





def AC_run():
    for step in range(nsteps):
        update_AC()
        # plot on the given interval
        if (step == 0) or (step % plot_interval == 0):
            print phi
            my_plot()



# Plotting routines
def my_plot():
    pass

def visualize_op():

    sumop = np.zeros(nx,ny)

    for j in range(ny):
        for i in range(nx):
            sum = 0.0
            for k in range(P):
                sum += phi(k,i,j)**2
            sumop[i][j] = sum



# Boilerplate
def main():
    AC_run()

if __name__ == "__main__":
    main()
