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
nsteps = 5000
plot_interval = 100

P = 10 # the number of order parameters to use
sumop = 0.0 # the sum of the order parameters (at each node point?)

phi = np.full(P*nxy,0.1)
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
KX, KY = np.meshgrid(kx,ky,sparse=True)

k2 = KX*KX + KY*KY

domain_size = []


def AC_update():
    global phi
    global df
    global df_FFT
    global phi_FFT

    # calculate df/dphi
    for i in range(1,P):
        sumj = np.zeros((nx,ny))
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
        phi_FFT[i] = (phi_FFT[i] - dt*L*df_FFT[i]) / (1.0 + dt*k2*L*kappa)
        phi[i][:][:] = np.fft.ifft2(np.squeeze(phi_FFT[i][:][:]))*nxy

# Plotting/Aux Calc routines
def average_grain_size(sumop):

    global domain_size

    sumop_FFT = np.fft.fft2(sumop)/nxy
    Sk = abs(sumop_FFT) * abs(sumop_FFT)
    kx2 = sum(sum(KX*KX*Sk))/sum(sum(Sk))
    ky2 = sum(sum(KY*KY*Sk))/sum(sum(Sk))
    Lx = 2*(np.pi)/kx2**.5
    Ly = 2*(np.pi)/ky2**.5
    domain_size.append(((Lx + Ly) / 2)**2)


def visualize_op(sumop,step):

    filename = 'op-step' + str(step) + '.png'
    plt.imshow(sumop)
    plt.colorbar()
    plt.clim(0,1)
    plt.savefig(filename)
    plt.clf()


def AC_run():
    for step in range(nsteps):
        AC_update()
        # print "Finished step:" + str(step)
        # plot on the given interval
        if (step == 0) or (step % plot_interval == 0):

            sumop = np.zeros((nx,ny))
            for i in range(nx):
                for j in range(ny):
                    op_sum = 0.0
                    for k in range(P):
                        op_sum += phi[k][i][j]**2
                    sumop[i][j] = op_sum

            average_grain_size(sumop)
            # for some reason, the grain size computation doesnt work on another thread.
            # do it sequentially for now
            try:
                phi_compute = copy.copy(sumop)
                # thread1 = Process(target=average_grain_size,args=(phi_compute,))
                # thread1.start()
                thread2 = Process(target=visualize_op,args=(phi_compute,step))
                thread2.start()
            except:
                print "Could not create plot thread."


def main():
    AC_run()
    global domain_size
    xdata = np.arange(0,nsteps,plot_interval)

    # print xdata.size
    # print len(domain_size)
    # plt.plot(xdata,domain_size,'ro')
    # plt.title('Average Domain Size')
    # plt.savefig("size.png")
    # plt.clf()

    # plot evolution of domain size at end of simulation
    domain_size= np.log(domain_size)
    xdata = np.log(xdata)
    plt.plot(xdata,domain_size,'ro')
    plt.title('Average Domain Size')
    plt.savefig("sizeln.png")
    plt.clf()

# Boilerplate
if __name__ == "__main__":
    main()
