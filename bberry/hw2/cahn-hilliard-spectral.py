# Do some global imports
from __future__ import division
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import thread
import copy


# Do some global defines and initialization
# Time step criteria dt is about (dx^4)
# This code produces spinodal Decompositions

K = 1.0   # controls interface width
dx = 1.0
dy = 1.0
dt = 1.0
M = 1.0

nx = 256
ny = 256
nxy = nx*ny
nsteps = 1500
plot_interval = 10

c = np.full(nxy,0.5)
for i in range(nxy):
    a = np.random.rand()
    if a < .5:
        c[i] -= 0.001*a
    else:
        c[i] += 0.001*a
c = c.reshape((nx,ny))


# c = 0.5 + .001*np.random.randn(nx,ny)
c0 = np.average(c)

L = []
total_c = []

total_c.append(sum(sum(c)))




# This code generates fourier-type meshing and indexing
# values go from 0 to pi and -pi to zero
kx = (((.5 + np.arange(nx)/nx) % 1) - .5)*(2*np.pi/dx)
ky = (((.5 + np.arange(ny)/ny) % 1) - .5)*(2*np.pi/dy)

KX, KY = np.meshgrid(kx,ky,sparse=True)

k2 = KX*KX + KY*KY
k4 = k2**2


# Update CH function
def update_CH():
    # first calculate df/dc
    global c
    df = c**3  - c

    # forwartd transfrom 'c' and 'df'

    df_FFT = np.fft.fft2(df) / nxy
    c_FFT = np.fft.fft2(c) / nxy


    # Update 'c' in fourier space, then invert

    c_FFT = (c_FFT - dt*k2*df_FFT) / (1.0 + dt*k4)
    c = np.real(np.fft.ifft2(c_FFT))*nxy


def domain_size():
        c0_FFT = np.fft.fft2(c - c0) / nxy
        Sk = abs(c0_FFT) * abs(c0_FFT)
        kx2 = sum(sum(KX*KX*Sk))/sum(sum(Sk))
        ky2 = sum(sum(KY*KY*Sk))/sum(sum(Sk))
        Lx = 2*(np.pi)/kx2**.5
        Ly = 2*(np.pi)/ky2**.5
        L.append((Lx + Ly) / 2)

# Run function
def CH_run():
    for step in range(nsteps):
        update_CH()

        print "Finished step:" + str(step)
        # plot on the given interval
        if (step == 0) or (step % plot_interval == 0):
            total_c.append(sum(sum(c)))

            try:
                c_plot = copy.copy(c)
                # c_plot = c
                thread.start_new_thread ( my_plot, (c_plot,step) )
                # my_plot(step)
            except:
                print "Could not create plot thread."

            domain_size()

# Plotting routines
def my_plot(c_plot,step):
    filename = 'chs-step' + str(step) + '.png'
    plt.imshow(c_plot)
    plt.colorbar()
    plt.clim(-1,1)
    plt.savefig(filename)
    plt.clf()

# Boilerplate
def main():
    CH_run()
    # Make an animation and cleanup the results
    # Currently has a problem. There should be a more pythonic way of doing
    # the animation.
    os.system('rm chs-animation.gif')
    os.system('convert -delay 100 -loop 0 chs-step* chs-animation.gif')
    # os.system('rm chs-step*')

    # L(t) plot code
    def fit_func(x,A,n):
        return A*(x**n)


    xdata = np.arange(0,nsteps,plot_interval)
    y = fit_func(xdata,1,1)
    popt, pcov = curve_fit(fit_func,xdata,L)
    x = np.arange(0,nsteps)
    y = fit_func(x,popt[0],popt[1])
    plt.plot(xdata,L,'ro',x,y,'b')
    plt.title('L(t)')
    plt.text(10,max(y)-2,'A=' + str(popt[0]))
    plt.text(10,max(y)-4,'n=' + str(popt[1]))
    plt.savefig("L-trend.png")
    plt.clf()

    print total_c
    plt.plot(total_c)
    plt.savefig("c_evolution.png")
    plt.clf()

if __name__ == "__main__":
    main()
