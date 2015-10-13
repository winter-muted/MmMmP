# Do some global imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


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
nsteps = 500
plot_interval = 100

c = 0.0 + 0.01*(np.random.rand(nx,ny) - 0.5)
# f1 = plt.figure(1)

# we can also initialize with a step function
# c = np.ones((nx,1))
# if i < nx/2 for i in range(nx):
#     c[i] = -1.0

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

    df_FFT = np.fft.fft2(df) / float(nxy)
    c_FFT = np.fft.fft2(c) / float(nxy)

    # Update 'c' in fourier space, then invert

    c_FFT = (c_FFT - dt*k2*df_FFT) / (1.0 + dt*k4)

    c = np.real(np.fft.ifft2(c_FFT))*nxy



# Run function
def CH_run():
    for step in range(nsteps):
        update_CH()
        # plot on the given interval
        if (step == 0) or (step % plot_interval == 0):
            run(data)
    ani = animation.FuncAnimation(fig, run, data_gen, blit=True, interval=10,
        repeat=False)

def run(data):


# Plotting routines
def my_plot():
    plt.imshow(c)
    plt.colorbar()
    plt.clim(-1,1)
    plt.show()

# Boilerplate
def main():
    CH_run()

if __name__ == "__main__":
    main()
