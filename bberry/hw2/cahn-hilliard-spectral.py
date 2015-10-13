# Do some global imports
import numpy as np
import matplotlib as plt

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
plot_invterval = 100

c = 0.0 + 0.01*(np.random.rand(ny,nx) - 0.5)
# we can also initialize with a step function
# c = np.ones((nx,1))
# if i < nx/2 for i in range(nx):
#     c[i] = -1.0

# here is some matlab code to fix indicies for fourier space
kx1 = mod( .5 + (0:(nx-1))/nx,1) - .5
kx = kx1*(2*pi/dx)

ky1 = mod( .5 + (0:(ny-1))/ny,1) - .5
ky = ky1*(2*pi/dy)

[KX,KY] = meshgrid(kx,ky)

k2 = KX.*KX + KY.*KY
k4 = k2**2








# Update CH function
def update_CH():
    # first calculate df/dc
    df = c.*c.*c.  - c

    # forwartd transfrom 'c' and 'df'

    df_FFT = fft2(df) / float(nxy)
    c_FFT = fft2(dfc) / float(nxy)

    # Update 'c' in fourier space, then invert

    c_FFT = (c_FFT - dt*k2.*df_FFT) ./ (1.0 + dt*k4)

    c = ifft2(c_FFT)*float(nxy)



# Run function
def CH_run():
    for step in range(nsteps):
        update_CH()
        # plot on the given interval
        if (step == 0) or (step % plot_interval == 0):
            print c
            my_plot()



# Plotting routines
def my_plot():
    pass

# Boilerplate
def main():
    CH_run()

if __name__ == "__main__":
    main()
