# Do some global imports
import numpy as np
import matplotlib as plt
from __future__ import division
# Do some global defines and initialization


kappa = 1.0   # controls interface width
dx = 1.0
dy = 1.0
dt = 1.0
L = 1.0

nx = 256
ny = 256
nxy = nx*ny
nsteps = 500
plot_invterval = 100


phi = 1.0 + 0.00*(np.random.rand(nx,ny) - 0.5)
phi(1:5,:) = -1.0
c = 0.5 + 0.4*(np.random.rand(nx,ny) - 0.5)

kx1 = mod( .5 + (0:(nx-1))/nx,1) - .5
kx = kx1*(2*pi/dx)

ky1 = mod( .5 + (0:(ny-1))/ny,1) - .5
ky = ky1*(2*pi/dy)

[KX,KY] = meshgrid(kx,ky)

k2 = KX.*KX + KY.*KY
k4 = k2**2

def AC_update():
    # calculate df/dphi
    df = phi**3 - phi + 0.01 * (1 - 2*phi**2 + phi**4) - c**2*(1 - c**2) + 0.5 * (0.5 - c)**2

    df_FFT = fft2(df)/nxy
    phi_FFT = fft2(phi)/nxy

    df_FFT = (df_FFT - dt*L*df_FFT) / (1.0 + dt*k2*L*kappa)
    phi = ifft2(phi_FFT)*nxy

def CH_update():
    A2 = 0.5*(phi + 1)
    dfdc = 2 * (1 - A2) * (4*c**3 - 6 * c**2 + 2*c) - 2.0 * A2 * (0.5 - c)

    dfdc_FFT = fft2(dfdc)/nxy
    c_FFT = fft2(c)/nxy





def eutectic_run():
    for step in range(nsteps):
        AC_update()
        CH_update()
        # plot on the given interval
        if (step == 0) or (step % plot_interval == 0):
            print phi
            my_plot()



# Plotting routines
def my_plot():
    pass

# Boilerplate
def main():
    eutectic_run()

if __name__ == "__main__":
    main()
