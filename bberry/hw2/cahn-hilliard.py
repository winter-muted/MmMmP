# Do some global imports
import numpy as np
import matplotlib as plt

# Do some global defines and initialization
# Time step criteria dt is about 1/(dx^4)

K = 1.0   # controls interface width
dx = 1.0
dt = 0.01
M = 1.0

nx = 51
nsteps = 1000
plot_invterval = 100

# c = 0.0 + 0.01*(np.random.rand(0,nx) - 0.5)
# we can also initialize with a step function
c = np.ones((nx,1))
if i < nx/2 for i in range(nx):
    c[i] = -1.0


mu = np.zeros((nx,1))
d2c = np.zeros((nx,1))
d2mu = np.zeros((nx,1))


# Update CH function
def update_CH():
    # Calculate the chemical potential

    # This is matlab syntax that needs to be pythonified
    #  multipling mu by a prefactor raises the energy barrier
    # and sharpens the interface
    # width of iterface:
    # l is proportional to sqrt(K/df_max)5
    d2c = (np.roll(c,-1) + np.roll(c,1) - 2.0*c)/(dx**2)
    mu = (c**3 - c) - 2.0*K*d2c

    # Update the concentration
    d2mu = (np.roll(mu,-1) + np.roll(mu,1) - 2.0*mu)/(dx**2)
    c = c + dt*M*d2mu



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


# Boilerplate
def main():
    CH_run()

if __name__ == "__main__":
    main()
