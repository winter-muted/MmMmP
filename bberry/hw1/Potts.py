# some global imports
import numpy as np
from random import random as rand
import matplotlib.pyplot as plt
import matplotlib.cm as cm

########## Parameters ###########
output_list = [100,200,300] # what time steps to output
ncycles = 500
Q = 10

nx = 50
ny = 50

J = 0.4
kT = 1.5
########## Initialization ###########
nxy = nx*ny
spin = np.zeros((nx,ny))

for i in xrange(nx):
    for j in xrange(ny):
        r = rand()
        spin[i][j] = np.random.randint(0,Q)

#########################################

def Attempt_Switch():
    # Pick a random site on the domain
    rand1 = np.random.randint(0,nx)
    rand2 = np.random.randint(0,ny)

    #  get energy of random site before swap
    E1 = Lattice_site_energy(rand1,rand2)


    # Select a random second nearest neighbor
    rand3 = periodic(np.random.randint(rand1-2,rand1+3))
    rand4 = periodic(np.random.randint(rand2-2,rand2+3))

    # Store the old spin and change the spin to that of the neighbor
    old_spin = spin[rand1][rand2]
    spin[rand1][rand2] = spin[rand3][rand4]

    E2 = Lattice_site_energy(rand1,rand2)
    dE = E2 - E1 # final - Initial

    # Decide if we want to keep the change, metropolis style
    if (dE > 0):
        if (np.random.rand() > np.exp(-dE/kT)):
            # reject the change
            spin[rand1][rand2] = old_spin

# Determine energy of node
def Lattice_site_energy(x_pos,y_pos):

    ie = periodic(x_pos + 1)
    iee = periodic(x_pos + 2)
    iw = periodic(x_pos -1)
    iww = periodic(x_pos - 2)
    jn = periodic(y_pos + 1)
    jnn = periodic(y_pos + 2)
    js = periodic(y_pos - 1)
    jss = periodic(y_pos - 2)

    Eij = 0.0

# Neighbors only contribute energy if the spin is different
    if (spin[x_pos][y_pos] != spin[ie][y_pos]):
        Eij += .5*J
    if (spin[x_pos][y_pos] != spin[iee][y_pos]):
        Eij += .5*J
    if (spin[x_pos][y_pos] != spin[iw][y_pos]):
        Eij += .5*J
    if (spin[x_pos][y_pos] != spin[iww][y_pos]):
        Eij += .5*J
    if (spin[x_pos][y_pos] != spin[x_pos][jn]):
        Eij += .5*J
    if (spin[x_pos][y_pos] != spin[x_pos][jnn]):
        Eij += .5*J
    if (spin[x_pos][y_pos] != spin[x_pos][js]):
        Eij += .5*J
    if (spin[x_pos][y_pos] != spin[x_pos][jss]):
        Eij += .5*J

    return Eij


# Fix up indicies to be periodic
def periodic(index):
    if index == nx:
        return 0
    if index == (nx + 1):
        return 1
    if index == -1:
        return (nx -1)
    if index == -2:
        return (nx -2)
    return index

# Program flow control
def potts():
    for cycle in xrange(ncycles):
        for i in xrange(nxy):
            Attempt_Switch()
        if int(cycle) in output_list:
            my_plot(cycle)


# Custom Plotting function
def my_plot(cycle):
    title = 'Time Step ' + str(cycle)
    filename = 'Potts' + str(cycle) + '_high'
    plt.imshow(spin,cmap=cm.Paired)
    plt.colorbar()
    plt.clim(0,10)
    plt.title(title)
    plt.savefig(filename)
    plt.clf()

# Boilerplate
def main():
    my_plot(0)
    potts()


if __name__ == "__main__":
    main()
