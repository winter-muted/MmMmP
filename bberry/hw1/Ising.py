# some global imports
import numpy as np
from random import random as rand
import matplotlib.pyplot as plt
import matplotlib.cm as cm

J = 0.4
HH = 0.5
kT = .1
nx = 50
ny = 50
nxy = nx*ny

filename = 'J04-HH5-kT1.png'
initial_biasing = .5  # what percentage of particles get assigned spin down at start
ncycles = 500

spin = np.zeros((nx,ny)) #Initalize an nx by ny array to zero
magn = []



def ising():
    for i in xrange(nx):
        for j in xrange(ny):
            if (rand() < initial_biasing):
                spin[i][j] = -1
            else:
                spin[i][j] = 1

    magn.append(np.sum(spin)/float(nxy))

    for cycle in xrange(ncycles):
        magn.append(np.sum(spin)/float(nxy))
        for i in xrange(nxy):
            Attempt_Switch()


    plt.imshow(spin,cmap=cm.summer)
    plt.colorbar()
    plt.clim(-1,1)
    plt.title('J=0.4,HH=0.5,kT=0.1')
    plt.savefig(filename)
    plt.clf()

def Attempt_Switch():
    rand1 = np.random.randint(0,nx)
    rand2 = np.random.randint(0,ny)


    #  get energy of random site before and after swap
    E1 = Lattice_site_energy(rand1,rand2)

    spin[rand1][rand2] = -1*spin[rand1][rand2]

    E2 = Lattice_site_energy(rand1,rand2)

    # Metropolis algo to decide acceptance

    dE = E2 - E1 # final - Initial

    if (dE > 0):
        if (np.random.rand() > np.exp(-dE/kT)):
            spin[rand1][rand2] = -1*spin[rand1][rand2]


def Lattice_site_energy(x_pos,y_pos):

    ie = x_pos + 1
    iw = x_pos - 1
    jn = y_pos + 1
    js = y_pos - 1

    if ie > nx-1:
        ie = 0
    if iw < 0:
        iw = nx -1
    if jn > ny-1:
        jn = 0
    if js < 0:
        js = ny -1

    So = spin[x_pos][y_pos]
    Se = spin[ie][y_pos]
    Sw = spin[iw][y_pos]
    Sn = spin[x_pos][jn]
    Ss = spin[x_pos][js]

    Ee = -.5*J*So*Se
    Ew = -.5*J*So*Sw
    En = -.5*J*So*Sn
    Es = -.5*J*So*Ss

    # External Field Energy

    Eo = -HH*So
    Eij = Ee + Ew + En + Es + Eo

    return Eij



def main():
    ising()


if __name__ == "__main__":
    main()
