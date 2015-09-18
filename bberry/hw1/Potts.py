# some global imports
import numpy as np
from random import random as rand
import matplotlib.pyplot as plt

J = 0.4
HH = 0.0
kT = 1.2
nx = 50
ny = 50
nxy = nx*ny

Q = 10
krondelt

spin = np.zeros((nx,ny))

for i in xrange(nx):
    for j in xrange(ny):
        r = rand()
        spin[i][j] = np.random.randint(0,Q-1) + 1

def potts():
    ncycles = 100

    for cycle in xrange(ncycles):
        for i in xrange(nxy):
            Attempt_Switch()
        # print cycle

    magn = np.sum(spin)/float(nxy)
    print magn
    plt.plot(spin)
    plt.savefig("Initial.png")
    plt.clf()

def Attempt_Switch():
    # generate a random int between "1" and nxy
    # map the integer into a 2-d index
    # aka pick a random lattice site
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
            # reject the change
            spin[rand1][rand2] = -1*spin[rand1][rand2]


def Lattice_site_energy(x_pos,y_pos):

    ie = x_pos + 1
    iw = x_pos - 1
    jn = y_pos + 1
    js = y_pos - 1

    if ie > nx-1:
        ie = 1
    if iw < 0:
        iw = nx -1
    if jn > ny-1:
        jn = 1
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
    potts()


if __name__ == "__main__":
    main()
