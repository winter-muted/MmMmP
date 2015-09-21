# some global imports
import numpy as np
from random import random as rand
import matplotlib.pyplot as plt
import matplotlib.cm as cm

J = -0.4
HH = 0.5
kT = 1.5
nx = 50
ny = 50
nxy = nx*ny
xA = 0.3

nA = xA*nxy

nAo = nA

ncycles = 500
filename = 'Alloy-A3-kT15-J4neg.png'

spin = np.zeros((nx,ny))

for i in xrange(nx):
    for j in xrange(ny):
        if (rand() < xA):
             spin[i][j] = 1
        else:
            spin[i][j] = -1


def ising():


    for cycle in xrange(ncycles):
        for i in xrange(nxy):
            Attempt_Switch()
        # print cycle

    magn = np.sum(spin)/float(nxy)

    plt.imshow(spin,cmap=cm.summer)
    plt.colorbar()
    plt.clim(-1,1)
    plt.title('Alloyed system, J=-0.4,XA=0.3,kT=1.5')
    plt.savefig(filename)
    plt.clf()


def Attempt_Switch():
    global nA
    # generate a random int between "1" and nxy
    # map the integer into a 2-d index
    # aka pick a random lattice site
    rand1 = np.random.randint(0,nx)
    rand2 = np.random.randint(0,ny)


    #  get energy of random site before and after swap
    E1 = Lattice_site_energy(rand1,rand2)

    spin[rand1][rand2] = -1*spin[rand1][rand2]
    nA += spin[rand1][rand2]

    E2 = Lattice_site_energy(rand1,rand2)

    # Metropolis algo to decide acceptance

    dE = E2 - E1 # final - Initial

    if (dE > 0):
        if (np.random.rand() > np.exp(-dE/kT)):
            # reject the change
            spin[rand1][rand2] = -1*spin[rand1][rand2]
            nA += spin[rand1][rand2]


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

    # Eo = -HH*So
    Eo = HH*pow((nA - nAo),2)
    Eij = Ee + Ew + En + Es + Eo

    return Eij



def main():
    ising()


if __name__ == "__main__":
    main()
