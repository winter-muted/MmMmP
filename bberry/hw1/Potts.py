# some global imports
import numpy as np
from random import random as rand
import matplotlib.pyplot as plt
import matplotlib.cm as cm

J = 0.4
HH = 0.0
kT = .8

nx = 50
ny = 50
nxy = nx*ny

ncycles = 500
Q = 10
filename = 'potts.png'

spin = np.zeros((nx,ny))

for i in xrange(nx):
    for j in xrange(ny):
        r = rand()
        spin[i][j] = np.random.randint(0,Q)

plt.imshow(spin,cmap=cm.Set3)
plt.colorbar()
plt.clim(0,10)
plt.title('Typical initial condition for 50/50 starting concentration')
plt.savefig('Potts.png')
plt.clf()



def potts():
    for cycle in xrange(ncycles):
        for i in xrange(nxy):
            Attempt_Switch()
        # print cycle

    magn = np.sum(spin)/float(nxy)

    plt.imshow(spin,cmap=cm.Set3)
    plt.colorbar()
    plt.clim(0,10)
    plt.title('Ending state for kT = 0.1')
    plt.savefig('Potts_final.png')
    plt.clf()


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



def Attempt_Switch():
    # generate a random int between "1" and nxy
    # map the integer into a 2-d index
    # aka pick a random lattice site
    rand1 = np.random.randint(0,nx)
    rand2 = np.random.randint(0,ny)


    #  get energy of random site before and after swap
    E1 = Lattice_site_energy(rand1,rand2)

    # Select a random second nearest neighbor
    rand3 = np.random.randint(rand1-2,rand1+2)
    rand4 = np.random.randint(rand2-2,rand2+2)

    rand3 = periodic(rand3)
    rand4 = periodic(rand4)

    E2 = Lattice_site_energy(rand3,rand4)



    dE = E2 - E1 # final - Initial

    if (dE > 0):
        if (np.random.rand() > np.exp(-dE/kT)):
            # reject the change
            spin[rand1][rand2] = -1*spin[rand1][rand2]


def Lattice_site_energy(x_pos,y_pos):
    index_list = ['ie','iee','iw','iww','jn','jnn','js','jss']
    ie = x_pos + 1
    iee = x_pos + 2
    iw = x_pos - 1
    iww = x_pos - 2
    jn = y_pos + 1
    jnn = y_pos + 2
    js = y_pos - 1
    jss = y_pos - 1

    ie = periodic(ie)
    iee = periodic(iee)
    iw = periodic(iw)
    iww = periodic(iww)
    jn = periodic(jn)
    jnn = periodic(jnn)
    js = periodic(js)
    jss = periodic(jss)


    So = spin[x_pos][y_pos]
    Se = spin[ie][y_pos] + spin[iee][y_pos]
    Sw = spin[iw][y_pos] + spin[iww][y_pos]
    Sn = spin[x_pos][jn] + spin[x_pos][jnn]
    Ss = spin[x_pos][js] + spin[x_pos][jss]

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
