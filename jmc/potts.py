# ising model for two state magnetic spin system
from numpy import *
from random import *
import matplotlib.pyplot as pp
import matplotlib.cm as cm

def ising(cycles):
	"""
	function that runs an ising model simulation for a binary state magnetic
	spin system.

	INPUTS:

	cycles = the number of monte carlo cycles to run the simulation for

	OUTPUTS:

	magn = a 1-D numpy array with length = cycles that stores the total
			magnetization of the two-state spin system at the end of each
			monte-carlo cycle
	"""

	# declare global variables
	global nxy

	# run simulation for the given number of monte carlo cycles (one cycle =
	# nx*ny randome monte carlo switch attempts)
	for cyc in range(cycles):
		cycle_name = cyc + 1
		print 'cycle %i completed' % cycle_name
		for i in range(nxy):
			attempt_switch()

	return 0;

def attempt_switch():
	'''
	This function attempts to make a switch for a randomely selected spin using
	the metropolis algorithm.
	'''

	# declare global variables
	global nx,ny,spin,kT,Q

	# select random lattice site (get the indices of the cite in the spin matrix)
	i = randint(0,nx-1)
	j = randint(0,ny-1)

	# calculate the energy of the lattice site
	e1 = lattice_site_energy(i,j)

	# make a trial switch to the randomely selected spin
	spin_old = spin[i,j]
	q_rand = randint(1,Q) 
	spin[i,j] = q_rand
	
	# calculate energy of switched spin
	e2 = lattice_site_energy(i,j)

	# use the metropolis algorithm to determine whether to accept switch or not
	# calculate change in energy
	de = e2 - e1

	if de > 0:
		r2 = random()
		if r2 > exp(-de/kT): # reject switch
			spin[i,j] = spin_old

	return 0;

def lattice_site_energy(i,j):
	"""
	function calculates energy of site spin(i,j)
	"""

	# declare global variables
	global spin,JJ,nx,ny,krondelt

	# get neighbor spin values
	ie = i + 1		# east
	iw = i - 1		# west
	jn = j + 1		# north
	js = j - 1		# south

	# implement periodic boundary conditions
	if ie > nx-1:
		ie = 0
	if iw < 0:
		iw = nx -1
	if jn > ny -1:
		jn = 0
	if js < 0:
		js = ny -1

	# assign neighbor spins
	So = spin[i,j]
	Se = spin[ie,j]
	Sw = spin[iw,j]
	Sn = spin[i,jn]
	Ss = spin[i,js]
	
	Snw = spin[iw,jn]
	Sne = spin[ie,jn]
	Ssw = spin[iw,js]
	Sse = spin[ie,js]

	# calculate neighbor interaction energies
	Ee = -0.5*JJ*(krondelt[So-1,Se-1] - 1)
	Ew = -0.5*JJ*(krondelt[So-1,Sw-1] - 1)
	En = -0.5*JJ*(krondelt[So-1,Sn-1] - 1)
	Es = -0.5*JJ*(krondelt[So-1,Ss-1] - 1)
	
	Enw = -0.5*JJ*(krondelt[So-1,Snw-1] - 1)
	Ene = -0.5*JJ*(krondelt[So-1,Sne-1] - 1)
	Esw = -0.5*JJ*(krondelt[So-1,Ssw-1] - 1)
	Ese = -0.5*JJ*(krondelt[So-1,Sse-1] - 1)

	# sum up neighbor interaction and field energies to get total site energy
	eij = Ee + Ew + En + Es + Enw + Ene + Esw + Ese

	return eij;

##############################################################################
# global variables
##############################################################################

nx = 50					# number of lattice nodes in x
ny = 50					# number of lattice nodes in y
nxy = nx*ny				# total number of lattice sites
JJ = 0.4				# spin energy (for the ising model

kT = 0.1				# boltzmann's constant times temperature
seed()					# seed the random number generator
Q = 10

krondelt = identity(Q)

# for ising model initialize spin lattice to random distributions
spin = zeros((nx,ny))
for j in range(ny):
	for i in range(nx):
		spin[i,j] = randint(1,Q)

###############################################################################
# main program
###############################################################################

cycles = 50

ising(cycles)

# plot final state after monte carlo cycles
pp.imshow(spin,cmap = cm.terrain)
pp.colorbar()
pp.caxis = ([0,Q-1])
pp.show()
###############################################################################
# end program
###############################################################################
