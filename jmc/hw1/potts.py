# Q-State Potts Model Code
from numpy import *
from random import *
import matplotlib.pyplot as pp
import matplotlib.cm as cm

def ising(cycles):
	"""
	main driver function for the Potts model
	"""

	# declare global variables
	global nxy

	# figure counter
	fcount = 1

	# path to directory that I will save figure files in
	path = '/home/joseph/Dropbox/graduate_school/fall_2015/multi-scale_materials_modeling/hw1/'

	# run simulation for the given number of monte carlo cycles (one cycle =
	# nx*ny randome monte carlo switch attempts)
	for cyc in range(cycles):
		cycle_name = cyc + 1
		print 'cycle %i completed' % cycle_name
		for i in range(nxy):
			attempt_switch()

		if cyc % 25 == 0:
			pp.imshow(spin,cmap = cm.terrain)

			if fcount == 1:
				pp.colorbar()
				pp.clim = (1,Q)

			pp.title('Cycle '+str(cyc))
			filename='hw1Cfig'+str(fcount)+'.eps'
			pp.savefig(path+filename)
			fcount = fcount + 1

	return 0;

def attempt_switch():
	'''
	This function attempts to make a switch for a randomely selected neighbor
	grain using the metropolis algorithm.
	'''

	# declare global variables
	global nx,ny,spin,kT,Q

	# select random lattice site (get the indices of the cite in the spin matrix)
	i = randint(0,nx-1)
	j = randint(0,ny-1)

	# calculate the energy of the lattice site
	e1 = lattice_site_energy(i,j)

	# make a random trial switch to one of the neighbor values
	spin_old = spin[i,j]			# store current state

	ii,jj = get_rand_neigh(i,j)	# randomely select a neighbor
	spin_new = spin[ii,jj]

	# make the trail switch
	spin[i,j] = spin_new

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
	function calculates energy of site i,j
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

	# first nearest neighbors
	So = spin[i,j]
	Se = spin[ie,j]
	Sw = spin[iw,j]
	Sn = spin[i,jn]
	Ss = spin[i,js]

	# second nearest neighbors
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

def get_rand_neigh(i,j):
	'''
	This function accepts the indices of a spin and returns a random neighbor
	from the eight different options.
	'''
	rn = randint(1,8)
	if rn == 1:				# N neighbor
		ii = i
		jj = j + 1
	elif rn == 2:			# NE neighbor
		ii = i + 1
		jj = j + 1
	elif rn == 3:			# E neighbor
		ii = i + 1
		jj = j
	elif rn == 4:			# SE neighbor
		ii = i + 1
		jj = j - 1
	elif rn == 5:			# S neighbor
		ii = i
		jj = j - 1
	elif rn == 6:			# SW neighbor
		ii = i - 1
		jj = j - 1
	elif rn == 7:			# W neighbor
		ii = i - 1
		jj = j
	else:					# NW neighbor
		ii = i - 1
		jj = j + 1

	# implement periodic boundary conditions
	if ii > nx-1:
		ii = 0
	if ii < 0:
		ii = nx -1
	if jj > ny -1:
		jj = 0
	if jj < 0:
		jj = ny -1

	return ii,jj;

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

krondelt = identity(Q)	# the kronecker delta function in matrix form

# for ising model initialize spin lattice to random distributions
spin = zeros((nx,ny))
for j in range(ny):
	for i in range(nx):
		spin[i,j] = randint(1,Q)

###############################################################################
# main program
###############################################################################

# figure counter
fcount = 1

# path to directory that I will save figure files in
path = '/home/joseph/Dropbox/graduate_school/fall_2015/multi-scale_materials_modeling/hw1/'
cycles = 75

ising(cycles)

# plot final state after monte carlo cycles
pp.imshow(spin,cmap = cm.terrain)
pp.title('Cycle '+str(cycles))
pp.savefig(path+'hw1Cfig4.eps')

###############################################################################
# end program
###############################################################################
