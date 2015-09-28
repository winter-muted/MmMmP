# ising model for binary alloy
from numpy import *
from random import *
import matplotlib.pyplot as pp
import matplotlib.cm as cm

def ising(cycles):
	"""
	main driver function for the ising model
	"""

	# declare global variables
	global nxy,spin,nx,ny,HH,JJ,nA

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
	global nx,ny,spin,nA,nAo

	# select random lattice site (get the indices of the cite in the spin matrix)
	i = randint(0,nx-1)
	j = randint(0,ny-1)

	# calculate the energy of the lattice site
	e1 = lattice_site_energy(i,j)

	# make a trial switch to the randomely selected spin
	spin[i,j] = -spin[i,j] # flip spin
	nA = nA + spin[i,j]

	# calculate energy of switched spin
	e2 = lattice_site_energy(i,j)

	# use the metropolis algorithm to determine whether to accept switch or not
	# calculate change in energy
	de = e2 - e1

	if de > 0:
		r2 = random()
		if r2 > exp(-de/kT): # reject switch
			spin[i,j] = -spin[i,j]
			nA = nA + spin[i,j]

	return 0;

def lattice_site_energy(i,j):
	"""
	function calculates energy of site spin(i,j)
	"""

	# declare global variables
	global spin,JJ,HH,nx,ny,nA,nAo

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

	# calculate neighbor interaction energies
	Ee = -0.5*JJ*So*Se
	Ew = -0.5*JJ*So*Sw
	En = -0.5*JJ*So*Sn
	Es = -0.5*JJ*So*Ss

	# External field energy
	Eo = HH*(nA - nAo)**2.0

	# sum up neighbor interaction and field energies to get total site energy
	eij = Ee + Ew + En + Es + Eo

	return eij;

##############################################################################
# global variables
##############################################################################

nx = 50					# number of lattice nodes in x
ny = 50					# number of lattice nodes in y
nxy = nx*ny				# total number of lattice sites
JJ = -0.4				# spin energy (for the ising model
HH = 0.05				# conserve xA parameter
kT = 1.5				# boltzmann's constant times temperature
xA = 0.3				# fraction of A atoms
nA = int(round(xA*nxy))		# number of A atoms in system
nAo = nA				# number of A atoms at the outset
seed()					# seed the random number generator

# for ising model initialize spin lattice to random distributions
spin = zeros((nx,ny))
for j in range(ny):
	for i in range(nx):
		ran = random()
		if ran < xA:
			spin[i,j] = 1
		else:
			spin[i,j] = -1

###############################################################################
# main program
###############################################################################

cycles = 500

ising(cycles)

# plot final state after monte carlo cycles
pp.imshow(spin,cmap = cm.summer)
pp.colorbar()
title = 'xA = '+str(xA)+', kT = '+str(kT)+', JJ = '+str(JJ)
pp.title(title)
pp.clim(-1.0,1.0)
path = '/home/joseph/Dropbox/graduate_school/fall_2015/multi-scale_materials_modeling/hw1/'
filename='hw1Bfig7.eps'
pp.savefig(path+filename)
###############################################################################
# end program
###############################################################################
