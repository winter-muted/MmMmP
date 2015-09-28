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
	# import numpy as n

	# declare global variables
	global nxy

	# declare array to store the total magnetization after each cycle
	magn = zeros(cycles)

	# run simulation for the given number of monte carlo cycles (one cycle =
	# nx*ny randome monte carlo switch attempts)
	for cyc in range(cycles):
		cycle_name = cyc + 1
		print 'cycle %i completed' % cycle_name
		for i in range(nxy):
			attempt_switch()

	magn = sum(sum(spin))/nxy

	return magn;

def attempt_switch():
	'''
	This function attempts to make a switch for a randomely selected spin using
	the metropolis algorithm.
	'''
	# import numpy as n
	# import random as r

	# declare global variables
	global nx,ny,spin

	# select random lattice site (get the indices of the cite in the spin matrix)
	i = randint(0,nx-1)
	j = randint(0,ny-1)

	# calculate the energy of the lattice site
	e1 = lattice_site_energy(i,j)

	# make a trial switch to the randomely selected spin
	spin[i,j] = -spin[i,j] # flip spin
	# calculate energy of switched spin
	e2 = lattice_site_energy(i,j)

	# use the metropolis algorithm to determine whether to accept switch or not
	# calculate change in energy
	de = e2 - e1

	if de > 0:
		r2 = random()
		if r2 > exp(-de/kT): # reject switch
			spin[i,j] = -spin[i,j]

	return 0;

def lattice_site_energy(i,j):
	"""
	function calculates energy of site spin(i,j)
	"""

	# declare global variables
	global spin,JJ,HH,nx,ny

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
	Eo = -HH*So

	# sum up neighbor interaction and field energies to get total site energy
	eij = Ee + Ew + En + Es + Eo

	return eij;

##############################################################################
# global variables
##############################################################################

nx = 50					# number of lattice nodes in x
ny = 50					# number of lattice nodes in y
nxy = nx*ny				# total number of lattice sites
JJ = 0.4				# spin energy (for the ising model
HH = 0.5				# external magnetic field (ising model)
kT = 0.1				# boltzmann's constant times temperature
seed()					# seed the random number generator

# for ising model initialize spin lattice to random distributions
spin = zeros((nx,ny))
for j in range(ny):
	for i in range(nx):
		ran = random()
		if ran < 0.4:
			spin[i,j] = 1
		else:
			spin[i,j] = -1

###############################################################################
# main program
###############################################################################

cycles = 500

magn = ising(cycles)

# plot final state after monte carlo cycles
pp.imshow(spin,cmap = cm.winter)
pp.colorbar()
title = 'HH = '+str(HH)+', kT = '+str(kT)+', JJ = '+str(JJ)
pp.title(title)
pp.clim(-1.0,1.0)
pp.show()
###############################################################################
# end program
###############################################################################
