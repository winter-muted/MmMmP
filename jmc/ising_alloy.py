# ising model for homework
from numpy import *
from random import *

# global variables
nx = 50					# number of lattice nodes in x
ny = 50					# number of lattice nodes in y
nxy = nx*ny				# total number of lattice sites
JJ = 0.4				# spin energy (for the ising model
HH = 0.05				# Energy Penalty constant for composition conservation
kT = 0.1				# boltzmann's constant times temperature

xA = 0.3
nA = int(xA*nxy)
nAo = nA

# for ising model initialize spin lattice to random distributions
seed()
spin = zeros((nx,ny))
for j in range(ny):
	for i in range(nx):
		ran = random()
		if ran < xA:
			spin[i,j] = 1
		else:
			spin[i,j] = -1

def ising(cycles):
		
	# run simulation
	
	for cyc in range[cycles]:
		for i in range[nxy]:
			attempt_switch()
			
	magn = sum(sum(spin))
	
	# plot
	
	
	return 0;

def

def attempt_switch(nA):
	
	# select random lattice site (get the indices of the cite in the spin matrix)
	
	# calculate the energy of the lattice site
	
	e1 = lattice_site_energy()
	
	# make a trial switch to the randomely selected spin
	spin[i,j] = -spin[i,j] # flip spin
	nA = nA + spin[i,j]
	
	# calculate energy of switched spin
	e2 = lattice_site_energy()
	
	# use the metropolis algorithm to determine whether to accept switch or not
	# calculate change in energy
	de = e2 - e1
	
	if de > 0:
		r2 = random()
		if r2 > exp(-de/kT): # reject switch
			spin[i,j] = -spin[i,j]
			nA = nA + spin[i,j]
	
	return 0;

def lattice_site_energy(i,j,nA,nAo):
	"""
	function calculates energy of site spin(i,j)
	"""
	
	# get neighbor spin values
	ie = i + 1		# east
	iw = i - 1		# west
	jn = j + 1		# north
	js = j - 1		# south
	
	# implement periodic boundary conditions
	if ie > nx-1:
		ei = 0
	if iw < 0:
		iw = nx -1
	if jn > ny:
		jn = ny -1
	if js < 0:
		js = 0
	
	# assign neighbor spins
	So = spin[i,j]
	Se = spin[ie,j]
	Sw = spin[iw,j]
	Sn = spin[i,jn]
	Ss = spin[i,js]
	
	# Sum neighbor interaction energies
	
	Ee = -0.5*JJ*So*Se
	Ew = -0.5*JJ*So*Sw
	En = -0.5*JJ*So*Sn
	Es = -0.5*JJ*So*Ss
	
	# External field energy
	
	#Eo = -HH*So
	Eo = HH*(nA -nAo)**2
	
	eij = Ee + Ew + En + Es + Eo
	
	return eij;
