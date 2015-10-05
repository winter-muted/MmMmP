from numpy import *
import matplotlib.pyplot as plt

#################################################
#		global variables		#
#################################################

# scalar parameters
kappa = 1.0			# 
dx = 1.0			# 
dy = 1.0			#
dt = 1.0			# 
L = 1.0				# 
nx = 256			# number of points in x
ny = 256			# number of points in y
nxy = nx*ny			# total number of points
nstr = 0.15			# noise strength
nsteps = 500			# 
skip = 10			# 
bias = 0.0			# free energy bias par

# arrays
phi = zeros((nx,ny))		# concentration array

kx1 = mod( 1.0/2 + arange(float(nx))/nx,1.0) -1.0/2
kx = kx1*(2.0*pi/dx)

ky1 = mod( 1.0/2 + arange(float(ny))/ny,1.0) -1.0/2
ky = ky1*(2.0*pi/dy)

KX,KY = meshgrid(kx,ky)

k2 = KX*KX + KY*KY

#################################################
#	initializations of phi			#
#################################################

# homogenous mixture with noise
noise = random.uniform(-nstr,nstr,(nx,ny))
phi = phi + noise

#################################################
#	Allen-Cahn Sim				#
#################################################

for step in range(nsteps):
	
	# calculate df in real space
	df = phi*phi*phi - phi + bias*(1 - 2*phi**2 + phi**4) 
	
	# get the fft of df and phi
	df_fft = fft.fft2(df)/nxy
	phi_fft = fft.fft2(phi)/nxy
	
	# update the CH solution
	phi_fft = (phi_fft - dt*L*k2*df_fft)/(1.0 + dt*kappa*L*k2)
	
	# transform phi_fft back to real space to get c again
	phi = real(fft.ifft2(phi_fft))*nxy
	
	print 'step %i completed.'%step 
	# plot phi
	
	#if step % skip == 0:
		

# plot final step
plt.imshow(phi)
plt.colorbar()
plt.clim(-1,1)
plt.show()
