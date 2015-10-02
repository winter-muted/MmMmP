from numpy import *
import matplotlib.pyplot as plt

#################################################
#				global variables			    #
#################################################

# scalar parameters
K = 1.0					#
dx = 1.0				#
dy = 1.0				#
dt = 1.0				#
M = 1.0					#
nx = 8				# number of points in x
ny = 8				# number of points in y
nxy = nx*ny				# total number of points
nstr = 0.15				# noise strength
nsteps = 1			#
skip = 10				#

# arrays
c = zeros((nx,ny))		# concentration array

kx1 = mod( 1/2 + arange(float(nx))/nx,1) -1/2
kx = kx1*(2*pi/dx)

print kx1

ky1 = mod( 1/2 + arange(float(ny))/ny,1) -1/2
ky = ky1*(2*pi/dy)

KX,KY = meshgrid(kx,ky)

k2 = KX*KX + KY*KY
k4 = k2*k2

#################################################
#			initializations of c				#
#################################################

# homogenous mixture with noise
noise = random.uniform(-nstr,nstr,(nx,ny))
c = c + noise

# plot inititial condition
#plt.imshow(c)
#plt.colorbar()
#plt.clim(-1,1)
#plt.show()

# step function

#################################################
#				run Cahn-Hilliard			    #
#################################################

for step in range(nsteps):

	# calculate df in real space
	df = c*c*c - c

	# get the fft of df and c
	df_fft = fft.fft2(df)/nxy
	c_fft = fft.fft2(c)/nxy

	# update the CH solution
	c_fft = (c_fft - dt*M*k2*df_fft)/(1.0 + dt*K*M*k4)

	# transform c_fft back to real space to get c again
	c = real(fft.ifft2(c_fft))*nxy

	# plot c

	#if step % skip == 0:


# plot final step
plt.imshow(c)
plt.colorbar()
plt.clim(-1,1)
plt.show()
