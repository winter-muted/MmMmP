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
nx = 256				# number of points in x
ny = 256				# number of points in y
nxy = nx*ny				# total number of points
nstr = 0.15				# noise strength
nsteps = 500			# 
skip = 50				# 

# arrays
c = zeros(nx,ny)		# concentration array

kx1 = mod( 1/2 + (0:(nx-1))/nx,1) -1/2
kx = =kx1*(2*pi/dx)

ky1 = mod( 1/2 + (0:(ny-1))/ny,1) -1/2
ky = ky1*(2*pi/dy)

KX,KY = meshgrid(kx,ky)

k2 = multiply(KX,KX) + multiply(KY,KY)
k4 = multiply(k2,k2)

#################################################
#			initializations of c				#
#################################################

# homogenous mixture with noise
noise = random.uniform(-nstr,nstr,(nx,ny))
c = c + noise

# step function

#################################################
#				run Cahn-Hilliard			    #
#################################################

for step in range(nsteps):
	
	# update the CH solution
	df = c*c*c - c
	df_fft = fft2(df)/nxy
	c_fft = fft2(c)/nxy
	
	c_fft = (c_fft - dt*M*k2*df_fft)/(1.0 + dt*K*M*k4) # figure out how to do element wise division
	c = ifft2(c_fft)*nxy
	
	# plot c
	
	#if step == 1 or step % skip == 0:
		
		
# plot final step
plt.imshow(c)
plt.colorbar()
plt.clim(-1,1)
plt.show()
