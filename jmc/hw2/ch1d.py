from numpy import *
import matplotlib.pyplot as plt

#################################################
#				global variables			    #
#################################################

# scalar parameters
K = 1.0					# 
dx = 1.0				# 
dt = 0.01				# must be ~ 1/dx^4
M = 1.0					# 
nx = 51					# number of points
nstr = 0.15				# noise strength
nsteps = 1			# 
skip = 1000				# 

# arrays
c = zeros(nx)			# 
mu = zeros(nx)			# 
d2c = zeros(nx)			# 
d2mu = zeros(nx)		# 
x = arange(0,nx,dx)	# array for plotting

#################################################
#			initializations of c				#
#################################################

# homogenous mixture with noise
noise = random.uniform(-nstr,nstr,nx)
c = c + noise

# step function

#################################################
#				run Cahn-Hilliard			    #
#################################################

for step in range(nsteps):
	
	# update the CH solution
	dc2 = (roll(c,-1) + roll(c,1) - 2.0*c)/(dx*dx)
	mu = multiply(c,multiply(c,c)) - c - 2.0*K*d2c
	d2mu = (roll(mu,-1) + roll(mu,1) - 2.0*mu)/(dx*dx)
	c = c + dt*M*d2mu
	
	
	# plot c
	
	#if step == 1 or step % skip == 0:
		#plt.plot(x,c,'b-*')
		#plt.show()
		
plt.plot(x,c,'b-*')
plt.ylim((-1,1))
plt.show()
