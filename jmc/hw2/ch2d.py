from numpy import *
import matplotlib.pyplot as plt
import os

#################################################
#	global variables			#
#################################################

# scalar parameters
K = 1.			# kappa interface parameter
dx = 1.0		# spacing in x-dir
dy = 1.0		# spacing in y-dir
dt = 1.0		# time step
M = 1.0			# mobility
nx = 128                # number of points in x
ny = 128                # number of points in y
nxy = nx*ny             # total number of points
nstr = 0.15             # noise strength
nsteps = 10000          # number of steps to run sim
skip = 500		# number of steps to skip between plots

# arrays
c = zeros((nx,ny))	# concentration array

kx1 = mod( 1.0/2 + arange(float(nx))/nx,1.0) -1.0/2
kx = kx1*(2.0*pi/dx)

ky1 = mod( 1.0/2 + arange(float(ny))/ny,1) -1.0/2
ky = ky1*(2.0*pi/dy)

KX,KY = meshgrid(kx,ky)

k2 = KX*KX + KY*KY
k4 = k2*k2

#################################################
#	initializations of c			#
#################################################

# homogenous mixture with noise
noise = random.uniform(-nstr,nstr,(nx,ny))
c = c + noise

#################################################
#		file io stuff   		#
#################################################

# name of directory to store output vtk files in
dirname = 'animationfiles'

# delete old vtk files from previous runs
if os.path.exists('./' + dirname):
    os.system('rm -r '+ dirname)

# create a directory for the animation files
script_dir = os.path.dirname(os.path.abspath(__file__))
dest_dir = os.path.join(script_dir, dirname)
try:
    os.makedirs(dest_dir)
except OSError:
    pass # already exists

#plot inititial condition
plt.imshow(c)
plt.title('step = 1')
plt.colorbar()
plt.clim(-1,1)
plt.savefig(dest_dir + '/ch2d_step_1.png')

# step function

#################################################
#           run Cahn-Hilliard                   #
#################################################

for step in range(nsteps):

    # calculate df in real space
    df = c*c*c - c

    # get the fft of df 
    df_fft = fft.fft2(df)/nxy
    c_fft = fft.fft2(c)/nxy

    # update the CH solution
    c_fft = (c_fft - dt*M*k2*df_fft)/(1.0 + dt*K*M*k4)

    # transform c_fft back to real space to get c again
    c = real(fft.ifft2(c_fft))*nxy

    # plot c

    if step % skip == 0 and step > 1:
        plt.imshow(c)
        filename='ch2d_step_'+str(step)+'.png'
        plt.savefig(dest_dir +'/'+ filename)

    print 'step = %i '%step


# plot final step
plt.plot(c)
filename='ch2d_step_'+str(step+1)+'.png'
plt.savefig(dest_dir +'/'+ filename)

# make movie from .png files
print('Making movie animation.mpg - this make take a while')
os.system("mencoder 'mf://ch2d_step_*.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.mpg")
