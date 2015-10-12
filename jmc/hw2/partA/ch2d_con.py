from numpy import *
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

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
nsteps = 15000          # number of steps to run sim
skip = 500		# number of steps to skip between plots
co = 0.25		# average concenctration of species 'a'
op = nsteps/skip + 1    # number of animation outpus


# arrays
sim_time = arange(0.0,nsteps+skip,skip)*dt
c = ones((nx,ny))	# concentration array
l = zeros(nsteps/skip+1)# average domain size array

kx1 = mod( 1.0/2 + arange(float(nx))/nx,1.0) -1.0/2
kx = kx1*(2.0*pi/dx)

ky1 = mod( 1.0/2 + arange(float(ny))/ny,1) -1.0/2
ky = ky1*(2.0*pi/dy)

KX,KY = meshgrid(kx,ky)

k2 = KX*KX + KY*KY
k4 = k2*k2

#################################################
#		initializations of c		#
#################################################

# homogenous mixture with noise
c = c*co		# initializing c right at co
noise = random.uniform(-nstr,nstr,(nx,ny))
c = c + noise

#################################################
#		file io stuff   		#
#################################################

# name of directory to store output vtk files in
dirname = 'ch2d_con_animationfiles'

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

# plot setup
f1 = plt.figure(1)
plt.imshow(c)
plt.colorbar()
plt.clim(0,1)

# step function

#################################################
#           run Cahn-Hilliard                   #
#################################################

cnt = 0;                    # animation counter

for step in range(nsteps):

    # calculate df in real space
    df = 2.0*c - 6.0*c**2 + 4.0*c**3

    # get the fft of df and c - co
    df_fft = fft.fft2(df)/nxy
    c_fft = fft.fft2(c)/nxy

    # update the CH solution
    c_fft = (c_fft - dt*M*k2*df_fft)/(1.0 + dt*K*M*k4)

    # transform c_fft back to real space to get c again
    c = real(fft.ifft2(c_fft))*nxy

    # plot c and caclulate average domain size

    if (step+1) % skip == 0 or step == 0 :
        #print out cnt value
        print "step = %i, cnt = %i"%(step+1,cnt)

        #plot
        plt.imshow(c)
	plt.title('step = '+str(step+1))
        filename=str(step+1).zfill(10)+'.png'
        plt.savefig(dest_dir +'/'+ filename)

        # caclulate average domain size
        co_fft = fft.fft2(c - co)/nxy
        Sk = abs(co_fft)*abs(co_fft)
        kx2 = sum(sum(KX*KX*Sk))/sum(sum(Sk))
        ky2 = sum(sum(KY*KY*Sk))/sum(sum(Sk))
        lx = 2.0*pi/kx2**0.5
        ly = 2.0*pi/ky2**0.5
        lstep = (lx + ly)/2.0
        l[cnt] = lstep
        cnt = cnt + 1

# plot the average domain size, l(t), as a function of simulation time and fit curve to it.

# curve fitting
def fitfunct(t,A,n):
    return A*t**n;
fitPar,fitCov = curve_fit(fitfunct,sim_time,l)
t = linspace(0.0,nsteps*dt)
lfit = fitPar[0]*t**fitPar[1]
par_string = '$A =$ '+str(fitPar[0])+', $n =$ '+str(fitPar[1])+', $C_A = $'+str(co)

# plotting
f2 = plt.figure(2)
ax = f2.add_subplot(111)
plt.plot(sim_time,l,'ko',t,lfit,'b')
plt.title('average domain size data with curve of the\nform $A*t^n$ fit to the data')
plt.xlabel('time in reduced units')
plt.ylabel('$l(t)$')
ax.annotate(par_string,xy=(400.0,15.0))
plt.savefig('adsplot_ch2d_con.png')

# play an animation of the simulation
os.system("animate -delay 30 ch2d_con_animationfiles/*.png & proc=$!")
