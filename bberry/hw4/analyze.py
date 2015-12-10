from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def fit_func(x,A,B):
    return A*x + B



filename = 'sigma33-small.txt'
strain_rate = .0001538*2 # .02 / 130 angstroms, engineering strain
energy_factor = 160.2 # eV/ang^3 to GPA
# strain_rate = 0
sim_steps = 30000
strain_interval = 20
plot_interval = 5

# elastic_range = 1350 # visually determined end of elastic deformation
elastic_range = 1050

data = np.genfromtxt(filename,delimiter=',',names=['x','y'])

strain_x = data['x']*(strain_rate*plot_interval/strain_interval)

fit_x = np.linspace(0,strain_x[elastic_range])
fit_data = data['y'][0:elastic_range]*energy_factor


params = curve_fit(fit_func,strain_x[0:elastic_range],fit_data)
[E,C] = params[0]
# #
fit_y = fit_x*E + C
print ("Young's Modulus = %f" % E)

plt.title('8x8x36 Z stress')
plt.xlabel('Strain')
plt.ylabel('sigma-zz (Pa)')
line1 = plt.plot(strain_x,data['y']*energy_factor)
line2 = plt.plot(strain_x[0:elastic_range],fit_data,fit_x,fit_y)
plt.legend(['Data','E=%f' %E])
plt.show()
plt.savefig('sigma33-small-fit.png')
