from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def fit_func(x,A,B):
    return A*x + B




filename = 'sigma33-small.txt'
strain_rate = .0001538 # .02 / 130 angstroms
# strain_rate = 0
sim_steps = 30000
strain_interval = 20
plot_interval = 5

elastic_range = 1000 # visually determined end of elastic deformation

data = np.genfromtxt(filename,delimiter=',',names=['x','y'])

strain_x = np.linspace(0,strain_rate*plot_interval/(strain_interval),sim_steps/plot_interval + 1)

fit_x = np.linspace(0,strain_rate*plot_interval/(strain_interval),elastic_range)
fit_data = data['y'][0:elastic_range]



params = curve_fit(fit_func,fit_x,fit_data)
[E,C] = params[0]
# #
fit_y = fit_x*E + C
#
# print("Slope = %f" % E)
# print(params[0])
print(fit_x.size)
print(fit_data.size)
print(strain_x.size)
print(data['y'].size)
# for i in range(0,20):
#     print(fit_x[i],strain_x[i])


plt.xlabel('Elongation Factor')
plt.ylabel('sigma-zz (Pa/angstrom)')
plt.plot(fit_x,fit_data)
plt.plot(strain_x,data['y'])
# # plt.legend()
plt.show()
# # plt.savefig('sigma33-small.png')
