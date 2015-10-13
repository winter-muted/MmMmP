# A plot tester
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm
# from mpl_toolkits.mplot3d import Axes3D
#
#
#
# nx = 10
# for nx in range(nx):
#     kx1 = (mod(.5 + (nx-1)/nx,1) - .5)
#     print kx1
#
# x = np.arange(-2, 2, 0.01)
# y = np.arange(-2, 2, 0.01)
# xx,yy = np.meshgrid(x,y,sparse=True)
# z = xx*np.exp((-xx**2)-(yy**2))
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(xx,yy,z,cmap=cm.coolwarm)
# plt.show()
"""
An animated image
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure()

def f(x, y):
    return np.sin(x) + np.cos(y)

x = np.linspace(0, 2 * np.pi, 120)
y = np.linspace(0, 2 * np.pi, 100).reshape(-1, 1)

im = plt.imshow(f(x, y), cmap=plt.get_cmap('jet'))

def updatefig(*args):
    global x,y
    x += np.pi / 15.
    y += np.pi / 20.
    im.set_array(f(x,y))
    return im,

ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)
plt.show()
