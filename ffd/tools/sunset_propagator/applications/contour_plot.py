#!/usr/bin/env python3
'''
======================
3D surface (color map)
======================

Demonstrates plotting a 3D surface colored with the coolwarm color map.
The surface is made opaque by using antialiased=False.

Also demonstrates using the LinearLocator and custom formatting for the
z axis tick labels.
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

# L = 129 #equal to half of the physical L + 1

L = 33
# L = 65

# data = np.loadtxt("fermisurface_omega0_bare.dat")
# data = np.loadtxt("realpart_omega0_bare.dat")
data = np.loadtxt("imagpart_omega0_BareThird.dat")
X = data[:, 0]
Y = data[:, 1]
Z = data[:, 2]
ZZ = Z
# ZZ = 1./(.01+Z*Z)


fig = plt.figure()
# ax = fig.gca(projection='3d')


# Make data.
# R = np.sqrt(X**2 + Y**2)

# Plot the surface.
surf=plt.pcolor(X.reshape(L,L), Y.reshape(L, L), ZZ.reshape(L, L))
# surf = ax.plot_surface(X.reshape(L, L), Y.reshape(L, L), Z.reshape(L, L), cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)



# surf = ax.pcolor(Z.reshape(L, L)
# Customize the z axis.
# ax.set_zlim(-.2, .2)
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
