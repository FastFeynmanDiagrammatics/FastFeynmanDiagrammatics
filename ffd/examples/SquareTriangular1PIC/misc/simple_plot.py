#!/usr/bin/env python3
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from mpl_toolkits import mplot3d

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

#import matplotlib.pyplot as plt
#import numpy as np
import os, sys
from math import sqrt

U = 8.
Omega = 0
Component = 'i' # i, r, A, 
order_minus = 1

if len(sys.argv) < 2:
    print("!!! You need more arguments !!!")
    exit()

print(sys.argv[1])
os.chdir(sys.argv[1])

if not os.path.exists("parameters"):
    print("!!!pameters file not found!!!!")
    exit()
else:
    with open("parameters/mu0") as ifile:
        mu0 = float(ifile.read())
    with open("parameters/Beta") as ifile:
        Beta = float(ifile.read())
    with open("parameters/geometry") as ifile:
        geometry = ifile.read()

def dispersion_relation(Kx, Ky):
    if geometry == "s":
        return -2*(cos(Kx)+cos(Ky)) - mu0
    elif geometry == "t":
        return -mu0
        
print("mu0= ", mu0)

den = 0.
err_den = 0.
sigma = {}
for fname in os.listdir():
    if fname.endswith(".series"):
        stem = fname.replace(".series", "")
        print(fname)
        o, val, err = np.loadtxt(fname, delimiter=' ', unpack=True)
        res = 0.
        stat_err = 0.
        last_order = o[len(val)-order_minus-1]
        trunc_err = (abs(val[len(val)-order_minus-1])+abs(err[len(val)-order_minus-1]))*pow(U, last_order)
        for j in range(len(val)-order_minus):
            res += val[j]*pow(U, o[j])
            stat_err += pow(abs(err[j]*pow(U, o[j])), 2)
        if order_minus > 0:
            j = len(val)-order_minus
#            res += .5*val[j]*pow(U, o[j])
            stat_err += pow(.5*abs(err[j]*pow(U, o[j])), 2)
        tot_err = sqrt(stat_err+pow(trunc_err, 2))
        print("result_error = ", res, tot_err)
        if stem == "density":
            den = res
            err_den = tot_err
        elif stem == "alpha":
            pass
        else:
            if not os.path.exists(stem+"/parameters"):
                print("parameters not found")
                exit()
            else:
                with open(stem+"/parameters/kx") as ifile:
                    kx = float(ifile.read())
                with open(stem+"/parameters/ky") as ifile:
                    ky = float(ifile.read())
                with open(stem+"/parameters/omega") as ifile:
                    omega = int(ifile.read())
                with open(stem+"/parameters/component") as ifile:
                    component = ifile.read()
                sigma[(kx, ky, omega, component)] = (res, tot_err)



                # plt.plot(x,y)
kx_ = np.array([])
ky_ = np.array([])
data = np.array([])
print("density =", den, err_den)
for (kx, ky, omega, component) in sigma:
    # print(kx, ky, omega, component, sigma[(kx, ky, omega, component)])
    if omega == Omega and component == Component:
        kx_ = np.append(kx_, kx)
        ky_ = np.append(ky_, ky)
        (val, err) = sigma[(kx, ky, omega, component)]
        data = np.append(data, -val)
        print(kx, ky, val, "+/-",err)

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter3D(kx_, ky_, sigma_, c=sigma_);
# plt.show()

# fig = plt.figure()
# ax = fig.gca(projection='3d')

# # Make data.
# # X = np.arange(-5, 5, 0.25)
# # Y = np.arange(-5, 5, 0.25)
# # X, Y = np.meshgrid(X, Y)
# # R = np.sqrt(X**2 + Y**2)
# # Z = np.sin(R)

# # Plot the surface.
# surf = ax.plot_surface(kx_, ky_, sigma_, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)

# # Customize the z axis.
# ax.set_zlim(-1.01, 1.01)
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))


# # Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)

#plt.show()


plt.scatter(kx_, ky_, c=data, cmap="viridis", s=100, marker="H")
if geometry == "t":
    plt.xlim(-0.1, 4.3)
    plt.ylim(-4.3, 0.1)
plt.colorbar()
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('Interesting Graph\nCheck it out')
# plt.legend()
plt.show()
