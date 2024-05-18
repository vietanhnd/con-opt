"""
Script that creates a 2D contour plot of the Himmelblau function

along with the constraints

f(x)=(x1^2 +x2 -11)^2 +(x1 +x2^2 -7)^2

s.t:

(x1+2)^2 -x2 >= 0
-4x1 +10x2 >= 0

-5<=x1<=5
-5<=x2<=5
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Define the Himmelblau function
def himmelblau(x1, x2):
    return (x1**2 + x2 - 11)**2 + (x1 + x2**2 - 7)**2

# Define the constraints
def constraint1(x1, x2):
    return (x1 + 2)**2 - x2

def constraint2(x1, x2):
    return -4*x1 + 10*x2

# Create the meshgrid
x1 = np.linspace(-6, 6, 100)
x2 = np.linspace(-6, 6, 100)

X1, X2 = np.meshgrid(x1, x2)
Y = himmelblau(X1, X2)

# Create the contour plot

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# 2D contour plot
ax.contour(X1, X2, Y, 50, cmap=cm.coolwarm)
ax.set_xlabel(r'$x_1$')
ax.set_ylabel(r'$x_2$')
ax.set_zlabel(r'$f(x)$')
ax.set_title('Himmelblau function')

plt.show()

# Add the constraints
