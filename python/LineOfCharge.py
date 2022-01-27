# LineOfCharge
# Assignment 2, Exercises 3-5
# Calculates the electrostatic potential and electric field due to short line 
# of charge in a 2D plane and produces vector, surface, and contour plots of
# the data
"""
Created on Mon Jan 24 15:33:58 2022

@author: masenpitts
"""

import numpy as np
import matplotlib.pyplot as plt

# Function that calculates the electric field and electrostatic potential
# at every point given a 2D array of values representing charges at points in
# space. 
# Parameters:
#   charge - A 2D array of values that represent the charge values and every
#           point in the region under consideration
#   d - The physical distance between mesh points (the function assumes that
#       mesh is rectangular and that the points are evenly spaces).
#   eX - The array used to store the values of the x-components of the electric
#       field at every point in the region.
#   eY - The array used to store the values of the y-components of the electric
#       field at every point in the region. 
#   V - The array used to store the values of the electrostatic potential at
#       every point in the region
def coulomb(charge, d, eX, eY, V):
    xM = charge.shape[0]
    yM = charge.shape[1]
    k = 8.9875517923e9
    
    for cX in range(xM):
        for cY in range(yM):
            if charge[cX,cY] != 0:
                kq = k*charge[cX,cY]
                for pX in range(xM):
                    dx = (pX-cX)*d
                    for pY in range (yM):
                        if charge[pX,pY] == 0:
                            dy = (pY-cY)*d
                            r = 1/np.sqrt(dx*dx + dy*dy)
                            r3 = r*r*r
                            kqR3 = kq*r3
                            eX[pX,pY] += kqR3*dx
                            eY[pX,pY] += kqR3*dy
                            V[pX,pY] += kq*r
                        


# Gives the number of mesh points extending in each direction
meshRes = 100

# Variables stored to initialize arrays and perform calculations
xMesh = meshRes + 1
yMesh = xMesh

# Physical boundaries of the region
xlim = 5
ylim = xlim

# Stores the physical distance between mesh points
d = xlim/meshRes

# Tuple used to represent the shape of the arrays
S = (xMesh, yMesh)

# Create lists representing the physical x and y positions of the 
# cross-section points
x = np.linspace(0, xlim, xMesh)
y = np.linspace(0, ylim, yMesh)

# Array storing the value of the charge at each location in Coulombs
charge = np.zeros(S)

# Arrays storing electric field component data and potential value data
# for each point in the 2D space
electricFieldX = np.full(S, 0, dtype=np.float64)
electricFieldY = np.full(S, 0, dtype=np.float64)
potential = np.full(S, 0, dtype=np.float64)

length = 3 # Length of the line of charge in meters
q = 1 # Total charge of the line in Coulombs

# Variables used to set the locations of the point charges
xCenter = int(xlim/(2*d))
yCenter = int(ylim/(2*d))
L = int(length/(2*d))

# Sets nonzero charge values on a line of points along the x-midline of
# the region
for X in range(xCenter-L, xCenter+L+1):
    charge[X, yCenter] = q/(2*L)

# Calculates the electric field and electrostatic potential at each mesh point
# in the region based on the given charge
coulomb(charge, d, electricFieldX, electricFieldY, potential)

Emag = np.sqrt(electricFieldX*electricFieldX + electricFieldY*electricFieldY)    

# Plots a 2D Vector Plot of the electric field in the region
plt.streamplot(x, y, electricFieldX.T, electricFieldY.T, color=Emag, cmap=plt.cm.plasma)
plt.title("Electric Field Surrounding Line")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.colorbar(label="E Field Magnitude (N/C)")
plt.show()

# Plots a 2D Color filled Contour plot of the potential in the region
plt.contourf(x, y, potential.T)
plt.title("Electrostatic Potential Surrounding Line")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.colorbar(label="Electric Potential (V)")
plt.show()

# Plots a 3D Surface Plot of the potential in the region
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
X, Y = np.meshgrid(x, y)
plt.title("Electrostatic Potential Surrounding Line")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
ax.plot_surface(X, Y, potential.T, cmap=plt.cm.plasma)
ax.set_zlabel("(V)")
plt.show()

