# LineOfCharge
# Assignment 2, Exercises 3-5
# Calculates the electrostatic potential and electric field due to short line 
# of charge in a 2D rectangular plane and produces vector, surface, and contour 
# plots of the data
"""
Created on Mon Jan 24 15:33:58 2022

@author: masenpitts
"""

import numpy as np
import matplotlib.pyplot as plt

# Function that calculates the electric field and electrostatic potential
# in SI units at every point given a 2D array of values representing charges at 
# points in space.
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
    xM = charge.shape[0] # Determines the shape of the "charge" array and uses
    yM = charge.shape[1] # it as the shape of the 2D region
    k = 8.9875517923e9 # Coulomb's Constant in N*m^2/C^2
    
    # The upper 2 for loops run through the "charge" array, checking for nonzero
    # values. If a nonzero charge value is detected the lower two for loops run 
    # through all elements of the eX, eY, and V arrays (which are assumed to have
    # the same shape as the "charge" array) and calculate the electric field and
    # electric potential at each point where there is no charge.
    for sourceX in range(xM):
        for sourceY in range(yM):
            # Check for nonzero charge value
            # Note: the variables kq, r, r3 and kqR3 are strategically placed
            # so as to avoid the unnecessary repetition of calculations and
            # optimize performance.
            if charge[sourceX,sourceY] != 0:
                kq = k*charge[sourceX,sourceY] 
                for pointX in range(xM):
                    # Physical x distance between source charge and test point
                    dx = (pointX-sourceX)*d
                    for pointY in range (yM):
                        # Check to make sure there is no charge at test point
                        if charge[pointX,pointY] == 0:
                            # Physical y distance between source charge and 
                            # test point
                            dy = (pointY-sourceY)*d
                            # Calculated for optimization
                            r = 1/np.sqrt(dx*dx + dy*dy)
                            r3 = r*r*r
                            kqR3 = kq*r3
                            # Store the calculated values for the electric field
                            # and electric potential at test point
                            eX[pointX,pointY] += kqR3*dx
                            eY[pointX,pointY] += kqR3*dy
                            V[pointX,pointY] += kq*r
                        


# Gives the number of mesh points extending in each direction
meshRes = 100

# Variables stored to initialize arrays and perform calculations
xMesh = meshRes + 1
yMesh = xMesh

# Physical boundaries of the region in meters
xMin = -2.5; xMax = 2.5
yMin = xMin; yMax = xMax

# Stores the total physical distance between the boundaries of the region
plotXRange = xMax - xMin; plotYRange = yMax - yMin
 
# Stores the physical distance between mesh points
d = plotXRange/meshRes

# Tuple used to represent the shape of the arrays
S = (xMesh, yMesh)

# Create lists representing the physical x and y positions of the 
# cross-section points

x = np.linspace(xMin, xMax, xMesh)
y = np.linspace(yMin, yMax, yMesh)

# Array storing the value of the charge at each location in Coulombs
charge = np.zeros(S)

# Arrays storing electric field component data and potential value data
# for each point in the 2D space
electricFieldX = np.full(S, 0, dtype=np.float64)
electricFieldY = np.full(S, 0, dtype=np.float64)
potential = np.full(S, 0, dtype=np.float64)

length = 3 # Length of the line of charge in meters
q = 1e-6 # Total charge of the line in Coulombs

# Variables used to set the locations of the point charges
xCenter = int(plotXRange/(2*d))
yCenter = int(plotYRange/(2*d))
L = int(length/(2*d))

# Sets nonzero charge values on a line of points along the x-midline of
# the region
for X in range(xCenter-L, xCenter+L+1):
    charge[X, yCenter] = q/(2*L)

# Calculates the electric field and electrostatic potential at each mesh point
# in the region based on the given charge
coulomb(charge, d, electricFieldX, electricFieldY, potential)

# Calculate magnitude of the electric field vectors to use on the streamplot colormap
Emag = np.sqrt(electricFieldX*electricFieldX + electricFieldY*electricFieldY)    

# Plots a 2D Vector Plot of the electric field in the region
vec = plt.figure()
plt.streamplot(x, y, electricFieldX.T, electricFieldY.T, color=Emag, cmap=plt.cm.magma, 
                density=1.3)
plt.title("Electric Field Around Line of Charge")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.colorbar(label="E Field Magnitude (N/C)")
plt.show()

# Plots a 2D Color filled Contour plot of the potential in the region
contour = plt.figure()
plt.contourf(x, y, potential.T)
plt.title("Electrostatic Potential Around Line of Charge")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.colorbar(label="Electric Potential (V)")
plt.show()

# Plots a 3D Surface Plot of the potential in the region
projection = plt.figure()
ax = projection.add_subplot(projection="3d")
X, Y = np.meshgrid(x, y)
plt.title("Electrostatic Potential Around Line of Charge")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
#ax.plot_surface(X, Y, potential.T, cmap=plt.cm.plasma)
ax.contour3D(X, Y, potential.T, 50, cmap=plt.cm.plasma)
ax.set_zlabel("Potential (V)")
plt.show()

