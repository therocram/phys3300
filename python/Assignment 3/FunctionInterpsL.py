# FunctionInterpsL
# Implementation of linear interpolation on some example functions
"""
Created on Fri Feb  4 20:26:44 2022

@author: masenpitts
"""

import matplotlib.pyplot as plt
import numpy as np
import LagrangeInterp as LI

# Get data and model parameters from console interface
funcSelect, xdata, fdata, x_value = LI.runConsole()

# Allows user to determine the order of the Lagrangian interpolation
n = int(input("Enter order of Lagrangian Interpolation: "))

# Perform Lagrangian interpolation at the point chosen by the user
f_LagInterp = LI.interp(xdata, fdata, n, x_value)

# Plot the true data along with the interpolated point
plt.plot(xdata, fdata, "bo", label="True Data")
plt.plot(x_value, f_LagInterp, "go", label="Target Point")

# Create a grid of values to approximate a continuous interpolatio
xgrid = np.arange(xdata[0], xdata[len(xdata)-1], 0.01)
finterp = [] # List of interpolated points
error = [] # List of error values at each point

# Perform Lagrangian interpolation and error calculation at each point on the 1D grid
for each in xgrid:
    interpf = LI.interp(xdata, fdata, n, each)
    actualf = LI.selectFunc(each, funcSelect)
    
    finterp.append(interpf)
    error.append(LI.percentDiff(interpf, actualf)/100)   
    
# Plot interpolated grid and error values
plt.plot(xgrid, finterp, color="red", label="Interpolation")
plt.plot(xgrid, error, color="orange", label="Uncertainty")
plt.legend()

# Print results for interpolated point directly chosen by user
LI.printResults(f_LagInterp, x_value, funcSelect)

plt.show()