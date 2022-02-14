# FunctionInterps
# Implementation of linear interpolation on some example functions
"""
Created on Wed Feb  2 14:37:49 2022

@author: masenpitts
"""

import matplotlib.pyplot as plt
import numpy as np
import LagrangeInterp as LI

# Get data and model parameters from console interface
funcSelect, xdata, fdata, x_value = LI.runConsole()

# Perform linear interpolation at the point chosen by the user
iBefore = LI.getIndex(xdata, x_value)
f_linInterp = LI.linearInterp(xdata, fdata, iBefore, x_value)

# Plot the true data along with the interpolated point
plt.plot(xdata, fdata, "bo", label="True Data")
plt.plot(x_value, f_linInterp, "go", label="Target Point")

# Create a grid of values to approximate a continuous interpolation
xgrid = np.arange(xdata[0], xdata[len(xdata)-1], 0.01)
finterp = [] # List of interpolated points
error = [] # List of error values at each point

# Perform linear interpolation and error calculation at each point on the 1D grid
for each in xgrid:
    eachI = LI.getIndex(xdata, each)
    interpf = LI.linearInterp(xdata, fdata, eachI, each)
    actualf = LI.selectFunc(each, funcSelect)
    
    finterp.append(interpf)
    error.append(LI.percentDiff(interpf, actualf)/100)   
    
# Plot interpolated grid and error values
plt.plot(xgrid, finterp, color="red", label="Interpolation")
plt.plot(xgrid, error, color="orange", label="Uncertainty")
plt.legend()

# Print results for interpolated point directly chosen by user
LI.printResults(f_linInterp, x_value, funcSelect)

plt.show()