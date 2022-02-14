# LagrangeInterp
# A file of helper functions used for linear
# and Lagrangian interpolations.
"""
Created on Mon Jan 31 14:45:39 2022

@author: masenpitts
"""

import numpy as np
import sys

# Interpolation Functions
##############################################################################

# Helper method that determines the index value of the closest
# true data point BEFORE the given x position where interpolation
# is to occur
# x - list of independent variable positions in true data
# x_pos - target x position where interpolation is to occur
def getIndex(x, x_pos):
    for index in range(len(x)):
        if x[index] > x_pos:
            return index-1
        
# Performs linear interpolation at the target position
# using the given true data.
# x - list of independent variable positions in true data
# f - list of depended variable values in true data
# i - index of closest true data point BEFORE target position
# x_pos - target x position where interpolation is to occur
def linearInterp(x, f, i, x_pos):
    if x[i+1] < x[i]:
        nextI = i-1
    else:
        nextI = i+1
    result = f[i] + (x_pos - x[i])/(x[nextI] - x[i])*(f[nextI] - f[i])
    return result

# Performs an nth order Lagrangian interpolation at the target position
# using the given true data.
# x - list of independent variable positions in true data
# f - list of depended variable values in true data
# n - gives the order of the interpolation
# i - index of closest true data point BEFORE target position
# x_pos - target x position where interpolation is to occur
def lagrangeInterp(x, f, n, i, x_pos):
    indexList = []  # List storing indices of true data points used in interpolation
    
    # A series of helper variables that assist in adding elements to the index list
    
    if x[i+1] < x[i]:
        forward = False
    else:
        forward = True
    countforward = 0
    countbackward = 0
    
    indexList.append(i)
    
    # Alternates between adding true data points located before and after the target
    # position to the index list.
    # After this loop the index list contains n + 1 elements.
    for each in range(n):
        if forward:
            indexList.append(i + countforward + 1)
            countforward += 1
        else:
            indexList.append(i - countbackward - 1)
            countbackward += 1
        forward = not forward
        
    # Helper variables used to compute interpolation
    result = 0
    product = 1
    
    # Performs nth order Lagrangian interpolation using the collected data points
    for j in indexList:
        for k in indexList:
            if k != j:
                product *= (x_pos - x[k])/(x[j] - x[k])
        result += f[j]*product
        product = 1
    return result

# Primary function. Determines the proximity of the point to be interpolated
# to the boundaries of the region and chooses to use linear or Lagrangian
# interpolation accordingly.
# x - list of independent variable positions in true data
# f - list of depended variable values in true data
# n - gives the order of the interpolation
# x_pos - target x position where interpolation is to occur
def interp(x, f, n, x_pos):
    iBefore = getIndex(x, x_pos) # Determines index of true data point before
                                 # target point
    
    # Determines whether Lagrangian interpolation is suitable by
    # considering the proximity of target point to the bounds of the data
    # set.
    if (iBefore < n/2) or (iBefore > len(x) - 1 - n/2) or (n==1):
        result = linearInterp(x, f, iBefore, x_pos)
    else:
        result = lagrangeInterp(x, f, n, iBefore, x_pos)
    return result

##############################################################################


# Console Functions
##############################################################################

# Allows the computation of function values at a position or list of positions
# from a given selection of functions.
# x - value the selected function is to be evaluated at
# selection - choice of function
def selectFunc(x, selection):
    if selection.strip() == "1":
        return np.sin(x**2)
    elif selection.strip() == "2":
        return np.exp(np.sin(x))
    elif selection.strip() == "3":
        return 0.2/((x - 3.2)**2 + 0.04)
    else:
        sys.exit("\nError: Please select an equation using its corresponding number")

# Calculates the absolute value of the percent difference between provided approximate
# and exact values.
def percentDiff(approx, exact):
    return 100*abs((approx-exact)/exact)

# Provides a console-based user interface for performing interpolations.
def runConsole():
    equationList = "1. f(x) = sin(x^2)\n2. f(x) = exp(sin(x))\n3. f(x) = 0.2/((x-3.2)^2 + 0.04)"
    
    print("Welcome to the Interpolator 1.4\nSelect an equation from the list below\n" + equationList)
    
    functionSelection = input("Select an equation from the list above: ")
    
    dataPointsSelection = int(input("Enter total number of desired starting points: "))
    
    if dataPointsSelection <= 0:
        sys.exit("Error: Value provided for starting points must be a nonzero positive integer")
    
    minRangeSelection = float(input("Enter Minimum of Data Range: "))
    maxRangeSelection = float(input("Enter Maximum of Data Range: "))
    
    if minRangeSelection > maxRangeSelection:
        sys.exit("Error: Enter Valid Data Range")
    
    interpSelection = float(input("Enter Location of Target Point: "))
        
    if (interpSelection < minRangeSelection) or (interpSelection > maxRangeSelection):
        sys.exit("Error: Ensure that target point is within provided data range")
    
    # Creates and returns the lists of true data based on the user input
    interval = (maxRangeSelection - minRangeSelection)/dataPointsSelection
    xdata = np.arange(minRangeSelection, maxRangeSelection + interval, interval)
    fdata = selectFunc(xdata, functionSelection)
    
    return functionSelection, xdata, fdata, interpSelection

# Prints the interpolation results and accuracy measurements
def printResults(f_interp, x_value, funcSelect):
    actual = selectFunc(x_value, funcSelect)
    
    print("\nInterpolated Value: " + str(f_interp))
    print("Actual Value: " + str(actual))
    print("Percent Difference: " + str(percentDiff(f_interp, actual)))
##############################################################################