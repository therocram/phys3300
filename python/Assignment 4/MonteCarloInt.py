# MonteCarloInt
"""
Created on Fri Feb 18 08:58:05 2022

@author: masenpitts
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
from random import random

# Used to obtain function values for f(x) = e^x
def getFunc(x):
    return np.exp(x)

# Calculates the absolute value of the percent difference between provided approximate
# and exact values.
def percentDiff(approx, exact):
    return 100*abs((approx-exact)/exact)

# Performs a Monte-Carlo Integration from x = a to x = b using n total points
# in the estimation.
# Uses the getFunc method to obtain values of f(x)
def monteCarlo(a, b, n):
    # Ensures that integration limits are valid
    if a > b:
        print("Error: Choose valid limits of integration")
        sys.exit()
        
    interval = b - a
    fN = 0 # Stores sum approximating average value of f(x) over the interval
    
    # Aprroximates the average value of f(x) over the interval using
    # a total of n randomly determined points along the interval
    for i in range(n):
        randX = interval*random() + a
        fN += getFunc(randX)
    
    return interval*fN/n

# Range over which the function is integrated
leftLim = 0
rightLim = 1

### Performs single estimation of integral value using Monte-Carlo integration
# with an N determined by the user, printing the result, the actual value of
# the integral, and the percent difference.
n = int(input("Enter total number of points to be used in estimation: "))

approx = monteCarlo(leftLim, rightLim, n)
actual = getFunc(rightLim) - getFunc(leftLim)

print("Result = {:.6f}".format(approx))
print("Actaul Value = {:.6f}".format(actual))
print("Percent Difference = {:.6f}".format(percentDiff(approx, actual)))
###

# Creates a list of N values frin
tableN = np.arange(100, 1100, 100)

### Creates a tab seperated table of the Monte-Carlo approximation value and 
# percent difference for each value of N in the tableN list
print("\n\nN\tResult\tPercent Diff")

tableLine = "{num:d}\t{result:.3f}\t{perDiff:.3f}"

for N in tableN:
    approx = monteCarlo(leftLim, rightLim, N)
    print(tableLine.format(num = N, result = approx,
                           perDiff = percentDiff(approx, actual)))
###    

### Creates a list of N values from nMin to nMax with the given step and is used
# to plot I (the value of the integral) as a function of N.
nMin = 1
nMax = 1000
step = 1

plotN = np.arange(nMin, nMax + step, step)
###

# Used to store the calculated value of I at each value in plotN
Idata = []

# Obtains Idata values using Monte-Carlo Integration
for N in plotN:
    Idata.append(monteCarlo(leftLim, rightLim, N))
    
# Creates a plot of the value of the estimated value of the integral vs. the
# value of N.
plt.plot(plotN, Idata, color="red")
plt.title("Value of Integral as a function of N")
plt.xlabel("N")
plt.ylabel("I")


