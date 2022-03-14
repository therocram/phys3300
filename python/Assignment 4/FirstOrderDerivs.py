# FirstOrderDerivs
"""
Created on Mon Feb 14 14:50:37 2022

@author: masenpitts
"""

import numpy as np
import matplotlib.pyplot as plt

# Approximates the derivative at a point in a given set of data
# using a single point or "forward difference" approximation
def forwardDiff(x, f, i, dx):
    return (f[i+1] - f[i])/dx

# Approximates the derivative at a point in a given set of data
# using three surrounding points or "centered difference" approximation
def centeredDiff(x, f, i, dx):
    return (f[i+1] - f[i-1])/(2*dx)

# Approximates the derivative at a point in a given set of data
# using five surrounding points.
def fivePoint(x, f, i, dx):
    return (f[i-2] - 8*f[i-1] + 8*f[i+1] - f[i+2])/(12*dx)

# Approximates the value of the integral over the given data sets using 
# Simpson's Rule. Integrates "f" with respect to "x"
def simpsonRule(x, f):
    integral = 0
    n = len(f) - 1
    for i in range(1, n, 2):
        dx = x[i+1] - x[i]
        integral += dx*(f[i-1] + 4.0*f[i] + f[i+1])/3.0
    if (n+1) % 2 == 0:
        return integral + dx*(5.0*f[n] + 8.0*f[n-1] - f[n-2])/12.0
    else:
        return integral
    
# Calculates the absolute value of the percent difference between provided approximate
# and exact values.
def percentDiff(approx, exact):
    return 100*abs((approx-exact)/exact)


#**** USE THIS VARIABLE TO MODIFY NUMBER OF DATA POINTS
nPoints = 51 # Total number of starting data points
#****


min = 0         # Determines domain over which data is plotted
max = 2*np.pi   # 

# Stores the distance of seperation between data points
interval = (max - min)/nPoints

# Creates two 1D grids of data values representing the function input and
# output values respectively.
x = np.arange(min, max, interval)
f = np.sin(x)

# This 1D grid stores the analytically determined derivative of the function
# data
trueDfDx = np.cos(x)

# Lists storing the numerically determined derivative values for each of the
# approximation types.
twoPointData = np.zeros(len(x)-1)
threePointData = np.zeros(len(x)-2)
fivePointData = np.zeros(len(x)-4)

# Lists storing the x positions corresponding to each derivative data set.
x2Point = [0]*len(twoPointData)
x3Point = [0]*len(threePointData)
x5Point = [0]*len(fivePointData)

# These loops determine the derivative at every point on each grid using the
# respective differentiation technique.
for i in range(len(twoPointData)):
    twoPointData[i] = forwardDiff(x, f, i, interval)
    x2Point[i] = x[i]

for j in range(len(threePointData)):
    threePointData[j] = centeredDiff(x, f, j + 1, interval)
    x3Point[j] = x[j + 1]
 
for k in range(len(fivePointData)):
    fivePointData[k] = fivePoint(x, f, k + 2, interval)
    x5Point[k] = x[k + 2]

# Plots sin(x), its analytical derivative cos(x), and each of the numerically
# determined derivative data sets. 
plt.plot(x, f, label="f(x) = sin(x)")
plt.plot(x, trueDfDx, label="f'(x) = cos(x)")
plt.plot(x2Point, twoPointData, label="Two Point")
plt.plot(x3Point, threePointData, label="Three Point")
plt.plot(x5Point, fivePointData, label="Five Point")
plt.title("Derivative Approximations for Sin(x) with " + str(nPoints) + " Data Points")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.show()

#############

# The actual value of the integral of cos(x) from x = 0 to x = 2Pi
actual = np.sin(2*np.pi) - np.sin(0)

# Integrates over each of the derivative data sets to approximate 
# the value of the integral using Simpson's Rule
integral2Point = simpsonRule(x2Point, twoPointData)
integral3Point = simpsonRule(x3Point, threePointData)
integral5Point = simpsonRule(x5Point, fivePointData)

# String giving the formatting for printing the approximations of the integral
# for each derivative type.
toString = "{name:s}:\nIntegral = {approx:.6f}\n"

# Prints the numerically calculated integral values for each
# derivative type.
print(toString.format(name = "One Point", approx = integral2Point))
print(toString.format(name = "Three Point", approx = integral3Point))
print(toString.format(name = "Five Point", approx = integral5Point))