# odeSolvers
# File of helper functions used to solve ODE's using
# the Euler, 4th order Runge-Kutta, and Leap Frog algorithms.

"""
Created on Sun Mar 13 16:45:35 2022

@author: masenpitts
"""
import sys
import numpy as np

# Each solver method, except where otherwise specified, has the following
# parameter configurations:
#   y - array representing the dynamical variable vector
#   g - array representing the generalized velocity vector
#   d - step size of the algorithm

def euler(y, g, d):
    for i in range(len(y)):
        y[i] += g[i](y)*d
    return y

def rungeKutta4(y, g, d):
    
    c1 = c2 = c3 = c4 = np.zeros(len(y))
    
    for i in range(len(y)):
        c1[i] = d * g[i](y)
        c2[i] = d * g[i](y + 0.5*c1)
        c3[i] = d * g[i](y + 0.5*c2)
        c4[i] = d * g[i](y + c3)
        
    return y + (c1 + 2*c2 + 2*c3 + c4)/6

# The Leap Frog method is only applicable to second order systems governed
# by an equation that can be written in the form x'' = f(x)
#
# The "y" lists is assumed to store the positions and
# velocities of a series of systems. That is, "y" is assumed to have the form
# y = [pos1, vel1, pos2, vel2, ...] 
#
# The "g" list is assumed to store the generalized velocities of the quantities
# stored in the "y" list. That is, "g" is assumed to have the form
# g = [vel1, accel1, vel2, accel2, ...] 
def leapFrog(y, g, d):
    # This loops works with values at 2 different index positions in each
    # array and so iterates in steps of two.
    for i in range(len(y), 2):
        oldg = g[i+1](y) # Stores old acceleration
        y[i] += g[i](y)*d + 0.5*oldg*d*d # Update positions
        newg = g[i+1](y) # Stores new acceleration
        y[i+1] += 0.5*(oldg + newg)*d # Update velocities
    
    return y
 
# This overhead function allows the user to easily select a solver method.   
def solver(y, g, d, selection):
    if selection == "1":
        return euler(y, g, d)
    elif selection == "2":
        return rungeKutta4(y, g, d)
    elif selection == "3":
        return leapFrog(y, g, d)
    else:
        print("Error: Enter valid selection for solver method")
        sys.exit()