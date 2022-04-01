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
    
    # Stores the parameters of the terms
    # up to the 4th order term of the series expansion of
    # the solution function
    c1 = np.zeros(len(y))
    c2 = np.zeros(len(y))
    c3 = np.zeros(len(y))
    c4 = np.zeros(len(y))
    
    ytemp = np.array(y)

    # Assign values to the parameters
    for i in range(len(y)):
        c1[i] = d * g[i](y)
        ytemp[i] = y[i] + 0.5*c1[i]
        c2[i] = d * g[i](ytemp)
        ytemp[i] = y[i] + 0.5*c2[i]
        c3[i] = d * g[i](ytemp)
        ytemp[i] = y[i] + c3[i]
        c4[i] = d * g[i](ytemp)
        
        # Approximate the solution value using a 4th order 
        # series expansion
        y[i] += (c1[i] + 2*c2[i] + 2*c3[i] + c4[i])/6.
        
    return y

# The Leap Frog method is only applicable to second order systems governed
# by an equation that can be written in the form x'' = f(x)
# Unlike the other subroutines, this function expects an array to store
# solution values "solutionList" and solves the system based on the shape
# of the array.
# This method is only for second order equations, so "y" and "g" are both
# expected to have two elements each corresponding to position and velocity
def leapFrog(y, g, dt, solutionList): 
    # Add initial conditions to solution list
    solutionList[0] = y
    
    # Advance velocity by half a time step
    y[1] += 0.5*dt*g[1](y)
    
    N = solutionList.shape[0]
    
    for t in range(1, N):
        y[0] += y[1]*dt
        y[1] += g[1](y)*dt
        
        # Store unsynced velocity
        unsyncVel = y[1]
        
        # Synchronize velocity using Euler and add it to output
        # array
        y[1] -= 0.5*dt*g[1](y)
        solutionList[t] = y
        
        # Restore unsynced veloctiy to continue calculations
        y[1] = unsyncVel
    
    return solutionList

# This overhead function allows the user to easily select a solver method.   
def solver(y, g, d, selection, solutionList = [], t = 0):
    # Time independent variations
    if t == 0: 
        if selection == "1":
            return euler(y, g, d)
        elif selection == "2":
            return rungeKutta4(y, g, d)
        elif selection == "3":
            if solutionList == []:
                print("Error: Provide list to store solution values for Leap Frog")
                sys.exit()
            return leapFrog(y, g, d, solutionList)
        else:
            print("Error: Enter valid selection for solver method")
            sys.exit()
    # Time dependent variations
    else:
        if selection == "2":
            return rungeKutta4t(y, g, d, t)
        else:
            print("Error: Enter valid selection for solver method")
            sys.exit()
        
# Time-dependent variation of rk4
def rungeKutta4t(y, g, d, t):
    
    # Stores the parameters of the terms
    # up to the 4th order term of the series expansion of
    # the solution function
    c1 = np.zeros(len(y))
    c2 = np.zeros(len(y))
    c3 = np.zeros(len(y))
    c4 = np.zeros(len(y))
    
    ytemp = np.array(y)

    # Assign values to the parameters
    for i in range(len(y)):
        c1[i] = d * g[i](y, t)
        ytemp[i] = y[i] + 0.5*c1[i]
        c2[i] = d * g[i](ytemp, t + 0.5*d)
        ytemp[i] = y[i] + 0.5*c2[i]
        c3[i] = d * g[i](ytemp, t + 0.5*d)
        ytemp[i] = y[i] + c3[i]
        c4[i] = d * g[i](ytemp, t + d)
        
        # Approximate the solution value using a 4th order 
        # series expansion
        y[i] += (c1[i] + 2*c2[i] + 2*c3[i] + c4[i])/6
        
    return y