# rectangularBox
# Implements the methods and functions of fakeTime.py to 
# solve for the ground state and first two excited states
# of a particle in a rectangular infinite square well.
#
# Created by Masen Pitts on 4/17/2022
# Last updated on 4/17/2022
"""
Created on Sun Apr 17 15:52:23 2022

@author: masenpitts
"""

import fakeTime as ft
import matplotlib.pyplot as plt

def V(X, Y):
    return 0

def grdGuess(X, Y):
    return X*(a-X)*Y*(b-Y)

def excGuess1(X, Y):
    x = X*(a-X)*((0.5*a) - X)
    y = Y*(b-Y)*((0.5*b) - Y)
    return x*y

def excGuess2(X, Y):
    x = X*(a-X)*((a/3.0) - X)*((2*a/3.0) - X)
    y = Y*(b-Y)*((b/3.0) - Y)*((2*b/3.0) - Y)
    return x*y

a = 3
b = 5

meshRes = 50
boundaryValue = 0
dtau = 0.0003
N = 500

grdState, excited1, excited2, energyList = ft.quantumSolver2D(0, 
                                             a, 0, b, meshRes, 
                                             V, grdGuess, excGuess1, 
                                             excGuess2, boundaryValue, 
                                             dtau, N)

grdPlot = plt.figure()
axG = grdPlot.add_subplot(projection="3d")
plt.title("Ground State Wavefunction")
plt.show()