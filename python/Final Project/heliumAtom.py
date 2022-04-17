# heliumAtom
# Implements the methods and functions of fakeTime.py to 
# solve for the ground state and first two excited states
# of the electrons in a simplified helium atom model
#
# Created by Masen Pitts on 4/10/2022
# Last updated on 4/17/2022
"""
Created on Sun Apr 10 18:07:14 2022

@author: masenpitts
"""

import fakeTime as ft
import numpy as np
import matplotlib as plt

def V(R1, R2):
    return -2/R1 -2/R2 + 1/np.fmax(R1,R2)

def grdGuess(R1, R2):
    return R1*np.exp(-R1)*R2*np.exp(-R2)

def excGuess1(R1, R2):
    hydR1 = (1-(1/2)*R1)*R1*np.exp(-R1/2)
    hydR2 = (1-(1/2)*R2)*R2*np.exp(-R2/2)
    return hydR1*hydR2

def excGuess2(R1, R2):
    hydR1 = (1-(2/3)*R1+(2/27)*np.power(R1,2))*R1*np.exp(-R1/3)
    hydR2 = (1-(2/3)*R2+(2/27)*np.power(R2,2))*R2*np.exp(-R2/3)
    return hydR1*hydR2

rmin = 0
rmax = 15

meshRes = 300
boundaryValue = 0
dtau = 0.0012
N = 5000

grdState, excited1, excited2, energyList = ft.quantumSolver2D(rmin, 
                                             rmax, rmin, rmax, meshRes, 
                                             V, grdGuess, excGuess1, 
                                             excGuess2, boundaryValue, 
                                             dtau, N)

# Debug simpson2D
# Start making plots to visualize stuff
# Next Steps:
#   -Improve accuracy w/ simpson rule, better guesses, larger region,
#    more grid points, etc.
#   -Are plotted wavefunctions reduced radial wavefunctions?



