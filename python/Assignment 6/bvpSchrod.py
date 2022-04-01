# bvpSchrod
# Solve the time independent Schrodinger equation (TISE) using boundary 
# value problem solving methods.
# Both the harmonic and anharmonic oscillators are modeled within a
# square infinite well of width 2*L (see parameters below)
"""
Created on Tue Mar 22 17:03:10 2022

@author: masenpitts
"""

import numpy as np
import matplotlib.pyplot as plt

# Relevant physical quantities
e = 1.6022e-19 # Conversion factor from eV to J
m = 9.1094e-31 # Mass of electron in kg
hb = 1.0546e-34 # Hbar in J*s

# Constant Parameters
V0 = 50*e # J
a = 1e-11 # m

L = 6*a # Half-width of infinite well at boundary
        # (Use 10*a for Harmonic osciallator, 6*a for Anharmonic)
N = 1000 # Total number of points in mesh grid
dx = L/N # Step size for x 

tolerance = e/10000000

# Potential energy function
def V(x):
    # Harmonic Oscillator: V0*(x*x)/(a*a)
    # Anharmonic Oscillator: V0*(x**4)/(a**4)
    return V0*(x**4)/(a**4)

# Implements a 4th Order Runge-Kutta integrator
def rk4(y, g, d, x, E):
    
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
        c1[i] = d * g[i](y, x, E)
        ytemp[i] = y[i] + 0.5*c1[i]
        c2[i] = d * g[i](ytemp, x + 0.5*d, E)
        ytemp[i] = y[i] + 0.5*c2[i]
        c3[i] = d * g[i](ytemp, x + 0.5*d, E)
        ytemp[i] = y[i] + c3[i]
        c4[i] = d * g[i](ytemp, x + d, E)
        
        # Approximate the solution value using a 4th order 
        # series expansion
        y[i] += (c1[i] + 2*c2[i] + 2*c3[i] + c4[i])/6.
        
    return y

# Components of generalized velocity vector of system
# Note: In the text g0 is phi = dpsi/dx and g1 is dphi/dx
def g0(y, x, E):
    return y[1]

def g1(y, x, E):
    return (-2*m/hb**2)*(E - V(x))*y[0]

# Solves the TISE for a given energy value and returns the solution
# value at the endpoint
def trial(E, g, d, xGrid):
    y = [0, 1] # Arbitrary inital conditions (but y[0] MUST be 0)
    
    for x in xGrid:
        y = rk4(y, g, d, x, E)
    
    return y[0]

# Uses the secant method to find an energy eigenvalue that solves the
# TISE and satisfies the boundary conditions within a specified 
# tolerance, given a starting range (E1, E2)
def findE(g, d, E1, E2, tolerance, xGrid):
    psiNext = trial(E1, g, d, xGrid)
    
    while abs(E1-E2) > tolerance:
        psi = psiNext
        psiNext = trial(E2, g, d, xGrid)
        E1, E2 = E2, E2 - psiNext*(E2-E1)/(psiNext-psi)
        
    return E2

# Builds a wavefunction on the given 1D grid of x values by solving the
# TISE with a given energy value.
def solveSchrod(g, d, E, xGrid):
    y = [0, 1] # Arbitrary inital conditions (but y[0] MUST be 0)
    psi = [] # Stores wavefunction values
    
    for x in xGrid:
        y = rk4(y, g, d, x, E)
        psi.append(y[0])
        
    return psi

# Uses Simpson's Rule to estimate the integral of a given list of
# function values over the given list of x values.
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

# Normalizes the given wavefunction over the gived 1D grid of
# x values. Here it is assumed that the wavefunction has decreased
# to zero at both ends of the grid, such that integrating over the
# grid is approximately the same as integrating the wavefunction over
# all of space.
def normalize(psi, xGrid):
    psiTemp = np.array(psi)
    probDens = np.power(psiTemp, 2)
    normConst = 1/np.sqrt(simpsonRule(xGrid, probDens))
    psiTemp = normConst*psiTemp
    return psiTemp

# Returns the probability of finding the particle corresponding
# to the given normalized probability density function over the given
# interval in position space, centered at xPos. In this case
# xList represents is the 1D grid over which the probability density
# is defined.
def probability(probDensNorm, xList, xPos, interval):
    xGrid = np.array(xList)
    probDens = np.array(probDensNorm)
    
    dx = xGrid[1] - xGrid[0] # Distance between grid points
    indexShift = int(0.5*interval/dx)
    targetIndex = -1
    
    # Finds the index postition of the grid point nearest to the 
    # target postion
    for i in range(len(xGrid)-1):
        if xGrid[i] == xPos:
            targetIndex = i
        elif xGrid[i] < xPos and xGrid[i+1] > xPos:
            avgX = (xGrid[i]+xGrid[i+1])/2
            if xPos > avgX:
                targetIndex = i+1
            else:
                targetIndex = i
    
    # Returns -1 if target index is too close to the endpoints
    # of the grid
    if targetIndex < indexShift or targetIndex > len(xGrid) - indexShift:
        return -1
    
    # Splices the given arrays to only include the portion of the grid
    # included in the interval of integration
    xRange = xGrid[targetIndex-indexShift: targetIndex+indexShift+1]
    probDensRange = probDens[targetIndex-indexShift: targetIndex+indexShift+1]
    
    # Uses Simpson Rule to calculate probability over interval
    prob = simpsonRule(xRange, probDensRange)
    return prob

# Guesses for Harmonic Oscillator:
# Ground: 100*e
# 1st Excited: 250*e
# 2nd Excited: 500*e
# 3rd Excited: 800*e
#
# Guesses for Anharmonic Oscillator:
# Ground: 100*e
# 1st Excited: 500*e
# 2nd Excited: 1000*e
# 3rd Excited: 2000*e

# Generalized velocity vector 
g = [g0, g1]

# Initial rough guess for ground state energy
Eguess = 100*e

# Grid of x values used to plot wavefunctions
xGrid = np.arange(-L, L, dx)

# Solves for ground state energy eigenvalue within specified
# tolerance
Eground = findE(g, dx, 0, Eguess, tolerance, xGrid)

# Builds ground state wavefunction using calculated eigenvalue 
psiGround = solveSchrod(g, dx, Eground, xGrid)

# Similarily solves for and builds the first three excited states
# of the system.
Eguess = 500*e
E1 = findE(g, dx, 0, Eguess, tolerance, xGrid)
psi1 = solveSchrod(g, dx, E1, xGrid)

Eguess = 1000*e
E2 = findE(g, dx, 0, Eguess, tolerance, xGrid)
psi2 = solveSchrod(g, dx, E2, xGrid)

Eguess = 2000*e
E3 = findE(g, dx, 0, Eguess, tolerance, xGrid)
psi3 = solveSchrod(g, dx, E3, xGrid)

# List storing the relevant wavefunctions of the system
wavefList = [psiGround, psi1, psi2, psi3]

# Normalizes each of the wavefunctions in the list
for i in range(len(wavefList)):
    wavefList[i] = normalize(wavefList[i], xGrid)

# Plots the 1D wavefunctions of the ground state and first three excited
# states of the system
wavefunctions = plt.figure()
plt.plot(xGrid, wavefList[0], color="blue", label="Ground State")
plt.plot(xGrid, wavefList[1], color="red", label="1st Excited")
plt.plot(xGrid, wavefList[2], color="green", label="2nd Excited")
plt.plot(xGrid, wavefList[3], color="orange", label="3rd Excited")
#plt.title("Wavefunctions for Quantum Harmonic Oscillator")
plt.title("Wavefunctions for Quantum Anharmonic Oscillator")
plt.xlabel("x (m)")
plt.ylabel("psi(x)")
plt.xlim(-L, L)
plt.legend()

# List storing values of calculated energy eigenvalues in eV
EList = np.array([Eground, E1, E2, E3])/e

# Plots an energy level diagram of the system
energyLevels = plt.figure()
for E in EList:
    plt.axhline(E, color="purple")
plt.ylabel("Energy (eV)")
plt.xticks([])
#plt.title("Energy Level Diagram for Quantum Harmonic Oscillator")
plt.title("Energy Level Diagram for Quantum Anharmonic Oscillator")

interval = a # Interval over which the probability distribution values
             # are calculated
probDistList = [[], [], [], []] # List that stores the probability
                                # distributions associated with each wavefunction
                                
xRange = [] # List storing the x values over which it is valid to
            # calculate probability distribution values

# Probability density of the ground state wavefunction
probDens = np.power(np.array(wavefList[0]), 2)

# Builds the probability distribution and corresponding 1D grid
# pf x values for the ground state wavefunction
for x in xGrid:
    prob = probability(probDens, xGrid, x, interval)
    if prob != -1:
        probDistList[0].append(prob)
        xRange.append(x)

# Similarily builds the probabilty distributions for each of the
# other wavefunctions in the system
for i in range(1, len(probDistList)):
    probDens = np.power(np.array(wavefList[i]), 2)
    
    for x in xRange:
        prob = probability(probDens, xGrid, x, interval)
        probDistList[i].append(prob)
    
# Plots the Probability distributions for the ground state and first
# three excited states of the system.
probabilityDists = plt.figure()
plt.plot(xRange, probDistList[0], color="blue", label="Ground State")
plt.plot(xRange, probDistList[1], color="red", label="1st Excited")
plt.plot(xRange, probDistList[2], color="green", label="2nd Excited")
plt.plot(xRange, probDistList[3], color="orange", label="3rd Excited")
#plt.title("Probability Distributions for Quantum Harmonic Oscillator")
plt.title("Probability Distributions for Quantum Anharmonic Oscillator")
plt.xlabel("x (m)")
plt.ylabel("P(x)")
plt.xlim(-(L-0.5*interval), L-0.5*interval)
plt.legend()

plt.show()
