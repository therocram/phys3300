# fakeTime
# Library of functions that implements a forward difference fake time
# algorithm to solve the time independent Schrodinger equation (TISE)
# over an arbitrary two dimensional region.
#
# Created by Masen Pitts on 4/10/2022
# Last updated on 4/21/2022

#**************************************************#
# These functions are entirely self-organizing. The only
# function you need to call to execute the full algorithm
# is quantumSolver2D, which has the following inputs and
# outputs:
# 
# quantumSolver2D(xmin, xmax, ymin, ymax, meshRes, potenE, grdGuess, 
#                 excGuess1, excGuess2, boundaryValue, dtau, N)
#
# Inputs:
# - xmin, xmax: minimum and maximum values of grid region in the first 
#   dimension.
#
# - ymin, ymax: minimum and maximum values of grid region in the second 
#   dimension.
#
# - meshRes: integer value that determines the total number of grid
#   points in the first dimension of the grid and solution arrays. The
#   total number of points in the second dimension is determined by
#   the function using ymax and the grid spacing in the first dimension.
#
# - potenE, grdGuess, excGuess1, excGuess2: Functions written by the user
#   that take (X, Y) as inputs, where X and Y are multi-dimensional arrays
#   of the same shape. These functions must be able to handle calculations
#   across entire arrays, so the use of NumPy functions is highly recommended.
#
#
# - boundaryValue: integer or float number that the setBoundaries
#   function will use to enforce the desired boundary condition.
#
# - dtau: step size for fake time variable.
#
# - N: total number of fake time steps the algorithm will run for.
#
#
# Outputs:
# - grdState, excited1, excited2: 2D numpy arrays that contain the
#   solution values for the ground state and first excited state
#   wavefunctions.
#
# - X, Y: 2D numpy arrays containing the values of the grid points.
#   these can be used to visualize the wavefunctions.
#
# - energyList: a one dimensional list containing the energy values
#   calculated by the algorithm. The values are sorted in the array
#   as [ground state, 1st excited, 2nd excited].
#
# For detailed descriptions of the helper functions see the "Functions"
# subsection of Section 2 of the report.
#**************************************************#

import numpy as np

def setBoundaries(array, boundaryValue):
    # Find the ends of the array in each dimension
    edgeY = array.shape[0] - 1
    edgeX = array.shape[1] - 1
    
    # Use a vectorized assignment operation to set boundary values
    array[0,:] = boundaryValue
    array[:,0] = boundaryValue
    array[edgeY,:] = boundaryValue
    array[:,edgeX] = boundaryValue
    
    return array
    
def integral2D(f2D, X, Y):
    d = X[0,1] - X[0,0] # Determine grid spacing
    temp = f2D*d*d
    return np.sum(temp)

def forwardStep(wavefunc, V, dtau, d):
    psiAvg = (1/4)*(wavefunc[left]+wavefunc[right]+wavefunc[up]+wavefunc[down])
    d2 = d*d
    firstTerm = (2/d2)*(psiAvg-wavefunc[nonBound])
    secondTerm = V*wavefunc[nonBound]
    wavefunc[nonBound] += dtau*(firstTerm - secondTerm)
    
    return wavefunc

def normalize(wavefunc, X, Y):
    waveSqr = np.power(wavefunc, 2)
    integral = integral2D(waveSqr, X, Y)
    normConst = 1/np.sqrt(integral)
    wavefunc *= normConst
    
    return wavefunc
    
def orthogonalize(wavefunc1, wavefunc2, X, Y):
    product = wavefunc1*wavefunc2
    innerProd = integral2D(product, X, Y)
    wavefunc1[nonBound] -= innerProd*wavefunc2[nonBound]
    
    return wavefunc1
    
def getEnergy(wavefunc, V, X, Y, d):
    psiAvg = np.zeros(wavefunc.shape)
    waveSqr = np.power(wavefunc, 2)
    
    psiAvg[nonBound] += (1/4)*(wavefunc[left]+wavefunc[right]+wavefunc[up]+wavefunc[down])
    
    firstTerm = -(2/(d*d))*(wavefunc*psiAvg - waveSqr)
    
    potenE = np.zeros(wavefunc.shape)
    potenE[nonBound] += V
    secondTerm = potenE*waveSqr
    
    integrand = firstTerm + secondTerm
    
    return integral2D(integrand, X, Y)

def groundState(wavefunc, X, Y, V, N, dtau, d):
    for tau in range(N):
        forwardStep(wavefunc, V, dtau, d)
        normalize(wavefunc, X, Y)
        
    return wavefunc
    
def firstExcited(wavefunc, grdState, X, Y, V, N, dtau, d):
    for tau in range(N):
        forwardStep(wavefunc, V, dtau, d)
        normalize(wavefunc, X, Y)
        orthogonalize(wavefunc, grdState, X, Y)
        normalize(wavefunc, X, Y)
        
    return wavefunc
    
def secondExcited(wavefunc, grdState, excited1, X, Y, V, N, dtau, d):
    for tau in range(N):
        forwardStep(wavefunc, V, dtau, d)
        normalize(wavefunc, X, Y)
        orthogonalize(wavefunc, grdState, X, Y)
        orthogonalize(wavefunc, excited1, X, Y)
        normalize(wavefunc, X, Y)
        
    return wavefunc

def quantumSolver2D(xmin, xmax, ymin, ymax, meshRes, potenE, grdGuess, excGuess1, 
                    excGuess2, boundaryValue, dtau, N):
    # Helper varibles used to create grids
    x = np.linspace(xmin, xmax, meshRes)
    d = x[1] - x[0]
    y = np.arange(ymin, ymax+d, d)
    meshX = meshRes
    meshY = len(y)
    
    # Tuples used for indexing
    global nonBound, left, right, up, down
    nonBound = (slice(1, meshY-1), slice(1, meshX-1))
    left = (slice(0, meshY-2), slice(1, meshX-1))
    right = (slice(2, meshY), slice(1, meshX-1))
    up = (slice(1, meshY-1), slice(2, meshX))
    down = (slice(1, meshY-1), slice(0, meshX-2))
    
    X, Y = np.meshgrid(x, y)
    
    # Initialize potential energy and wavefunction arrays
    V = potenE(X[nonBound], Y[nonBound])
    grdState = grdGuess(X, Y)
    excited1 = excGuess1(X, Y)
    excited2 = excGuess2(X, Y)

    # Enforce boundary conditions
    setBoundaries(grdState, boundaryValue)
    setBoundaries(excited1, boundaryValue)
    setBoundaries(excited2, boundaryValue)

    normalize(grdState, X, Y)
    normalize(excited1, X, Y)
    normalize(excited2, X, Y)

    groundState(grdState, X, Y, V, N, dtau, d)
    Egrd = getEnergy(grdState, V, X, Y, d)
    
    firstExcited(excited1, grdState, X, Y, V, N, dtau, d)
    E1 = getEnergy(excited1, V, X, Y, d)
    
    secondExcited(excited2, grdState, excited1, X, Y, V, N, dtau, d)
    E2 = getEnergy(excited2, V, X, Y, d)
    
    energyList = [Egrd, E1, E2]
    print("System Energy Eigenvalues in Natural Units:")
    for E in energyList:
        print("\t" + str(E))
    
    return grdState, excited1, excited2, X, Y, energyList