# fakeTime
# Library of functions that implements a forward difference fake time
# algorithm to solve the time independent Schrodinger equation (TISE)
# over an arbitrary two dimensional region.
#
# Created by Masen Pitts on 4/10/2022
# Last updated on 4/17/2022
"""
Created on Sun Apr 10 17:26:37 2022

@author: masenpitts
"""

import numpy as np
#import matplotlib as plt

def setBoundaries(array, boundaryValue):
    # Find the ends of the array in each dimension
    edgeX = array.shape[0] - 1
    edgeY = array.shape[1] - 1
    
    # Use a vectorized assignment operation to set boundary values
    array[0,:] = boundaryValue
    array[:,0] = boundaryValue
    array[edgeX,:] = boundaryValue
    array[:,edgeY] = boundaryValue

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
    
def simpson2D(f2D, X, Y):
    '''
    intList = []
    
    for y in range(Y.shape[0]):
        intList.append(simpsonRule(f2D[y], X[y]))
        
    return simpsonRule(intList, Y[:,0])
    '''
    d = X[0,1] - X[0,0]
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
    integral = simpson2D(waveSqr, X, Y)
    #print(integral)
    normConst = 1/np.sqrt(integral)
    wavefunc *= normConst
    
    return wavefunc
    
def orthogonalize(wavefunc1, wavefunc2, X, Y):
    product = wavefunc1*wavefunc2
    innerProd = simpson2D(product, X, Y)
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
    
    return simpson2D(integrand, X, Y)

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
'''
def plotData(grdState, excited1, excited2, X, Y):
    grdPlot = plt.figure()
    axG = grdPlot.add_subplot(projection="3d")
    plt.title("Ground State Wavefunction")
    plt.show()
    return
'''
def quantumSolver2D(xmin, xmax, ymin, ymax, meshRes, potenE, grdGuess, excGuess1, 
                    excGuess2, boundaryValue, dtau, N):
    x = np.linspace(xmin, xmax, meshRes)
    d = x[1] - x[0]
    y = np.arange(ymin, ymax+d, d)
    meshX = meshRes
    meshY = len(y)
    global nonBound, left, right, up, down
    nonBound = (slice(1, meshY-1), slice(1, meshX-1))
    left = (slice(0, meshY-2), slice(1, meshX-1))
    right = (slice(2, meshY), slice(1, meshX-1))
    up = (slice(1, meshY-1), slice(2, meshX))
    down = (slice(1, meshY-1), slice(0, meshX-2))
    
    X, Y = np.meshgrid(x, y)
    
    V = potenE(X[nonBound], Y[nonBound])
    grdState = grdGuess(X, Y)
    excited1 = excGuess1(X, Y)
    excited2 = excGuess2(X, Y)

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
    
    #plotData(grdState, excited1, excited2, X, Y)
    
    energyList = [Egrd, E1, E2]
    print("System Energy Eigenvalues in Hartree Atomic Units:")
    for E in energyList:
        print("\t" + str(E))
    
    return grdState, excited1, excited2, energyList