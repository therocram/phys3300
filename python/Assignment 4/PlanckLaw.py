# PlanckLaw
"""
Created on Wed Feb 16 15:37:16 2022

@author: masenpitts
"""

import numpy as np
import matplotlib.pyplot as plt

# Approximates the derivative at a point in a given set of data
# using five surrounding points or "centered difference" approximation
def fivePoint(x, f, i):
    dx = x[i+1] - x[i]
    return (f[i-2] - 8*f[i-1] + 8*f[i+1] - f[i+2])/(12*dx)

# Returns spectral radiance in units of J/(m^2*s)
# - wavelen: Wavelength in m
# - temp: Temperature in K
def planckSpec(wavlen, temp):
    h = 6.63e-34 # Planck's Constant in J*s
    c = 3.00e8 # Speed of light in m/s
    k = 1.381e-23 # Boltzmann's Constant in J/K
    
    # This calculation will be performed thousands of times so some helper
    # variables have been created to increase efficiency.
    hc = h*c
    l5 = np.power(wavlen, 5)
    
    numerator = 2*hc*c
    denominator = np.multiply(l5, np.exp(hc/(wavlen*k*temp)) - 1)
    
    return np.divide(numerator, denominator)

# Finds the wavelength at which the spectral radiance is maximum over a grid
# of wavelength values and a temperature value using the properties of first
# order derivatives.
# - wavelen: Wavelength in m
# - temp: Temperature in K
def findMax(wavlen, temp):
    # Calculates corresponding grid of spectral radiance values
    spectralRad = planckSpec(wavlen, temp)
    
    # A five point approximation is used to calculate the derivative,
    # so the first two points and last two points in the wavlen grid
    # cannot be used.
    derivRange = len(wavlen) - 4
    
    # Store the values of the maximum spectral radiance and wavelength at
    # which said maximum occurs respectively
    radMax = 0
    maxWavLen = wavlen[0]
    
    # Searches for a local maximum, which occurs when dI/dlamda = 0, and
    # determines the maximum spectral radiance and corresponding wavelength
    # out of those maxima.
    for l in range(derivRange):
        j = l + 2
        deriv = fivePoint(wavlen, spectralRad, j)
        # It is unlikely that our approximate derivative will ever be exactly
        # zero so we simply look for very small values.
        if abs(deriv) < 0.00001:
            if spectralRad[j] > radMax:
                radMax = spectralRad[j]
                maxWavLen = wavlen[j]
    return maxWavLen

lambdaMin = 0.01 # Maximum and minimum wavelengths of analysis range
lambdaMax = 20 # in micrometers
dlambda = 0.01 # Interval between wavelength grid points in micrometers

# 1D Grid of wavelength values over which the spectral radiance will be calculated
lambdaGrid = np.arange(lambdaMin, lambdaMax, dlambda)

tempMin = 150  # Maximum and minimum temperature values of analysis range
tempMax = 3000 # in kelvin
dT = 1         # Interval between temperature grid points in kelvin

# Scales the wavelength grid to be in micrometers
waveLenRange = np.multiply(lambdaGrid, 1e-6)

# 1D Grid of temeprature values over which we will perform our analysis
tempGrid = np.arange(tempMin, tempMax, dT)

# List storing the wavelength at which max spectral radiance occurs for each
# temperature value
maxLambdas = []

# Finds the wavelength at which max spectral radiance occurs for each
# point in tempGrid
for T in tempGrid:
    maxLambdas.append(findMax(waveLenRange, T))

# Plots the wavelength at which max spectral radiance occurs as a function
# of temperature
plt.plot(tempGrid, np.divide(maxLambdas, 1e-6), color = "black")
plt.xlabel("Temperature (K)")
plt.ylabel("Wavelength (micrometers)")
plt.title("Wavelength at Max Spectral Radiance vs. Temperature",
          fontsize=13)
    