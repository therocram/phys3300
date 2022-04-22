# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 20:59:03 2022

@author: masenpitts
"""

import numpy as np

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
    intList = []
    
    for y in range(Y.shape[0]):
        intList.append(simpsonRule(f2D[y], X[y]))
        
    return simpsonRule(intList, Y[:,0])

def integral2D(f2D, X, Y):
    d = X[0,1] - X[0,0]
    temp = f2D*d*d
    return np.sum(temp)

mesh = 100
r1 = np.linspace(0,1,mesh)
r2 = np.linspace(0,1,2*mesh)

R1, R2 = np.meshgrid(r1, r2)

func2D = ((np.sin(2*R1*np.pi))**2)*((np.sin(2*R2*np.pi))**2)

print(integral2D(func2D, R1, R2))
print(simpson2D(func2D, R1, R2))

print(R1.shape, R2.shape)

