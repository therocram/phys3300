# IntegralApprox
"""
Created on Wed Feb 16 14:38:58 2022

@author: masenpitts
"""

import numpy as np
import matplotlib.pyplot as plt

def trapRule(x, f):
    integral = 0
    n = len(f) - 1
    for i in range(n):
        dx = x[i+1] - x[i]
        integral += (dx/2.0)*(f[i] + f[i+1])
    return integral

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
    
    
xMin = 0
xMax = np.pi
dx = 0.1

actual = -(np.cos(xMax) - np.cos(xMin))

x = np.arange(xMin, xMax, dx)
f = np.sin(x)

plt.plot(x, f, label="f(x)")
plt.legend()
plt.show()

result1 = trapRule(x, f)
result2 = simpsonRule(x, f)

error1 = 100*abs(actual - result1)/actual
error2 = 100*abs(actual - result2)/actual

text = "{name:s} Rule = {approx:.6f} error = {error:.6f}%"

print(text.format(name = "Trapezoid", approx = result1, error = error1))
print(text.format(name = "Simpson", approx = result2, error = error2))

