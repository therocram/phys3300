# InterpolationTemplate
"""
Created on Fri Jan 28 14:49:53 2022

@author: masenpitts
"""

import matplotlib.pyplot as plt
import numpy as np

def getData():
    f = [2, 3, 19, 122, 370, 416, 2]
    x = [1, 2, 2.7, 3.05, 3.7, 5, 6]
    return f, x

def getInput():
    x_value = input("Enter position of target value: ")
    return float(x_value)
    
def getIndex(x, x_value):
    for i in range(len(x)):
        if x[i] > x_value:
            return i-1, i
        
def interp(x, f, i, x_value):
    result = f[i] + (x_value - x[i])/(x[i+1] - x[i])*(f[i+1] - f[i])
    return result

f, x = getData()

x_value = getInput()
i, i1 = getIndex(x, x_value)
result = interp(x, f, i, x_value)

x_grid = np.arange(1,6,0.01)
f_interp = []

for each in x_grid:
    index, index1 = getIndex(x, each)
    f_interp.append(interp(x, f, index, each))

print(result)

plt.plot(x, f, "bo")
plt.plot(x_grid, f_interp)
plt.show()