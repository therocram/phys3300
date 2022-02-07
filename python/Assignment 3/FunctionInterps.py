# FunctionInterps
"""
Created on Wed Feb  2 14:37:49 2022

@author: masenpitts
"""

import matplotlib.pyplot as plt
import LagrangeInterp as LI

funcSelect, xdata, fdata, x_value = LI.runConsole()

iBefore = LI.getIndex(xdata, x_value)
f_linInterp = LI.linearInterp(xdata, fdata, iBefore, x_value)

plt.plot(xdata, fdata, "bo")

LI.printResults(f_linInterp, x_value, funcSelect)

'''
min = 0
max = 6
n = 4
interval = 0.1

f, x = getSin(min, max, interval)

x_grid = np.arange(min, max, 0.01)
f_interp = []

x_value = float(input("Enter X-value for Error Approximation: "))
f_value = LI.interp(x, f, n, x_value)

for each in x_grid:
    f_interp.append(LI.interp(x, f, n, each))
    
print(percentDiff(f_value, myFunc(x_value))) 
    
plt.plot(x, f, "bo")
plt.plot(x_grid, f_interp, "r")
plt.plot(x_value, f_value, "go")
plt.show()
'''