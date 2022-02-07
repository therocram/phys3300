# FunctionInterpsL
"""
Created on Fri Feb  4 20:26:44 2022

@author: masenpitts
"""

import matplotlib.pyplot as plt
import LagrangeInterp as LI

funcSelect, xdata, fdata, x_value = LI.runConsole()

n = int(input("Enter order of Lagrangian Interpolation: "))

f_linInterp = LI.interp(xdata, fdata, n, x_value)

plt.plot(xdata, fdata, "bo")

LI.printResults(f_linInterp, x_value, funcSelect)