# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 14:30:13 2022

@author: masenpitts
"""

import numpy as np
import matplotlib.pyplot as plt

def g0(y):
    return y[1]

def g1(x):
    k = 1000
    m = 10
    return -(k/m)*y[0]

def g2(y):
    return -y[0]

g = [g0, g1, g2]
y = [0, 1, 1]
N = 1000
dt = 0.001
timeList = np.zeros(N)
yList = np.zeros((N, len(y)))

for t in range(N):
    timeList[t] = t
    for i in range(len(g)):
       y[i] += g[i](y)*dt
    yList[t] = y

plt.plot(timeList, yList[:,1])