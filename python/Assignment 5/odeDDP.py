# odeDDP
# Solving the damped driven pendulum (DDP) using ODE solver methods
"""
Created on Tue Mar 15 13:15:37 2022

@author: masenpitts
"""

import numpy as np
import odeSolves as os
import matplotlib.pyplot as plt

# Simulation Parameters
q = 0.5
b = 1.15
omega0 = 2/3.

# Generalized velocity vector components
def g0(y, t):
    return y[1]

def g1(y, t):
    return b*np.cos(omega0*t) - np.sin(y[0]) - q*y[1]

# Dynamic Variable Vector w/ initial conditions
y = [0, 2]

# Generalized velocity vector
g = [g0, g1]

N = 10000 # Total number of time steps in simulation
n = 1000  # Total number of data points plotted
nt = int(N/n) # Determines how frequently points are plotted

dt = 0.1 # Time step of simulation

# Initializes solution list
solution = np.zeros((n, len(y)))
solution[0] = y

# Counting variable that allows the program to solve the problem
# at N time steps but only plot n out of those points.
c = 1

for j in range(1, N):
    y = os.solver(y, g, dt, "2", t = j*dt)
    
    # Shifts an angle by 2pi if it would go out of the plot bounds
    if y[0] > np.pi:
        y[0] -= 2*np.pi
    elif y[0] < -np.pi:
        y[0] += 2*np.pi
    
    if c == nt:
        solution[int(j/c)] = y
        c = 0
    
    c += 1

# Plots n points of solution
plt.plot(solution[:, 0], solution[:, 1], ".", color="black")
plt.xlim(-np.pi, np.pi)
plt.ylim(-3, 3)
plt.axhline(0, color="black")
plt.axvline(0, color="black")
plt.xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
plt.yticks(np.arange(-3, 4, 1))
plt.title("Phase Diagram of DDP with q = " + str(q) + ", b = " + str(b) 
          + ", w0 = {:.2f}".format(omega0))
plt.xlabel("Theta")
plt.ylabel("Omega")
plt.show()