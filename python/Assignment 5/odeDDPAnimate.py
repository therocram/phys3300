# odeDDPAnimate
# Solving the damped driven pendulum (DDP) using ODE solver methods
# and animates the solution plot
"""
Created on Sat Mar 19 11:53:59 2022

@author: masenpitts
"""

import numpy as np
import odeSolves as os
import matplotlib.pyplot as plt
import matplotlib.animation as animation

%matplotlib qt 

# Simulation Parameters
q = 0.5
b = 0.9
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

# Animate a plot of n points of the solution

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim = (-np.pi, np.pi), ylim = (-3, 3),
                     xticks = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi], yticks = np.arange(-3, 4, 1))
ax.axhline(0, color="black")
ax.axvline(0, color="black")
ax.set_title("Phase Diagram of DDP with q = " + str(q) + ", b = " + str(b) 
          + ", w0 = {:.2f}".format(omega0))
ax.set_xlabel("Theta")
ax.set_ylabel("Omega")

dataPlot, = ax.plot([], [], ".", color="black")
time = "time = {:.1f}s"
timeDisplay = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def start():
    dataPlot.set_data([], [])
    timeDisplay.set_text('')
    return dataPlot, timeDisplay

def animate(i):
    dataPlot.set_data(solution[:i+1, 0], solution[:i+1, 1])
    timeDisplay.set_text(time.format(i*nt*dt))
    
    return dataPlot, timeDisplay

ani = animation.FuncAnimation(fig, animate, np.arange(1, n), interval = 5,
                              blit=True, init_func=start)

plt.show()