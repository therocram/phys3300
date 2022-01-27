# SinglePendulum
# Assignment 2, Exercise 1
# Creates a 2D Line Plot of the Angular Position of simple pendulum
# vs. Time
"""
Created on Fri Jan 14 15:27:43 2022

@author: masenpitts
"""

import matplotlib.pyplot as plt
import numpy as np

phi = np.pi/2 # Initial angular position. Keeps track of last angular
              # position (rad)
omega = 0 # Angular velocity (rads/s)
g = 9.81 # Acceleration of gravity (m/s^2)
runtime = 30 # Duration of simulation (s)
dt = 0.01 # Time step (s)
steps = int(runtime/dt) # Total number of steps in simulation. Also determines
                        # the total number of points on the graph
                        
length = 3.4 # Length of pendulum (m)

angPosition = np.zeros(steps) # Array storing angular positions of pendulum at
                              # each simulation step
                              
# Implements Euler Algorithm for the calculated number of steps
for t in range(steps):
    angPosition[t] = phi
    alpha = -(g/length)*np.sin(phi) # Angular acceleration (rad/s^2)
    omega += alpha*dt
    phi += omega*dt

plt.plot(range(steps), angPosition, color="Red")
plt.title("Angular Position of Simple Pendulum from Vertical vs. Time")
plt.xlim(0, steps)
plt.xticks(range(0, steps+1, int(steps/10)), range(0, runtime+1, int(runtime/10)))
plt.xlabel("Time (s)")
plt.ylabel("Angular Position from Vertical (rads)")

plt.show()



