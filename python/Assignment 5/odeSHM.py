# odeSHM
# Solving the simple harmonic oscillator (SHM) using
# several ODE solver methods
"""
Created on Sun Mar 13 22:23:22 2022

@author: masenpitts
"""

import numpy as np
import odeSolves as os
import matplotlib.pyplot as plt

# Generalized velocity vector components
def g0(y):
    return y[1]

def g1(y):
    return -y[0]

# Calculates the total mechanical energy in the system
def mechE(y):
    return 0.5*y[1]*y[1] - 0.5*y[0]*y[0]

# Initial conditions
yinit = [1, 1]

# Dynamic Variable Vector
y = yinit[:]

# Generalized velocity vector
g = [g0, g1]

# Allows the user to determine the time step size
# using an interactive console.
print("Welcome to the SHM Solver 1.4\n") 
dt = float(input("Enter time step: "))

# Time Step and time array
start = 0
stop = 50
timeList = np.arange(start, stop, dt)
N = len(timeList)

# Stores solution data for each solver method
eulerSol = np.zeros((N, len(y)))
rk4Sol = np.zeros((N, len(y)))
leapFrogSol = np.zeros((N, len(y)))

# Stores the values of total mechanical energy at each time step
# for each method
#eulerE = np.zeros(N)
#rk4E = np.zeros(N)
#leapFrogE = np.zeros(N)

# Solves the system for the length of time provided using each 
# solver method by running Leap Frog method once
# and by running every other solver method once at each time step

os.solver(y, g, dt, "3", leapFrogSol)

y = yinit[:] # Reset system to initial conditions
eulerSol[0] = y
for t in range(1, N):
    y = os.solver(y, g, dt, "1")
    eulerSol[t] = y
    
y = yinit[:]
rk4Sol[0] = y
for t in range(1, N):
    y = os.solver(y, g, dt, "2")
    rk4Sol[t] = y

# Use to calculate and plot the total mechanical energy in the system
# over time for each solver method. Allows one to measure the stability 
# of each solver method.
#for state in range(N):
#    eulerE[state] = mechE(eulerSol[state])
#    rk4E[state] = mechE(rk4Sol[state])
#    leapFrogE[state] = mechE(leapFrogSol[state])


plt.rcParams["figure.figsize"] = (5, 10)

# Plots phase diagrams (position vs. velocity in this case) of
# each solution method for comparison 
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = True)
fig.suptitle("Phase Diagrams for Simple\n Harmonic Oscillator with dt = " + str(dt), fontsize=15)

ax1.plot(eulerSol[:, 1], eulerSol[:, 0], "o", color="blue")
ax1.set_title("Euler")
ax1.set_ylabel("Position")
ax1.set_xlabel("Velocity")
#ax1.set_aspect("equal")
ax1.grid()

ax2.plot(rk4Sol[:, 1], rk4Sol[:, 0], "o", color="red")
ax2.set_title("Runge-Kutta 4th Order")
ax2.set_ylabel("Position")
ax2.set_xlabel("Velocity")
#ax2.set_aspect("equal")
ax2.grid()

ax3.plot(leapFrogSol[:, 1], leapFrogSol[:, 0], "o", color="green")
ax3.set_title("Leap Frog")
ax3.set_ylabel("Position")
ax3.set_xlabel("Velocity")
#ax3.set_aspect("equal")
ax3.grid()

plt.rcParams.update(plt.rcParamsDefault)

# Plots postion and veclocity of system over time for
# each solution method for comparison 
fig2, (bx1, bx2) = plt.subplots(2, sharex = True, sharey = True)
fig2.suptitle("Simple Harmonic Oscillator with dt = " + str(dt), fontsize=15)

bx1.plot(timeList, eulerSol[:, 0], label="Euler", color="blue")
bx1.plot(timeList, rk4Sol[:, 0], label="RK4", color="red")
bx1.plot(timeList, leapFrogSol[:, 0], label="Leap Frog", color="green")
bx1.legend()
bx1.set_ylim(-1.5, 1.5)
bx1.set_xlabel("Time")
bx1.set_ylabel("Position")

bx2.plot(timeList, eulerSol[:, 1], label="Euler", color="blue")
bx2.plot(timeList, rk4Sol[:, 1], label="RK4", color="red")
bx2.plot(timeList, leapFrogSol[:, 1], label="Leap Frog", color="green")
bx2.legend()
bx2.set_ylim(-1.5, 1.5)
bx2.set_xlabel("Time")
bx2.set_ylabel("Velocity")

# Used to plot energy in solution over time for each solution method
#energyFig = plt.figure()
#plt.plot(timeList, eulerE, color="blue")
#plt.plot(timeList, rk4E, color="red")
#plt.plot(timeList, leapFrogE, color="green")

plt.show()