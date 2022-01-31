# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 07:55:20 2022

@author: masenpitts
"""

from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

%matplotlib qt 

# The following code was taken from an example animation code page on the 
# official MatPlotLib website
# URL: https://matplotlib.org/2.0.2/examples/animation/double_pendulum_animated.html
# Everything within the lines of pounds was taken from this website. No additional
# comments were added. Everything after this block is either transformed from the
# borrowed code or originally created.
###########################################################################

"""
===========================
The double pendulum problem
===========================

This animation illustrates the double pendulum problem.
"""

# Double pendulum formula translated from the C code at
# http://www.physics.usyd.edu.au/~wheat/dpend_html/solve_dpend.c

G = 9.8  # acceleration due to gravity, in m/s^2
L1 = 1.0  # length of pendulum 1 in m
L2 = 1.0  # length of pendulum 2 in m
M1 = 1.0  # mass of pendulum 1 in kg
M2 = 1.0  # mass of pendulum 2 in kg


def derivs(state, t):

    dydx = np.zeros_like(state)
    dydx[0] = state[1]

    del_ = state[2] - state[0]
    den1 = (M1 + M2)*L1 - M2*L1*cos(del_)*cos(del_)
    dydx[1] = (M2*L1*state[1]*state[1]*sin(del_)*cos(del_) +
               M2*G*sin(state[2])*cos(del_) +
               M2*L2*state[3]*state[3]*sin(del_) -
               (M1 + M2)*G*sin(state[0]))/den1

    dydx[2] = state[3]

    den2 = (L2/L1)*den1
    dydx[3] = (-M2*L2*state[3]*state[3]*sin(del_)*cos(del_) +
               (M1 + M2)*G*sin(state[0])*cos(del_) -
               (M1 + M2)*L1*state[1]*state[1]*sin(del_) -
               (M1 + M2)*G*sin(state[2]))/den2

    return dydx

# create a time array from 0..100 sampled at 0.05 second steps
dt = 0.05
t = np.arange(0.0, 20, dt)

# th1 and th2 are the initial angles (degrees)
# w10 and w20 are the initial angular velocities (degrees per second)
th1 = 120.0
w1 = 0.0
th2 = -10.0
w2 = 0.0

# initial state
state = np.radians([th1, w1, th2, w2])

# integrate your ODE using scipy.integrate.
y = integrate.odeint(derivs, state, t)

x1 = L1*sin(y[:, 0])
y1 = -L1*cos(y[:, 0])

x2 = L2*sin(y[:, 2]) + x1
y2 = -L2*cos(y[:, 2]) + y1

###########################################################################

time_template = "time = %.1fs"  # Text template for displaying the time in
                                # seconds on both plots during the animations

fig1 = plt.figure("Pendulum 1") # Creates the figure with the animation for the
                                # upper pendulum

oMax1 = max(y[:,1]) # Determines the limits of the vertical axis for the upper
                    # pendulum's Poincare Plot

# This loop finds values of angular position outside of the graph limits and
# their values by 2Pi to fit on the graph. This has no physical affect on the
# system since the sinusoidal functions used to solve the differential equation
# were 2Pi periodic.
# My reasoning was that having the edges of the plots at -Pi and Pi would allow
# us to determine more easily when each pendulum went over the top.
for theta in range(len(y[:, 0])):
    if y[theta, 0] > np.pi:
        y[theta, 0] -= 2*np.pi
    if y[theta, 0] < -np.pi:
        y[theta, 0] += 2*np.pi

# Subplot used to format the pendulum 1 graph
ax1 = fig1.add_subplot(111, autoscale_on=False, xlim=(-np.pi, np.pi), 
                       ylim=(-oMax1, oMax1))
ax1.grid() # Creates a grid for the pendulum 1 graph

# Formatting for the Pendulum 1 graph
plt.title("Poincare Diagram of Upper Pendulum")
plt.xlabel("Angular Position (radians)")
plt.ylabel("Angular Velcocity (rads/s)")

# A text object and data list continually updated by the animation functions
# that plots the data points on the Poincare Diagram and displays the current
# time in seconds
time_text1 = ax1.text(0.05, 0.9, "", transform=ax1.transAxes)
pend1, = ax1.plot([], [], "o")

# Initialization function for the Pendulum 1 animation
def init1():
    # Switches current figure to the Pendulum 1 graph
    plt.figure("Pendulum 1")
    
    # Initializes graph elements with default settings
    pend1.set_data([], [])
    time_text1.set_text("")
    
    return pend1, time_text1


def animate1(i):
    plt.figure("Pendulum 1")
    
    # Updates Pendulum 1 graph with all data points from the start of the
    # animation up until the current animation step.
    pend1.set_data(y[:i+1, 0], y[:i+1, 1])
    
    # Updates the text displaying the animation time
    time_text1.set_text(time_template % (i*dt))
    
    return pend1, time_text1

# Animation function for the Upper Pendulum. Runs simultaneously with other
# animation
ani1 = animation.FuncAnimation(fig1, animate1, np.arange(1, len(y)),
                              interval=25, blit=True, init_func=init1)

# The variables and functions below are analogous to those defined for
# Pendulum 1, but are used to animate the Poincare Diagram for the Lower
# Pendulum.
fig2 = plt.figure("Pendulum 2")

oMax2 = max(y[:, 3])

for theta in range(len(y[:, 2])):
    if y[theta, 2] > np.pi:
        y[theta, 2] -= 2*np.pi
    if y[theta, 2] < -np.pi:
        y[theta, 2] += 2*np.pi

ax2 = fig2.add_subplot(111, xlim=(-np.pi, np.pi), 
                       ylim=(-oMax2, oMax2))
ax2.grid()

plt.title("Poincare Diagram of Lower Pendulum")
plt.xlabel("Angular Position (radians)")
plt.ylabel("Angular Velcocity (rads/s)")

pend2, = ax2.plot([], [], "o", color="red")
time_text2 = ax2.text(0.05, 0.9, '', transform=ax2.transAxes)

def init2():
    plt.figure("Pendulum 2")
    pend2.set_data([], [])
    
    time_text2.set_text("")
    
    return pend2, time_text2


def animate2(i):
    plt.figure("Pendulum 2")
    pend2.set_data(y[:i+1, 2], y[:i+1, 3])
    
    time_text2.set_text(time_template % (i*dt))
    
    return pend2, time_text2

ani2 = animation.FuncAnimation(fig2, animate2, np.arange(1, len(y)),
                              interval=25, blit=True, init_func=init2)

plt.show()
