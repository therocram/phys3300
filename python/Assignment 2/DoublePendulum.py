# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 07:55:20 2022

@author: masenpitts
"""

"""
===========================
The double pendulum problem
===========================

This animation illustrates the double pendulum problem.
"""

# Double pendulum formula translated from the C code at
# http://www.physics.usyd.edu.au/~wheat/dpend_html/solve_dpend.c

from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

%matplotlib qt 

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


'''
fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
'''

fig1 = plt.figure()
time_template = "time = %.1fs"

oMax = max(y[:,3])

for omega in range(len(y[:, 2])):
    if y[omega, 2] > np.pi:
        y[omega, 2] -= 2*np.pi
    if y[omega, 2] < -np.pi:
        y[omega, 2] += 2*np.pi

ax1 = fig1.add_subplot(111, autoscale_on=False, xlim=(-np.pi, np.pi), 
                       ylim=(-2*np.pi, 2*np.pi))
ax1.grid()
time_text1 = ax1.text(0.05, 0.9, "", transform=ax1.transAxes)
pend1, = ax1.plot([], [], "o-")

ax2 = fig1.add_subplot(111, xlim=(-2*np.pi, 2*np.pi), 
                       ylim=(-2*np.pi, 2*np.pi))
ax2.grid()
pend2, = ax2.plot([], [], "o-", color="red")
time_text2 = ax2.text(0.05, 0.9, '', transform=ax2.transAxes)

def init():
    pend1.set_data([], [])
    pend2.set_data([], [])
    
    time_text1.set_text("")
    time_text2.set_text("")
    
    return pend1, pend2, time_text1, time_text2


def animate(i):
    pend1.set_data(y[:i+1, 0], y[:i+1, 1])
    pend2.set_data(y[:i+1, 2], y[:i+1, 3])
    
    time_text1.set_text(time_template % (i*dt))
    time_text2.set_text(time_template % (i*dt))
    
    return pend1, pend2, time_text1, time_text2

#ani = animation.FuncAnimation(fig, animate, np.arange(1, len(y)),
#                              interval=25, blit=True, init_func=init)

ani1 = animation.FuncAnimation(fig1, animate, np.arange(1, len(y)),
                              interval=25, blit=True, init_func=init)

#print(y[:,0], y[:,1])
#plt.plot(y[:,0], y[:,1])

# ani.save('double_pendulum.mp4', fps=15)
plt.show()
