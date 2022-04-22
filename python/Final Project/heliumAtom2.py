# heliumAtom2
# Implements the methods and functions of fakeTime.py to 
# solve for the ground state and first two excited states
# of the electrons in a simplified helium atom model
#
# Created by Masen Pitts on 4/20/2022
# Last updated on 4/21/2022

#**************************************************#
# The parameters that the user will most commonly modify are:
#
# - rmin: minimum range of the configuration space region. this
#   should always be kept at zero.
#
# - rmax: maximum range of the configuration space region in Bohr radii.
#   If you only care about the ground state this can be set as low as 6
#   or 7 (be aware that you will get very inaccurate estimates for the
#   excited states if you do this). A range of 16-18 seems to be ok for
#   the excited states. Pushing it much further will start giving innaccurate
#   ground state values if meshRes is too small.
#
# - meshRes: this gives the total number of grid points in both dimensions
#   since the configuration space region is square. The total number of
#   grid point is the SQUARE of this value (e.g. a meshRes of 100 corresponds)
#   to 10,000 total grid points. Be aware that significantly increasing
#   meshRes will likely make the algorithm unstable if dtau is not lowered
#   by an appropriate amount.
#
# - boundaryValue: this sets the boundary condition for the system. This
#   should always be kept at zero for this system.
#
# - dtau: this determines the fake time step for the algorithm. This will
#   almost always be much less than 0.01.
#
# - N: this determines the total number of algorithm steps. The larger this
#   value is the more accurate the energy eignvalues will be and the longer
#   the program will run.
#
# The V should be kept the same since changing it will change the system 
# that you are working with.
#
# It is highly recommended that grdGuess, excGuess1, and excGuess2 be kept
# the same, though feel free to try out other guess functions if you think
# they will do better.
#
# The Matplotlib 3D surface plots done at the bottom of the file are
# completely optional. Feel free to visualize and analyze the data
# however you want to. 
# The following variables returned by the algorithm contain relevant data:
#
# - grdState, excited1, excited2: numpy arrays that contain the
#   solution values for the ground state and first excited state
#   reduced wavefunctions.
#
# - R1, R2: 2D these are the configuration space grids that should be
#   used for plotting the wavefunctions (see the surface plots below as 
#   an example).
#
# - energyList: this contains all of the energy eigenvalues sorted from
#   lowest to highest.
#
# See Section 3 of the report for a detailed analysis of the results of
# this program's use.
#**************************************************#

import fakeTime as ft
import numpy as np
import matplotlib.pyplot as plt

def V(R1, R2):
    return -2/R1 -2/R2 + 1/np.fmax(R1,R2)

def grdGuess(R1, R2):
    return R1*np.exp(-R1)*R2*np.exp(-R2)

def excGuess1(R1, R2):
    # Antisymmetric Guess
    psi12 = R1*np.exp(-R1)*(1-(1/2)*R2)*R2*np.exp(-R2/2)
    psi21 = R2*np.exp(-R2)*(1-(1/2)*R1)*R1*np.exp(-R1/2)
    
    return psi12 - psi21

def excGuess2(R1, R2):
    # Symmetric Guess
    psi12 = R1*np.exp(-R1)*(1-(1/2)*R2)*R2*np.exp(-R2/2)
    psi21 = R2*np.exp(-R2)*(1-(1/2)*R1)*R1*np.exp(-R1/2)
    
    return psi12 + psi21
    
    

rmin = 0
rmax = 17

meshRes = 400
boundaryValue = 0
dtau = 0.0008
N = 6000

grdState, excited1, excited2, R1, R2, energyList = ft.quantumSolver2D(rmin, 
                                                rmax, rmin, rmax, meshRes, 
                                                V, grdGuess, excGuess1, 
                                                excGuess2, boundaryValue, 
                                                dtau, N)

# Next Steps:
#   -Improve accuracy w/ simpson rule, reduced radial -> radial
#   -Use density plots instead of 3D surface plots
#   -Take a closer look at the reduced radial wavefunction vs.
#    radial wavefunction thing (done)
#   -Experiment with asymmetric guess functions. How long do we need
#    to let them run before they reproduce symmetry/antip-symmetry?

nonBound = (slice(1, R1.shape[0]-1), slice(1, R1.shape[1]-1))

grdState[nonBound] /= (R1[nonBound]*R2[nonBound])
excited1[nonBound] /= (R1[nonBound]*R2[nonBound])
excited2[nonBound] /= (R1[nonBound]*R2[nonBound])

ft.normalize(grdState, R1, R2)
ft.normalize(excited1, R1, R2)
ft.normalize(excited2, R1, R2)

grd = grdState[nonBound]
exc1 = excited1[nonBound]
exc2 = excited2[nonBound]

r1 = R1[nonBound]
r2 = R2[nonBound]

grdPlot = plt.figure()
axG = grdPlot.add_subplot(projection="3d")
plt.title("Ground State Wavefunction")
axG.set_xlabel(r"$r_1$")
axG.set_ylabel(r"$r_2$")
axG.set_zlabel(r"$\psi (r_1,r_2)$")
axG.plot_surface(r1, r2, grd, cmap=plt.cm.plasma)
axG.set_xlim3d(rmin, rmax)
axG.set_ylim3d(rmin, rmax)

exc1Plot = plt.figure()
ax1 = exc1Plot.add_subplot(projection="3d")
plt.title("First Excited State Wavefunction")
ax1.set_xlabel(r"$r_1$")
ax1.set_ylabel(r"$r_2$")
ax1.set_zlabel(r"$\psi (r_1,r_2)$")
ax1.plot_surface(r1, r2, exc1, cmap=plt.cm.plasma)
ax1.set_xlim3d(rmin, rmax)
ax1.set_ylim3d(rmin, rmax)

exc2Plot = plt.figure()
ax2 = exc2Plot.add_subplot(projection="3d")
plt.title("Second Excited State Wavefunction")
ax2.set_xlabel(r"$r_1$")
ax2.set_ylabel(r"$r_2$")
ax2.set_zlabel(r"$\psi (r_1,r_2)$")
ax2.plot_surface(r1, r2, exc2, cmap=plt.cm.plasma)
ax2.set_xlim3d(rmin, rmax)
ax2.set_ylim3d(rmin, rmax)
plt.show()