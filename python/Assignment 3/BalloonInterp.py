# BalloonInterp
"""
Created on Wed Feb  9 08:13:49 2022

@author: masenpitts
"""

import matplotlib.pyplot as plt
import numpy as np
import LagrangeInterp as LI

# Write filepath to data here
filePath = "C:/Users/masenpitts/Downloads/balloon_sounding.txt"

# Open file to be read by program
balloonFile = open(filePath, "r")

# Stores information relating pressure data to temperature data
# [[pressure], [temperature]]
pressureTemp = [[], []]

# Stores information relating pressure data to wind data
# [[pressure], [wind direction], [wind speed]]
pressureWind = [[], [], []]

# Variables used to store the maximum and minimum pressure values at which
# there exists valid information for temperature and wind data respectively.
# Used to dermine intervals on which interpolation will be performed
tPMin = 100000
tPMax = 0
wPMin = 100000
wPMax = 0

# These variables are used to make sure that at least one element gets added
# to the data arrays in the pressureTemp and pressureWind lists
countT = 0
countW = 0

# Extracts valid data from the given file and adds it to the appropriate arrays.
# Also keeps track of maximum and minimum pressure values for temperature and
# wind data.
for eachLine in balloonFile:
    # Ignore any lines that start with "#"
    if eachLine[0] != "#":
        # Finds the column elements and puts them into an array
        thisLine = eachLine.split()
        
        # Filters out invalid temperature data
        if thisLine[2] != "99999":
            P = float(thisLine[0])
            # Prevents the addition of identical data points to the data lists
            if countT == 0 or pressureTemp[0][len(pressureTemp[0]) - 1] != P:
                pressureTemp[0].append(P)
                pressureTemp[1].append(float(thisLine[2]))
                
                # Keeps track of min and max pressure values for temp data
                if P < tPMin:
                    tPMin = P
                if P > tPMax:
                    tPMax = P
          
            countT += 1
        
        # Filters out invalid wind data
        if thisLine[4] != "99999":
            P = float(thisLine[0])
            if countW == 0 or pressureWind[0][len(pressureWind[0]) - 1] != P:
                pressureWind[0].append(P)
                pressureWind[1].append(float(thisLine[4]))
                pressureWind[2].append(float(thisLine[5]))
                
                # Keeps track of min and max pressure values for wind data
                if P < wPMin:
                    wPMin = P
                if P > wPMax:
                    wPMax = P
           
            countW += 1

# Orders of Lagrangian Interpolations for temperature, wind direction, and
# wind speed respectively.
nT = 4
nWD = 2
nWS = 2

# In the text file the data entries are listed in descending order. The
# interpolation program assumes that the independent variable list is in 
# ascending order so the collected lists of data need to be reversed.
pressureTemp[0].reverse()
pressureTemp[1].reverse()
pressureWind[0].reverse()
pressureWind[1].reverse()
pressureWind[2].reverse()

# Arrays of pressure values at which interpolation will be performed and the
# results plotted.
PTgrid = np.arange(tPMin, tPMax, 1)
PWgrid = np.arange(wPMin, wPMax, 1)

# Stores information relating to the linear and Lagrangian interpolations
# of the collected data sets
# [[temperature], [wind direction], [wind speed]]
linInterps = [[], [], []]
lagranInterps = [[], [], []]

# Interpolates the temperature value at each point on the pressure-temperature
# grid
for P in PTgrid:
    linT = LI.interp(pressureTemp[0], pressureTemp[1], 1, P)
    lagranT = LI.interp(pressureTemp[0], pressureTemp[1], nT, P)
    
    linInterps[0].append(linT)
    lagranInterps[0].append(lagranT)
 
# Interpolates the values of wind direction and speed at each point on the 
# pressure-wind grid
for P in PWgrid:
    linWD = LI.interp(pressureWind[0], pressureWind[1], 1, P)
    lagranWD = LI.interp(pressureWind[0], pressureWind[1], nWD, P)
    
    linInterps[1].append(linWD)
    lagranInterps[1].append(lagranWD)
    
    linWS = LI.interp(pressureWind[0], pressureWind[2], 1, P)
    lagranWS = LI.interp(pressureWind[0], pressureWind[2], nWS, P)
    
    linInterps[2].append(linWS)
    lagranInterps[2].append(lagranWS)
    
# Creates a plot of pressure vs. temperature including the actual data as well as
# the linear and Lagrangian interpolation results.
PTfig = plt.figure("Pressure vs. Temperature")
plt.plot(pressureTemp[1], pressureTemp[0], "go", label="Data")
plt.plot(linInterps[0], PTgrid, color="blue", label="Linear")
plt.plot(lagranInterps[0], PTgrid, color="red", 
             label=("Lagrangian (Order " + str(nT) + ")"))
plt.title("Pressure vs. Temperature")
plt.xlabel("Temperature (1/10 degree Celsius)")
plt.ylabel("Pressure (1/10 mBar)")
plt.ylim(9000, 0)
plt.legend()

# Plot of pressure vs. wind direction with same factors plotted
PWDfig = plt.figure("Pressure vs. Wind Direction")
plt.plot(pressureWind[1], pressureWind[0], "go", label="Data")
plt.plot(linInterps[1], PWgrid, color="blue", label="Linear")
plt.plot(lagranInterps[1], PWgrid, color="red", 
             label=("Lagrangian (Order " + str(nWD) + ")"))
plt.title("Pressure vs. Wind Direction")
plt.xlabel("Wind Direction (degrees)")
plt.ylabel("Pressure (1/10 mBar)")
plt.ylim(9000, 0)
plt.legend()

# Plot of pressure vs. wind speed with same factors plotted
PWSfig = plt.figure("Pressure vs. Wind Speed")
plt.plot(pressureWind[2], pressureWind[0], "go", label="Data")
plt.plot(linInterps[2], PWgrid, color="blue", label="Linear")
plt.plot(lagranInterps[2], PWgrid, color="red", 
             label=("Lagrangian (Order " + str(nWS) + ")"))
plt.title("Pressure vs. Wind Speed")
plt.xlabel("Wind Speed (1/10 m/s)")
plt.ylabel("Pressure (1/10 mBar)")
plt.ylim(9000, 0)
plt.legend()

plt.show()