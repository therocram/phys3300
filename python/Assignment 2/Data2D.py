# Data2D
# Assignment 2, Exercise 2
# Creates a bar plot of 2D time series data read from a file
# using pandas and matplotlib
"""
Created on Fri Jan 21 12:32:21 2022

@author: masenpitts
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Reads data from csv file
oilData = pd.read_csv("C:/Users/masenpitts/Downloads/Imports Crude Oil.csv")

# Array storing error bar values (assuming an uncertainty of 5%)
error = np.multiply(oilData["oil"], 0.05)

# Create bar plot of data
plt.bar(oilData["month"], oilData["oil"], yerr=error, color="green")
# Configure figure settings
plt.title("Monthly Imports of all Grades of Crude Oil from World to U.S")
plt.xlabel("Month (Year Month - XXXX XX)")
plt.ylabel("Crude Oil Imported (thousands of barrels)")
plt.show()