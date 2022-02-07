# Correlation2D
# Assignment 2, Exercise 2
"""
Created on Sat Jan 29 23:28:27 2022

@author: masenpitts
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

deathsData = pd.read_csv("C:/Users/masenpitts/Documents/drownCorr.csv")

drownData = deathsData["drown"]
gastricData = deathsData["gastric"]

drownErr = np.multiply(0.05, drownData)
gastricErr = np.multiply(0.05, gastricData)

#plt.plot(drownData, gastricData, "o")
plt.errorbar(drownData, gastricData, color="purple", yerr=gastricErr, xerr=drownErr, fmt="o",
             ecolor="green")
plt.title("Drownings caused by an accident involving a fishing boat\n vs. Deaths caused by inhalation of gastric contents")
plt.xlabel("Drownings")
plt.ylabel("Inhalation of Gastric Contents")
#Drownings caused by an accident involving a fishing boat
#Deaths caused by inhalation of gastric contents