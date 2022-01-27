# VisualizationPractice
# Assignment 2, Exercise 1
"""
Created on Fri Jan 14 14:26:18 2022

@author: Masen Pitts
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

data = pd.read_csv("C:/Users/masenpitts/Downloads/tips.csv")

display(data.head(20))
print("\n\n\n")

#plt.scatter(data['day'], data['tip'], c=data['size'], s=data['total_bill'])

#plt.plot(data['tip'], label='Tip')
#plt.plot(data['size'], label='Size')

#plt.bar(data['day'], data['tip'])

#plt.hist(data['total_bill'])

#sns.lineplot(x='sex', y='total_bill', data=data)

#sns.scatterplot(x='day', y='tip', data=data, hue='sex')

#sns.lineplot(x='day', y='tip', data=data)

sns.histplot(x='total_bill', data=data, kde=True, hue='sex')

plt.title('Bill Totals')
plt.xlabel('Days')
plt.ylabel('Tip')

#plt.colorbar()
#plt.legend()

print(data)

plt.show()

