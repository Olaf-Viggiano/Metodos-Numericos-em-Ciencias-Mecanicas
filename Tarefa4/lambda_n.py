#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 17:31:30 2025

@author: olafviggiano
"""

import numpy as np

x=np.arange(0,50,0.0001)
Bi=200*0.005/2.8

y=np.abs(Bi/np.tan(x)-x)

from scipy.signal import find_peaks
import matplotlib.pyplot as plt

# Find minima by inverting the signal
minima_indices, _ = find_peaks(-y)

print("Indices of minima:", minima_indices)
print("Values of minima:", y[minima_indices])

# Plot to visualize
plt.plot(y, label='Signal')
plt.plot(minima_indices, y[minima_indices], 'ro', label='Minima')
plt.ylim(-0.5,10)
plt.legend()
plt.show()

for i in minima_indices:
    print(x[i])
    