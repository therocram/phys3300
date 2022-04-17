# -*- coding: utf-8 -*-
"""
Created on Sat Apr  9 21:20:20 2022

@author: masenpitts
"""

import numpy as np

def add(X, Y):
    return X + Y

def multiply(X, Y):
    return X*Y

mesh = 10
r1 = np.linspace(0,1,mesh)
r2 = np.linspace(0,1,2*mesh)

R1, R2 = np.meshgrid(r1, r2)
meshX = len(r1)
meshY = len(r2)

R1 = np.ones(R1.shape)

nonBound = (slice(1, meshY-1), slice(1, meshX-1))
left = (slice(0, meshY-2), slice(1, meshX-1))
right = (slice(2, meshY), slice(1, meshX-1))
up = (slice(1, meshY-1), slice(2, meshX))
down = (slice(1, meshY-1), slice(0, meshX-2))

a = (1/4)*(R1[left]+R1[right]+R1[up]+R1[down])

print(a, R1)