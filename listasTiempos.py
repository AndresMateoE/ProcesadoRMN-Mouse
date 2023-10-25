# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 15:21:04 2023

@author: Kea500
"""
import numpy as np
 
Tmin = 0.1
Tmax = 3000
N = 16

amax = np.log10(Tmax)
amin = np.log10(Tmin)
astep = (amax-amin)/(N-1)
t1Axis = []

for i in range(N):
    t1Axis.append(10**(amin+astep*(i)))

print(t1Axis)