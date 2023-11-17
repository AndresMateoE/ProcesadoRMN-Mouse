# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 15:21:04 2023

@author: Kea500
"""
import numpy as np
import core_IO_Andres as IO
import core_plot_am as graph
 
Tmin = 0.1
Tmax = 5000
N = 32

amax = np.log10(Tmax)
amin = np.log10(Tmin)
astep = (amax-amin)/(N-1)
t1Axis = []

for i in range(N):
    t1Axis.append(10**(amin+astep*(i)))

print(t1Axis)
print(np.sum(t1Axis))

# =============================================================================
# Lineal = np.around(np.linspace(Tmin,Tmax, N), 2)
# print(Lineal)
# =============================================================================

# =============================================================================
# pwd_dode = ('H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Dodecano/')
# Z_dode, Daxis , T1axis =  IO.read_T1D(pwd_dode)
# 
# 
# filas_a_mantener = [0, 4, 9, 13, 18, 22, 27, 31]  
# columnas_a_mantener = [0, 4, 9, 13, 18, 22, 27, 31]  
# 
# Z_reshape = Z_dode[filas_a_mantener][:, columnas_a_mantener]
# 
# print(Z_reshape)
# =============================================================================

