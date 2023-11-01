# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 16:25:23 2023

@author: Lenovo
"""

# Mapa T1-D Mouse
import core_IO_Andres as IO
import core_plot_am as graph
import pandas as pd
import numpy as np


#Carpeta con las mediciones

pwd_hept = ('G:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Heptano/')
pwd_dode = ('G:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Dodecano/')
pwd = ('G:/Unidades compartidas/TF-Andres/Reshapes/Dodecano_Heptano/32/')
#Leemos los archivos necesarios:
param = IO.read_acq(pwd_dode)
Z_dode, T1axis , Daxis =  IO.read_T1D(pwd_dode)
Z_hept, T1axis , Daxis =  IO.read_T1D(pwd_hept)

Z = Z_dode + Z_hept

# =============================================================================
# posicionesT1 = [0, 4, 9, 13, 18, 22, 27, 31]
# posicionesD = [0, 1, 2, 4, 7, 11, 19, 31]
# 
# Z_reshape = np.array(Z[posicionesT1][:, posicionesD])
# Daxis = np.array([Daxis[i] for i in posicionesD])
# T1axis = np.array([T1axis[i] for i in posicionesT1])
# 
# Z = Z_reshape
# =============================================================================

nT1, nD = len(T1axis), len(Daxis)


alpha = 0.01
T1min, T1max = 1/2, 4
Dmin, Dmax = -1, 1


T1xx = 10
Dxx = 10

tEcho = param['echoTime']

S0, T1, D, K1, K2 = IO.initKernelT1D(nT1, nD, T1axis, Daxis, 
                                     T1min, T1max, Dmin, Dmax)     
print(f'Starting NLI: Alpha = {alpha}.')
S, iter = IO.NLI_FISTA_2D(K1, K2, Z, alpha, S0)
if iter < 100000:
    print('Inversion ready!')
else:
    print('Warning!')
    print('Maximum number of iterations reached!')
    print('Try modifying T2Range and/or alpha settings.')


#Guardamos los datos procesados

Transformada = pd.DataFrame(S)
VectorT1 = pd.DataFrame(T1)
VectorD = pd.DataFrame(D)
Transformada.to_csv(pwd+"Transformada.txt", index=False, header=False)
VectorT1.to_csv(pwd+"VectorT1.txt", index=False, header=False)
VectorD.to_csv(pwd+"VectorD.txt", index=False, header=False)


#Graficos

graph.MapaT1D(T1axis, Daxis, Z, T1, D, S, pwd, alpha, T1min, T1max, Dmin, Dmax, T1xx, Dxx)
graph.PlotT1D_D(Daxis, D, S, pwd, alpha, Dmin, Dmax)
graph.PlotT1D_T1(T1axis, T1, S, pwd, alpha, T1min, T1max)
#graph.PlotT1D_T1_NoNorm(T1axis, T1, S, pwd, alpha, T1min, T1max)
#graph.PlotT1D_D_NoNorm(Daxis, D, S, pwd, alpha, Dmin, Dmax)
# CODIGO CHEVA
