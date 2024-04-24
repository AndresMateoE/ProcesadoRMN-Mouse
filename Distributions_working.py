# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 20:03:21 2024

@author: Lenovo
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.stats import moment
import core_plot_am as graph
import core_IO_Andres as IO

Dxx = 10

pwdDode = 'H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Dodecano/'
Dif_dode, S_dode, Dax_dode, T1ax_dode = IO.read_dif(pwdDode, Dxx)

pwdHept = 'H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Heptano/'
Dif_hept, S_hept, Dax_hept, T1ax_hept = IO.read_dif(pwdHept, Dxx)

pwdPet = 'H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Petroleo/'
Dif_pet, S_pet, Dax_pet, T1ax_pet = IO.read_dif(pwdPet, Dxx)

Dmin, Dmax, Dxx = -1, 1, 10

graph.plot_difusion(Dax_dode, Dif_dode, Dmin, Dmax, "Difusión Dodecano", "Red")
graph.plot_difusion(Dax_hept, Dif_hept, Dmin, Dmax, "Difusión Heptano", "Green")

# Son 3, 6, 15 y 30 nm

Deff_n5 = np.array([0.25, 2.02, 3.5, 4.0])
Deff_n7 = np.array([0.11, 1.05, 1.98, 2.16])
Deff_c6 = np.array([0.042, 0.51, 1.02, 1.05])

def D_ajuste(N,a,b):
    return a* N**(-b-0.7)

# =============================================================================
# params, covariance = curve_fit(D_ajuste, (5,7), (Deff_n5[0],Deff_n7[0]))
# a_fit, b_fit = params
# print(a_fit,b_fit)
# =============================================================================
    
# Resultados Ajuste
v = 0.7
A, B = 317.80, 1.77  # Con mediciones Heptano y Dodecano
# A, B = 273, 1.62     # Los del paper 
# A, B = 76.22, 1.13   # 30 nm
# A, B = 53.39, 0.993  # 15 nm
# A, B = 46.19, 1.24   # 6 nm
# A, B = 12.69, 1.74   # 3 nm

D_med_dode = 1 / np.mean(Dif_dode**(1.42))
N_med_dode = (A**(1.42)*D_med_dode)**(v/(v+B))
DiNi_dode = A * N_med_dode**(-B)
N_x_dode = (DiNi_dode / (Dax_dode))**(-1.42)
N_y_dode = (DiNi_dode / (Dif_dode))**(-1.42)

D_med_hept = 1 / np.mean(Dif_hept**(1.42))
N_med_hept = (A**(1.42)*D_med_hept)**(v/(v+B))
DiNi_hept = A * N_med_hept**(-B)
N_x_hept = (DiNi_hept / (Dax_hept))**(-1.42)
N_y_hept = (DiNi_hept / (Dif_hept))**(-1.42)



graph.plot_cadena(N_x_dode, N_y_dode, 30, 'Dodecano', 'red')
graph.plot_cadena(N_x_hept, N_y_hept, 100, 'Heptano', 'green')


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
