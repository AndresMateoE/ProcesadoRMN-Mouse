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

pwdDode = 'H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Dodecano/'

S_dode = pd.read_csv(pwdDode+"Transformada.txt", header=None).to_numpy()
T1_dode = pd.read_csv(pwdDode+"VectorT1.txt", header=None).to_numpy()
Dax_dode = pd.read_csv(pwdDode+"VectorD.txt", header=None).to_numpy()

pwdHept = 'H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Heptano/'

S_hept = pd.read_csv(pwdHept+"Transformada.txt", header=None).to_numpy()
T1_hept = pd.read_csv(pwdHept+"VectorT1.txt", header=None).to_numpy()
Dax_hept = pd.read_csv(pwdHept+"VectorD.txt", header=None).to_numpy()

pwdPet = 'H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Petroleo/'
 
S_oil = pd.read_csv(pwdPet+"Transformada.txt", header=None).to_numpy()
T1_oil = pd.read_csv(pwdPet+"VectorT1.txt", header=None).to_numpy()
Dax_oil = pd.read_csv(pwdPet+"VectorD.txt", header=None).to_numpy()

Dmin = -1
Dmax = 1
Dxx = 10
for i in range(0,Dxx):
    S_dode[i,:] = 0
    S_dode[-i,:] = 0
    S_hept[i,:] = 0
    S_hept[-i,:] = 0
    S_oil[i,:] = 0
    S_oil[-i,:] = 0
    
Dif_hept = np.sum(S_hept, axis=0)
#Dif_hept = Dif_hept / np.max(Dif_hept)

Dif_dode = np.sum(S_dode, axis=0)
#Dif_dode = Dif_dode / np.max(Dif_dode)

Dif_oil = np.sum(S_oil, axis=0)
#Dif_oil = Dif_oil / np.max(Dif_oil)


graph.plot_difusion(Dax_dode, Dif_dode, Dmin, Dmax, "Difusión Dodecano", "Red")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
