# -*- coding: utf-8 -*-
"""
Distribuciones de Difusión a largo de cadena

"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
from scipy.optimize import minimize
from scipy.optimize import curve_fit

# Primero tengo que cargar los datos, vamos a poder trabajar con Heptano y Dodecano

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

# Recortemos los bordes
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


# Ahora tengo que sumar en T1 para obtener las distribuciones de Difusión

Dif_hept = np.sum(S_hept, axis=0)
#Dif_hept = Dif_hept / np.max(Dif_hept)

Dif_dode = np.sum(S_dode, axis=0)
#Dif_dode = Dif_dode / np.max(Dif_dode)

Dif_oil = np.sum(S_oil, axis=0)


fig, ax = plt.subplots(dpi=600)
ax.plot(Dax_dode, Dif_dode, label = 'Distrib.', color = 'red')
ax.set_title("Difusión Dodecano")
ax.set_xlabel(r'$D$ [10$^{-9}$ m$^{2}$/s]')
ax.set_xscale('log')
#ax.set_ylim(-0.02, 1.2)
ax.set_xlim(10.0**Dmin, 10.0**Dmax)
plt.show()

fig, ax = plt.subplots(dpi=600)
ax.plot(Dax_hept, Dif_hept, label = 'Distrib.', color = 'Green')
ax.set_title("Difusión Heptano")
ax.set_xlabel(r'$D$ [10$^{-9}$ m$^{2}$/s]')
ax.set_xscale('log')
#ax.set_ylim(-0.02, 1.2)
ax.set_xlim(10.0**Dmin, 10.0**Dmax)
plt.show()

fig, ax = plt.subplots(dpi=600)
ax.plot(Dax_oil, Dif_oil, label = 'Distrib.', color = 'purple')
ax.set_title("Difusión Petroleo")
ax.set_xlabel(r'$D$ [10$^{-9}$ m$^{2}$/s]')
ax.set_xscale('log')
#ax.set_ylim(-0.02, 1.2)
ax.set_xlim(10.0**Dmin, 10.0**Dmax)
plt.show()

# Ahora tengo que hacer los ajustes
# Defino las funciones



def D_ajuste(N,a,b):
    return a * N**(-b-0.7)

#result = D_ajuste(12,a,b)
#type(result)

params, covariance = curve_fit(D_ajuste, (7,12), (2.55,0.67))
a_fit, b_fit = params
print(a_fit,b_fit)

A, B = 317.80, 1.77

# Entonces ahora podemos usar esto para armar la otra distribución

def Largo_medio(D,a,b,v):
    return (a**(1/v)*D**(-1/v))**(v/(v+b))
















