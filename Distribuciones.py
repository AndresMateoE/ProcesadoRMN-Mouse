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
from scipy.stats import moment
import core_plot_am as graph

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
#Dif_oil = Dif_oil / np.max(Dif_oil)

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
# Solo tengo Heptano y Dodecano, asi que con esos dos puntos hago el ajuste para conocer A y B

def D_ajuste(N,a,b):
    return a* N**(-b-0.7)

# =============================================================================
# params, covariance = curve_fit(D_ajuste, (7,12), (2.55,0.67))
# a_fit, b_fit = params
# print(a_fit,b_fit)
# =============================================================================

# =============================================================================
# fig, ax = plt.subplots(dpi=600)
# ax.scatter((7,12), (2.55,0.67), label = 'P/ajuste.', color = 'purple', lw=4)
# ax.plot(np.linspace(0, 40),D_ajuste(np.linspace(0, 40),1,b_fit), color = 'red', zorder=-2)
# #ax.set_title("Difusión Petroleo")
# #ax.set_xlabel(r'$D$ [10$^{-9}$ m$^{2}$/s]')
# #ax.set_xscale('log')
# ax.set_ylim(-3, 10)
# ax.set_xlim(5, 40)
# plt.show()
# =============================================================================


A, B = 317.80, 1.77
#A, B = 1, 1.77
v = 0.7
# Entonces ahora podemos usar esto para armar la otra distribución
# Primero podemos calcular N media como sugiere el paper y utilizando la distribución de petroleo



D_med = 1 / np.mean(Dif_oil**(1.42))
N_med = (A**(1/v)*D_med)**(v/(v+B))

# Por último deberia encontrar la distribución de largos de cadena
# Con todo esto tenemos la relacion DiNi

DiNi = A * N_med**(-B)

N_x = (DiNi / (Dax_oil))**(-1/v)
N_y = (DiNi / (Dif_oil))**(-1/v)

N_y = N_y / np.max(N_y)


fig, ax = plt.subplots(dpi=600)
ax.plot(N_x, N_y, label = 'Distrib.', color = 'orange')
ax.set_title("Distribución largo de cadena Petroleo")
ax.set_xlabel('Largo de cadena')
#ax.set_ylim(-0.02, 1.2)
ax.set_xlim(0, 300)
plt.show()

####################################################################

D_med_dode = 1 / np.mean(Dif_dode**(1.42))
#N_med_dode = (A**(1/v) * D_med_dode**(-1/v))**(v/(v+B))
N_med_dode = (A**(1/v) * D_med_dode)**(v/(v+B))

# Por último deberia encontrar la distribución de largos de cadena
# Con todo esto tenemos la relacion DiNi

DiNi_dode = A * N_med_dode**(-B)

N_x_dode = (DiNi_dode / (Dax_dode))**(-1/v)
N_y_dode = (DiNi_dode / (Dif_dode))**(-1/v)

N_y_dode = N_y_dode / np.max(N_y_dode)


fig, ax = plt.subplots(dpi=600)
ax.plot(N_x, N_y, label = 'Distrib.', color = 'red')
ax.set_title("Distribución largo de cadena Dodecano")
ax.set_xlabel('Largo de cadena')
#ax.set_ylim(-0.02, 1.2)
ax.set_xlim(0, 20)
plt.show()

#################################################################################

D_med_hept = 1 / np.mean(Dif_hept**(1.42))
N_med_hept = (A**(1/v) * D_med_hept)**(v/(v+B))

# Por último deberia encontrar la distribución de largos de cadena
# Con todo esto tenemos la relacion DiNi

DiNi_hept = A * N_med_hept**(-B)

N_x_hept = (DiNi_hept / (Dax_hept))**(-1/v)
N_y_hept = (DiNi_hept / (Dif_hept))**(-1/v)

N_y_hept = N_y_hept / np.max(N_y_hept)


fig, ax = plt.subplots(dpi=600)
ax.plot(N_x_hept, N_y_hept, label = 'Distrib.', color = 'green')
ax.set_title("Distribución largo de cadena Heptano")
ax.set_xlabel('Largo de cadena')
#ax.set_ylim(-0.02, 1.2)
#ax.set_xlim(0, 20)
plt.show()

graph.plot_cadena(N_x_hept, N_y_hept, 30, 'Heptano', 'green')

# def Largo_medio(D,a,b,v):
    # return (a**(1/v)*D**(-1/v))**(v/(v+b))
















