# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 14:08:29 2023

Codigo para graficar varias proyecciones

@author: Andres
"""

import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import find_peaks

#Elijo los datos a procesar

#(Daxis, D, S, pwd, alpha, Dmin, Dmax)

pwd1659 = 'H:/Unidades compartidas/TF-Andres/Mediciones/Dodecano_Bentheimer/Dodecano_Bentheimer/231002/Mojado_Agua/1659/231002_T1D/1/'

S_1659 = pd.read_csv(pwd1659+"Transformada.txt", header=None).to_numpy()
T1_1659 = pd.read_csv(pwd1659+"VectorT1.txt", header=None).to_numpy()
D_1659 = pd.read_csv(pwd1659+"VectorD.txt", header=None).to_numpy()

pwd1544 = 'H:/Unidades compartidas/TF-Andres/Mediciones/Dodecano_Bentheimer/Dodecano_Bentheimer/231002/Mojado_Agua/1544/231002_T1D/1/'

S_1544 = pd.read_csv(pwd1544+"Transformada.txt", header=None).to_numpy()
T1_1544 = pd.read_csv(pwd1544+"VectorT1.txt", header=None).to_numpy()
D_1544 = pd.read_csv(pwd1544+"VectorD.txt", header=None).to_numpy()

pwd1358 = 'H:/Unidades compartidas/TF-Andres/Mediciones/Dodecano_Bentheimer/Dodecano_Bentheimer/231002/Mojado_Agua/1358/231002_T1D/1/'

S_1358 = pd.read_csv(pwd1358+"Transformada.txt", header=None).to_numpy()
T1_1358 = pd.read_csv(pwd1358+"VectorT1.txt", header=None).to_numpy()
D_1358 = pd.read_csv(pwd1358+"VectorD.txt", header=None).to_numpy()

pwd1152 = 'H:/Unidades compartidas/TF-Andres/Mediciones/Dodecano_Bentheimer/Dodecano_Bentheimer/231002/Mojado_Agua/1152/231002_T1D/1/'

S_1152 = pd.read_csv(pwd1152+"Transformada.txt", header=None).to_numpy()
T1_1152 = pd.read_csv(pwd1152+"VectorT1.txt", header=None).to_numpy()
D_1152 = pd.read_csv(pwd1152+"VectorD.txt", header=None).to_numpy()

Dmin , Dmax = -1, 1

projD1659 = np.sum(S_1659, axis=0)
projD1544 = np.sum(S_1544, axis=0)
projD1358 = np.sum(S_1358, axis=0)
projD1152 = np.sum(S_1152, axis=0)


plt.rcParams.update({'font.size': 8})

fig, axs = plt.subplots(4,1)
fig.set_size_inches(8/2.54, 20/2.54)
#fig.set_size_inches(10/2.54, 10/2.54)    
# Projected T2 distribution

############# 1152
peaks2, _ = find_peaks(projD1152, height=0.005, distance = 5)
peaks2x, peaks2y = D_1152[peaks2], projD1152[peaks2]
# grafico los peaks

axs[0].plot(D_1152, projD1152, label = 'Distrib.', color = 'teal')
for i in range(len(peaks2x)):
    axs[0].plot(peaks2x[i], peaks2y[i] * 1.05, lw = 0.1, marker=2, 
                  color='black')
    axs[0].annotate(f'{peaks2x[i,0]:.2f}', 
                      xy = (peaks2x[i], peaks2y[i] * 1.07), 
                      fontsize=10, ha = 'center')
axs[0].set_xlabel(r'$D$ [10$^{-9}$ m$^{2}$/s]')
axs[0].set_xscale('log')
axs[0].set_ylim(-0.001, 1.2*np.max(projD1152))
axs[0].set_xlim(10.0**Dmin, 10.0**Dmax)
#calculo comulativos
cumD = np.cumsum(projD1152)
cumD /= cumD[-1]
ax = axs[0].twinx()
ax.plot(D_1152, cumD, label = 'Cumul.', color = 'coral')
ax.set_ylim(-0.02, 1.2)
#ax.set_aspect('equal', adjustable='datalim')    

############## 1358
peaks2, _ = find_peaks(projD1358, height=0.005, distance = 5)
peaks2x, peaks2y = D_1358[peaks2], projD1358[peaks2]
# grafico los peaks
axs[1].plot(D_1358, projD1358, label = 'Distrib.', color = 'teal')
for i in range(len(peaks2x)):
    axs[1].plot(peaks2x[i], peaks2y[i] * 1.05, lw = 0.1, marker=2, 
                  color='black')
    axs[1].annotate(f'{peaks2x[i,0]:.2f}', 
                      xy = (peaks2x[i], peaks2y[i] * 1.07), 
                      fontsize=10, ha = 'center')
axs[1].set_xlabel(r'$D$ [10$^{-9}$ m$^{2}$/s]')
axs[1].set_xscale('log')
axs[1].set_ylim(-0.001, 1.2*np.max(projD1358))
axs[1].set_xlim(10.0**Dmin, 10.0**Dmax)
#calculo comulativos
cumD = np.cumsum(projD1358)
cumD /= cumD[-1]
ax = axs[1].twinx()
ax.plot(D_1358, cumD, label = 'Cumul.', color = 'coral')
ax.set_ylim(-0.02, 1.2)

############## 1544
peaks2, _ = find_peaks(projD1544, height=0.005, distance = 5)
peaks2x, peaks2y = D_1544[peaks2], projD1544[peaks2]
# grafico los peaks
axs[2].plot(D_1544, projD1544, label = 'Distrib.', color = 'teal')
for i in range(len(peaks2x)):
    axs[2].plot(peaks2x[i], peaks2y[i] * 1.05, lw = 0.1, marker=2, 
                  color='black')
    axs[2].annotate(f'{peaks2x[i,0]:.2f}', 
                      xy = (peaks2x[i], peaks2y[i] * 1.07), 
                      fontsize=10, ha = 'center')
axs[2].set_xlabel(r'$D$ [10$^{-9}$ m$^{2}$/s]')
axs[2].set_xscale('log')
axs[2].set_ylim(-0.001, 1.2*np.max(projD1544))
axs[2].set_xlim(10.0**Dmin, 10.0**Dmax)
#calculo comulativos
cumD = np.cumsum(projD1544)
cumD /= cumD[-1]
ax = axs[2].twinx()
ax.plot(D_1544, cumD, label = 'Cumul.', color = 'coral')
ax.set_ylim(-0.02, 1.2)

############## 1659
peaks2, _ = find_peaks(projD1659, height=0.005, distance = 5)
peaks2x, peaks2y = D_1659[peaks2], projD1659[peaks2]
# grafico los peaks
axs[3].plot(D_1659, projD1659, label = 'Distrib.', color = 'teal')
for i in range(len(peaks2x)):
    axs[3].plot(peaks2x[i], peaks2y[i] * 1.05, lw = 0.1, marker=2, 
                  color='black')
    axs[3].annotate(f'{peaks2x[i,0]:.2f}', 
                      xy = (peaks2x[i], peaks2y[i] * 1.07), 
                      fontsize=10, ha = 'center')
axs[3].set_xlabel(r'$D$ [10$^{-9}$ m$^{2}$/s]')
axs[3].set_xscale('log')
axs[3].set_ylim(-0.001, 1.2*np.max(projD1659))
axs[3].set_xlim(10.0**Dmin, 10.0**Dmax)
#calculo comulativos
cumD = np.cumsum(projD1659)
cumD /= cumD[-1]
ax = axs[3].twinx()
ax.plot(D_1659, cumD, label = 'Cumul.', color = 'coral')
ax.set_ylim(-0.02, 1.2)

pwdd = ('G:/Unidades compartidas/TF-Andres/Mediciones/Dodecano_Bentheimer/Dodecano_Bentheimer/231002/Mojado_Agua/GraficosJuntos/')
plt.savefig(pwdd+"Diff_1D_NoNorm_Mojado", dpi=300)
plt.show()