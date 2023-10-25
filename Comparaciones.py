# -*- coding: utf-8 -*-
"""
Comparacion de mapas de distintas dimensiones

@author: Lenovo
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt




pwdHept8 = 'H:/Unidades compartidas/TF-Andres/Presentaciones Martes/231024/Tiempos_Heptano/8x8/1/'

S_H8 = pd.read_csv(pwdHept8+"Transformada.txt", header=None).to_numpy()
T1_H8 = pd.read_csv(pwdHept8+"VectorT1.txt", header=None).to_numpy()
D_H8 = pd.read_csv(pwdHept8+"VectorD.txt", header=None).to_numpy()

pwdHept16 = 'H:/Unidades compartidas/TF-Andres/Presentaciones Martes/231024/Tiempos_Heptano/16x16/1/'

S_H16 = pd.read_csv(pwdHept16+"Transformada.txt", header=None).to_numpy()
T1_H16 = pd.read_csv(pwdHept16+"VectorT1.txt", header=None).to_numpy()
D_H16 = pd.read_csv(pwdHept16+"VectorD.txt", header=None).to_numpy()

pwdDode16 = 'H:/Unidades compartidas/TF-Andres/Presentaciones Martes/231024/Tiempos_Dodecano/16x16/1/'

S_D16 = pd.read_csv(pwdDode16+"Transformada.txt", header=None).to_numpy()
T1_D16 = pd.read_csv(pwdDode16+"VectorT1.txt", header=None).to_numpy()
D_D16 = pd.read_csv(pwdDode16+"VectorD.txt", header=None).to_numpy()

pwdDode8 = 'H:/Unidades compartidas/TF-Andres/Presentaciones Martes/231024/Tiempos_Dodecano/8x8/1/'

S_D8 = pd.read_csv(pwdDode8+"Transformada.txt", header=None).to_numpy()
T1_D8 = pd.read_csv(pwdDode8+"VectorT1.txt", header=None).to_numpy()
D_D8 = pd.read_csv(pwdDode8+"VectorD.txt", header=None).to_numpy()


T1min, T1max = 2, 4
Dmin, Dmax = -1, 1

maxi = np.max([T1min, Dmin])
mini = np.min([T1max, Dmax])

fig, ax = plt.subplots()
fig.set_size_inches(10/2.54, 10/2.54)
ax.plot([10.0**mini, 10.0**maxi], [10.0**mini, 10.0**maxi], 
                  color='black', ls='-', alpha=0.7, zorder=-2, 
                  label = r'$T_1$ = $T_2$')
#Eliminamos los bordes
T1xx = 20
Dxx = 10
for i in range(0,T1xx):
    S_D8[i,:] = 0
    S_D8[-i,:] = 0
    S_H8[i,:] = 0
    S_H8[-i,:] = 0
    S_D16[i,:] = 0
    S_D16[-i,:] = 0
    S_H16[i,:] = 0
    S_H16[-i,:] = 0
    
for i  in range(0,Dxx):
    S_D8[:,i] = 0
    S_D8[:,-i] = 0
    S_D16[:,i] = 0
    S_D16[:,-i] = 0
    S_H8[:,i] = 0
    S_H8[:,-i] = 0
    S_H16[:,i] = 0
    S_H16[:,-i] = 0

#Filtramos
# =============================================================================
# A_dode = S_dode
# valor_umbral = 0.05 * np.max(S_dode)
# A_dode[A_dode < valor_umbral] = 0
# 
# A_agua = S_agua
# valor_umbral = 0.05 * np.max(S_agua)
# A_agua[A_agua < valor_umbral] = 0
# 
# A_hept = S_hept
# valor_umbral = 0.05 * np.max(S_hept)
# A_hept[A_hept < valor_umbral] = 0
# =============================================================================


pwdd = 'H:/Unidades compartidas/TF-Andres/Presentaciones Martes/231024/'

ax.contour(T1_D16[:,0], D_D16[:,0], S_D16.T, 7, cmap= 'Blues_r')
ax.contour(T1_H16[:,0], D_H16[:,0], S_H16.T, 7, cmap= 'Reds_r')
#ax.contour(T1_hept[:,0], D_hept[:,0], S_hept.T, 7, cmap= 'Reds_r')
ax.set_xlabel(r'$T_1$ [ms]')
ax.set_ylabel(r'$D$ [10$^{-9}$ m$^{2}$/s]', fontsize=10)
ax.set_xlim(10.0**T1min, 10.0**T1max)
ax.set_ylim(10.0**Dmin, 10.0**Dmax)
ax.set_xscale('log')
ax.set_yscale('log')
plt.tick_params(axis='both', which='major', labelsize=10)
plt.savefig(pwdd+"Mapa_T1D_H16_D16", dpi=300)
plt.show()

fig, ax = plt.subplots()
fig.set_size_inches(10/2.54, 10/2.54)
ax.plot([10.0**mini, 10.0**maxi], [10.0**mini, 10.0**maxi], 
                  color='black', ls='-', alpha=0.7, zorder=-2, 
                  label = r'$T_1$ = $T_2$')

ax.contour(T1_D8[:,0], D_D8[:,0], S_D8.T, 7, cmap= 'Blues_r')
ax.contour(T1_H8[:,0], D_H8[:,0], S_H8.T, 7, cmap= 'Reds_r')
#ax.contour(T1_hept[:,0], D_hept[:,0], S_hept.T, 7, cmap= 'Reds_r')
ax.set_xlabel(r'$T_1$ [ms]')
ax.set_ylabel(r'$D$ [10$^{-9}$ m$^{2}$/s]', fontsize=10)
ax.set_xlim(10.0**T1min, 10.0**T1max)
ax.set_ylim(10.0**Dmin, 10.0**Dmax)
ax.set_xscale('log')
ax.set_yscale('log')
plt.tick_params(axis='both', which='major', labelsize=10)
plt.savefig(pwdd+"Mapa_T1D_H8_D8", dpi=300)
plt.show()