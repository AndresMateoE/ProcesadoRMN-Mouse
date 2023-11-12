# -*- coding: utf-8 -*-
"""
Editor de gráficos

La idea es ejecutar cada codigo de mapa para que quede guardado el valor de S.
Luego lo modificas a gusto en esta sección.

"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl

########
# DT2
########

# =============================================================================
# pwd = 'H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Agua_5000/'
# 
# S = pd.read_csv(pwd+"Transformada.txt", header=None).to_numpy()
# T2 = pd.read_csv(pwd+"VectorT2.txt", header=None).to_numpy()
# D = pd.read_csv(pwd+"VectorD.txt", header=None).to_numpy()
# 
# Dmin, Dmax = -1, 1
# T2min, T2max = 1, 4
# 
# fig, ax = plt.subplots(dpi=600)
# fig.set_size_inches(10/2.54, 10/2.54)
# maxi = np.max([Dmin, T2min])
# mini = np.min([Dmax, T2max])
# ax.plot([10.0**mini, 10.0**maxi], [10.0**mini, 10.0**maxi], 
#                   color='black', ls='-', alpha=0.7, zorder=-2)
# 
# A = S
# valor_umbral = 0.05* np.max(S)
# A[A < valor_umbral] = 0
# alpha = 0.001
# 
# T1xx = 10
# Dxx = 20
# for i in range(0,Dxx):
#     S[i,:] = 0
#     S[-i,:] = 0
# 
# for i  in range(0,T1xx):
#     S[:,i] = 0
#     S[:,-i] = 0
# 
# pwdd = 'H:/Unidades compartidas/TF-Andres/Graficos_Finales/'
# 
# bounds = np.linspace(0, np.max(S))
# cmap = mpl.cm.Blues_r
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
# N=8
# 
# ax.set_title(rf'$\alpha$ = {alpha}', fontsize=10)
# ax.contour(T2[:,0], D[:,0], A, N, vmax=0.004, cmap= cmap)
# 
# colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#               ax=ax, spacing='uniform', pad=0.01, ticks=np.round(np.linspace(0, np.max(S), N), decimals=5))
#               #ax=ax, spacing='uniform', pad=0.01, ticks=np.linspace(0, np.max(S), N))
# colorbar.ax.tick_params(labelsize=5)
# colorbar.minorticks_off()
# 
# ax.set_xlabel(r'$T_2$ [ms]', fontsize=10)
# ax.set_ylabel(r'$D$ [10$^{-9}$ m$^{2}$/s]', fontsize=10)
# ax.set_ylim(10.0**Dmin, 10.0**Dmax)
# ax.set_xlim(10.0**T2min, 10.0**T2max)
# ax.set_xscale('log')
# ax.set_yscale('log')
# plt.tick_params(axis='both', which='major', labelsize=10)
# #plt.savefig(pwdd+"Mapa_DT2_Agua", dpi=600)
# plt.show()
# =============================================================================

########
# T1D
########

pwd = 'H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Agua_3000/'

S = pd.read_csv(pwd+"Transformada.txt", header=None).to_numpy()
T1 = pd.read_csv(pwd+"VectorT1.txt", header=None).to_numpy()
D = pd.read_csv(pwd+"VectorD.txt", header=None).to_numpy()


T1min, T1max = 1, 4
Dmin, Dmax = -1, 1

maxi = np.max([T1min, Dmin])
mini = np.min([T1max, Dmax])

fig, ax = plt.subplots(dpi=300)
fig.set_size_inches(10/2.54, 10/2.54)
ax.plot([10.0**mini, 10.0**maxi], [10.0**mini, 10.0**maxi], 
                  color='black', ls='-', alpha=0.7, zorder=-2, 
                  label = r'$T_1$ = $T_2$')
alpha = 0.0001
A = S
valor_umbral = 0.05*np.max(S)
A[A < valor_umbral] = 0

T1xx = 10
Dxx = 10
for i in range(0,Dxx):
    S[i,:] = 0
    S[-i,:] = 0

for i  in range(0,T1xx):
    S[:,i] = 0
    S[:,-i] = 0

pwdd = 'H:/Unidades compartidas/TF-Andres/Graficos_Finales/'

bounds = np.linspace(0, np.max(S))
cmap = mpl.cm.Blues_r
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
N=8

ax.set_title(rf'$\alpha$ = {alpha}', fontsize=10)
contour = ax.contour(T1[:,0], D[:,0], A.T, N, cmap=cmap)

colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
              ax=ax, spacing='uniform', pad=0.01, ticks=np.around(np.linspace(0, np.max(S), N), decimals=4))
colorbar.ax.tick_params(labelsize=5)
colorbar.minorticks_off()

ax.set_xlabel(r'$T_1$ [ms]', fontsize=10)
ax.set_ylabel(r'$D$ [10$^{-9}$ m$^{2}$/s]', fontsize=10)
ax.set_xlim(10.0**T1min, 10.0**T1max)
ax.set_ylim(10.0**Dmin, 10.0**Dmax)
ax.set_xscale('log')
ax.set_yscale('log')
plt.tick_params(axis='both', which='major', labelsize=10)
plt.savefig(pwdd+"Mapa_T1D_Agua_3000", dpi=600)
plt.show()
    

############
#   Mapas dobles DT2
############

# =============================================================================
# pwdAgua = 'H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/DT2_Todos/Agua/'
# 
# S_agua = pd.read_csv(pwdAgua+"Transformada.txt", header=None).to_numpy()
# T2_agua = pd.read_csv(pwdAgua+"VectorT2.txt", header=None).to_numpy()
# D_agua = pd.read_csv(pwdAgua+"VectorD.txt", header=None).to_numpy()
# 
# pwdDode = 'H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/DT2_Todos/Dodecano/'
# 
# S_dode = pd.read_csv(pwdDode+"Transformada.txt", header=None).to_numpy()
# T2_dode = pd.read_csv(pwdDode+"VectorT2.txt", header=None).to_numpy()
# D_dode = pd.read_csv(pwdDode+"VectorD.txt", header=None).to_numpy()
# 
# pwdHept = 'H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/DT2_Todos/Heptano/'
# 
# S_hept = pd.read_csv(pwdHept+"Transformada.txt", header=None).to_numpy()
# T2_hept = pd.read_csv(pwdHept+"VectorT2.txt", header=None).to_numpy()
# D_hept = pd.read_csv(pwdHept+"VectorD.txt", header=None).to_numpy()
# 
# 
# T2min, T2max = 1, 3
# Dmin, Dmax = -1, 1
# 
# maxi = np.max([T2min, Dmin])
# mini = np.min([T2max, Dmax])
# 
# fig, ax = plt.subplots()
# fig.set_size_inches(10/2.54, 10/2.54)
# ax.plot([10.0**mini, 10.0**maxi], [10.0**mini, 10.0**maxi], 
#                   color='black', ls='-', alpha=0.7, zorder=-2, 
#                   label = r'$T_1$ = $T_2$')
# 
# #Eliminamos los bordes
# T2xx = Dxx = 5
# for i in range(0,Dxx):
#     S_dode[i,:] = 0
#     S_dode[-i,:] = 0
#     S_hept[i,:] = 0
#     S_hept[-i,:] = 0
#     S_agua[i,:] = 0
#     S_agua[-i,:] = 0
# 
# for i  in range(0,T2xx):
#     S_dode[:,i] = 0
#     S_dode[:,-i] = 0
#     S_hept[:,i] = 0
#     S_hept[:,-i] = 0
#     S_agua[:,i] = 0
#     S_agua[:,-i] = 0
# 
# 
# #Filtramos
# A_dode = S_dode
# valor_umbral = 0.1 * np.max(S_dode)
# A_dode[A_dode < valor_umbral] = 0
# 
# A_agua = S_agua
# valor_umbral = 0.1 * np.max(S_agua)
# A_agua[A_agua < valor_umbral] = 0
# 
# A_hept = S_hept
# valor_umbral = 0.1 * np.max(S_hept)
# A_hept[A_hept < valor_umbral] = 0
# 
# # =============================================================================
# # S_todo = S_dode + S_hept + S_agua
# # 
# # A_todo = S_todo
# # valor_umbral = 0.001 * np.max(S_todo)
# # A_todo[A_todo < valor_umbral] = 0
# # =============================================================================
# 
# #print(np.max(S_agua))
# #print(np.max(S_hept))
# #print(np.max(S_dode))
# 
# 
# maxi = np.max([Dmin, T2min])
# mini = np.min([Dmax, T2max])
# ax.plot([10.0**mini, 10.0**maxi], [10.0**mini, 10.0**maxi], 
#                   color='black', ls='-', alpha=0.7, zorder=-2)
# 
# 
# 
# pwdd = 'H:/Unidades compartidas/TF-Andres/Presentaciones Martes/231024/'
# 
# ax.contour(T2_agua[:,0], D_agua[:,0], A_agua, 10, vmax=0.004, cmap= 'Blues_r')
# ax.contour(T2_dode[:,0], D_dode[:,0], A_dode, 10, vmax=0.004, cmap= 'Greens_r')
# ax.contour(T2_hept[:,0], D_hept[:,0], A_hept, 10, vmax=0.004, cmap= 'Reds_r')
# ax.set_xlabel(r'$T_2$ [ms]')
# ax.set_ylabel(r'$D$ [10$^{-9}$ m$^{2}$/s]', fontsize=10)
# ax.set_ylim(10.0**Dmin, 10.0**Dmax)
# ax.set_xlim(10.0**T2min, 10.0**T2max)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.yaxis.set_minor_formatter(plt.NullFormatter())
# plt.tick_params(axis='both', which='major', labelsize=10)
# plt.savefig(pwdd+"Mapa_DT2_Separado", dpi=300)
# plt.show()
# =============================================================================


############
#   Mapas dobles T1D
############

# =============================================================================
# pwdAgua = 'G:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Agua_5000/'
# 
# S_agua = pd.read_csv(pwdAgua+"Transformada.txt", header=None).to_numpy()
# T1_agua = pd.read_csv(pwdAgua+"VectorT1.txt", header=None).to_numpy()
# D_agua = pd.read_csv(pwdAgua+"VectorD.txt", header=None).to_numpy()
# 
# pwdDode = 'G:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Dodecano/'
# 
# S_dode = pd.read_csv(pwdDode+"Transformada.txt", header=None).to_numpy()
# T1_dode = pd.read_csv(pwdDode+"VectorT1.txt", header=None).to_numpy()
# D_dode = pd.read_csv(pwdDode+"VectorD.txt", header=None).to_numpy()
# 
# pwdHept = 'G:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Heptano/'
# 
# S_hept = pd.read_csv(pwdHept+"Transformada.txt", header=None).to_numpy()
# T1_hept = pd.read_csv(pwdHept+"VectorT1.txt", header=None).to_numpy()
# D_hept = pd.read_csv(pwdHept+"VectorD.txt", header=None).to_numpy()
# 
# T1min, T1max = 0, 4
# Dmin, Dmax = -1, 1
# 
# maxi = np.max([T1min, Dmin])
# mini = np.min([T1max, Dmax])
# 
# fig, ax = plt.subplots(dpi=300)
# fig.set_size_inches(10/2.54, 10/2.54)
# ax.plot([10.0**mini, 10.0**maxi], [10.0**mini, 10.0**maxi], 
#                   color='black', ls='-', alpha=0.7, zorder=-2, 
#                   label = r'$T_1$ = $T_2$')
# #Eliminamos los bordes
# T1xx = Dxx = 10
# for i in range(0,Dxx):
#     S_dode[i,:] = 0
#     S_dode[-i,:] = 0
#     S_hept[i,:] = 0
#     S_hept[-i,:] = 0
#     S_agua[i,:] = 0
#     S_agua[-i,:] = 0
# 
# for i  in range(0,T1xx):
#     S_dode[:,i] = 0
#     S_dode[:,-i] = 0
#     S_hept[:,i] = 0
#     S_hept[:,-i] = 0
#     S_agua[:,i] = 0
#     S_agua[:,-i] = 0
# 
# 
# #Filtramos
# A_dode = S_dode
# valor_umbral = 0.1 * np.max(S_dode)
# A_dode[A_dode < valor_umbral] = 0
# 
# A_agua = S_agua
# valor_umbral = 0.1 * np.max(S_agua)
# A_agua[A_agua < valor_umbral] = 0
# 
# A_hept = S_hept
# valor_umbral = 0.1 * np.max(S_hept)
# A_hept[A_hept < valor_umbral] = 0
# 
# S = S_hept + S_dode + S_agua
# 
# 
# pwdd = 'H:/Unidades compartidas/TF-Andres/Presentaciones Martes/231024/'
# 
# #ax.contour(T1_dode[:,0], D_dode[:,0], A_dode.T, 7, cmap= 'Blues_r')
# #ax.contour(T1_agua[:,0], D_agua[:,0], A_agua.T, 7, cmap= 'Greens_r')
# #ax.contour(T1_hept[:,0], D_hept[:,0], A_hept.T, 7, cmap= 'Reds_r')
# ax.contour(T1_hept[:,0], D_hept[:,0], S.T, 7, cmap= 'Reds_r')
# ax.set_xlabel(r'$T_1$ [ms]')
# ax.set_ylabel(r'$D$ [10$^{-9}$ m$^{2}$/s]', fontsize=10)
# ax.set_xlim(10.0**T1min, 10.0**T1max)
# ax.set_ylim(10.0**Dmin, 10.0**Dmax)
# ax.set_xscale('log')
# ax.set_yscale('log')
# plt.tick_params(axis='both', which='major', labelsize=10)
# #plt.savefig(pwdd+"Mapa_T1D_AHD", dpi=300)
# plt.show()
# 
# =============================================================================





