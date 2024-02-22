# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 17:16:54 2023

@author: Andres
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

pwd = 'G:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/DT2_Todos/Dodecano/'
pwdd = 'G:/Unidades compartidas/TF-Andres/Graficos_Finales/'

S = pd.read_csv(pwd+"data2D.dat", header=None,delim_whitespace=True).to_numpy()
T2 = pd.read_csv(pwd+"T2Axis.dat", header=None,delim_whitespace=True).to_numpy()
D = pd.read_csv(pwd+"DiffAxis.dat", header=None,delim_whitespace=True).to_numpy()

fig = plt.figure(dpi=600)
ax = fig.add_subplot(111, projection='3d')
ax.set_title('$D$-$T_2$', fontsize=10)
surf = ax.plot_surface(T2, D.T, S.T, cmap='inferno')
fig.colorbar(surf)
ax.view_init(elev=35, azim=50)
ax.set_xlabel(r'$Tiempo$ [ms]', fontsize=10)
ax.set_ylabel(r'b-value [$10^{9}$ $s/m^{2}$]', fontsize=10)
ax.yaxis.set_tick_params(labelsize=10)
ax.xaxis.set_tick_params(labelsize=10)
ax.zaxis.set_tick_params(labelsize=10)
#ax.set_zlabel('Eje Z')
#plt.savefig(pwdd+"Datos_Surf_DT2", dpi=600)
plt.show()