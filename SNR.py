# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:30:44 2024

@author: Lenovo
"""

import numpy as np
import core_IO_Andres as IO
import matplotlib.pyplot as plt

pwdT2 = ('H:/Unidades compartidas/TF-Andres/2024-Dif-T1/Mediciones/DT2_Todos/Agua/')
pwdT1 = ('H:/Unidades compartidas/TF-Andres/2024-Dif-T1/Mediciones/T1D_Todos/Agua_5000/')

zT2, __, T2axis = IO.read_DT2(pwdT2)
zT1, __, T1axis = IO.read_T1D(pwdT1)


T2 = np.sum(zT2, axis=0)
T1 = np.sum(zT1, axis=1)

SNRT1 = T1[-1] / T1[0]
n = len(T2)
SNRT2 = sum(T2[:5]) / sum(T2[n-5:n])

print(SNRT2)
print(SNRT1)


# =============================================================================
# fig, axs = plt.subplots()
# axs.plot(T2axis, T2, label = 'Distrib.', color = 'teal')
# axs.set_xlabel(r'$Tiempo$ [ms]')
# plt.show()
# 
# fig, axs = plt.subplots()
# axs.plot(T1axis, T1, label = 'Distrib.', color = 'teal')
# axs.set_xlabel(r'$Tiempo$ [ms]')
# plt.show()
# =============================================================================
