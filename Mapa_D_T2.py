# Mapa T1-D Mouse

import core_IO_Andres as IO
import core_plot_am as graph
import pandas as pd

#Carpeta con las mediciones
pwd = ('H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/DT2_Todos/Dodecano/')
#Leemos los archivos necesarios:
param = IO.read_acq(pwd)
Z, Daxis, T2axis = IO.read_DT2(pwd)

nT2, nD = len(T2axis), len(Daxis)

alpha = 0.0001
Dmin, Dmax = -1, 1
T2min, T2max = 1, 4
T2xx = Dxx = 5

tEcho = param['echoTime']
S0, D, T2, K1, K2 = IO.initKernelDT2(nD, nT2, Daxis, T2axis, 
                                     Dmin, Dmax, T2min, T2max)

print(f'Starting NLI: Alpha = {alpha}.')
S, iter = IO.NLI_FISTA_2D(K1, K2, Z, alpha, S0)
if iter < 100000:
    print('Inversion ready!')
else:
    print('Warning!')
    print('Maximum number of iterations reached!')
    print('Try modifying T2Range and/or alpha settings.')



#Guardado de datos ya procesados, nos interesa S0,S,T2 y D

Transformada = pd.DataFrame(S)
VectorT2 = pd.DataFrame(T2)
VectorD = pd.DataFrame(D)
Transformada.to_csv(pwd+"Transformada.txt", index=False, header=False)
VectorT2.to_csv(pwd+"VectorT2.txt", index=False, header=False)
VectorD.to_csv(pwd+"VectorD.txt", index=False, header=False)

#Graficos

graph.MapaDT2(Daxis, T2axis, Z, D, T2, S.T, pwd, alpha, Dmin, Dmax, T2min, T2max, T2xx, Dxx)
graph.PlotDT2_D(Daxis, D, S, pwd, alpha, Dmin, Dmax)
graph.PlotDT2_T2(T2axis, T2, S, pwd, alpha, T2min, T2max)




### Agregar codigo de guardado de datos procesados
#IO.writeSRCPMG(T1, T2, S, root)

### Codigo Cheva
# =============================================================================
# MLap_D, MLap_T2 = IO.fitLapMag_DT2(Daxis, T2axis, D, T2, S)
# params='hola'   
# turri.SRCPMG(Daxis, T2axis, Z, D, T2, S, MLap_D, MLap_T2, pwd, 
#              alpha, Dmin, Dmax, T2min, T2max, params, tEcho)
# =============================================================================


