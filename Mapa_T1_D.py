# Mapa T1-D Mouse
import core_IO_Andres as IO
import core_plot_am as graph
import pandas as pd


#Carpeta con las mediciones

# =============================================================================
# pwd_agua = ('G:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Agua_5000/')
# pwd_dode = ('G:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Heptano/')
# pwd_hept = ('G:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Dodecano/')
# =============================================================================

pwd = ('H:/Unidades compartidas/TF-Andres/Mediciones/Mediciones finales/T1D_Todos/Dodecano/')

#Leemos los archivos necesarios:

param = IO.read_acq(pwd)
# =============================================================================
# Z_agua, T1axis, Daxis = IO.read_T1D(pwd_agua)
# Z_hept, _, _ = IO.read_T1D(pwd_hept)
# Z_dode, _, _ = IO.read_T1D(pwd_dode)
# Z = Z_dode + Z_agua + Z_hept
# =============================================================================

Z, T1axis, Daxis = IO.read_T1D(pwd)


nT1, nD = len(T1axis), len(Daxis)


alpha = 0.0001
T1min, T1max = 1/2, 5
Dmin, Dmax = -1, 1


T1xx = 10
Dxx = 10

tEcho = param['echoTime']

S0, T1, D, K1, K2 = IO.initKernelT1D(nT1, nD, T1axis, Daxis, 
                                     T1min, T1max, Dmin, Dmax)     
print(f'Starting NLI: Alpha = {alpha}.')
S, iter = IO.NLI_FISTA_2D(K1, K2, Z, alpha, S0)
if iter < 100000:
    print('Inversion ready!')
else:
    print('Warning!')
    print('Maximum number of iterations reached!')
    print('Try modifying T2Range and/or alpha settings.')


#Guardamos los datos procesados

Transformada = pd.DataFrame(S)
VectorT1 = pd.DataFrame(T1)
VectorD = pd.DataFrame(D)
Transformada.to_csv(pwd+"Transformada.txt", index=False, header=False)
VectorT1.to_csv(pwd+"VectorT1.txt", index=False, header=False)
VectorD.to_csv(pwd+"VectorD.txt", index=False, header=False)



#Graficos

graph.MapaT1D(T1axis, Daxis, Z, T1, D, S, pwd, alpha, T1min, T1max, Dmin, Dmax, T1xx, Dxx)
graph.PlotT1D_D(Daxis, D, S, pwd, alpha, Dmin, Dmax)
graph.PlotT1D_T1(T1axis, T1, S, pwd, alpha, T1min, T1max)
#graph.PlotT1D_T1_NoNorm(T1axis, T1, S, pwd, alpha, T1min, T1max)
#graph.PlotT1D_D_NoNorm(Daxis, D, S, pwd, alpha, Dmin, Dmax)
# CODIGO CHEVA

# =============================================================================
# MLap_D, MLap_T1 = IO.fitLapMag_DT1(Daxis, Taxis, D, T2, S)
# params='hola'   
# print('Plotting CPMG processed data...')
# turri.SRCPMG(T1axis, Daxis, Z, T1, D, S, MLap_T1, MLap_D, pwd, 
#              alpha, T1min, T1max, Dmin, Dmax, params, tEcho)
# =============================================================================



