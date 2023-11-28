# Voy copiando las cosas que hagan falta sin orden (en principio).
#El primer objetivo son los mapas T1-D

import numpy as np
import pandas as pd
import scipy.fft as FT
from scipy.optimize import curve_fit



###############################################################################
########################### COSAS VARIAS
###############################################################################


def read_acq(pwd):
    
    param = pd.read_csv(pwd+'acqu.par', header = None, delim_whitespace = True).to_numpy()
    
    parametros = {}
    for parametro in param:
        nombre_parametro = parametro[0]
        valor_parametro = parametro[2]
        parametros[nombre_parametro] = valor_parametro
        
    
    return parametros

def read_T1D(pwd):
    data = pd.read_csv(pwd+'Data_2D_T1D.dat', header = None,  delim_whitespace=True).to_numpy()
    T1axis = pd.read_csv(pwd+'TsatAxis.dat', header = None, delim_whitespace = True).to_numpy()
    Daxis = pd.read_csv(pwd+'DiffAxis.dat', header = None, delim_whitespace = True).to_numpy()

    return data,T1axis,Daxis


def read_DT2(pwd):
    data = pd.read_csv(pwd+'data2D.dat', header = None,  delim_whitespace=True).to_numpy()
    T2axis = pd.read_csv(pwd+'T2Axis.dat', header = None, delim_whitespace = True).to_numpy()
    Daxis = pd.read_csv(pwd+'DiffAxis.dat', header = None, delim_whitespace = True).to_numpy()

    return data,Daxis,T2axis
###############################################################################
######################## Para Mapas T1-D 
###############################################################################

def initKernelT1D(nP1, nP2, T1axis, Daxis, T1min, T1max, Dmin, Dmax):
    '''
    Initialize variables for Laplace transform.
    '''


    nBinx = nBiny = 200
    S0 = np.ones((nBinx, nBiny))
    T1 = np.logspace(T1min, T1max, nBinx)
    D  = np.logspace(Dmin, Dmax, nBiny) 


    K1 = 1 - np.exp(- T1axis / T1)
    K2 = np.exp(- Daxis * D)

    return S0, T1, D, K1, K2

def initKernelT1D_c(nP1, nP2, T1axis, Daxis, T1min, T1max, Dmin, Dmax):
    '''
    Initialize variables for Laplace transform.
    '''


    nBinx = nBiny = 200
    S0 = np.ones((nBinx, nBiny))
    T1 = np.logspace(T1min, T1max, nBinx)
    D  = np.logspace(Dmin, Dmax, nBiny) 


    K1 = np.exp(- T1axis / T1)
    K2 = np.exp(- Daxis * D)

    return S0, T1, D, K1, K2


def initKernelDT2(nP1, nP2, Daxis, T2axis, Dmin, Dmax, T2min, T2max):
    '''
    Initialize variables for Laplace transform.
    '''

    nBinx = nBiny = 150
    S0 = np.ones((nBinx, nBiny))
    D  = np.logspace(Dmin, Dmax, nBinx) 
    T2 = np.logspace(T2min, T2max, nBiny)


    K1 = np.exp(- Daxis * D)
    '''
    Defino K2 de difusion T2 seria D
    '''
    
    K2 = np.exp(- T2axis / T2) 

    return S0, D, T2, K1, K2



def NLI_FISTA_2D(K1, K2, Z, alpha, S):
    '''
    Numeric Laplace inversion, based on FISTA.
    '''

    K1TK1 = K1.T @ K1
    K2TK2 = K2.T @ K2
    K1TZK2 = K1.T @ Z @ K2
    ZZT = np.trace(Z @ Z.T)

    invL = 1 / (np.trace(K1TK1) * np.trace(K2TK2) + alpha)
    factor = 1 - alpha * invL

    Y = S
    tstep = 1
    lastRes = np.inf

    for iter in range(100000):
        term2 = K1TZK2 - K1TK1 @ Y @ K2TK2
        Snew = factor * Y + invL * term2
        Snew[Snew<0] = 0

        tnew = 0.5 * (1 + np.sqrt(1 + 4 * tstep**2))
        tRatio = (tstep - 1) / tnew
        Y = Snew + tRatio * (Snew - S)
        tstep = tnew
        S = Snew

        if iter % 1000 == 0:
            TikhTerm = alpha * np.linalg.norm(S)**2
            ObjFunc = ZZT - 2 * np.trace(S.T @ K1TZK2) + np.trace(S.T @ K1TK1 @ S @ K2TK2) + TikhTerm

            Res = np.abs(ObjFunc - lastRes) / ObjFunc
            lastRes = ObjFunc
            print(f'\t# It = {iter} >>> Residue = {Res:.6f}')

            if Res < 1E-5:
                break

    return S, iter


def fitLapMag_T1D(T1axis, Daxis, T1, D, S):
    '''
    Fits decay from T1 and T2 distributions.
    '''

    print(f'\tFitting T1 projection from 2D-Laplace in time domain...')

    t1 = range(len(T1axis))
    d1 = range(len(T1))
    S1 = np.sum(S, axis=1)
    M1 = []

    for i in t1:
        m1 = 0
        for j in d1:
            m1 += S1[j] * (1 - np.exp(-T1axis[i] / T1[j]))
        M1.append(m1[0])

    print(f'\tFitting T2 projection from 2D-Laplace in time domain...')

    t2 = range(len(Daxis))
    d2 = range(len(D))
    S2 = np.sum(S, axis=0)
    M2 = []

    for i in t2:
        m2 = 0
        for j in d2:
            m2 += S2[j] * np.exp(-Daxis[i] * D[j])
        M2.append(m2[0])

    return np.array(M1), np.array(M2)

def fitLapMag_DT2(Daxis, T2axis, D, T2, S):
    '''
    Fits decay from T1 and T2 distributions.
    '''

    print(f'\tFitting T1 projection from 2D-Laplace in time domain...')

    t1 = range(len(Daxis))
    d1 = range(len(D))
    S1 = np.sum(S, axis=1)
    M1 = []

    for i in t1:
        m1 = 0
        for j in d1:
            m1 += S1[j] * (1 - np.exp(-Daxis[i] * D[j]))
        M1.append(m1[0])

    print(f'\tFitting T2 projection from 2D-Laplace in time domain...')

    t2 = range(len(T2axis))
    d2 = range(len(T2))
    S2 = np.sum(S, axis=0)
    M2 = []

    for i in t2:
        m2 = 0
        for j in d2:
            m2 += S2[j] * np.exp(-T2axis[i] / T2[j])
        M2.append(m2[0])

    return np.array(M1), np.array(M2)