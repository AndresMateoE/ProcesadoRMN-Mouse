
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import matplotlib as mpl

# PARAMETROS MATPLOTLIB ---------------------------------

plt.rcParams["font.weight"] = "normal"
plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "Verdana"

plt.rcParams["axes.labelweight"] = "normal"
plt.rcParams["axes.linewidth"] = 1

plt.rcParams['xtick.major.size'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 1.5
plt.rcParams['xtick.minor.width'] = 1.5
plt.rcParams['ytick.major.size'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 1.5
plt.rcParams['ytick.minor.width'] = 1.5

plt.rcParams["legend.loc"] = 'best'
plt.rcParams["legend.frameon"] = True
plt.rcParams["legend.fancybox"] = True
plt.rcParams["legend.shadow"] = True
plt.rcParams["legend.fontsize"] = 5
plt.rcParams["legend.edgecolor"] = 'black'

plt.rcParams["lines.linewidth"] = 1.5
plt.rcParams["lines.markersize"] = 1
plt.rcParams["lines.linestyle"] = '-'

plt.rcParams["figure.autolayout"] = True

#-----------------------------------------------------


def MapaT1D(T1axis, Daxis, Z, T1, D, S, pwd, 
            alpha, T1min, T1max, Dmin, Dmax, T1xx, Dxx):

    #Eliminamos las curvas mas chicas del mapa
    A = S
    valor_umbral = 0.1 * np.max(S)
    A[A < valor_umbral] = 0
    #Eliminamos los bordes
    for i in range(0,Dxx):
        S[i,:] = 0
        S[-i,:] = 0

    for i  in range(0,T1xx):
        S[:,i] = 0
        S[:,-i] = 0

    #Graficamos
    maxi = np.max([T1min, Dmin])
    mini = np.min([T1max, Dmax])
    
    fig, ax = plt.subplots(dpi=600)
    fig.set_size_inches(10/2.54, 10/2.54)
    
# =============================================================================
#     ax.plot([10.0**mini, 10.0**maxi], [10.0**mini, 10.0**maxi], 
#                       color='black', ls='-', alpha=0.7, zorder=-2, 
#                       label = r'$T_1$ = $T_2$')
# =============================================================================
    bounds = np.linspace(0, np.max(S))
    cmap = mpl.cm.magma
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    N=8
    
    ax.set_title(rf'$\alpha$ = {alpha}', fontsize=10)
    contour = ax.contour(T1, D, A.T, N, cmap=cmap)
    x1, y1 = 856.87, 0.7
    x2, y2 = 1380, 2.57
    x = (np.linspace(10, 10000, 1000))
    y = (0.000000007)* x**(np.log10(y2/y1)/np.log10(x2/x1))
    ax.plot(x,y, linestyle='dashed', color='green', lw=1, zorder=-4)
    ax.plot([10,100000],[2.18,2.18], linestyle='dashed', color='blue', lw=1, zorder=-4)

    ax.set_xlabel(r'$T_1$ [ms]', fontsize=10)
    ax.set_ylabel(r'$D$ [10$^{-9}$ m$^{2}$/s]', fontsize=10)
    ax.set_xlim(10.0**T1min, 10.0**T1max)
    ax.set_ylim(10.0**Dmin, 10.0**Dmax)
    ax.set_xscale('log')
    ax.set_yscale('log')
#    ax.set_aspect('equal', adjustable='datalim')
    plt.tick_params(axis='both', which='major', labelsize=10)
    colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=ax, spacing='uniform', pad=0.01, ticks=np.around(np.linspace(0, np.max(S), N), decimals=4))
    colorbar.ax.tick_params(labelsize=5)
    colorbar.minorticks_off()
    plt.savefig(pwd+"Mapa_T1D", dpi=600)
    plt.show()
    


def MapaDT2(Daxis, T2axis, Z, D, T2, S, pwd, 
            alpha, Dmin, Dmax, T2min, T2max, T2xx, Dxx):
    #Eliminamos las primeras curvas
    A = S
    valor_umbral = 0.05 * np.max(S)
    A[A < valor_umbral] = 0
    #Eliminamos bordes
    for i in range(0,Dxx):
        S[i,:] = 0
        S[-i,:] = 0

    for i  in range(0,T2xx):
        S[:,i] = 0
        S[:,-i] = 0
    #Graficamos
    
    fig, ax = plt.subplots(dpi=300)
    fig.set_size_inches(10/2.54, 10/2.54)
    maxi = np.max([Dmin, T2min])
    mini = np.min([Dmax, T2max])
    ax.plot([10.0**mini, 10.0**maxi], [10.0**mini, 10.0**maxi], 
                      color='black', ls='-', alpha=0.2, zorder=-2)

    ax.set_title(rf'$\alpha$ = {alpha}', fontsize=10)
    ax.contour(T2, D, A.T, 100, cmap= 'inferno')
    ax.set_xlabel(r'$T_2$ [ms]', fontsize=10)
    ax.set_ylabel(r'$D$ [10$^{-9}$ m$^{2}$/s]', fontsize=10)
    ax.set_ylim(10.0**Dmin, 10.0**Dmax)
    ax.set_xlim(10.0**T2min, 10.0**T2max)
    ax.set_xscale('log')
    ax.set_yscale('log')
  #  ax.set_aspect('equal', adjustable='datalim')
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.savefig(pwd+"Mapa_DT2", dpi=600)
    plt.show()
    



def PlotT1D_D(Daxis, D, S, pwd, alpha, Dmin, Dmax):

    fig, axs = plt.subplots()
    #fig.set_size_inches(10/2.54, 10/2.54)    
    # Projected T2 distribution
    projD = np.sum(S, axis=0)
    projD = projD / np.max(projD)
    #Calculo los peaks
    peaks2, _ = find_peaks(projD, height=0.025, distance = 5)
    peaks2x, peaks2y = D[peaks2], projD[peaks2]
    # grafico los peaks
    axs.plot(D, projD, label = 'Distrib.', color = 'teal')
    for i in range(len(peaks2x)):
        axs.plot(peaks2x[i], peaks2y[i] + 0.05, lw = 0.2, marker=2, 
                      color='black')
        axs.annotate(f'{peaks2x[i]:.2f}', 
                          xy = (peaks2x[i], peaks2y[i] + 0.07), 
                          fontsize=10, ha = 'center')
    axs.set_xlabel(r'$D$ [10$^{-9}$ m$^{2}$/s]')
    axs.set_xscale('log')
    axs.set_ylim(-0.02, 1.2)
    axs.set_xlim(10.0**Dmin, 10.0**Dmax)
    #calculo comulativos
    cumD = np.cumsum(projD)
    cumD /= cumD[-1]
    ax = axs.twinx()
    ax.plot(D, cumD, label = 'Cumul.', color = 'coral')
    ax.set_ylim(-0.02, 1.2)
    #ax.set_aspect('equal', adjustable='datalim')    
    plt.savefig(pwd+"Diff_1D", dpi=600)
    plt.show()
    
def PlotT1D_D_NoNorm(Daxis, D, S, pwd, alpha, Dmin, Dmax):

    fig, axs = plt.subplots()
    #fig.set_size_inches(10/2.54, 10/2.54)    
    # Projected T2 distribution
    projD = np.sum(S, axis=0)
    #projD = projD / np.max(projD)
    #Calculo los peaks
    peaks2, _ = find_peaks(projD, height=0.005, distance = 5)
    peaks2x, peaks2y = D[peaks2], projD[peaks2]
    # grafico los peaks
    axs.plot(D, projD, label = 'Distrib.', color = 'teal')
    for i in range(len(peaks2x)):
        axs.plot(peaks2x[i], peaks2y[i] * 1.05, lw = 0.2, marker=2, 
                      color='black')
        axs.annotate(f'{peaks2x[i]:.2f}', 
                          xy = (peaks2x[i], peaks2y[i] * 1.07), 
                          fontsize=10, ha = 'center')
    axs.set_xlabel(r'$D$ [10$^{-9}$ m$^{2}$/s]')
    axs.set_xscale('log')
    axs.set_ylim(-0.001, 1.2*np.max(projD))
    axs.set_xlim(10.0**Dmin, 10.0**Dmax)
    #calculo comulativos
    cumD = np.cumsum(projD)
    cumD /= cumD[-1]
    ax = axs.twinx()
    ax.plot(D, cumD, label = 'Cumul.', color = 'coral')
    ax.set_ylim(-0.02, 1.2)
    #ax.set_aspect('equal', adjustable='datalim')    
    plt.savefig(pwd+"Diff_1D_NoNorm", dpi=600)
    plt.show()

def PlotT1D_T1(T1axis,T1, S, pwd, alpha, T1min, T1max):
    
    fig, axs = plt.subplots()
    #fig.set_size_inches(10/2.54, 10/2.54)
    projT1 = np.sum(S, axis=1)
    projT1 = projT1 / np.max(projT1)
    peaks1, _ = find_peaks(projT1, height=0.025, distance = 5)
    peaks1x, peaks1y = T1[peaks1], projT1[peaks1]
    
    
    axs.plot(T1, projT1, label = 'Distrib.', color = 'teal')
    for i in range(len(peaks1x)):
        axs.plot(peaks1x[i], peaks1y[i] + 0.05, lw = 0, marker=11, 
                      color='black')
        axs.annotate(f'{peaks1x[i]:.2f}', 
                          xy = (peaks1x[i], peaks1y[i] + 0.07), 
                          fontsize=10, ha = 'center')
    axs.set_xlabel(r'$T_1$ [ms]')
    axs.set_ylabel(r'Intensidad de Señal', fontsize=14)
    axs.set_xscale('log')
    axs.set_ylim(-0.02, 1.2)
    axs.set_xlim(10.0**T1min, 10.0**T1max)
    
    cumT1 = np.cumsum(projT1)
    cumT1 /= cumT1[-1]
    ax = axs.twinx()
    ax.plot(T1, cumT1, label = 'Cumul.', color = 'coral')
    ax.set_ylabel(r'Señal Acumulada', fontsize=14)
    ax.set_ylim(-0.02, 1.2)
    #ax.set_aspect('equal', adjustable='datalim')
    plt.savefig(pwd+"T1_1D", dpi=300)
    plt.show()
    
    
def PlotT1D_T1_NoNorm(T1axis,T1, S, pwd, alpha, T1min, T1max):
    
    fig, axs = plt.subplots()
    #fig.set_size_inches(10/2.54, 10/2.54)
    projT1 = np.sum(S, axis=1)
    #projT1 = projT1 / np.max(projT1)
    peaks1, _ = find_peaks(projT1, height=0.005, distance = 5)
    peaks1x, peaks1y = T1[peaks1], projT1[peaks1]
    
    
    axs.plot(T1, projT1, label = 'Distrib.', color = 'teal')
    for i in range(len(peaks1x)):
        axs.plot(peaks1x[i], peaks1y[i]*1.05, lw = 0, marker=11, 
                      color='black')
        axs.annotate(f'{peaks1x[i]:.2f}', 
                          xy = (peaks1x[i], peaks1y[i]*1.07), 
                          fontsize=10, ha = 'center')
    axs.set_xlabel(r'$T_1$ [ms]')
    axs.set_xscale('log')
    axs.set_ylim(-0.001, 1.2*np.max(projT1))
    axs.set_xlim(10.0**T1min, 10.0**T1max)
    
    cumT1 = np.cumsum(projT1)
    cumT1 /= cumT1[-1]
    ax = axs.twinx()
    ax.plot(T1, cumT1, label = 'Cumul.', color = 'coral')
    ax.set_ylim(-0.02, 1.2)
    #ax.set_aspect('equal', adjustable='datalim')
    plt.savefig(pwd+"T1_1D_NoNorm", dpi=600)
    plt.show()

def PlotDT2_T2(T2axis,T2, S, pwd, alpha, T2min, T2max):
    
    fig, axs = plt.subplots()
    #fig.set_size_inches(10/2.54, 10/2.54)
    projT2 = np.sum(S, axis=0)
    projT2 = projT2 / np.max(projT2)
    peaks1, _ = find_peaks(projT2, height=0.025, distance = 5)
    peaks1x, peaks1y = T2[peaks1], projT2[peaks1]
    
    
    axs.plot(T2, projT2, label = 'Distrib.', color = 'teal')
    for i in range(len(peaks1x)):
        axs.plot(peaks1x[i], peaks1y[i] + 0.05, lw = 0, marker=11, 
                      color='black')
        axs.annotate(f'{peaks1x[i]:.2f}', 
                          xy = (peaks1x[i], peaks1y[i] + 0.07), 
                          fontsize=10, ha = 'center')
    axs.set_xlabel(r'$T_2$ [ms]')
    axs.set_xscale('log')
    axs.set_ylim(-0.02, 1.2)
    axs.set_xlim(10.0**T2min, 10.0**T2max)
    
    cumT2 = np.cumsum(projT2)
    cumT2 /= cumT2[-1]
    ax = axs.twinx()
    ax.plot(T2, cumT2, label = 'Cumul.', color = 'coral')
    ax.set_ylim(-0.02, 1.2)
    #ax.set_aspect('equal', adjustable='datalim')
    plt.savefig(pwd+"T2_1D", dpi=600)
    plt.show()

def PlotDT2_T2_NoNorm(T2axis,T2, S, pwd, alpha, T2min, T2max):
    
    fig, axs = plt.subplots()
    #fig.set_size_inches(10/2.54, 10/2.54)
    projT2 = np.sum(S, axis=0)
    #projT2 = projT2 / np.max(projT2)
    peaks1, _ = find_peaks(projT2, height=0.005, distance = 5)
    peaks1x, peaks1y = T2[peaks1], projT2[peaks1]
    
    
    axs.plot(T2, projT2, label = 'Distrib.', color = 'teal')
    for i in range(len(peaks1x)):
        axs.plot(peaks1x[i], peaks1y[i] * 1.05, lw = 0, marker=11, 
                      color='black')
        axs.annotate(f'{peaks1x[i]:.2f}', 
                          xy = (peaks1x[i], peaks1y[i] * 1.07), 
                          fontsize=10, ha = 'center')
    axs.set_xlabel(r'$T_2$ [ms]')
    axs.set_xscale('log')
    axs.set_ylim(-0.001, 1.2*np.max(projT2))
    axs.set_xlim(10.0**T2min, 10.0**T2max)
    
    cumT2 = np.cumsum(projT2)
    cumT2 /= cumT2[-1]
    ax = axs.twinx()
    ax.plot(T2, cumT2, label = 'Cumul.', color = 'coral')
    ax.set_ylim(-0.02, 1.2)
    #ax.set_aspect('equal', adjustable='datalim')
    plt.savefig(pwd+"T2_1D_NoNorm", dpi=600)
    plt.show()

    
def PlotDT2_D(Daxis, D, S, pwd, alpha, Dmin, Dmax):

    fig, axs = plt.subplots()
    #fig.set_size_inches(10/2.54, 10/2.54)
    # Projected T2 distribution
    projD = np.sum(S, axis=1)
    projD = projD / np.max(projD)
    #Calculo los peaks
    peaks2, _ = find_peaks(projD, height=0.025, distance = 5)
    peaks2x, peaks2y = D[peaks2], projD[peaks2]
    # grafico los peaks
    axs.plot(D, projD, label = 'Distrib.', color = 'teal')
    for i in range(len(peaks2x)):
        axs.plot(peaks2x[i], peaks2y[i] + 0.05, lw = 0.2, marker=2, 
                      color='black')
        axs.annotate(f'{peaks2x[i]:.2f}', 
                          xy = (peaks2x[i], peaks2y[i] + 0.07), 
                          fontsize=10, ha = 'center')
    axs.set_xlabel(r'$D$ [10$^{-9}$ m$^{2}$/s]')
    axs.set_xscale('log')
    axs.set_ylim(-0.02, 1.2)
    axs.set_xlim(10.0**Dmin, 10.0**Dmax)
    #calculo comulativos
    cumD = np.cumsum(projD)
    cumD /= cumD[-1]
    ax = axs.twinx()
    ax.plot(D, cumD, label = 'Cumul.', color = 'coral')
    ax.set_ylim(-0.02, 1.2)
    #ax.set_aspect('equal', adjustable='datalim')
    plt.savefig(pwd+"Diff_1D", dpi=600)
    plt.show()