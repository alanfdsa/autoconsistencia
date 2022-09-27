# Esse arquivo é parte do projeto de iniciação ciêntifica "Magnetismo e efeito magnetocalórico no composto GdNi4Al" da autoria de 
# Alan Fillipe de Souzaz Almeida (2022)

from math import tanh, sinh, log
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib 

kb = 0.086    # em meV/K
mub = 0.05788 # em meV/T
R = 8.316     # Gas constant in J/(mol.K)

#####################################################################################
# Bext = Campo magnético externo aplicado, tem que ser uma lista
# n = Número de íons magnéticos por célular unitária
# L = Momento angular orbital total
# S = Momento angular de spin total
# Tc = Temperatura crítica do material
# Ti = Temperatura inicial
# Tf = Temperatura final

def CalculoMxT(Bext,n,L,S,Tc,Ti,Tf):
    J = L+S # Momento angular total
    gj = 1 + (J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1)) # Fator de Landé
    Ms = n * gj * J * mub # Magnetização de saturação
    plambda = (3*kb*Tc)/(n*(mub**2)*((gj)**2)*J*(J+1)) # Intensidade do campo molecular
    Ti = Ti + 0.001 # Temperatura inicial não nula pra não dar divisão por zero
    Mo = Ms # Chute da magnetização inicial

    fig, axis = plt.subplots(ncols=1,figsize=(13,10))
    for Bext in Bext:
        listaT = []
        listaM = []
        for T in np.arange(Ti,Tf,0.1):
            y = (gj * mub * J * (plambda * Mo+Bext)) / (kb * T)
            M = Ms * (((2 * J + 1) / (2 * J)) * (1 / (tanh(((2 * J + 1) / (2 * J)) * y))) - (1 / (2 * J)) * (
                        1 / tanh(y / (2 * J))))
            while abs(M - Mo) > 10e-9:
                Mo = M
                y = (gj * mub * J * (plambda * Mo + Bext)) / (kb * T)
                M = Ms * (((2 * J + 1) / (2 * J)) * (1 / (tanh(((2 * J + 1) / (2 * J)) * y))) - (1 / (2 * J)) * (
                            1 / tanh(y / (2 * J))))
            listaM.append(M/mub)
            listaT.append(T)
        axis.plot(listaT,listaM,label="B$_{ext}$ = "+str(Bext)+" T")

    axis.set_xlabel("T (K)",fontsize=25)
    axis.set_ylabel("M/$\mu_B$",fontsize=25)
    axis.set_xlim(Ti,Tf)
    axis.set_ylim(0,)
    axis.legend(loc='center', bbox_to_anchor=(0.85, 0.83),fontsize=20)
    plot = plt.tight_layout()
    
#####################################################################################
# T = Campo magnético externo aplicado, tem que ser uma lista
# n = Número de íons magnéticos por célular unitária
# L = Momento angular orbital total
# S = Momento angular de spin total
# Tc = Temperatura crítica do material
# Bi = Temperatura inicial
# Bf = Temperatura final

def CalculoMxB(T,n,L,S,Tc,Bi,Bf):
    J = L+S # Momento angular total
    gj = 1 + (J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1)) # Fator de Landé
    Ms = n * gj * J * mub # Magnetização de saturação
    plambda = (3*kb*Tc)/(n*(mub**2)*((gj)**2)*J*(J+1)) # Intensidade do campo molecular
    Bi = Bi + 0.001 # Campo magnético inicial não nula pra não dar divisão por zero
    Mo = Ms # Chute da magnetização inicial
    
    fig, axis = plt.subplots(ncols=1,figsize=(13,10))
    for T in T:
      listaBext = []
      listaM = []
      for Bext in np.arange(Bi, Bf,0.001):
          y = (gj * mub * J * (plambda * Mo+Bext)) / (kb * T)
          M = Ms * (((2 * J + 1) / (2 * J)) * (1 / (tanh(((2 * J + 1) / (2 * J)) * y))) - (1 / (2 * J)) * (
                    1 / tanh(y / (2 * J))))
          while abs(M - Mo) > 10e-8:
              Mo = M
              y = (gj * mub * J * (plambda * Mo + Bext)) / (kb * T)
              M  = Ms * (((2 * J + 1) / (2 * J)) * (1 / (tanh(((2 * J + 1) / (2 * J)) * y))) - (1 / (2 * J)) * (
                        1 / tanh(y / (2 * J))))
          listaM.append(M/mub)
          listaBext.append(Bext)
      axis.plot(listaBext,listaM,label="T = "+str(T)+" K")


    axis.set_xlabel("B$_{ext}$ (T)",fontsize=25)
    axis.set_ylabel("M/$\mu_B$",fontsize=25)
    axis.set_xlim(Bi,Bf)
    axis.set_ylim(0,)
    axis.legend(loc='center', bbox_to_anchor=(0.8, 0.2),fontsize=20)

    plt.tight_layout()