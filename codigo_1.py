#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 16:02:06 2025

@author: olafviggiano
"""


import numpy as np
import matplotlib.pyplot as plt

def RK_method(callback,y0,t,St,Re): # Função para resolver as EDOs em problemas de valor inicial
                                    # aplicando o método de Runge-Kutta de quarta ordem 
                                    # callback é o nome da função da EDO a ser resolvida
                                    # y0 é a condição inicial, St é o número de Stokes e Re o número de Reynolds
                                    # t é o vetor temporal
                                    
    y = np.zeros(len(t))            # definição da variável de saída
    y[0] = y0
    for i in range(0,len(t)-1):     # Para todo os instantes de tempo, calcula-se os fatores e determina-se o próximo y
        h = t[i+1]-t[i]
        F1 = h*callback(y[i],t[i],St,Re)
        F2 = h*callback((y[i]+F1/2),(t[i]+h/2),St,Re)
        F3 = h*callback((y[i]+F2/2),(t[i]+h/2),St,Re)
        F4 = h*callback((y[i]+F3),(t[i]+h),St,Re)
        y[i+1] = y[i] + 1/6*(F1 + 2*F2 + 2*F3 + F4)
    return y


def Fl(vz,t,St,Re):                 # Função da primeira EDO, desconsiderando os efeitos inerciais
    return (-vz+1)/St


    


St=2 
Re=0        
h=0.1       # Passo de tempo do método de Runge-Kutta
t_fin=10    # Tempo final

t = np.arange(0,t_fin+h,h)  # Vetor temporal para aplicar no Runge-Kutta
y0=0                        # Condição inicial 

t_ana = np.linspace(0,t_fin,10000)  #Vetor temporal para plotar a solução analítica


#%%
#Primeiro Item do enunciado: Solução da primeira EDO para diferentes Stokes

All_Stokes=[0.1,0.5,1,2,5,10]   # Lista contendo diferentes Stokes

plt.figure(num=1, dpi=300, figsize=(7, 5))
plt.rcParams["font.family"] = "serif"
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['savefig.facecolor']='white'


for St in All_Stokes:
    y_RK = RK_method(Fl,y0,t,St,Re)     # Solução numérica da primeira EDO
    y_ana=1-np.exp(-t_ana/St)           # Solução analítica da primeira EDO
    plt.plot(t_ana,y_ana,c="black")     # Plota em linhas pretas as soluções analíticas para cada Stokes
    plt.plot(t,y_RK,"o", label = "St="+str(St),markersize=2.5)  # Plota em pontos coloridos as soluções numéricas para cada Stokes
    


plt.title("Solução da sedimentação desconsiderando inércia (Re=0)")
plt.legend(loc=4)
plt.grid(True, color='lightgray')
plt.xlabel("$t$")
plt.ylabel("$y$")

plt.show()

#%%
#Segundo Item do enunciado: Efeito do passo de tempo na solução numérica
#Para esse caso, vamos considerar St=5 fixo e um tempo final maior t_fin=20
 
St=5
t_fin=20

All_h=[1,2,4,6,8,10]  # Lista dos passos de tempo do método de Runge-Kutta

plt.figure(num=1, dpi=300, figsize=(7, 5))
plt.rcParams["font.family"] = "serif"
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['savefig.facecolor']='white'

t_ana = np.linspace(0,t_fin,10000)  #Vetor temporal para plotar a solução analítica
y_ana=1-np.exp(-t_ana/St)           # Solução analítica da primeira EDO

for h in All_h:
    t = np.arange(0,t_fin+h,h)
    y_RK = RK_method(Fl,y0,t,St,Re)     # Solução numérica da primeira EDO
    
    plt.plot(t,y_RK,"-", label = "h="+str(h),markersize=2.5)  # Plota em pontos coloridos as soluções numéricas para cada Stokes
plt.plot(t_ana,y_ana,c="black")     # Plota em linhas pretas as soluções analíticas para cada Stokes
plt.title("Efeito do passo de tempo na solução numérica")
plt.legend(loc=4)
plt.grid(True, color='lightgray')
plt.xlabel("$t$")
plt.ylabel("$y$")
plt.xlim([-1,20])
plt.show()

#%%
#Itens 3,4 e 5: Considere agora o termo de inércia, adimensionalize a equação e 
# compare a solução numérica com a solução analítica apresentada no artigo
# Vamos considerar St=1 fixo, para análise somente da influência de Re


def Fl_inerc(vz,t,St,Re):           # Função da segunda EDO, considerando os efeitos inerciais
    return (-vz-(3/8)*Re*vz**2+1)/St


St=1

h=0.1       # Passo de tempo do método de Runge-Kutta
t_fin=10    # Tempo final

t = np.arange(0,t_fin+h,h)  # Vetor temporal para aplicar no Runge-Kutta'
t_ana = np.linspace(0,t_fin,10000)  #Vetor temporal para plotar a solução analítica


All_Reynolds=[0.001,0.1,0.5,1,2,5]   # Lista contendo diferentes Reynolds

plt.figure(num=1, dpi=300, figsize=(7, 5))
plt.rcParams["font.family"] = "serif"
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['savefig.facecolor']='white'

for Re in All_Reynolds:
    y_RK=RK_method(Fl_inerc, y0, t, St,Re)  #Solução numérica da segunda EDO
    vt=(-8+np.sqrt(64+96*Re))/(6*Re)    # Velocidade terminal da esfera, considerando efeitos de inércia
    eps=(3/8)*Re
    A=St
    y_ana=vt*(1+((1+2*eps*vt)/((-1-eps*vt)*np.exp((1+2*eps*vt)*t_ana/A)-(eps*vt))))     #Solução analítica apresentada no artigo
    plt.plot(t_ana,y_ana,c="black")     # Plota em linhas pretas as soluções analíticas para cada Reynolds
    plt.plot(t,y_RK,"o", label = "Re="+str(Re),markersize=2.5)  # Plota em pontos coloridos as soluções numéricas para cada Reynolds
    
plt.title("Solução da sedimentação considerando inércia")
plt.legend(loc=4)
plt.grid(True, color='lightgray')
plt.xlabel("$t$")
plt.ylabel("$y$")

plt.show()



#No limite Re->0, podemos comparar as duas soluções:
Re=0.0001

plt.figure(num=1, dpi=300, figsize=(7, 5))
plt.rcParams["font.family"] = "serif"
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['savefig.facecolor']='white'
y_RK=RK_method(Fl_inerc, y0, t, St,Re)  #Solução numérica da segunda EDO
y_ana=1-np.exp(-t_ana/St)           # Solução analítica da primeira EDO
plt.plot(t_ana,y_ana,c="black",label="Solução analítica para Re=0")     # Plota em linhas pretas as soluções analíticas para cada Reynolds
plt.plot(t,y_RK,"o", label = "Solução numérica para Re-->0",markersize=2.5)  # Plota em pontos coloridos as soluções numéricas para cada Reynolds
    
plt.title("Limite em que Re-->0")
plt.legend(loc=4)
plt.grid(True, color='lightgray')
plt.xlabel("$t$")
plt.ylabel("$y$")

plt.show()




