#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 16:02:06 2025

@author: olafviggiano
"""


import numpy as np
import matplotlib.pyplot as plt

def RK_method(callback,y0,t,St,Re):
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(0,len(t)-1):
        h = t[i+1]-t[i]
        F1 = h*callback(y[i],t[i],St,Re)
        F2 = h*callback((y[i]+F1/2),(t[i]+h/2),St,Re)
        F3 = h*callback((y[i]+F2/2),(t[i]+h/2),St,Re)
        F4 = h*callback((y[i]+F3),(t[i]+h),St,Re)
        y[i+1] = y[i] + 1/6*(F1 + 2*F2 + 2*F3 + F4)
    return y


def Fl(vz,t,St,Re):
    return (-vz+1)/St

def Fl_inerc(vz,t,St,Re):
    return (-vz-(3/8)*Re*vz**2+1)/St
    


St=2
Re=0
h=0.1
t_fin=10


t = np.arange(0,t_fin+h,h)

h = t[2]-t[1]
y0=0

t_ana = np.linspace(0,t_fin,10000)

y_ana=1-np.exp(-t_ana/St)
y_RK = RK_method(Fl,y0,t,St,Re)


plt.figure(num=1, dpi=300, figsize=(7, 5))
plt.rcParams["font.family"] = "serif"
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['savefig.facecolor']='white'
plt.plot(t_ana,y_ana,'r-',label = "Analytical Solution")
plt.plot(t,y_RK,'ko', label = "Runge Kutta",markersize=2.5)
plt.legend(loc=4)
plt.grid(True, color='lightgray')
plt.xlabel("$t$ (-)")
plt.ylabel("$y$ (-)")
# plt.ylim([-2,2])
# plt.text(s='Euler\'s Method Vs. Runge Kutta 4th Order Method', x=1, y=1.03, fontsize=15, ha='center', va='center')
# plt.text(s=f'n = {points}, h = {h:.3f}', x=1, y=8000, fontsize=11, ha='center', va='center')
# plt.savefig('euler_rk_output.png', transparent=False, bbox_inches="tight")
plt.show()

#%%

Re=0.5

y_RK=RK_method(Fl_inerc, y0, t, St,Re)

vt=(-8+np.sqrt(64+96*Re))/(6*Re)
eps=(3/8)*Re
A=St

y_ana=vt*(1+((1+2*eps*vt)/((-1-eps*vt)*np.exp((1+2*eps*vt)*t_ana/A)-(eps*vt))))

# y_ana = vt - (2*eps*vt + 1) / (eps + (eps*vt + 1)/vt * np.exp((2*eps*vt + 1)/A * t_ana))


plt.figure(num=1, dpi=300, figsize=(7, 5))
plt.rcParams["font.family"] = "serif"
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['savefig.facecolor']='white'
plt.plot(t_ana,y_ana,'r-',label = "Analytical Solution")
plt.plot(t,y_RK,'ko', label = "Runge Kutta",markersize=2.5)
plt.legend(loc=4)
plt.grid(True, color='lightgray')
plt.xlabel("$t$ (-)")
plt.ylabel("$y$ (-)")
# plt.ylim([-2,2])
# plt.text(s='Euler\'s Method Vs. Runge Kutta 4th Order Method', x=1, y=1.03, fontsize=15, ha='center', va='center')
# plt.text(s=f'n = {points}, h = {h:.3f}', x=1, y=8000, fontsize=11, ha='center', va='center')
# plt.savefig('euler_rk_output.png', transparent=False, bbox_inches="tight")
plt.show()



