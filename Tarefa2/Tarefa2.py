#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 27 18:58:15 2025

@author: olafviggiano

---------------------------------------------------

COMENTARIOS SOBRE OS RESULTADOS NO FINAL DO ARQUIVO

---------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt

#Primeiramente vamos definir as funções do problema

def u(x,y,Ca_in,tau,k0,E,R):
    return (-Ca_in+x+(tau*k0*np.exp((-E/R)/y))*x)#/Ca_in

def v(x,y,DH,V,k0,E,R,rho,cp,Qdot,UA,Tc,Te):
    return (DH*V*k0*np.exp((-E/R)/y)*x-y*(rho*cp*Qdot+UA)+Te*(rho*cp*Qdot)+(Tc*UA))#/(Tc*(rho*cp*Qdot+UA))

def dudx(y,tau,k0,E,R):
    return (+1+tau*k0*np.exp((-E/R)/y))#/Ca_in

def dudy(x,y,tau,k0,E,R):
    return (-tau*k0*(E/R)*np.exp((-E/R)/y)*x/y**2)#/Ca_in
    
def dvdx(y,DH,V,k0,E,R):
    return (DH*V*k0*np.exp((-E/R)/y))#/(Tc*(rho*cp*Qdot+UA))

def dvdy(x,y,DH,V,k0,E,R,rho,cp,Qdot,UA):
    return (-DH*V*k0*(E/R)*np.exp((-E/R)/y)*(x/y**2)-(rho*cp*Qdot+UA))#/(Tc*(rho*cp*Qdot+UA))

def Det_J(J11,J12,J21,J22):
    return J11*J22-(J21*J12)


#%%
#Agora, define-se as variáveis relevantes do problema

V = 1.0         #m3
tau = 100         #s
Qdot=V/tau
k0 = 1.0 * 10**6 #  s−1
E = 8.0 * 10**4 #J/mol
R = 8.314 #J/mol.K
DH = 4.0 * 10**4 #J/mol
rho = 1000 #kg/m3
cp = 4180 #J/kg.K
UA = 2.0 * 10**4 #W/K
Tc = 300 #K

#Agora, determinamos as condições de entrada iniciais.

Ca_in = 2000    #mol/m3
Te = 330# K

erro=1
tol=10**-5

#Agora, vamos implementar um método para resolver as equações do problema

x=Ca_in    #Primeiro, definimos as variáveis iterativas iguais as condições de entrada
y=Te            
it=0

while erro>tol:     #fazendo iterações até que o erro seja menor do que a tolerância estabelecida
    it=it+1
    # print("Iteração= ",it)
    # print("Erro=",erro)
    #Agora, utilizando as funções definidas no topo do programa, calculamos os termos de cada iteração

    ui=u(x, y, Ca_in, tau, k0, E, R)
    vi=v(x, y, DH, V, k0, E, R, rho, cp, Qdot, UA, Tc, Te)
    
    dui_dx=dudx(y, tau, k0, E, R)
    dui_dy=dudy(x, y, tau, k0, E, R)
    dvi_dx=dvdx(y, DH, V, k0, E, R)
    dvi_dy=dvdy(x, y, DH, V, k0, E, R, rho, cp, Qdot, UA)
    
    DetJ=Det_J(dui_dx, dui_dy , dvi_dx , dvi_dy )
    
    x_new= x - 1/DetJ*(ui*dvi_dy-vi*dui_dy)
    y_new= y -1/DetJ*(vi*dui_dx-ui*dvi_dx)
    
    Rx=x-x_new # Calculamos a diferença entre a iteração atual e a anterior, para analizar a convergência
    Ry=y-y_new
    
    erro=np.max([np.abs(Rx),np.abs(Ry)])
    
    x=np.copy(x_new)
    y=np.copy(y_new)
    

print("Primeira parte do problema")
print()



print("Iterações= ",it)
print("Erro=",erro)
print()

print("Concentração na entrada=",Ca_in,"mol/m^3 ")
print("Temperatura na entrada=",Te,"°K")
print()
print("Concentração na saída=",x,"mol/m^3 ")
print("Temperatura na saída=",y,"°K")
print("--------------------------------------------------")


#%%

#Agora, modificamos os valores das variáveis fixas:
    
V = 1.0         #m3
tau = 50         #s
Qdot=V/tau
k0 = 1.0e8 #  s−1
E = 7.0e4 #J/mol
R = 8.314 #J/mol.K
DH = 2.0e5 #J/mol
rho = 1e3 #kg/m3
cp = 4180 #J/kg.K
UA = 1.5e3 #W/K
Tc = 280 #K

#Agora, determinamos as condições de entrada iniciais.

Ca_in = 600   #mol/m3
Te = 305# K

erro=1
tol=10**-8

#Agora, vamos implementar um método para resolver as equações do problema

x=Ca_in    #Primeiro, definimos as variáveis iterativas iguais as condições de entrada
y=Te            
it=0

while erro>tol:     #fazendo iterações até que o erro seja menor do que a tolerância estabelecida
    it=it+1
    # print("Iteração= ",it)
    # print("Erro=",erro)
    #Agora, utilizando as funções definidas no topo do programa, calculamos os termos de cada iteração

    ui=u(x, y, Ca_in, tau, k0, E, R)
    vi=v(x, y, DH, V, k0, E, R, rho, cp, Qdot, UA, Tc, Te)
    
    dui_dx=dudx(y, tau, k0, E, R)
    dui_dy=dudy(x, y, tau, k0, E, R)
    dvi_dx=dvdx(y, DH, V, k0, E, R)
    dvi_dy=dvdy(x, y, DH, V, k0, E, R, rho, cp, Qdot, UA)
    
    DetJ=Det_J(dui_dx, dui_dy , dvi_dx , dvi_dy )
    
    x_new= x - 1/DetJ*(ui*dvi_dy-vi*dui_dy)
    y_new= y -1/DetJ*(vi*dui_dx-ui*dvi_dx)
    
    Rx=x-x_new # Calculamos a diferença entre a iteração atual e a anterior, para analizar a convergência
    Ry=y-y_new
    
    erro=np.max([np.abs(Rx),np.abs(Ry)])
    
    x=np.copy(x_new)
    y=np.copy(y_new)
    
    if it>1000:
        break
    
print("Segunda parte do problema")
print()
print("Iterações= ",it)
print("Erro=",erro)
print()
print("Concentração na entrada=",Ca_in,"mol/m^3 ")
print("Temperatura na entrada=",Te,"°K")
print()
print("Concentração na saída=",x,"mol/m^3 ")
print("Temperatura na saída=",y,"°K")
print("--------------------------------------------------")
print()

Ca_in = 1200   #mol/m3
Te = 350# K

erro=1
tol=10**-8

#Agora, vamos implementar um método para resolver as equações do problema

x=Ca_in    #Primeiro, definimos as variáveis iterativas iguais as condições de entrada
y=Te            
it=0

while erro>tol:     #fazendo iterações até que o erro seja menor do que a tolerância estabelecida
    it=it+1
    # print("Iteração= ",it)
    # print("Erro=",erro)
    #Agora, utilizando as funções definidas no topo do programa, calculamos os termos de cada iteração

    ui=u(x, y, Ca_in, tau, k0, E, R)
    vi=v(x, y, DH, V, k0, E, R, rho, cp, Qdot, UA, Tc, Te)
    
    dui_dx=dudx(y, tau, k0, E, R)
    dui_dy=dudy(x, y, tau, k0, E, R)
    dvi_dx=dvdx(y, DH, V, k0, E, R)
    dvi_dy=dvdy(x, y, DH, V, k0, E, R, rho, cp, Qdot, UA)
    
    DetJ=Det_J(dui_dx, dui_dy , dvi_dx , dvi_dy )
    
    x_new= x - 1/DetJ*(ui*dvi_dy-vi*dui_dy)
    y_new= y -1/DetJ*(vi*dui_dx-ui*dvi_dx)
    
    Rx=x-x_new # Calculamos a diferença entre a iteração atual e a anterior, para analizar a convergência
    Ry=y-y_new
    
    erro=np.max([np.abs(Rx),np.abs(Ry)])
    
    x=np.copy(x_new)
    y=np.copy(y_new)
    
    if it>1000:
        break
    

print("Iterações= ",it)
print("Erro=",erro)
print()
print("Concentração na entrada=",Ca_in,"mol/m^3 ")
print("Temperatura na entrada=",Te,"°K")
print()
print("Concentração na saída=",x,"mol/m^3 ")
print("Temperatura na saída=",y,"°K")
print("--------------------------------------------------")
print()


Ca_in = 2000   #mol/m3
Te = 420# K

erro=1
tol=10**-8

#Agora, vamos implementar um método para resolver as equações do problema

x=Ca_in    #Primeiro, definimos as variáveis iterativas iguais as condições de entrada
y=Te            
it=0

while erro>tol:     #fazendo iterações até que o erro seja menor do que a tolerância estabelecida
    it=it+1
    # print("Iteração= ",it)
    # print("Erro=",erro)
    #Agora, utilizando as funções definidas no topo do programa, calculamos os termos de cada iteração

    ui=u(x, y, Ca_in, tau, k0, E, R)
    vi=v(x, y, DH, V, k0, E, R, rho, cp, Qdot, UA, Tc, Te)
    
    dui_dx=dudx(y, tau, k0, E, R)
    dui_dy=dudy(x, y, tau, k0, E, R)
    dvi_dx=dvdx(y, DH, V, k0, E, R)
    dvi_dy=dvdy(x, y, DH, V, k0, E, R, rho, cp, Qdot, UA)
    
    DetJ=Det_J(dui_dx, dui_dy , dvi_dx , dvi_dy )
    
    x_new= x - 1/DetJ*(ui*dvi_dy-vi*dui_dy)
    y_new= y -1/DetJ*(vi*dui_dx-ui*dvi_dx)
    
    Rx=x-x_new # Calculamos a diferença entre a iteração atual e a anterior, para analizar a convergência
    Ry=y-y_new
    
    erro=np.max([np.abs(Rx),np.abs(Ry)])
    
    x=np.copy(x_new)
    y=np.copy(y_new)
    
    if it>1000:
        break
    

print("Iterações= ",it)
print("Erro=",erro)
print()
print("Concentração na entrada=",Ca_in,"mol/m^3 ")
print("Temperatura na entrada=",Te,"°K")
print()
print("Concentração na saída=",x,"mol/m^3 ")
print("Temperatura na saída=",y,"°K")
print("--------------------------------------------------")
print()


#%%

V = 1.0         #m3
tau = 50         #s
Qdot=V/tau
k0 = 1.0e8 #  s−1
E = 7.0e4 #J/mol
R = 8.314 #J/mol.K
DH = 2.0e5 #J/mol
rho = 1e3 #kg/m3
cp = 4180 #J/kg.K
UA = 1.5e3 #W/K
Tc = 280 #K

N=100

All_Ca=np.linspace(500, 3000,N)

All_T=np.linspace(250, 400,N)

Ca_out=np.zeros((N,N))
Ca_out_norm=np.zeros((N,N))
T_out=np.zeros((N,N))



for i,Ca_in in enumerate(All_Ca):
    for j,Te in enumerate(All_T):
        
        erro=1
        tol=10**-8

        #Agora, vamos implementar um método para resolver as equações do problema

        x=Ca_in    #Primeiro, definimos as variáveis iterativas iguais as condições de entrada
        y=Te            
        it=0

        while erro>tol:     #fazendo iterações até que o erro seja menor do que a tolerância estabelecida
            it=it+1
            # print("Iteração= ",it)
            # print("Erro=",erro)
            #Agora, utilizando as funções definidas no topo do programa, calculamos os termos de cada iteração

            ui=u(x, y, Ca_in, tau, k0, E, R)
            vi=v(x, y, DH, V, k0, E, R, rho, cp, Qdot, UA, Tc, Te)
            
            dui_dx=dudx(y, tau, k0, E, R)
            dui_dy=dudy(x, y, tau, k0, E, R)
            dvi_dx=dvdx(y, DH, V, k0, E, R)
            dvi_dy=dvdy(x, y, DH, V, k0, E, R, rho, cp, Qdot, UA)
            
            DetJ=Det_J(dui_dx, dui_dy , dvi_dx , dvi_dy )
            
            x_new= x - 1/DetJ*(ui*dvi_dy-vi*dui_dy)
            y_new= y -1/DetJ*(vi*dui_dx-ui*dvi_dx)
            
            Rx=x-x_new # Calculamos a diferença entre a iteração atual e a anterior, para analizar a convergência
            Ry=y-y_new
            
            erro=np.max([np.abs(Rx),np.abs(Ry)])
            
            x=np.copy(x_new)
            y=np.copy(y_new)
            
            if it>1000:
                break
            

        # print("Iteração= ",it)
        # print("Erro=",erro)
        Ca_out[i,j]=x
        Ca_out_norm[i,j]=x/All_Ca[i]
        T_out[i,j]=y
        
#%%




x_vec=All_Ca
y_vec=All_T

plt.figure(dpi=300) 
plt.title("Concentração")
# plt.pcolormesh(x_vec,y_vec,Temp[:,:,i],cmap='plasma',vmin=0, vmax=1)
plt.contourf(x_vec,y_vec,np.transpose(Ca_out[:,:]),levels=10,cmap='viridis')
plt.colorbar()
plt.show()

plt.figure(dpi=300) 
plt.title("Concentração normalizada")
# plt.pcolormesh(x_vec,y_vec,Temp[:,:,i],cmap='plasma',vmin=0, vmax=1)
plt.contourf(x_vec,y_vec,np.transpose(Ca_out_norm[:,:]),levels=10,cmap='viridis')
plt.colorbar()
plt.show()


plt.figure(dpi=300) 
plt.title("Temperatura")
# plt.pcolormesh(x_vec,y_vec,Temp[:,:,i],cmap='plasma',vmin=0, vmax=1)
plt.contourf(x_vec,y_vec,np.transpose(T_out[:,:]),levels=10,cmap='plasma',vmin=250,vmax=500)
plt.colorbar()
plt.show()

"""
Comentários sobre os resultados do programa

O esperado era que se formassem três regiões nos gráficos de concentração e 
temperatura. Essa visualização fica mais fácil a partir do gráfico de 
concentração normalizada, que é calculada dividindo a concentração final pela
concentração inicial. Assim, quando a concentração normalizada é 1, nenhuma 
reação ocorreu, e quando ela é zero, a concentração de saída é muito pequena
em relação a concentração de entrada, indicando que a maior parte da espécie A
reagiu. No gráfico de concentração normalizada observa-se três regiões, uma em
que a concentração normalizada é sempre 1, uma em que é sempre zero, e uma 
terceira região de transição, em que apenas uma parcela da espécie A reagiu no
reator. 

Esse mesmo efeito pode ser observado no gráfico de temperatura, em que quando
a temperatura de entrada é muito pequena, a temperatura de saída se mantêm a 
mesma. Para maiores temperaturas de entrada, observa-se um aumento da temperatura
de saída, indicando que a reação está ocorrendo e liberando calor. Tambem pode-se
observar que para maiores temperaturas de entrada, a temperatura de saída
depende da concentração, de modo que para maiores concentrações iniciais 
observa-se uma maior temperatura de saída. 

"""

































