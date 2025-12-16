# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

# Para resolver o problema, é necessário assumir alguns valores iniciais de concentração
# e volume dos reatores, que não foram especificados. 
# Vamos então, por motivos de simplificação, assumir que o volume dos reatores é
# 1 unidade de volume, e a concentração inicial de cada reator é igual a 1
# Como possuimos um sistema de equações acopladas, para resolver cada passo de
# tempo a partir de um método de Runge-Kutta de 4a ordem, precisamos calcular
# os fatores do método para cada equação simultaneamente, ou seja, primeiro
# calculamos todos os k1s, depois todos os k2s, e assim por diante, para obter
# o próximo passo de tempo para todas as equações simultaneamente

#Primeiro vamos definir os valores de concentrações iniciais e as vazoes Q.

Q01=5
c01=10

Q03=8
c03=20

Q15=3
Q12=3
Q31=1

Q23=1
Q25=1
Q24=1

Q34=8

Q54=2
Q44=11

Q55=2

c1_t=[1]
c2_t=[1]
c3_t=[1]
c4_t=[1]
c5_t=[1]

#Agora, vamos definir 5 funções para cada uma das equações diferenciais. 

def dc1dt (c1,c2,c3,c4,c5):
    Q01=5
    c01=10
    Q31=1
    Q12=3
    Q15=3
    return Q01*c01+Q31*c3-Q12*c1-Q15*c1

def dc2dt (c1,c2,c3,c4,c5):
    Q12=3
    Q23=1
    Q25=1
    Q24=1
    
    return Q12*c1-Q23*c2-Q24*c2-Q25*c2

def dc3dt (c1,c2,c3,c4,c5):
    Q03=8
    c03=20
    Q23=1
    Q31=1
    Q34=8
    
    return Q03*c03+Q23*c2-Q31*c3-Q34*c3

def dc4dt (c1,c2,c3,c4,c5):
    Q24=1
    Q34=8
    Q54=2
    Q44=11
    
    return Q24*c2+Q34*c3+Q54*c5-Q44*c4

def dc5dt (c1,c2,c3,c4,c5):
    Q55=2
    Q54=2
    Q15=3
    Q25=1
    
    return Q15*c1 +Q25*c2-Q54*c5-Q55*c5

#Agora, vamos implementar o método de Runge-Kutta de quarta ordem, calculando
#de forma simultanea todas as equações. 


h=0.01  #Passo de tempo
time_vec=np.arange(0,3,h)

c1=c1_t[0]
c2=c2_t[0]
c3=c3_t[0]
c4=c4_t[0]
c5=c5_t[0]

#%%

for i,t in enumerate(time_vec[1:]):

#primeiro calculamos todos os k1

    k1_c1=dc1dt(c1,c2,c3,c4,c5)
    k1_c2=dc2dt(c1,c2,c3,c4,c5) 
    k1_c3=dc3dt(c1,c2,c3,c4,c5)
    k1_c4=dc4dt(c1,c2,c3,c4,c5)
    k1_c5=dc5dt(c1,c2,c3,c4,c5)
    
    #Agora, com base nos k1s, calculamos todos os k2
    
    k2_c1=dc1dt(c1+1/2*k1_c1*h, c2+1/2*k1_c2*h,c3+1/2*k1_c3*h,c4+1/2*k1_c4*h,c5+1/2*k1_c5*h)
    k2_c2=dc2dt(c1+1/2*k1_c1*h, c2+1/2*k1_c2*h,c3+1/2*k1_c3*h,c4+1/2*k1_c4*h,c5+1/2*k1_c5*h)
    k2_c3=dc3dt(c1+1/2*k1_c1*h, c2+1/2*k1_c2*h,c3+1/2*k1_c3*h,c4+1/2*k1_c4*h,c5+1/2*k1_c5*h)
    k2_c4=dc4dt(c1+1/2*k1_c1*h, c2+1/2*k1_c2*h,c3+1/2*k1_c3*h,c4+1/2*k1_c4*h,c5+1/2*k1_c5*h)
    k2_c5=dc5dt(c1+1/2*k1_c1*h, c2+1/2*k1_c2*h,c3+1/2*k1_c3*h,c4+1/2*k1_c4*h,c5+1/2*k1_c5*h)
    
    #Agora, os k3
    
    k3_c1=dc1dt(c1+1/2*k2_c1*h, c2+1/2*k2_c2*h,c3+1/2*k2_c3*h,c4+1/2*k2_c4*h,c5+1/2*k2_c5*h)
    k3_c2=dc2dt(c1+1/2*k2_c1*h, c2+1/2*k2_c2*h,c3+1/2*k2_c3*h,c4+1/2*k2_c4*h,c5+1/2*k2_c5*h)
    k3_c3=dc3dt(c1+1/2*k2_c1*h, c2+1/2*k2_c2*h,c3+1/2*k2_c3*h,c4+1/2*k2_c4*h,c5+1/2*k2_c5*h)
    k3_c4=dc4dt(c1+1/2*k2_c1*h, c2+1/2*k2_c2*h,c3+1/2*k2_c3*h,c4+1/2*k2_c4*h,c5+1/2*k2_c5*h)
    k3_c5=dc5dt(c1+1/2*k2_c1*h, c2+1/2*k2_c2*h,c3+1/2*k2_c3*h,c4+1/2*k2_c4*h,c5+1/2*k2_c5*h)
    
    #Por fim, os k4
    
    k4_c1=dc1dt(c1+k3_c1*h, c2+k3_c2*h,c3+k3_c3*h,c4+k3_c4*h,c5+k3_c5*h)
    k4_c2=dc2dt(c1+k3_c1*h, c2+k3_c2*h,c3+k3_c3*h,c4+k3_c4*h,c5+k3_c5*h)
    k4_c3=dc3dt(c1+k3_c1*h, c2+k3_c2*h,c3+k3_c3*h,c4+k3_c4*h,c5+k3_c5*h)
    k4_c4=dc4dt(c1+k3_c1*h, c2+k3_c2*h,c3+k3_c3*h,c4+k3_c4*h,c5+k3_c5*h)
    k4_c5=dc5dt(c1+k3_c1*h, c2+k3_c2*h,c3+k3_c3*h,c4+k3_c4*h,c5+k3_c5*h)
    
    #Agora, calculamos as concentrações para o próximo passo de tempo, e salvamos
    #em uma lista para avaliar a evolução temporal das concentrações
        
    c1_new=c1+1/6*(k1_c1+2*k2_c1+2*k3_c1+k4_c1)*h
    c1_t.append(c1_new)
    c1=c1_new
        
    c2_new=c2+1/6*(k1_c2+2*k2_c2+2*k3_c2+k4_c2)*h
    c2_t.append(c2_new)
    c2=c2_new
    
    c3_new=c3+1/6*(k1_c3+2*k2_c3+2*k3_c3+k4_c3)*h
    c3_t.append(c3_new)
    c3=c3_new
    
    c4_new=c4+1/6*(k1_c4+2*k2_c4+2*k3_c4+k4_c4)*h
    c4_t.append(c4_new)
    c4=c4_new
    
    c5_new=c5+1/6*(k1_c5+2*k2_c5+2*k3_c5+k4_c5)*h
    c5_t.append(c5_new)
    c5=c5_new

#%%

plt.figure(dpi=300)
plt.title("Evolução temporal das concentrações em cada reator")
plt.plot(time_vec,c1_t,'-',label="c1")
plt.plot(time_vec,c2_t,'-',label="c2",lw=2,color="black")
plt.plot(time_vec,c3_t,'-',label="c3")
plt.plot(time_vec,c4_t,'-',label="c4")
plt.plot(time_vec,c5_t,'-.',label="c5")
plt.legend()
plt.xlabel("tempo(s)")
plt.ylabel("concentração")
plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    