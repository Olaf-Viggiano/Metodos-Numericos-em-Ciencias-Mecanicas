#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  4 17:02:12 2025

@author: olafviggiano
"""

import numpy as np
import matplotlib.pyplot as plt





def p1(x,c):       #Função para polinômio definido previamente
    return (x-1)*(x-3)*(x-5)*(x-7)*(x-9)

def pn(x,c):    #Função para calcular um polinômio de ordem n 
                #com coeficientes dados pelo usuário
    p_out=0
    for i,coef in enumerate(c):
        p_out=p_out+(x**i*coef)
        
    return p_out

#%%

tol=1e-5

print("------------------------------------------------------------------------------------------")
print("Programa para solução de raízes de polinômios utilizando os métodos da secante e de Muller")
print("------------------------------------------------------------------------------------------")
print()
print("Vamos encontrar as raízes do polinômio proposto: f(x)=(x-1)*(x-3)*(x-5)*(x-7)*(x-9)")
# print("Escolha o polinômio que você deseja encontrar as raízes:")
# print("[0]:polinômio formado por coeficientes dados pelo usuário")
# print("[1]:")
print()

# choice=int(input("Digite o número da opção: "))
# while (choice!=0 and choice!=1):
#     choice=int(input("Digite uma opção válida: "))
    
# print()
    
c=[]

# if choice==0:
#     func=pn
#     print("Vamos definir os coeficientes do polinômio desejado")
#     order=int(input("Digite a ordem do polinômio: "))
#     print()
    
#     for i in range(order+1):
#         c.append(float(input(f"Digite o coeficiente do termo de {i}a ordem: ")))
#     print("Coeficentes do polinômio: ",c)
    
    
# if choice==1:
#     func=p1
#     order=5

func=p1
order=5

print()
    
print(f"Escolha {order} valores de chute inicial, para encontrar as {order} raízes do polinômio")
xi=[]
for i in range(order):
    xi.append(float(input(f"{i+1}° chute: ")))
    
print("Valores de chute iniciais: ", xi)


print("Vamos calcular as raízes utilizando os métodos de Muller e da Secante para comparar a convergência dos dois métodos") 
print()
input("Pressione enter para continuar")   
print("------------------------------------------------------------------------------------------")
print("Calculando raízes...")




#%%
dx=0.05

Max_iter=30


x_out_sec=np.zeros((Max_iter,order))
x_out_mul=np.zeros((Max_iter,order))
x_out_sec[0,:]=xi
x_out_mul[0,:]=xi



print()
print("Método da Secante")
print()

for i,x0 in enumerate(xi):
    x1=x0+dx
    dif=2*tol
    k=0
    while dif>tol and abs(func(x0,c))>1e-15 and k<Max_iter-1:
        x_next=x1-(func(x1,c)*(x1-x0))/(func(x1,c)-func(x0,c))
        dif=np.abs(x_next-x1)
        x0=x1
        x1=x_next
        k=k+1
        x_out_sec[k,i]=x1
    
    print("Raiz encontrada=",x1)
    
# print(x_out_sec)    


print()
print("Método de Muller")
print()

for i,x0 in enumerate(xi):
    x1=x0+dx
    x2=x1+dx
    k=0
    dif=2*tol
    while dif>tol and abs(func(x0,c))>1e-15 and k<Max_iter:
        
        h0=x1-x0
        h1=x2-x1
        delta_0=(func(x1,c)-func(x0,c))/h0
        delta_1=(func(x2,c)-func(x1,c))/h1
        
        a=(delta_1-delta_0)/(h1+h0)
        b=a*h1+delta_1
        c1=func(x2,c)
        
        disc = np.sqrt(b**2 - 4*a*c1)
        denom = b + disc if abs(b + disc) > abs(b - disc) else b - disc 
        
        #Na equação da parabola, temos duas raízes, referentes ao sinal 
        #positivo ou negativo da raiz da equação. A linha acima escolhe a maior
        #das duas opções para evitar a divisão por um número muito pequeno,
        #que estava levando o método a divergir nas raízes 3 e 7 do polinomio
        #proposto. 
        
        x3 = x2 - (2*c1) / denom
        # x3=x2-(2*c)/(b+np.sqrt(b**2-(4*a*c)))
        
        dif=np.abs(x3-x0)
        # print(dif)
        # print("x3",x3)
        # print("Erro",x3-x2)
        # print("iteração",k)
        x0=x3
        x1=x0+dx
        x2=x1+dx
        x_out_mul[k,i]=x0
        k=k+1
        
        # print(x0)
    print("Raiz encontrada=",x2)
# print(x_out_mul)
    
    
#%%

#Agora, plotando os erros relativos para comparar a convergência dos dois métodos

erro_sec=np.zeros((Max_iter-1,order))
erro_mul=np.zeros((Max_iter-1,order))

for j in range(order):
    for i in range(Max_iter-1):
        if x_out_sec[i+1,j]!=0:
            erro_sec[i,j]=abs((x_out_sec[i+1,j]-x_out_sec[i,j]))
        if x_out_mul[i+1,j]!=0:    
            erro_mul[i,j]=abs(x_out_mul[i+1,j]-x_out_mul[i,j])
        
raiz=[1,3,5,7,9]
        
for i in range(order): #para cada raiz, plotar um gráfico comparando
    raiz_sec=x_out_sec[np.max(np.nonzero(x_out_sec[:,i]))-1,i]
    raiz_mul=x_out_mul[np.max(np.nonzero(x_out_mul[:,i]))-1,i]
    plt.figure(dpi=100)
    plt.title(f"Erro relativo por iteração. Raiz={raiz[i]}")
    # plt.figtext(0.5, -0.08, f"Raiz encontrada método da secante={raiz_sec}", 
    #         wrap=True, horizontalalignment='center', fontsize=10)
    # plt.figtext(0.5, -0.03, f"Raiz encontrada método de Muller={raiz_mul}", 
    #         wrap=True, horizontalalignment='center', fontsize=10)
    plt.plot(erro_sec[erro_sec[:,i] != 0,i],"o-",label="Método da Secante")
    plt.plot(erro_mul[erro_mul[:,i] != 0 ,i],"o-",label="Método de Muller")
    plt.legend()
    plt.ylabel("Erro relativo")
    plt.xlabel("Iteração")
    plt.yscale("log")


plt.show()



























    

