#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 16:38:00 2025

@author: olafviggiano
"""

#Considerando o exemplo dado em sala para uma função na forma F(x,y)=1/2 x•A•x-b•x
#Assim, podemos calcular o gradiente de f como gradf=A•x-b

import numpy as np
import matplotlib.pyplot as plt
import random


N=2

A = np.array([[4, -1],
               [ -1, 3]])

b = np.array([1, 2])

x_k = np.array([2, 1])

plt.figure(dpi=300)


x = np.arange(-0.5, 1.5, 0.1)
y = np.arange(0, 2, 0.1)

X, Y = np.meshgrid(x, y)


path_AM=[x_k.copy()]

print("Método do aclive máximo")
print()


#Método do Aclive máximo
for i in range(100):

    grad_f=np.dot(A,x_k)-b
    if np.max(np.abs(grad_f))<1e-9:
        break
    
    alpha_k=np.dot(np.transpose(grad_f),grad_f)/(np.dot(np.dot(np.transpose(grad_f),A),grad_f))
    
    x_k=x_k-alpha_k*grad_f
    print(i,x_k)
    path_AM=np.vstack((path_AM,x_k))

#Método dos gradientes conjugados - Fletcher Reeves
print()
print("Método dos gradientes conjugados")
print()


A = np.array([[4, -1],
               [ -1, 3]])

b = np.array([1, 2])

x_k = np.array([2, 1])

# plt.figure(dpi=300)
# plt.contour(X, Y, F, levels=50, cmap="plasma")
# plt.colorbar()
# plt.plot(x_k[0],x_k[1],'o')

path_GC = [x_k.copy()]
r_k = b - A @ x_k
p_k = r_k.copy()

for i in range(100):
    alpha_k = (r_k.T @ r_k) / (p_k.T @ A @ p_k)

    x_k = x_k + alpha_k * p_k
    path_GC.append(x_k.copy())
    print(i, x_k)
    # plt.plot(x_k[0],x_k[1],'o')
    r_next = r_k - alpha_k * (A @ p_k)

    if np.max(np.abs(r_next)) < 1e-9:
        break

    beta_k = (r_next.T @ r_next) / (r_k.T @ r_k)
    p_k = r_next + beta_k * p_k
    r_k = r_next

path_GC = np.array(path_GC)




x = np.linspace(-1, 3, 100)
y = np.linspace(-1, 3, 100)
X, Y = np.meshgrid(x, y)
Z = 0.5*(A[0,0]*X**2 + 2*A[0,1]*X*Y + A[1,1]*Y**2) - (b[0]*X + b[1]*Y)

ax = plt.figure(dpi=300).add_subplot(projection='3d')
plt.title("Visualização 3D da função quadrática")
ax.plot_surface(X,Y, np.transpose(Z), cmap="plasma", lw=0.5, rstride=1, cstride=1,)


plt.show()

plt.figure(figsize=(6,5), dpi=150)

plt.contour(X, Y, Z, levels=40, colors="black",linewidths=1)

plt.plot(path_AM[:,0], path_AM[:,1], "-o", color="black", markersize=5,label="Aclive Máximo")
plt.plot(path_GC[:,0], path_GC[:,1], "-o", color="red", markersize=3,lw=1,label="Gradiente Conjugado")


plt.legend()
plt.title("Comparação do caminho de cada método")
plt.xlabel("x")
plt.ylabel("y")
plt.gca().set_aspect("equal", adjustable="box")
plt.show()

#%%
#Método de Newton: Este método pode ser aplicado para qualquer tipo de função,
#Diferentemente do método do gradiente conjungado mostrado em sala, que pode
#ser aplicado somente para funções quadráticas. No caso do método de Newton,
#vamos utilizar uma função qualquer, e compara-lo com o método de Levenberg-Marquardt


#Primeiro definimos uma funcão f(x,y), o seu gradiente e a sua hessiana

def f(xy, a, b, c, d):
    x, y = xy
    return 3 - a*x**2 - b*y**2 + c*x**4 + d*y**4

def grad_f(xy, a, b, c, d):
    x, y = xy
    return np.array([
        -2*a*x + 4*c*x**3,
        -2*b*y + 4*d*y**3])

def Hess_f(xy, a, b, c, d):
    x, y = xy
    return np.array([
        [-2*a + 12*c*x**2,           0.0],
        [0.0,               -2*b + 12*d*y**2]])

def Newton(grad, hess, x0, a, b, c, d):

    x = x0
    path_Newt=x.copy()
    for k in range(100):
        grad = grad_f(x, a, b, c, d)
        Hess = Hess_f(x, a, b, c, d)
        x_new = x -np.linalg.inv(Hess) @ grad

        x = x_new
        path_Newt=np.vstack((path_Newt,x))
        # stop if gradient is tiny
        if np.max(np.abs(grad)) < 1e-9:
            return x, np.array(path_Newt)

    print("Reached maximum iterations.")
    return x, np.array(path_Newt)

a=1
b=1
c=0.1
d=0.1

print()
print("Método de Newton")


x0 = np.array([1, 1.5])
x_f,path_Newt = Newton(grad_f, Hess_f, x0,a,b,c,d)

print(path_Newt)
print()

#Método de Levenberg-Marquardt

def Lev_Marq(grad, hess, x0, a, b, c, d):
    
    x = x0
    alpha=0.01
    decay=0.5
    path_Newt=x.copy()
    for k in range(100):
        grad = grad_f(x, a, b, c, d)
        Hess = Hess_f(x, a, b, c, d)
        
        if alpha>1:
            x_new=x-grad/alpha
        else:
            Hess_mod=Hess+alpha*np.eye(2)
            x_new = x -np.linalg.inv(Hess_mod) @ grad

        x = x_new
        path_Newt=np.vstack((path_Newt,x))
        # stop if gradient is tiny
        if np.max(np.abs(grad)) < 1e-9:
            return x, np.array(path_Newt)
        alpha=alpha*decay

    print("Reached maximum iterations.")
    return x, np.array(path_Newt)

print("Método de Levenberg-Marquardt")

x_f2,path_LM=Lev_Marq(grad_f, Hess_f, x0, a, b, c, d)
print(path_LM)


x = np.linspace(-4, 4, 100)
y = np.linspace(-4, 4, 100)


X, Y = np.meshgrid(x, y)

F = 3 - a*X**2 - b*Y**2 + c*X**4 + d*Y**4


ax = plt.figure(dpi=300).add_subplot(projection='3d')
ax.plot_surface(X,Y, np.transpose(F), cmap="plasma", lw=0.5, rstride=1, cstride=1,)
plt.title("Visualização 3D da função")

plt.show()

plt.figure(dpi=300)
plt.contour(X, Y, F, levels=20, colors="black",linewidths=0.5)
plt.xlim(-4,4)
plt.ylim(-4,4)

plt.plot(path_Newt[:,0], path_Newt[:,1], "-o",ms=4,lw=1, color="blue",label="Newton")

plt.plot(path_Newt[-1,0], path_Newt[-1,1], "*", markersize=7,color="blue")

plt.plot(path_LM[:,0], path_LM[:,1], "-o",ms=2,lw=0.5, color="red",label="Lev-Marq")

plt.plot(path_LM[-1,0], path_LM[-1,1], "*", markersize=7,color="red")
plt.gca().set_aspect("equal", adjustable="box")
plt.title("Comparação dos caminhos de cada método")
plt.legend()
plt.show()




























