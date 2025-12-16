"""
Created on Wed Oct 29 17:31:30 2025

@author: olafviggiano
"""
import numpy as np
import matplotlib.pyplot as plt


def thomas_algorithm(e, f, g, d):
    #Algoritmo de Thomas para solução de sistemas lineares tridiagonais

    n = len(f)
    
    e = np.array(e, dtype=float)
    f = np.array(f, dtype=float)
    g = np.array(g, dtype=float)
    d = np.array(d, dtype=float)
    
    # Forward elimination
    for i in range(1, n):
        w = e[i] / f[i - 1]
        f[i] = f[i] - w * g[i - 1]
        d[i] = d[i] - w * d[i - 1]
    
    # Back substitution
    x = np.zeros(n)
    x[-1] = d[-1] / f[-1]
    
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - g[i] * x[i + 1]) / f[i]
    
    return x

def T_ana(x,Fo):
    #Função para cálculo do Theta analítico
    lambd_n=[0.5643,3.2510,6.3395,9.4625,12.5947,15.7307,18.8685,22.0074,25.1469,28.2870]
    Theta=0
    for n in range(10):
        Theta=Theta+((4*np.sin(lambd_n[n]))/(2*lambd_n[n]+np.sin(2*lambd_n[n]))*np.cos(lambd_n[n]*x)*np.exp(-(lambd_n[n]**2)*Fo))
    return Theta


#Definição dos parâmetros do problema, valores retirados da internet
L=0.5 #Comprimento do reator
N=100 #numero de pontos da malha
dt=0.1 #Passo de tempo

dx=L/(N-1) #tamanho da malha

T_inf=300 #Temperatura do fluido de resfriamento
T0=900 #Temperatura inicial do reator

k=2.8
rho=10970
Cp=270
h=200

alpha=k/(rho*Cp)


x=np.linspace(0,L,N) #Vetor x da malha

# for ti in t:
#     T=T_ana(x, ti)
#     plt.plot(x,T, label='tempo='+str(ti))
#     plt.legend()
#     plt.show()
    
    

g=np.zeros(N)   #Definição dos vetores que representam o sistema linear tridiagonal
f=np.zeros(N)
e=np.zeros(N)




Fo=alpha*dt/(dx**2) #Cálculo do Número de Fourier numérico
Bi=h*dx/k           #Cálculo do Número de Biot numérico

All_Fo_ana=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5] #lista para variar o Fourier analítico

#Sem Geração

plt.figure(dpi=300)
for Fo_ana in All_Fo_ana:
    
    #Primeiro definimos o sistema linear tridiagonal (g,f,e)
    g[:]=-Fo
    g[0]=-2*Fo
    
    e[:]=-Fo
    e[-1]=-2*Fo
    
    f[:]=(1+2*Fo)
    f[-1]=(1+2*Fo+2*Fo*Bi)
    
    #Agora definimos um vetor para salvar a temperatura em cada ponto da malha
    T=np.zeros_like(x)
    #Definimos a malha inteira no momento incial com temperatura inicial
    T[:]=T0
    
    q=0 #Sem geração neste caso
    A=q/(rho*Cp)
    
    tempo=0
    t_fin=Fo_ana*L**2/alpha #Tempo final da simulação
    
    #loop para calcular a temperatura para cada passo de tempo
    while tempo<t_fin:
        tempo=tempo+dt
        d=T+A
        d[-1]=T[-1]+Fo*(2*Bi*T_inf+q*dx**2/k)
        T=thomas_algorithm(e, f, g, d)
        
    #Plot da temperatura para cada número de Fourier numérico vs analítico
    plt.plot(x/L,(T-T_inf)/(T0-T_inf), label='Fo='+str(Fo_ana))
    plt.plot(x/L,T_ana(x/L, Fo_ana),"--",c="black",lw=1)

plt.title("$\Theta$ vs x/L")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.xlabel("x/L")
plt.ylabel("$\Theta$")
plt.show()


#%%
#Com Geração


Fo_ana=1 #Considerando Fo=1
q_vec=np.linspace(0,10e5,50) #Vetor de q para ver o limite assintotico q-->0



T_center=np.zeros_like(q_vec) #Vetor para armazenar a temperatura no centro

for i,q in enumerate(q_vec):
    
    g[:]=-Fo
    g[0]=-2*Fo
    
    e[:]=-Fo
    e[-1]=-2*Fo
    
    f[:]=(1+2*Fo)
    f[-1]=(1+2*Fo+2*Fo*Bi)
    
    T=np.zeros_like(x)
    T[:]=T0
       
    
    A=q/(rho*Cp)

    tempo=0
    t_fin=Fo_ana*L**2/alpha 
    
    while tempo<t_fin:
        tempo=tempo+dt
        d=T+A
        d[-1]=T[-1]+Fo*(2*Bi*T_inf+q*dx**2/k)
        T=thomas_algorithm(e, f, g, d)
    
    T_center[i]=T[0]
    print("q=",q)
    
#Vetor para armazenar a temperatura analítica no centro
T_center_ana=np.zeros_like(T_center)
T_center_ana[:]=T_ana(0,Fo_ana)

plt.figure(dpi=300)
plt.title("$\Theta$ no centro vs $\dot q$")
plt.plot(q_vec,(T_center-T_inf)/(T0-T_inf),label="$\Theta(\dot q)$")
plt.plot(q_vec,T_center_ana,"--",c="black",label="Limite assintótico $\dot q =0$")
plt.legend(loc='upper left')
plt.tight_layout()
plt.xlabel("$\dot q$")
plt.ylabel("$\Theta_{centro}$")
plt.show()


#%%

#Agora, vamos analisar a diferença entre os resultados com e sem geração
# variando o número de Fourier, ou seja, o tempo adimensional do problema

Fo_ana=20 #Vamos analisar até o regime permanente, atingido em Fo~15

T_center=[]
T_center_ana=[]
Fo_vec=[]

g[:]=-Fo
g[0]=-2*Fo

e[:]=-Fo
e[-1]=-2*Fo

f[:]=(1+2*Fo)
f[-1]=(1+2*Fo+2*Fo*Bi)

T=np.zeros_like(x)
T[:]=T0

q=10e5
A=q/(rho*Cp)

tempo=0
t_fin=Fo_ana*L**2/alpha 

while tempo<t_fin:
    tempo=tempo+dt
    d=T+A
    d[-1]=T[-1]+Fo*(2*Bi*T_inf+q*dx**2/k)
    T=thomas_algorithm(e, f, g, d)
    T_center.append((T[0]-T_inf)/(T0-T_inf))
    T_center_ana.append(T_ana(0,tempo*alpha/L**2))
    Fo_vec.append(tempo*alpha/L**2)

# T_center[i]=T[0]
# T_center_ana[i]=T_ana(0,Fo_ana)
print("Fo=",Fo_ana)

plt.figure(dpi=300)
plt.title("$\Theta$ no centro vs Fourier")
plt.plot(Fo_vec,T_center,label="$\Theta$ com geração")
plt.plot(Fo_vec,T_center_ana,"--",c="black",label="$\Theta$ sem geração")
plt.legend(loc='upper right')
plt.tight_layout()
plt.xlabel("$Fo$")
plt.ylabel("$\Theta_{centro}$")
plt.show()
    
    
    
    
    
    





