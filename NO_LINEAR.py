import numpy as np


def TDMA(a,b,c,d):
    
    N = len(a)
    cp = np.zeros(N,dtype='float64') # store tranformed c or c'
    dp = np.zeros(N,dtype='float64') # store transformed d or d'
    X = np.zeros(N,dtype='float64') # store unknown coefficients
    cp[0] = c[0]/b[0]  
    dp[0] = d[0]/b[0]
    for i in np.arange(1,(N),1):
        dnum = b[i] - a[i]*cp[i-1]
        cp[i] = c[i]/dnum
        dp[i] = (d[i]-a[i]*dp[i-1])/dnum
    
    # Perform Back Substitution
    X[(N-1)] = dp[N-1]  # Obtain last xn 

    for i in np.arange((N-2),-1,-1):  # use x[i+1] to obtain x[i]
        X[i] = (dp[i]) - (cp[i])*(X[i+1])
    
    return(X)

L=0.04                          # Longitud de pared [m]
nodos=3                         # Número de nodos
dx=L/(nodos-1)                  # Delta x [m]
k=28                            # Conductividad térmica [W/mK]
alfa=12.5*10**-6                # Difusividad térmica [m^2/s]
TL=0                            # Temperatura izquierda fija [°C]
e=5*10**6                       # Generación [W/m^3]
Tamb=30                         # Temperatura del medio [°C]
h=45                            # Coeficiente de convección [W/m2K]
dt=15                           # Intervalo de tiempo [s]
f=(alfa*dt)/(dx**2)             # Número de Fourier
interv=10                       # Número de intervalos de tiempo
t=np.ones(interv)               # Malla temporal
t[0]=0                          # Tiempo inicial [s]

for i in range(1, interv):
    t[i]=dt*i

T_inicial=np.ones(nodos)*200       # Temperatura incial de pared
T_inicial[0]=TL
 
Ti=np.ones(nodos)*200       # Temperatura en el tiempo i+1
T_res=[]

# INICIALIZAR VECTORES TDMA
a=np.zeros(nodos)               # Lower
b=np.zeros(nodos)               # Main
c=np.zeros(nodos)               # Upper
d=np.zeros(nodos)               # Right hand

# CICLO DE ITERACIONES PARA ESTADO TRANSITORIO
for j in range(0, len(t)):
    Ti[0]=TL
    
    # Ensamblar [a]
    a[0]=0
    a[-1]=2*f

    for i in range(1, nodos-1):
        a[i]=f

    # Ensamblar [c]
    c[0]=0
    c[-1]=0

    for i in range(1, nodos-1):
        c[i]=f
        
    # Ensamblar [b]
    b[0]=1
    b[-1]=-(1+2*f+2*f*(h*dx/k))

    for i in range(1, nodos-1):
        b[i]=-(1+2*f)
        
    #Ensamblar [d]
    d[0]=TL
    d[-1]=-Ti[-1]-((2*f*h*dx*Tamb)/k)-(f*e*dx**2/k)
    for i in range(1, nodos-1):
        d[i]=-Ti[i]-(f*e*dx**2/k)
    
    # LLAMAR A LA FUNCIÓN TDMA
    T = TDMA(a,b,c,d)
    T_res.append(T)
    Ti=np.copy(T)
    
import pandas as pd

T_res=np.array(T_res)
T_res=np.vstack([T_inicial, T_res])
DF = pd.DataFrame(T_res, columns = ['T nodo 0 (°C)','T nodo 1 (°C)','T nodo 2 (°C)'])
pd.options.display.float_format="{:,.2f}".format
print(DF)
