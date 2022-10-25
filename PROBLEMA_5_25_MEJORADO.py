
import numpy as np
import matplotlib.pyplot as plt

""" Problema 5.25, Transferencia de calor, Cengel"""


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

#DATOS DEL PROBLEMA
TL=200
nodos=6
L=0.05
dx=L/(nodos-1)
malla=np.linspace(0,L,nodos)
k=180
h=25
T_amb=25
Ts=290 #[K]
e=0.9
espesor=0.01
sigma=5.67*10**-8

# Datos trigonométricos

tan_theta=(espesor/2)/L
theta=np.arctan(tan_theta)
cos_theta=np.cos(theta)
sin_theta=np.sin(theta)

# Constantes

C=(h*dx**2)/(k*L*sin_theta)
D=(e*sigma*dx**2)/(k*L*sin_theta)
CT1=k*tan_theta
CT2=(h*dx)/cos_theta
CT3=(e*sigma*dx)/cos_theta

# Variables
A=np.zeros(nodos)
A[0]=(1+(dx/(2*L)))
for i in range(1,nodos):
    A[i]=(1-(i-0.5)*(dx/L))
    
B=np.zeros(nodos)
B[0]=(1-(dx/(2*L)))
for i in range(1,nodos):
    B[i]=(1-(i+0.5)*(dx/L))
    
# INICIALIZAR VECTORES
a=np.zeros(nodos)  # LOWER
b=np.zeros(nodos)  # MAIN
c=np.zeros(nodos)  # UPPER
d=np.zeros(nodos)  # RIGHT HAND

# ensamblar [a] 
a[0]=0
a[-1]=CT1

for i in range(1,nodos-1):
    a[i]=A[i]
    
# ensamblar [c]

c[0]=0
c[-1]=0

for i in range(1, nodos-1):
    c[i]=B[i]
    
# ensamblar [b]

b[0]=1
b[-1]=-(CT1+CT2)

for i in range(1, nodos-1):
    b[i]=-(A[i]+B[i]+C)
    



# VECTOR DE TEMPERATURAS SUPUESTAS PARA LA PRIMERA ITERACION
iter=0
tol=0.0000001
error=100

T_sup=np.ones(nodos)*100


while error>tol:
    T_sup[0]=TL

    #print(T_sup)  # con este vector, se calcula el primer vector [d]

    d_init=np.zeros(nodos)
    d_init[0]=T_sup[0]
    d_init[-1]=-CT3*(Ts**4-(T_sup[-1]+273)**4)-CT2*T_amb
    for i in range(1, nodos-1):
        d_init[i]=D*(Ts**4-(T_sup[i]+273)**4)-C*T_amb

    d=np.copy(d_init)

    T = TDMA(a,b,c,d)
    error=np.amax(np.abs(np.array(T_sup-T)))
    iter=iter+1
    T_sup=np.copy(T)
    
    
print(T)
print(iter)


