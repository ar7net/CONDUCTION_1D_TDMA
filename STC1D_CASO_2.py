
import numpy as np
import matplotlib.pyplot as plt

""" CASO 2: Formulación TDMA para conducción unidimensional en estado estable
    con generación (e), pared plana, temperaturas fijas en los extremos"""


def TDMA(a,b,c,d):
    
    N = len(a)
    cp = np.zeros(N,dtype='float64') # store tranformed c or c'
    dp = np.zeros(N,dtype='float64') # store transformed d or d'
    X = np.zeros(N,dtype='float64') # store unknown coefficients
    
    # Perform Forward Sweep
    # Equation 1 indexed as 0 in python
    cp[0] = c[0]/b[0]  
    dp[0] = d[0]/b[0]
    # Equation 2, ..., N (indexed 1 - N-1 in Python)
    for i in np.arange(1,(N),1):
        dnum = b[i] - a[i]*cp[i-1]
        cp[i] = c[i]/dnum
        dp[i] = (d[i]-a[i]*dp[i-1])/dnum
    
    # Perform Back Substitution
    X[(N-1)] = dp[N-1]  # Obtain last xn 

    for i in np.arange((N-2),-1,-1):  # use x[i+1] to obtain x[i]
        X[i] = (dp[i]) - (cp[i])*(X[i+1])
    
    return(X)

#conducción 1D estado estable. Temperaturas fijas en extremos, sin generación
TL=0
TR=100
nodos=100
L=0.1
dx=L/(nodos-1)
malla=np.linspace(0,L,nodos)
e=5000000
k=48

a=np.zeros(nodos)
b=np.zeros(nodos)
d=np.zeros(nodos)

# ensamblar [a] y [c]
a[0]=0
a[-1]=0

for i in range(1,nodos-1):
    a[i]=1./(dx**2)

c=np.copy(a)


# ensamblar [b] diagonal principal

b[0]=1
b[-1]=1

for i in range(1, nodos-1):
    b[i]=-2./(dx**2)


# ensamblar vector de términos del lado derecho
d[0]=TL
d[-1]=TR

for i in range(1, nodos-1):
    d[i]=-e/k

# LLamar a la función TDMA
T = TDMA(a,b,c,d)
print(T)

plt.plot(malla,T)
plt.show()

