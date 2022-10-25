import numpy as np


TL=200
L=0.05
n=3
dx=L/(n-1)
e=0.9
sigma=5.67*10**-8
Ts=290
Tamb=25
k=25
h=28
diam=0.005

area=(np.pi*diam**2)/4
p=np.pi*diam

C1=k*area/dx
C2=h*(p*dx)
C3=e*sigma*(p*dx)

print(C1, C2, C3)