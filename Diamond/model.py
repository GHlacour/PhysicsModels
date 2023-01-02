import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit

def einstein(T,TE):
    R=8.31
    if TE<1:
       TE=1
    T[0]=1
#    Cv=3*R*(TE**2/T**2)*np.exp(TE/T)/((np.exp(TE/T)-1)**2)
    Cv=3*R*(TE**2/T**2)*np.exp(TE/T)/((np.exp(TE/T)-1)**2)

    Cv[0]=0
    return Cv

def linear(T,TL):
    Cv=T/TL
    return Cv

def debyefunction(x):
    y=x**4*np.exp(x)/((np.exp(x)-1)**2)
    return y

def debye(T,TD):
    R=8.31
    xD=TD/T
    Cv=np.zeros(T.size)
    for i in range(1,T.size):
       # Cv[i]=3*R*(T[i]**3/TD**3)*integrate.quad(lambda x: debyefunction(x),0,xD[i])
       C=integrate.quad(lambda x: debyefunction(x),0,xD[i])
       Cv[i]=9*R*(T[i]**3/TD**3)*C[0]
    return Cv


data=np.loadtxt('Data.dat',skiprows=0)

T=data[:,0]
Cv=data[:,1]*4.184 # Conversion to Joule/mol/K

dCv=debye(T,1800)
eCv=einstein(T,1400)

# Find the best Debye Temperature 
popt,covs=curve_fit(debye,T,Cv,p0=1000)
dCv=debye(T,popt)
print('Debye Temperature: '+str(popt))

# Find the best Einstein Temperature
popt,covs=curve_fit(einstein,T,Cv,p0=1400)
eCv=einstein(T,popt)
print('Einstein Temperature: '+str(popt))

plt.plot(T,Cv,'k')
plt.plot(T,dCv,'b-')
plt.plot(T,eCv,'r:')
plt.show()
