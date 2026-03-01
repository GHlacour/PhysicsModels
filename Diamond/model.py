import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import curve_fit

def einstein(T,TE):
    R=8.31
#    if TE<1:
#       TE=1
#    T[0]=1
#    Cv=3*R*(TE**2/T**2)*np.exp(TE/T)/((np.exp(TE/T)-1)**2)
    Cv=3*R*(TE**2/T**2)*np.exp(-TE/T)/((np.exp(-TE/T)-1)**2)

    return Cv

def linear(T,TL):
    Cv=T/TL
    return Cv

def debyefunction(x):
    if x<0:
       x=0
#    y=x**4*np.exp(x)/((np.exp(x)-1)**2)
    y=x**4*np.exp(-x)/((np.exp(-x)-1)**2)
    return y

def debye(T,TD):
    R=8.31
    xD=TD/T
    Cv=np.zeros(T.size)
    for i in range(1,T.size):
       # Cv[i]=3*R*(T[i]**3/TD**3)*integrate.quad(lambda x: debyefunction(x),0,xD[i])
       C=integrate.quad(lambda x: debyefunction(x),0,xD[i])
       Cv[i]=9*R*(T[i]**3/TD**3)*C[0]
#       if (Cv[i]==0):
#           Cv[i]=0.0001
    return Cv

def de_function(x,TD):
    z=debye(np.array([x,x]),TD)/np.array([x,x])
    y=z[1]
    if x==0.0:
       y=0       
    if y!=y:
       y=0
#    print(y)
    return y

def debye_entropy(T,TD):
    S=integrate.quad(lambda x: de_function(x,TD),0,T)
    return S

def make_fit(T,Cv,compound):
   # Find the best Debye Temperature 
   popt,covs=curve_fit(debye,T,Cv,p0=1000)
   dCv=debye(T,popt)
   print('Debye Temperature: '+str(popt))

   # Find the best Einstein Temperature
   popt,covs=curve_fit(einstein,T,Cv,p0=1400)
   eCv=einstein(T,popt)
   print('Einstein Temperature: '+str(popt))

   plt.title(compound)
   plt.plot(T,Cv,'k',label='Exp.')
   plt.plot(T,dCv,'b-',label='Debye')
   plt.plot(T,eCv,'r:',label='Einstein')
   plt.xlabel('Temperature (K)')
   plt.ylabel('$C_V$ (J/mol K)')
   plt.legend()
   plt.show()

data=np.loadtxt('diamond.dat',skiprows=0)
T=data[:,0]
Cv=data[:,1]*4.184 # Conversion to Joule/mol/K
make_fit(T,Cv,'Diamond')

data=np.loadtxt('graphite.dat',skiprows=0)
T=data[:,0]
Cv=data[:,1]*1 # Keep in Joule/mol/K
make_fit(T,Cv,'Graphite')

print(debye(np.array([300.0]),1850.7))
print('Diamond  '+str(debye_entropy(298.0,1850.7)))
print('Graphite '+str(debye_entropy(298.0,1600)))
