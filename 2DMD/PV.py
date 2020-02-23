import numpy as np
import matplotlib.pyplot as plt

N=900
R=8.3144598e-3 # J/mol/K Gas constant
T=300 # K
V=np.array([1/0.1, 1/0.4, 1/0.7, 1/1.0,1/1.5,1/3,1/10])
P1=np.array([0.25, 1.05,2.0,2.99,4.83,12.22,216.79])
P2=np.array([0.26,0.97,1.69,2.26,3.36,6.1,44.21])
P3=np.array([0.26,1.06,1.77,2.55,3.77,7.43,25.34])
P=R*T/V

plt.plot(V,P,'r',label='Ideal')
plt.plot(V,P1,'b',label='Repulsive')
plt.plot(V,P2,'g',label='Attractive')
plt.plot(V,P3,'k',label='Non-interacting')
plt.xlabel('Molar volume nm$^{2}$')
plt.ylabel('Pressure (kJ/mol/nm$^2$)')
plt.xlim(0.1,1)
plt.legend()
plt.show()
plt.semilogy(V,P,'r',label='Ideal')
plt.semilogy(V,P1,'b',label='Repulsive')
plt.semilogy(V,P2,'g',label='Attractive')
plt.semilogy(V,P3,'k',label='Non-interacting')
plt.xlabel('Molar volume nm$^{2}$')
plt.ylabel('Pressure (kJ/mol/nm$^2$)')
#plt.xlim(0.1,1)
plt.legend()
plt.show()
