import numpy as np
from numba import jit
import matplotlib.pyplot as plt
# These are subroutines for molecular dynamics

# Put atoms on a grid
def roster(N,r):
  q=np.zeros((N,2))
  k=1
  NN=int(np.sqrt(N))
  for i in range(NN):
    for j in range(NN):
      k=i*NN+j
      if k<N:
        q[k,0]=(j+0.5)/NN*r
        q[k,1]=(i+0.5)/NN*r
  return q

# Find potential energy
@jit(nopython=True)
def find_energy(q,R,N,Aw,Aa,B):
  E=0
  # Find interaction energy between praticles
  for ii in range(N):
    for jj in range(ii+1,N):
      dx=q[ii,0]-q[jj,0]
      dy=q[ii,1]-q[jj,1]
      dist=np.sqrt(dx**2+dy**2)
      E=E+Aa*np.exp(-dist/B)

  # Find wall-particle energy
  for ii in range(N):
    E=E+Aw*np.exp(-q[ii,0]/B)
    E=E+Aw*np.exp(-(R-q[ii,0])/B)
    E=E+Aw*np.exp(-q[ii,1]/B)
    E=E+Aw*np.exp(-(R-q[ii,1])/B)
  return E

# Scale momenta to match temperature
#@jit(nopython=True)
def scale_velocity(N,p,m,Temp):
  R=8.3144598e-3 # J/mol/K Gas constant
  Ekin=find_kinetic(p,m)
  Ekin0=2*N*R*Temp/2;
  # Scale momenta
  p=p*np.sqrt(Ekin0/Ekin)
  return p

# Calculate the kinetic energy
#@jit(nopython=True)
def find_kinetic(p,m):
  E=(np.dot(p[:,0],p[:,0])+np.dot(p[:,1],p[:,1]))/2/m
  return E

# Calculate the forces
@jit(nopython=True)
def find_forces(q,r,N,Aw,Aa,B):
  F=np.zeros((N,2))
  P=0
  # Forces between particles
  for ii in range(N):
    for jj in range(N):
      if ii!=jj:
        dx=q[ii,0]-q[jj,0]
        dy=q[ii,1]-q[jj,1]
        dist=np.sqrt(dx**2+dy**2)
        F[ii,0]=F[ii,0]+Aa*np.exp(-dist/B)*dx/dist/B
        F[ii,1]=F[ii,1]+Aa*np.exp(-dist/B)*dy/dist/B

  # Forces with walls
  for ii in range(N):
    F[ii,0]=F[ii,0]+Aw*np.exp(-q[ii,0]/B)/B
    F[ii,0]=F[ii,0]-Aw*np.exp(-(r-q[ii,0])/B)/B
    F[ii,1]=F[ii,1]+Aw*np.exp(-q[ii,1]/B)/B
    F[ii,1]=F[ii,1]-Aw*np.exp(-(r-q[ii,1])/B)/B
    P=P+Aw*(np.exp(-q[ii,0]/B)+Aw*np.exp(-(r-q[ii,0])/B)+np.exp(-q[ii,1]/B)+Aw*np.exp(-(r-q[ii,1])/B))
  return (F,P)

# Plot miscrostate
def plot_state(q,r):
  plt.figure(1)
  plt.clf()
  # Make particle 0 blue
  plt.plot(q[0,0],q[0,1],'bo')
  # Make all other particles red
  plt.plot(q[1:,0],q[1:,1],'ro')
  plt.axis([0,r,0,r])
  plt.gca().set_aspect('equal',adjustable='box')
  plt.title('Simulation Ensemble',fontsize=18)
  plt.xlabel('x-axis (nm)',fontsize=18)
  plt.ylabel('y-axis (nm)',fontsize=18)
  plt.draw()
  plt.pause(0.0001)

# Plot physical variables
def plot_variables(t,Epot,Ekin,P,R,V,N):
  plt.figure(2)
  plt.plot(t,Epot,t,Ekin,t,Epot+Ekin)
  plt.title('Energies',fontsize=18)
  plt.xlabel('Time (ps)',fontsize=18)
  plt.ylabel('Energy (kJ/mol)',fontsize=18)
  plt.legend(('Potential','Kinetic','Total'),fontsize=18)
  plt.draw()
  plt.figure(3)
  plt.subplot(2,2,1)
  plt.plot(t,P)
  plt.xlabel('Time (ps)',fontsize=18)
  plt.ylabel('Pressure (kJ/mol/nm$^2$)',fontsize=18)
  plt.subplot(2,2,2)
  T=Ekin/R
  plt.plot(t,T)
  plt.xlabel('Time (ps)',fontsize=18)
  plt.ylabel('Temperature (K)',fontsize=18)
  plt.subplot(2,2,3)
  PVNRT=P*V/N/R/T
  plt.plot(t,PVNRT) 
  plt.xlabel('Time (ps)',fontsize=18)
  plt.ylabel('PVNRT',fontsize=18)
  plt.subplot(2,2,4)
  plt.plot(t,Ekin+Epot)
  plt.xlabel('Time (ps)',fontsize=18)
  plt.ylabel('Energy (kJ/mol)',fontsize=18)
  plt.tight_layout()
  plt.draw()
  # Exclude initial timesteps from further analysis
  start=50 
  # Make histograms skipping #start first points
  plt.figure(4)
  plt.subplot(2,2,1)
  plt.hist(P[start:],100)
  plt.xlabel('Pressure (kJ/mol/nm$^2$)',fontsize=18)
  plt.subplot(2,2,2)
  plt.hist(T[start:],100)
  plt.xlabel('Temperature (K)',fontsize=18)
  plt.subplot(2,2,3)
  plt.hist(PVNRT[start:],100)
  plt.xlabel('PVNRT',fontsize=18)
  plt.subplot(2,2,4)
  E=Ekin+Epot
  plt.hist(E[start:],100)
  plt.xlabel('Energy (kJ/mol)',fontsize=18)
  plt.tight_layout()
  plt.show()

# Subroutine for plotting Maxwell-Boltzmann like diagram (2D not 3D!)
def velocity_histogram(ps):
  # Exclude initial timesteps from further analysis
  start=50
  p1=ps[:,0]*ps[:,0]+ps[:,1]*ps[:,1]
  # print(ps.size)
  p=np.sqrt(p1)
  plt.figure(2)
  plt.hist(p[start:],100)
  plt.xlabel('Momentum (u nm/ps)',fontsize=18)
  plt.ylabel('Occurences',fontsize=18)
  plt.show()
