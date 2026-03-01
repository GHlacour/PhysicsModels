import numpy as np
from numba import jit
from numba import njit
import matplotlib.pyplot as plt
# These are subroutines for molecular dynamics

# Put atoms on an initial grid
@njit
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

# Find distance PBC
@njit
def find_distance(q2,q1,box):
  dx=q2[0]-q1[0]
  dy=q2[1]-q1[1]
  if dx>=box/2.0:
    dx=dx-box
  elif dx<=-box/2.0:
    dx=dx+box
  if dy>=box/2.0:
    dy=dy-box
  elif dy<=-box/2.0:
    dy=dy+box
  dist=np.sqrt(dx**2+dy**2)
  return [dx,dy,dist]

# Find new positions and apply PBC, Euler
@njit
def update_positions(q,p,dt,m,F,r):
  q=q+p*dt/m+F/m*dt**2/2.0
  # Move particles into box if outside
  q=q-r*(q>=r)
  q=q+r*(q<0)
  return q

# Find new positions and apply PBC, Verlet
@njit
def update_positions_verlet(q,q0,dt,m,F,r):
  q1=q
  q=2*q-q0+F/m*dt**2
  # Move particles into box if outside
  q0=q1-r*(q>=r)
  q0=q1+r*(q<0)
  q=q-r*(q>=r)
  q=q+r*(q<0)
  return [q,q0]

# Update momentum
@njit
def update_momentum(p,Fold,F,dt):
  p=p+(Fold+F)/2.0*dt
  return p

# Find potential energy
@njit
def find_energy(q,r,N,Aa,B):
  E=0
  # Find interaction energy between particles
  for ii in range(N):
    # We only need to loop over pairs of different particles
    for jj in range(ii+1,N):
      # Find the distance with PBC
      [dx,dy,dist]=find_distance(q[ii,:],q[jj,:],r)
      if dist<r/2:
        E=E+Aa*np.exp(-dist/B)
  return E

# Scale momenta to match temperature
@njit
def scale_velocity(N,p,m,Temp):
  R=8.3144598e-3 # J/mol/K Gas constant
  Ekin=find_kinetic(p,m)
  Ekin0=2*N*R*Temp/2
# Scale momenta
  p=p*np.sqrt(Ekin0/Ekin)
  return p

# Calculate the kinetic energy
@njit
def find_kinetic(p, m):
    E = 0.0
    for i in range(p.shape[0]):
        E += p[i, 0]**2 + p[i, 1]**2
    return E / (2 * m)

# Calculate the forces and pressure
@njit
def find_forces(q,r,N,Aa,B):
  F=np.zeros((N,2))
  P=0
  # Forces between particles
  for ii in range(N):
    for jj in range(ii+1,N):
      if ii!=jj:
        # Find the distance with PBC
        [dx,dy,dist]=find_distance(q[ii,:],q[jj,:],r)
        # Find force on first particle
        if dist<r/2:
          F0=Aa*np.exp(-dist/B)/dist/B
          Fx=F0*dx
          Fy=F0*dy
          F[ii,0]=F[ii,0]+Fx
          F[ii,1]=F[ii,1]+Fy
          # Use Newtons law to find force on the other particle
          F[jj,0]=F[jj,0]-Fx
          F[jj,1]=F[jj,1]-Fy
          # Pressure deviation from ideal as given on page 647 of
          # Pathria and Beale Edition 3
          P+=(Fx*dx+Fy*dy)

  return (F,P)

# Plot miscrostate
def plot_state(q,r):
  plt.figure(1,figsize=(10,10))
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
  plt.pause(0.0000001)

# Plot physical variables
def plot_variables(t,Epot,Ekin,P,R,V,N):
  plt.figure(2,figsize=(10,10))
  plt.plot(t,Epot,t,Ekin,t,Epot+Ekin)
  plt.title('Energies',fontsize=18)
  plt.xlabel('Time (ps)',fontsize=18)
  plt.ylabel('Energy (kJ/mol)',fontsize=18)
  plt.legend(('Potential','Kinetic','Total'),fontsize=18)
  plt.draw()
  plt.figure(3,figsize=(10,10))
  plt.subplot(2,2,1)
  plt.plot(t,P)
  plt.xlabel('Time (ps)',fontsize=18)
  plt.ylabel('Pressure (kJ/mol/nm$^2$)',fontsize=18)
  plt.subplot(2,2,2)
  T=Ekin/R/N
  plt.plot(t,T)
  plt.xlabel('Time (ps)',fontsize=18)
  plt.ylabel('Temperature (K)',fontsize=18)
  plt.subplot(2,2,3)
  PVNRT=P*V/N/R/T
  plt.plot(t,PVNRT) 
  plt.xlabel('Time (ps)',fontsize=18)
  plt.ylabel('PV/NRT',fontsize=18)
  plt.subplot(2,2,4)
  plt.plot(t,Ekin+Epot)
  plt.xlabel('Time (ps)',fontsize=18)
  plt.ylabel('Energy (kJ/mol)',fontsize=18)
  plt.tight_layout()
  plt.draw()
  # Exclude initial timesteps from further analysis
  start=50 
  # Make histograms skipping #start first points
  plt.figure(4,figsize=(10,10))
  plt.subplot(2,2,1)
  plt.hist(P[start:],100)
  plt.xlabel('Pressure (kJ/mol/nm$^2$)',fontsize=18)
  mytext='<P>=%.2f $\pm$ %.2f kJ/mol/nm$^2$' %(np.mean(P[start:]),np.std(P[start:])) 
  plt.title(mytext)
  plt.subplot(2,2,2)
  plt.hist(T[start:],100)
  plt.xlabel('Temperature (K)',fontsize=18)
  mytext='<T>=%.2f $\pm$ %.2f K' %(np.mean(T[start:]),np.std(T[start:]))
  plt.title(mytext)
  plt.subplot(2,2,3)
  plt.hist(PVNRT[start:],100)
  plt.xlabel('PV/NRT',fontsize=18)
  mytext='<PV/NRT>=%.2f $\pm$ %.2f' %(np.mean(PVNRT[start:]),np.std(PVNRT[start:]))
  plt.title(mytext)
  plt.subplot(2,2,4)
  E=Ekin+Epot
  plt.hist(E[start:],100)
  plt.xlabel('Energy (kJ/mol)',fontsize=18)
  mytext='<E>=%.2f $\pm$ %.2f kJ/mol' %(np.mean(E[start:]),np.std(E[start:]))
  plt.title(mytext)
  plt.tight_layout()
  plt.show()

# Subroutine for plotting Maxwell-Boltzmann like diagram (2D not 3D!)
def velocity_histogram(ps):
  # Exclude initial timesteps from further analysis
  start=50
  p1=ps[:,0]*ps[:,0]+ps[:,1]*ps[:,1]
  # print(ps.size)
  p=np.sqrt(p1)
  plt.figure(2,figsize=(10,10))
  plt.hist(p[start:],100)
  plt.xlabel('Momentum (u nm/ps)',fontsize=18)
  plt.ylabel('Occurences',fontsize=18)
  mytext='<p>=%.2f $\pm$ %.2f u nm/ps' %(np.mean(p[start:]),np.std(p[start:]))
  plt.title(mytext)
  plt.show()
