import numpy as np
import matplotlib.pyplot as plt
import subroutines_v3 as sr
import time as dotime

# Main Program
# N number of particles
# r length of one side in nm
# Temp temperature in Kelvin
# times number of snapshots to calculate
# fixed_T 1 is fixed 0 is fluctuating
# fast number of snapshots between plot on screen 100 give nice animation
#   0 is no plots at all during simulation
def dynamics2D(N,r,Temp,times,fixed_T,fast):
  start = dotime.process_time()
  print('The number density is ' + str(N/r**2) + ' 1/nm^2.')
  m=18; # Mass in a.m.u.
  R=8.3144598e-3 # Gas constant kJ/mol/K
  dt=0.005 # Timestep in ps
  Aw=180000 # kJ/mol Particle-wall repulsion Aw*exp(-dq/B)
  Aa=180000 # kJ/mol Particle-Particle repulsion Aa*exp(-dq/B)
  B=0.025 # nm 
  C=0.0024 # kJ/mol/nm^6
  N=int(np.floor(np.sqrt(N))**2) # Make number of particels a square number
  q=sr.roster(N,r) 
  Epot=sr.find_energy(q,r,N,Aw,Aa,B,C)
  # Initialize momenta
  p=np.random.normal(loc=0,scale=1,size=(N,2))
  #  print(p)
  # Scale momenta to get right kinetci energy
  p=sr.scale_velocity(N,p,m,Temp)
  # Initialize arrays
  Ekin=np.zeros(times)
  Epot=np.zeros(times)
  T=np.zeros(times)
  t=np.arange(0,times)*dt
  P=np.zeros(times)
  Epot[0]=sr.find_energy(q,r,N,Aw,Aa,B,C)
  Ekin[0]=sr.find_kinetic(p,m)
  ps=p
  # Initial step
  (F,Pt)=sr.find_forces(q,r,N,Aw,Aa,B,C)
  q=q+p*dt/m+F/m*dt**2/2  
  # Loop over all times
  for time in range(1,times):
    if (fast>0 and time%fast==0):
      print('Timesteps passed ' + str(time) + '.')
    #  Update Positions
    q=q+p/m*dt+F/m*dt**2/2
    # Find new forces
    Fold=F
    [F,Pt]=sr.find_forces(q,r,N,Aw,Aa,B,C)
    # Update momentum
    p=p+(Fold+F)/2*dt
    if (time % 10==0):
      ps=np.append(ps,p,axis=0)
    if (fixed_T==1):
       x=1
       p=sr.scale_velocity(N,p,m,Temp)
    # Plot microstate
    if (fast>0 and time%fast==0):
      sr.plot_state(q,r)
    P[time]=Pt
    Epot[time]=sr.find_energy(q,r,N,Aw,Aa,B,C)
    Ekin[time]=sr.find_kinetic(p,m)
  print('Time spend ' + str(dotime.process_time() - start)+' s')
  # Plot final microstate
  sr.plot_state(q,r)
  # Plot State variables
  sr.plot_variables(t,Epot,Ekin,P,R,r**2,N)
  sr.velocity_histogram(ps)

# Settings
# Number of atoms, size of box, temperature, length of trajectory, fixed_T
dynamics2D(900,30,300,1000,1,0)
#dynamics2D(900,30,310,1000,1,0)
#dynamics2D(900,30*np.sqrt(1.0/0.1),300,10000,1,0)
#dynamics2D(900,30*np.sqrt(1.0/0.4),300,10000,1,0)
#dynamics2D(900,30*np.sqrt(1.0/0.7),300,1000,1,0)
#dynamics2D(900,30*np.sqrt(1.0/10),300,10000,1,0)
#dynamics2D(900,30*np.sqrt(1.0/3),300,10000,1,0)
#dynamics2D(900,30*np.sqrt(1.0/1.5),300,10000,1,0)
