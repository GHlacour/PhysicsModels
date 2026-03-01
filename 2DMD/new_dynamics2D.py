import numpy as np
import matplotlib.pyplot as plt
import new_subroutines as sr
import time as dotime
from numba import jit
from numba import njit

# Molecular Dynamics Program
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
  R=8.3144598e-3 # kJ/mol/K Gas constant
  dt=0.005 # Timestep in ps
  # Aa=180000 # kJ/mol Particle-Particle repulsion Aa*exp(-dq/B)
  # http://sites.science.oregonstate.edu/~hetheriw/astro/rt/info/water/water_models.html
  #  epsilon=0.650 # kJ/mol
  #sigma=0.355 # nm
  # gamma=12.75
  #Aa=epsilon/(1-6/gamma)*6/gamma
  #print(Aa)
  #B=sigma/gamma
  #print(B)
  #B=0.025 # nm 
  #B=0.03 # nm
  Aa=180000 # kJ/mol Particle-Particle repulsion Aa*exp(-dq/B)
  B=0.025 # nm

  N=int(np.floor(np.sqrt(N))**2) # Make number of particles a square number
  q=sr.roster(N,r) 
  Epot=sr.find_energy(q,r,N,Aa,B)
  # Initialize momenta
  p=np.random.normal(loc=0,scale=1,size=(N,2))
  # Scale momenta to get right kinetci energy
  p=sr.scale_velocity(N,p,m,Temp)
  # Initialize arrays
  Ekin=np.zeros(times)
  Epot=np.zeros(times)
  T=np.zeros(times)
  t=np.arange(0,times)*dt
  P=np.zeros(times)
  Epot[0]=sr.find_energy(q,r,N,Aa,B)
  Ekin[0]=sr.find_kinetic(p,m)
  ps=p
  # Initial step
  (F,Pt)=sr.find_forces(q,r,N,Aa,B)
  P[0]=Ekin[0]/(r**2) +Pt/2/N
  q0=q
  q=sr.update_positions(q,p,dt,m,F,r)
  # Loop over all times
  for time in range(1,times):
    if (fast>0 and time%fast==0):
      print('Timesteps passed ' + str(time) + '.')
    #  Update Positions
    # q=sr.update_positions(q,p,dt,m,F,r)
    [q,q0]=sr.update_positions_verlet(q,q0,dt,m,F,r)
    # Find new forces
    Fold=F
    [F,Pt]=sr.find_forces(q,r,N,Aa,B)
    # Update momentum
    p=sr.update_momentum(p,Fold,F,dt)
    # Store the momentum every 10 steps
    if (time % 10==0):
      ps=np.append(ps,p,axis=0)
    if (fixed_T==1):
       p=sr.scale_velocity(N,p,m,Temp)
    # Plot microstate
    if (fast>0 and time%fast==0):
      sr.plot_state(q,r)
    # Add energies in array
    Epot[time]=sr.find_energy(q,r,N,Aa,B)
    Ekin[time]=sr.find_kinetic(p,m)
    # Add pressure in array
    P[time]=Ekin[time]/(r**2) +Pt/2/N
  print('Time spend ' + str(dotime.process_time() - start)+' s')
  # Plot final microstate
  sr.plot_state(q,r)
  # Plot State variables
  sr.plot_variables(t,Epot,Ekin,P,R,r**2,N)
  sr.velocity_histogram(ps)

# Settings
# Number of atoms, size of box (nm), temperature (K), length of trajectory (frames), fixed T (0/1),
# steps between snapshots to visualize (0=do not visualize)
#dynamics2D(1,30,300,20000,0,100) # Single particle demonstration
#dynamics2D(81,30,300,20000,0,100) # Small box demonstration
#dynamics2D(81,30,300,50000,1,100) # We fix temperature
#dynamics2D(81,8,300,10000,1,100) # Larger density by smaller box
dynamics2D(900,30,300,10000,1,0) # Larger density by more particles

#dynamics2D(900,30,310,10000,1,100) # Larger density different temperature
#dynamics2D(900,30*np.sqrt(1.0/0.1),300,10000,1,0) # Runs with different density
#dynamics2D(900,30*np.sqrt(1.0/0.4),300,10000,1,0)
#dynamics2D(900,30*np.sqrt(1.0/0.7),300,10000,1,0)
#dynamics2D(900,30*np.sqrt(1.0/10),300,10000,1,0)
#dynamics2D(900,30*np.sqrt(1.0/3),300,10000,1,0)
#dynamics2D(900,30*np.sqrt(1.0/1.5),300,10000,1,0)
