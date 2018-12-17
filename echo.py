from vpython import *
import numpy as np

class grid:
  def __init__(self,N,x,y,z):
    self.spins =[]
    self.total=[]
    for i in range(N):
      spin=arrow()
      spin.axis=vector(x[i],y[i],z[i])
      spin.pos=vector(0,0,0)
      spin.color=color.yellow
      self.spins.append(spin)
    
N=10
dt=0.01
x=np.ones((N))
y=np.zeros((N))
z=np.zeros((N))
myspins=grid(N,x,y,z)
v=np.random.lognormal(size=N)
sphere(pos=vector(0,0,0),radius=1,color=color.blue,opacity=0.3)
t=0
for it in range(1000):
    rate(60)
    x=np.cos(t*v)
    y=np.sin(t*v)
    for i in range(N):
      myspins.spins[i].axis=vector(x[i],y[i],z[i])
    t=t+dt
    
# Reverse Direction
v0=t*v
v2=np.random.lognormal(size=N)
T2=0.2
Tc=1
v=-(v*exp(-T2/Tc)+v2*(1-exp(-T2/Tc)))
t=0
for it in range(1000):
    rate(60)
    x=np.cos(v0+t*v)
    y=np.sin(v0+t*v)
    for i in range(N):
      myspins.spins[i].axis=vector(x[i],y[i],z[i])
    t=t+dt
