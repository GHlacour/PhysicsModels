from vpython import *

class water:
  def __init__(self,x,y,z):
    self.pos=x

    v1=cos(109.5)*z+sin(109.5)*y+x
    v2=-cos(109.5)*z+sin(109.5)*y+x
    self.oxygen=sphere(pos=self.pos,radius=0.5,color=color.red)
    self.hydrogen1=sphere(pos=v1,radius=0.2,color=color.white)
    self.hydrogen2=sphere(pos=v2,radius=0.2,color=color.white)
    self.bond1=helix(pos=self.pos,radius=0.1,axis=v1-self.pos,thickness=0.05)
    self.bond2=helix(pos=self.pos,radius=0.1,axis=v2-self.pos,thickness=0.05)

  def update(self,x,y,z,l1,l2):
    v1=l1*(cos(109.5)*z+sin(109.5)*y)+x
    v2=l2*(-cos(109.5)*z+sin(109.5)*y)+x
    self.hydrogen1.pos=v1
    self.hydrogen2.pos=v2
    self.bond1.axis=v1-self.pos
    self.bond2.axis=v2-self.pos

x=vector(0,0,0)
y=vector(1,0,0)
z=vector(0,1,0)
firstwater=water(x,y,z)
t=0
dt=0.01
while True:
    rate(60)
    t=t+dt
    l1=sin(10*t)*0.1+1
    l2=sin(10*t)*0.1+1
    firstwater.update(x,y,z,l1,l2)
    
