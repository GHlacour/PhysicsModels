from vpython import *

class photon:
  def __init__(self,x,v,z):
    self.pos=x
    self.axis=v
    self.pol=z

    self.myarrow=arrow(pos=self.pos,axis=self.axis,color=color.yellow)
    self.wave=curve(self.pos+0.125*self.axis)
    for ix in range(21):
      wpos=self.pos
      wpos=wpos+(ix+5)/40*self.axis
      wpos=wpos+self.pol*sin(ix*2*pi/10)*exp(-(ix-10)**2/40)
      self.wave.append(wpos)
    self.wave.color=color.yellow
    self.wave.emissive=True

  def update(self,x):
    self.myarrow.pos=x
    self.wave.origin=x

x=vector(0,0,0)
v=vector(1,0,0)
z=vector(0,1,0)
firstphoton=photon(x,v,z)
sphere(pos=x,radius=0.1)
dt=0.01
while True:
    rate(60)
    x=x+dt*v
    firstphoton.update(x)
    print(x)
