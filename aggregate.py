from vpython import *

class photon:
  def __init__(self,x,v,z):
    self.pos=x
    self.axis=v
    self.pol=z

    self.myarrow=arrow(pos=self.pos,axis=self.axis,color=color.red)
    self.wave=curve(self.pos+0.125*self.axis)
    for ix in range(21):
      wpos=self.pos
      wpos=wpos+(ix+5)/40*self.axis
      wpos=wpos+self.pol*sin(ix*2*pi/10)*exp(-(ix-10)**2/40)
      self.wave.append(wpos)
    self.wave.color=color.red
    self.wave.emissive=True

  def update(self,x):
    self.myarrow.pos=x
    self.wave.origin=x-self.pos

x=vector(0,3,0)
v=vector(0,-1,0)
z=vector(1,0,0)
scene=canvas(background=color.white)
scene.camera.pos=vector(4,3,1)
#scene.center=vector(2,-2,0)
#sphere(pos=x,radius=0.1)
box(pos=vector(0,-2,0),height=0.1,length=5,width=5,texture=textures.wood)
cylinder(pos=vector(-2,-1,0),axis=vector(4,0,0),color=color.green,texture=textures.rough)
for j in range(9):
  x=vector(0,3,0)
  v=vector(0,-1,0)
  z=vector(cos(pi*j/8),0,sin(pi*j/8))
  firstphoton=photon(x,v,z)
  dt=0.01
  for i in range(400):
      rate(60)
      x=x+dt*v
      firstphoton.update(x)

