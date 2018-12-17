from vpython import *

# Bruce Sherwood

class light:
        
    def __init__(self):
        # Define physical objects
        self.photon = []
        
        # Create physical objects
        make_photon(self, start, end, visible)

    # Make light 
    def make_spring(self, start, end, visible):
        spring = helix()
        spring.pos = start.pos
        spring.axis = end.pos-start.pos
        spring.visible = visible
        spring.thickness = 0.05
        spring.radius = 0.5*atom_radius
        spring.length = spacing
        spring.start = start
        spring.end = end
        spring.color = color.orange
        self.springs.append(spring)

c = crystal(N, atom_radius, spacing, 0.1*spacing*sqrt(k/m))

while True:
    rate(60)
    # Do Physics
    if photon.visible:
            atom.pos = atom.pos + atom.momentum/m*dt
    for spring in c.springs:
        spring.axis = spring.end.pos - spring.start.pos
        L = mag(spring.axis)
        spring.axis = spring.axis.norm()
        spring.pos = spring.start.pos+0.5*atom_radius*spring.axis
        Ls = L-atom_radius
        spring.length = Ls
        Fdt = spring.axis * (k*dt * (1-spacing/L))
        if spring.start.visible:
            spring.start.momentum = spring.start.momentum + Fdt
        if spring.end.visible:
            spring.end.momentum = spring.end.momentum - Fdt


