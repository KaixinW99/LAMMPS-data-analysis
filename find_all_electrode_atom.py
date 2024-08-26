import os, sys
os.environ['OVITO_GUI_MODE'] = '1'
import ovito
import numpy as np
import matplotlib.pyplot as plt

data_file = sys.argv[1]
data_name = os.path.basename(data_file)
data_name_without_extension = os.path.splitext(data_name)[0]

p_type = int(input("The type number of electrode (Default: 4):") or 4) #particle type

# Load the trajectory file using Ovito
pipeline = ovito.io.import_file(data_file)
data = pipeline.compute()
Info_particles = data.particles
particle_mask = (Info_particles.particle_type == p_type)

f = open("Pt_pos.dat","w")
f.write("The electrode atom position (x y z)\n")
for pos in Info_particles.positions[particle_mask]:
    f.write("%.3f %.3f %.3f\n"%(pos[0],pos[1],pos[2]))
f.close()

#plt.scatter(Info_particles.positions[particle_mask][:,0],Info_particles.positions[particle_mask][:,1])
#plt.show()
