import numpy as np
import matplotlib.pyplot as plt
import ovito
import sys
import os

data_file = sys.argv[1]
data_name = os.path.basename(data_file)
data_name_without_extension = os.path.splitext(data_name)[0]
# Load the trajectory file using Ovito
pipeline = ovito.io.import_file(data_file)

# Get the basic setting
data = pipeline.compute()
box_x, box_y, box_z = [r[n] for n,r in enumerate(data.cell)]
cutoff = int(input("The cutoff of the size of cluster(Angstrom) (Default: 10):") or 10)
p_type = int(input("The type number of atom (Default: 3):") or 3) #particle type

# Select the Particle Type
pipeline.modifiers.append(ovito.modifiers.SelectTypeModifier(property = "Particle Type", types = {p_type}))

# Cluster Analysis
pipeline.modifiers.append(ovito.modifiers.ClusterAnalysisModifier(
    cutoff=cutoff,
    sort_by_size=True,
    compute_com=True,
    only_selected=True,
    cluster_coloring=True,
    unwrap_particles=True,
))

# Move it to center
def move_to_center(frame,data):
    data.cell_.pbc = (True, True, False)
    ctable = data.tables["clusters"]
    cluster_coms = np.array(ctable["Center of Mass"])
    cluster_com = cluster_coms[0]
    data.particles_.positions_ -= cluster_com - np.array([box_x/2,box_y/2,cluster_com[2]])


pipeline.modifiers.append(move_to_center)
#ovito.io.export_file(pipeline, "output_trial2.lammpstrj", "lammps/dump", multiple_frames=True, columns=["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"])

yz_data, num_bins_y, num_bins_z = [], int(box_y), int(box_z)
num_frames = pipeline.source.num_frames

start = int(input("Start Frame (Default: total/2):") or num_frames//2)
end = int(input("End Frame (Default: total):") or num_frames)


# Map the data on the yoz plane
for i in range(start,end):
    data = pipeline.compute(i)
    Info_particles = data.particles
    particle_mask = (Info_particles.particle_type == p_type)
    yz_data += list(Info_particles.positions[particle_mask][:,1:])

bin_size = int(input("Bin size(Angstrom) (Default: 2):") or 2) # angstrom
y_edges = np.arange(0,num_bins_y+1,bin_size)
z_edges = np.arange(0,num_bins_z+1,bin_size)
hist2d, y_edges, z_edges = np.histogram2d(np.array(yz_data)[:,0], np.array(yz_data)[:,1], bins=(y_edges,z_edges))
mid_y_index = (num_bins_y-1)//(2*bin_size)
# Histogram does not follow Cartesian convention (see Notes),
# therefore transpose H for visualization purposes.
hist2d = hist2d.T/(bin_size)**2/(end-start)/box_x
y_edges, z_edges = (y_edges[1:]+y_edges[:-1])/2,(z_edges[1:]+z_edges[:-1])/2
print("!!!Please check the density level (half of the surface excess)!!!")

# Plot the histogram as a 2D map
fig, ax =plt.subplots()
ax.plot(z_edges, hist2d[:,mid_y_index])
# Set labels for the x and y axes with LaTeX formatting
ax.set_xlabel(r'z($\mathring{A}$)')
ax.set_ylabel(r'Density')
pic_name = input("The name of density distribution on z direction (E.x. %s_zden.png):"%(data_name_without_extension)) or "%s_zden.png"%(data_name_without_extension)
plt.savefig(pic_name,dpi=300,bbox_inches='tight')
plt.show()

fig, ax =plt.subplots()
# Plot the heatmap by imshow()
Y, Z = np.meshgrid(y_edges,z_edges)
heatmap  = ax.imshow(hist2d, cmap='plasma', origin='lower',extent=[y_edges[0], y_edges[-1], z_edges[0], z_edges[-1]])
cbar = fig.colorbar(heatmap)

# Plot the contour based on the density level
level = float(input("The density level (half of the surface excess):") or 0.0038)
contours = ax.contour(Y,Z,hist2d, levels=[level], colors="white")
#ax.set(xlim=[box_y/2,box_y/4*3],ylim=[0,box_y/4])
ax.set(xlim=[box_y/2,box_y],ylim=[0,box_y/2])
ax.set_xlabel(r'r($\mathring{A}$)')
ax.set_ylabel(r'z($\mathring{A}$)')
pic_name = input("The name of heatmap without fitting curve (E.x. %s_ConAng.png):"%(data_name_without_extension)) or "%s_ConAng.png"%(data_name_without_extension)
plt.savefig(pic_name,dpi=300,bbox_inches='tight')
print("!!!Please check the baseline of the bubble!!!")
plt.show()

# get the vertices of contour
paths = contours.collections[0].get_paths()[0]
vertices = paths.vertices
y, z = vertices[:,0], vertices[:,1]

# Set the baseline of nanobubble
baseline = float(input("Baseline of the bubble(Angstrom):") or 10.71)
y, z = y[z>baseline], z[z>baseline]

# Start to fit the curve
from scipy.optimize import least_squares
# Define the circle equation
def circle_equation(params, x, y):
    xc, yc, r = params
    return (x - xc)**2 + (y - yc)**2 - r**2

# Fit the circle equation to the data
params_guess = [np.mean(y), np.mean(z), 1.0]
res = least_squares(circle_equation, params_guess, args=(y,z))

# Get the fitted circle parameters
yc_fit, zc_fit, r_fit = res.x
print("Center: ({:.3f}, {:.3f}), Radius: {:.3f}".format(yc_fit, zc_fit, r_fit))

# generate data points on the circle
theta = np.linspace(0, 2*np.pi, 100)
y_c = yc_fit + r_fit * np.cos(theta)
z_c = zc_fit + r_fit * np.sin(theta)

# Calculate the contact angle
contactangle = 90 - np.degrees(np.arcsin((baseline-zc_fit)/r_fit))
print("Contact angle: {:.1f}".format(contactangle))

# plot the headmap and contour again
fig, ax =plt.subplots()
heatmap  = ax.imshow(hist2d, cmap='plasma', origin='lower',extent=[y_edges[0], y_edges[-1], z_edges[0], z_edges[-1]])
contours = ax.contour(Y,Z,hist2d, levels=[level], colors="white")
cbar = fig.colorbar(heatmap)

# plot the baseline
base = np.linspace(0,y_edges,100)
ax.plot(base,np.full_like(base,baseline))

# plot the fitting curve
ax.plot(y_c[z_c>0],z_c[z_c>0],"pink",label="Circle")
#ax.set(xlim=[box_y/2,box_y/4*3],ylim=[0,box_y/4])
ax.set(xlim=[box_y/2,box_y],ylim=[0,box_y/2])
ax.set_xlabel(r'r($\mathring{A}$)')
ax.set_ylabel(r'z($\mathring{A}$)')
ax.set_title(r'Radius: %.3f$\mathring{A}$, Contact Angle: %.1f$^{\circ}$'%(r_fit,contactangle))
pic_name = input("The name of heatmap with fitting curve (E.x. %s):"%(pic_name)) or pic_name
plt.savefig(pic_name,dpi=300,bbox_inches='tight')
plt.show()

#ovito.io.export_file(pipeline, "output_trial2.lammpstrj", "lammps/dump", multiple_frames=True, columns=["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"])
