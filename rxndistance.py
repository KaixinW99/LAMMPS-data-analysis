import numpy as np
import matplotlib.pyplot as plt
import os, sys
import re

#xcenter, ycenter, zcenter = 45, 45, 0
#def distance(x,y,z): 
#    return ((x-xcenter)**2+(y-ycenter)**2+(z-zcenter)**2)**0.5

data_file = sys.argv[2]
data_name = os.path.basename(data_file)
data_name_without_extension = os.path.splitext(data_name)[0]
electrode_atom_data = sys.argv[1]
pot = re.search(r'(\d+)', data_name)[0]

print()
print('%sns'%open(data_file).readline().strip())
print()
rxns = np.genfromtxt(data_file, unpack=False, skip_header=1)
electrode = np.genfromtxt(electrode_atom_data, unpack=False, skip_header=1)

distances = []
for rxn in rxns:
    distance = np.inf
    for atom in electrode:
        distance = min(distance, (sum((rxn-atom)**2))**0.5)
    distances.append(distance)
    
hist, bin_edges = np.histogram(distances, bins=20)
bin_centers = (bin_edges[:-1] + bin_edges[1:])/2

fig, ax = plt.subplots()
ax.plot(bin_centers, hist, "o-")
ax.set_xlabel("Distance (Angstrom)")
ax.set_ylabel("Reaction Frequency")
ax.set_xlim(0,10)
ax.set_title("%smV Reaction Frequency %.1f"%(pot,sum(hist[bin_centers<=5])/sum(hist)*100)+'%'+' below 5A')
plt.savefig("%s_rxnfreq.png"%data_name_without_extension, dpi=300, bbox_inches='tight')
plt.show()