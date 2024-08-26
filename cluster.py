### Kaixin Wang ###
import os, sys
os.environ['OVITO_GUI_MODE'] = '1'
import ovito
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

data_file = sys.argv[1]
data_name = os.path.basename(data_file)
data_name_without_extension = os.path.splitext(data_name)[0]
# Load the trajectory file using Ovito
pipeline = ovito.io.import_file(data_file)
clusterlist=[]

timestep = float(input("Timestep(fs) (Default: 5):") or 5)
cutoff = int(input("The cutoff of the size of cluster(Angstrom) (Default: 6):") or 6)
p_type = int(input("The type number of atom (Default: 3):") or 3) #particle type
start_t = float(input("Start time(ns) (Default: 0):") or 0)
end_t = float(input("End time(ns) (Default: 25):") or 25)

pipeline.modifiers.append(ovito.modifiers.ExpressionSelectionModifier(expression = r"ParticleType==%i"%p_type))
pipeline.modifiers.append(ovito.modifiers.ClusterAnalysisModifier(
    cutoff=cutoff, 
    sort_by_size=True, 
    compute_com=True,
    only_selected=True,
    cluster_coloring=True,
    unwrap_particles=True))
for i in tqdm(range(pipeline.source.num_frames),desc = 'calculate'):
    #export_file(pipeline, 'clusters.txt', 'txt/table', key='clusters')
    #export_file(pipeline, "output.xyz", "xyz", columns =["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z","Cluster"])
    data = pipeline.compute(i)
    #print(data.tables)
    cluster_table = data.tables['clusters']
    #print(cluster_table)
    try: 
        clusterlist.append(cluster_table['Cluster Size'][0])
    except IndexError:
        clusterlist.append(0)

t, cluster = zip(*enumerate(clusterlist))
t = np.array(t)*20000*timestep/1000000
cluster = np.array(cluster)
fit_size = cluster[(t>=start_t)&(t<=end_t)]
ave_size = np.average(fit_size)
err_size = np.std(fit_size,ddof=1)/((len(fit_size))**0.5)
std_size = np.std(fit_size,ddof=1)
#print(fit_size)
#print(err_size*((len(fit_size))**0.5))
f = open("%s_ClusterSize.dat"%data_name_without_extension,"w")
for i in range(len(t)):
    f.write('%.3f %d\n'%(t[i],cluster[i]))
f.close()

fig, ax =plt.subplots()
ax.plot(t,cluster,"-")
ax.set_xlabel(r"Time(ns)")
ax.set_ylabel(r"Nanobubble Size")
ax.set_title(r"Average NB Size: %.3f$\pm$%.3f; Time range: %d-%dns"%(ave_size,std_size,start_t,end_t))
print("std" + str(std_size))
plt.savefig("%s_ClusterSize.png"%data_name_without_extension,dpi=300,bbox_inches='tight')
plt.show()

