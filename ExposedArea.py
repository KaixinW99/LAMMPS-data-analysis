### Kaixin Wang ###
import numpy as np
import os
os.environ['OVITO_GUI_MODE'] = '1'
import matplotlib.pyplot as plt
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from tqdm import tqdm
import argparse
particletype = 4

parser = argparse.ArgumentParser(description='Process lammpstrj to get the Exposed Area.')
parser.add_argument('lammpstrj', metavar='path', type=str, nargs=1,
                    help='path to the lammpstrj file')
args = parser.parse_args()
lammpstrj_path = args.lammpstrj[0]
lammpstrj_name = os.path.basename(lammpstrj_path)
lammpstrj_name_without_extension = os.path.splitext(lammpstrj_name)[0]

reactive=[]
pipeline = import_file(lammpstrj_path)
pipeline.modifiers.append(ExpressionSelectionModifier(expression = "ParticleType==%i"%particletype))
pipeline.modifiers.append(ComputePropertyModifier(
    output_property = 'reactive_region',
    neighbor_expressions = ["ParticleType==2"],
    only_selected = True,
    cutoff_radius = 5
))
for i in tqdm(range(pipeline.source.num_frames),desc = 'calculate'):
    data = pipeline.compute(i)
    sel=data.particles.selection[...]
    num=np.count_nonzero(data.particles["reactive_region"][sel==1])
    reactive.append(num)

f = open("%s_ExposedArea_ptype%i.dat"%(lammpstrj_name_without_extension,particletype),"w")
for d,i in tqdm(enumerate(reactive),desc = 'write'):
    f.write('%.3f %.10f\n'%(d/10,i)) # 112 pt atom for square; 116 pt atom for disk; 216 pt atom for hemibot; 109 pt atom for hemisuf
f.close()

timestep = float(input("Timestep(ns) (Default: 5):") or 5)
t, rxn = zip(*enumerate(reactive))
t = np.array(t)*20000*timestep/1000000
rxn = np.array(rxn)
start_t = float(input("Start time(ns) (Default: 50):") or 50)
end_t = float(input("End time(ns) (Default: 100):") or 100)

fit_rxn = rxn[(t>=start_t)&(t<=end_t)]
ave_rxn = np.average(fit_rxn)
err_rxn = np.std(fit_rxn,ddof=1) #/((len(fit_rxn))**0.5)
fig, ax =plt.subplots()
ax.plot(t,rxn,"-")
ax.set_xlabel("Time(ns)")
ax.set_ylabel("Number of Exposed Atom")
ax.set_title(r"Average Exposed Atoms: %.3f$\pm$%.3f; Time range: %d-%dns"%(ave_rxn,err_rxn,start_t,end_t))
print("Std:" + str(err_rxn))
plt.savefig("%s_ExposedArea_ptype%i.png"%(lammpstrj_name_without_extension,particletype),dpi=300,bbox_inches='tight')
plt.show()