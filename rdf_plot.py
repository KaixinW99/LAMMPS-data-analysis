import os, sys
os.environ['OVITO_GUI_MODE'] = '1'
import ovito
import numpy as np
import matplotlib.pyplot as plt

data_file = sys.argv[1]
data_name = os.path.basename(data_file)
data_name_without_extension = os.path.splitext(data_name)[0]

r, rdf = np.genfromtxt(data_file,skip_header=2,usecols=[0,7],unpack=True)
fig,ax=plt.subplots(figsize=[4,3.5])
#ax.plot(td100,ed100,label=r"small hemispheric NP",ls="-",lw=1)
ax.plot(r,rdf,label=r"Pt-mW",ls="-")
ax.set(xlim=[0,10])#title="Radial distribution function")
# Set the legend font size
ax.legend(loc=0,fontsize=12)
# Set the tick label size
ax.tick_params(axis='both', which='both',labeltop=False,labelright=False,right=True,top=True,labelsize=10, direction="in")
# Set the axis label size
ax.set_xlabel(r"r($\AA$)", fontsize=12)
ax.set_ylabel('RDF', fontsize=12)
fig.tight_layout()
plt.savefig('PtmW_rdf_wo_bubble.png',dpi=300,bbox_inches='tight')
plt.show()