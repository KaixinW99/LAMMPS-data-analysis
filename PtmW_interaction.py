import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as SciSig
import os, sys

data_file = sys.argv[1]
data_name = os.path.basename(data_file)
data_name_without_extension = os.path.splitext(data_name)[0]

E = float(input("Enter the electric field strength in mV (default: 500): ") or 500)

F, R, T= 96485.3321, 8.3145, 298
r, rdf = np.genfromtxt(data_file,skip_header=2,usecols=[0,7],unpack=True) # from 0.05 to 9.95, step 0.1
proximity = np.exp(-np.arange(0.05,10,0.1))
#* potential term does not depend on change of distance
potential = np.exp(0.5*F/R/T*E/1000)

fig,ax=plt.subplots(figsize=[4,3.5])
interaction=rdf*proximity*potential
ax.plot(r,interaction,label=r"Pt-mW",ls="-",color="#0000a7")
max_index = np.argmax(interaction)
#ax.scatter(r[max_index],interaction[max_index],color="#c1272d")
plt.text(r[max_index], interaction[max_index]+0.0003, f'r={r[max_index]} at Max Point', color='#c1272d')
#ax.set(xlim=[0,10], ylim=[0,0.0040])
# Set the legend font size
ax.legend(loc=0,fontsize=12)
# Set the tick label size
ax.tick_params(axis='both', which='both',labeltop=False,labelright=False,right=True,top=True,labelsize=10, direction="in")
# Set the axis label size
ax.set_xlabel(r"r($\mathring{A}$)", fontsize=12)
ax.set_ylabel(r"RDF$\times$Reaction", fontsize=12)
ax.set_title("Pt-mW Interaction at %d mV"%E)
fig.tight_layout()
plt.savefig('Pt_mW_interaction_%dmV'%E,dpi=300,bbox_inches='tight')
plt.show()

#print(SciSig.convolve(rdf,proximity, mode="same"))


