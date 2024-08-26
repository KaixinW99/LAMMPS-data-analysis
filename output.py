import os
import re
import sys
import numpy as np
from numpy.polynomial import polynomial as P
import matplotlib.pyplot as plt

def output(FILENAME):
    global xlo, ylo, zlo, xhi, yhi, zhi
    file=open(FILENAME, "r")
    lines=file.readlines()
    file.close()
    c_react = 0 ### the counter of the reaction 
    xreact,yreact,zreact=[],[],[] ### the position of total reaction happened
    nreact = [] ### the accumulated number of total reaction happened
    datalist = [] ### all the computed data for each step
    for index,line in enumerate(lines):
        if re.search("Created orthogonal box",line):
            scope=[pos.split(" ") for pos in line.strip("Created orthogonal box = ( ) \n").split(") to (")]
            xlo, ylo, zlo = [float(item) for item in scope[0]]
            xhi, yhi, zhi = [float(item) for item in scope[1]]
            print("-"*36)
            print("The box range is")
            print("(xlo, xhi) = (%f, %f)"%(xlo,xhi))
            print("(ylo, yhi) = (%f, %f)"%(ylo,yhi))
            print("(zlo, zhi) = (%f, %f)"%(zlo,zhi))
        elif re.match(r"(\d+) (\w+)",line):
            print("-"*36)
            numoftype=line.strip("\n").split(" ")
            name=numoftype[-1]
            num =int(numoftype[0])
            print("Group: %s, Number: %d"%(name,num))
        elif re.search("Setting up Verlet run",line):
            break
    for Nindex, line in enumerate(lines[(index+5):]):
        if re.search("Step",line):
            print("-"*36)
            c_items = [item for item in line.strip("\n").split(" ") if item!=""]
            print("All computed items: (start from 0 index)")
            print(*c_items,"nreact","xreact","yreact","zreact",sep=",")
        elif re.search("freaction:",line):
            flist=re.findall(r"(\d+)\.?(\d+)?",line)
            xreact.append(float(flist[0][0]+"."+flist[0][1]))
            yreact.append(float(flist[1][0]+"."+flist[1][1]))
            zreact.append(float(flist[2][0]+"."+flist[2][1]))
            c_react+=1
        elif re.search("breaction:",line):
            c_react-=1
        elif re.search("Loop time",line):
            break
        elif re.search("ERROR",line):
            print("ERROR!!!: "+line)
            break
        else:
            all_items = [float(item) for item in line.strip(" \n").split(" ") if item!=""]
            datalist.append(all_items)
            # add some time check point for reaction
            xreact.append(["Time",all_items[0]])
            yreact.append(["Time",all_items[0]])
            zreact.append(["Time",all_items[0]])
            nreact.append(c_react)
    return *list(zip(*datalist)), nreact, xreact, yreact, zreact

def current(t,q,start_t,end_t):
    k,cov = P.polyfit(t[(t>start_t)&(t<end_t)],q[(t>start_t)&(t<end_t)],deg=1,full=True)
    #print(cov)
    #print("Current (nA):",k[1]*(10**9))
    return k[1]*(10**9) # Current: nA

data_file = sys.argv[1]
data_name = os.path.basename(data_file)
data_name_without_extension = os.path.splitext(data_name)[0]

Step,Temp,c_mTemp,TotEng,Press,nreact,x,y,z = output(data_file)

timestep = float(input("Timestep(angstrom) (Default: 5):") or 5)
time = np.array(Step)*timestep/1000000 # unit: ns
q = np.array(nreact)*(1.6e-10) # unit: nC
#xoy = np.column_stack((x,y))

# Output the charge Data
f = open("%s_charge.dat"%data_name_without_extension,"w")
for i in range(len(time)):
    f.write("%.3f %.8E \n"%(time[i],q[i]))
f.close()

interval = 1 # ns
# Plot the current vs time
fig, ax =plt.subplots()
xrange = np.arange(0,time[-1]+1,interval)
cur_list = []
for i in range(len(xrange)-1):
    cur_list.append(current(time,q,xrange[i],xrange[i+1]))
cur_list = np.array(cur_list)
ax.plot(xrange[:-1], cur_list,"+-")

f = open("%s_current_%dns.dat"%(data_name_without_extension,interval),"w")
for i in range(len(cur_list)):
    f.write("%.3f %.8E \n"%(xrange[i],cur_list[i]))
f.close()

f = open("%s_rate_%dns.dat"%(data_name_without_extension,interval),"w")
for i in range(len(cur_list)):
    f.write("%.3f %.3f \n"%(xrange[i],cur_list[i]/10**9/1.6e-10)) # /ns
f.close()

start_ave = float(input("The start time for averaging the current(ns) (Default: 50):") or 50)
end_ave = float(input("The end time for averaging the current(ns) (Default: 100):") or time[-1])
fit_cur = cur_list[(xrange[:-1]>=start_ave)&(xrange[:-1]<=end_ave)]
ave_cur = np.average(fit_cur)
std_cur = np.std(fit_cur,ddof=1)/(20**0.5)
error_cur = np.std(fit_cur,ddof=1)/((len(fit_cur))**0.5)
ax.set_title(r"Current %.3f$\pm$%.3fnA; Time range: %d-%dns"%(ave_cur,std_cur,start_ave,end_ave))
plt.savefig("%s_current_%dns.png"%(data_name_without_extension,interval),dpi=300,bbox_inches='tight')
plt.show()

# Rxn Map
def reaction_pos(pos):
    rxn_pos=[]
    time_chkpoint = time[-1]
    for item in pos[::-1]:
        if isinstance(item,list):
            time_chkpoint = item[1]*timestep/1000000
        else:
            rxn_pos.append([time_chkpoint,item])
    return rxn_pos

bin_size = int(input("Bin size for reaction map(Angstrom) (Default: 1):") or 1) # angstrom
x_edges = np.arange(0,180+1,bin_size)
y_edges = np.arange(0,180+1,bin_size)
start_map = float(input("The start time for reaction map (Default: 50):") or 50)
end_map = float(input("The end time for reaction map (Default: 100):") or time[-1])
x = [pos for time, pos in reaction_pos(x) if start_map<=time<=end_map]
y = [pos for time, pos in reaction_pos(y) if start_map<=time<=end_map]
z = [pos for time, pos in reaction_pos(z) if start_map<=time<=end_map]

f = open("%s_3drxnMap.dat"%data_name_without_extension,"w")
f.write("Reaction Map From %s and %s\n"%(start_map,end_map))
for i in range(len(x)):
    f.write("%.3f %.3f %.3f\n"%(x[i],y[i],z[i]))
f.close()

hist2d, x_edges, y_edges = np.histogram2d(x, y, bins=(x_edges,y_edges))
x_edges, y_edges = (x_edges[1:]+x_edges[:-1])/2,(y_edges[1:]+y_edges[:-1])/2

fig, ax =plt.subplots()
# Plot the heatmap by imshow()
heatmap  = ax.imshow(hist2d.T, cmap='plasma', origin='lower',extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]])
cbar = fig.colorbar(heatmap)
plt.savefig("%s_rxnMap.png"%data_name_without_extension,dpi=300,bbox_inches='tight')
plt.show()