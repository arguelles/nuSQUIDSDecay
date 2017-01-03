#! /usr/bin/python
import matplotlib as mpl

import numpy as np
import matplotlib.pyplot as plt
import sys


from pylab import rcParams
rcParams['figure.figsize'] = 10, 8
plt.rc('font',family='Times New Roman')
plt.rc('font', size=16)

colors=['#FF0000','#FF7F00','#FFFF00','#00FF00','#0000FF','#4B0082','#9400D3']
sheepcolors=['#000000']

nusquids_data=[]
nusheep_data=[]
tolerances=[]
sheeptol=[]

tolvals=[6,8,10,12,14,17,20]
for x in tolvals:
	fname="output/earth_dcy_e"+str(x)+"_10e2lt_regen.npz"
	nusquids_data.append(np.load(fname))
	tolerances.append("Tol: "+str(x))

sheep_tolvals=[6]


nusheep_data=[]
nusheep_data.append(np.load("energy_regen_numode_matter.npz"))
sheeptol.append("Tol: 6")


#for x in sheep_tolvals:
#	fname="energy_regen_e"+str(x)+"_matter.npz"
#	nusheep_data.append(np.load(fname))
#	sheeptol.append("Tol: "+str(x))

#-----------------------------------------------------------------#
#Plotting data

"""
# Plotting true sterile-muon oscillation probability.
trueE = np.logspace(2,5,10000)

earth_diameter = 64573027604478.22 #Earth diameter in inverse eV

def p_mumu(l,E,phase):
	dm2_42 = 0.999999994375 
	p = np.cos(phase+dm2_42*l/(4.0*E))**2
	return p

atm=-0.05825
pvec = np.vectorize(p_mumu)
trueP = pvec(earth_diameter,np.multiply(trueE,1e9),0)
shiftP = pvec(earth_diameter*(1.0+atm),np.multiply(trueE,1e9),0)

print "---------------------"
print "MAXIMUM ABS DEVIATION"
print np.max(np.abs(np.subtract(regen_flux[1,1,:],nusheep_data['_num_vec'])))
print "---------------------"
"""
fig,ax = plt.subplots()

handles_1 = []
n=0
hand, = ax.plot([],[],color="black",marker="",lw=2,ls='-',label="nuSQUIDS")
handles_1.append(hand)
hand, = ax.plot([],[],color="black",marker="",lw=2,ls='--',label="nuSHEEP: 1e-6")
handles_1.append(hand)

handles_2 = []

for n,tol in enumerate(tolvals):
	hand, = ax.plot((nusquids_data[n])['e'],((nusquids_data[n])['flux'])[0,1,:],color=colors[n],marker="",lw=2,ls='-',label=tolerances[n])
	handles_2.append(hand)
	
for n,tol in enumerate(sheep_tolvals):
	hand, = ax.plot(np.multiply(1e-9,(nusheep_data[n])['_E']),(nusheep_data[n])['_num_vec'],color=sheepcolors[n],marker="",lw=2,ls='--')
	handles_2.append(hand)
plt.xscale('log')
#plt.yscale('log')
#ax.set_ylim([0,1e-10])

ax.set_xlim([1e2,1e4])
ax.set_ylim([-0.1,0.7])
ax.set_xlabel("Neutrino Energy [GeV]")
ax.set_ylabel(r'$P_{\mu\rightarrow\mu}$')
plt.axhline(color='black')

# Create a legend for the first line.
leg1 = plt.legend(handles=handles_1, loc="upper right")

# Add the legend manually to the current Axes.
ax = plt.gca().add_artist(leg1)

# Create another legend for the second line.
plt.legend( handles=handles_2, loc="lower left")



#ax.legend(loc='lower left', shadow=False,fancybox=True)

plt.show()
