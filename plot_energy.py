#! /usr/bin/python
import matplotlib as mpl

import numpy as np
import matplotlib.pyplot as plt
import sys


from pylab import rcParams
rcParams['figure.figsize'] = 10, 8
plt.rc('font',family='Times New Roman')
plt.rc('font', size=16)


#regen_data=np.load("output/regen_orderflip.npz")
#regen_data=np.load("output/dcy_noregen_vaccuum_2flv.npz")
regen_regdata=np.load("output/earth_dcy_e20_10e2lt_regen.npz")
regen_regdatae6=np.load("output/earth_dcy_e6_10e2lt_regen.npz")
regen_regdatae7=np.load("output/earth_dcy_e7_10e2lt_regen.npz")
regen_data=np.load("output/earth_dcy_e20_10e2lt.npz")
#regen_datafast=np.load("output/dcy_noregen_vacuum_fastdcy_massfix.npz")
#noregen_data=np.load("output/dcy14_noregen_softail.npz")

#nusheep_data=np.load("energy_noregen_freshtest_nusheep_nomatter.npz")
nusheep_regdata=np.load("energy_regen_earthmatter_matter.npz")
nusheep_data=np.load("energy_noregen_earthmatter_matter.npz")
#nusheep_datafast=np.load("energy_noregen_round2_fastdcy_nomatter.npz")
#nusheep_data_noacmix=np.load("energy_noregen_freshtest_nusheep_noactivemix_matter.npz")

regen_energy=regen_regdata['e']
noregen_energy=regen_data['e']
#regen_energyfast=regen_datafast['e']
regen_flux=regen_regdata['flux']
regene6_flux=regen_regdatae6['flux']
regene7_flux=regen_regdatae7['flux']
noregen_flux=regen_data['flux']
#regen_fluxfast=regen_datafast['flux']
#noregen_energy=noregen_data['e']
#noregen_flux=noregen_data['flux']

pi=3.141592
#-----------------------------------------------------------------#
#Plotting data

print regen_energy
print nusheep_data['_E']

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
#ax.plot(regen_energy,regen_flux[1,1,:],color="red",marker="o",lw=2,ls='-',label="Nusquids Muon")
#ax.plot(regen_energy,regen_flux[1,1,:],color="red",marker="o",lw=2,ls='-',label="NuSQUIDS")
#ax.plot(np.multiply(1e-9,nusheep_data['_E']),nusheep_data['_num_vec'],color="blue",marker="o",lw=2,ls='-',label="Nusheep Muon")
ax.plot(np.multiply(1e-9,nusheep_regdata['_E']),nusheep_regdata['_num_vec'],color="blue",marker="",lw=2,ls='-',label="NuSHEEP Regen")
ax.plot(np.multiply(1e-9,nusheep_data['_E']),nusheep_data['_num_vec'],color="blue",marker="",lw=2,ls='--',label="NuSHEEP")
#ax.plot(trueE,trueP,color="black",lw=2,ls='-',label="True Muon")
#ax.plot(trueE,shiftP,color="black",lw=2,ls='--',label="True Muon L-Shifted")
#ax.plot(np.multiply(1e-9,nusheep_data_noacmix['_E']),nusheep_data_noacmix['_num_vec'],color="purple",marker="o",lw=2,ls='-',label="Nusheep Muon Noactive")
ax.plot(regen_energy,regen_flux[1,1,:],color="red",marker="",lw=2,ls='-',label="NuSQUIDS Regen")
ax.plot(regen_energy,regene6_flux[1,1,:],color="green",marker="",lw=2,ls='-',label="NuSQUIDS Regen")
ax.plot(regen_energy,regene7_flux[1,1,:],color="purple",marker="",lw=2,ls='-',label="NuSQUIDS Regen")
ax.plot(noregen_energy,noregen_flux[1,1,:],color="red",marker="",lw=2,ls='--',label="NuSQUIDS")
#ax.plot(noregen_energy,noregen_flux[1,2,:],color="blue",lw=2,marker="o",ls='-',label="Without Regeneration")
plt.xscale('log')
#plt.yscale('log')
#ax.set_ylim([0,1e-10])

ax.set_xlim([1e2,1e4])
ax.set_ylim([0.0,0.6])
ax.set_xlabel("Neutrino Energy [GeV]")
ax.set_ylabel(r'$P_{\mu\rightarrow\mu}$')
ax.legend(loc='lower left', shadow=False,fancybox=True)

"""
axes=[]

xax_vector=valvec

axes.append(plt.subplot2grid((4,4),(1,0),rowspan=3,colspan=4))
if (numcomp==False):
	axes[0].plot(xax_vector,pert_vec0,color="red",label="Standard Vacuum Solution",lw=2,ls='-')
#	axes[0].plot(xax_vector,pert_vec1,color="green",label="",lw=2,ls='--')
#	axes[0].plot(xax_vector,pert_vec2,color="green",label="",lw=2,ls='-.')
	axes[0].plot(xax_vector,O1_vec,color="green",label="Perturbation O1",lw=2,ls='-')
	axes[0].plot(xax_vector,O2_vec,color="blue",label="Perturbation O2",lw=2,ls='-')
#	axes[0].plot(xax_vector,diag_vec,color="black",label="",lw=0,ls='--')
#	axes[0].plot(xax_vector,diag_vec,color="black",label="",lw=0,ls='-.')
axes[0].plot(xax_vector,num_vec,color="black",label="Numerical Solution",lw=2,ls='-')

#if (numcomp==True):
#axes[0].plot(xax_vector,diag_vec,color="orange",label="Diagonal Solution",lw=2,ls='--')

axes[0].set_xlim([xax_vector[0],xax_vector[-1]])
axes[0].set_ylim([0.45,1.1])
axes[0].set_xlabel("Neutrino Energy [TeV]",labelpad=0)
axes[0].set_ylabel("Muon Neutrino Survival Fraction")
axes[0].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))

axes.append(plt.subplot2grid((4,4),(0,0),colspan=4,sharex=axes[0]))

if (numcomp==False):
	axes[1].plot(xax_vector,np.divide(np.subtract(pert_vec0,num_vec),num_vec),color="red",lw=2,ls='-')
#	axes[1].plot(xax_vector,np.divide(np.subtract(pert_vec1,num_vec),num_vec),color="green",lw=2,ls='--')
#	axes[1].plot(xax_vector,np.divide(np.subtract(pert_vec2,num_vec),num_vec),color="green",lw=2,ls='-.')
	axes[1].plot(xax_vector,np.divide(np.subtract(O1_vec,num_vec),num_vec),color="green",lw=2,ls='-')
#axes[1].plot(xax_vector,np.divide(np.subtract(diag_vec,num_vec),num_vec),color="pink",lw=2,ls='-')
	axes[1].plot(xax_vector,np.divide(np.subtract(O2_vec,num_vec),num_vec),color="blue",lw=2,ls='-')
axes[1].plot(xax_vector,np.divide(np.subtract(num_vec,num_vec),num_vec),color="black",lw=2,ls='-')

#if (numcomp==True):
#axes[1].plot(xax_vector,np.divide(np.subtract(num_vec,num_vec),num_vec),color="orange",lw=2,ls='--')

#axes[1].set_yscale('log')
axes[1].set_xlim([xax_vector[0],xax_vector[-1]])
axes[1].set_ylabel("Relative Residual")
axes[1].set_ylim([-0.2,0.2])
axes[1].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))
plt.setp(axes[1].get_xticklabels(), visible=False)

for ax in axes:
	ax.set_xscale('log')



legends=[]
frames=[]

single_axis=[]
single_axis.append(axes[0])

for ax in enumerate(single_axis):

	legends.append(ax[1].legend(loc='lower right', shadow=False,fancybox=True))

	# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
	frames.append(legends[ax[0]].get_frame())
	frames[ax[0]].set_facecolor('1.0')

	# Set the fontsize
	for label in legends[ax[0]].get_texts():
	    label.set_fontsize(16)
	
	for label in legends[ax[0]].get_lines():
	    label.set_linewidth(2.0)  # the legend line width

for ax in axes:
	ax.tick_params(axis='x', which='major', width=1,length=5)
	ax.tick_params(axis='x', which='minor', width=1,length=3)
	ax.tick_params(axis='y', which='major', width=1,length=5)
	ax.tick_params(axis='y', which='minor', width=1,length=3)


# Get the Bbox
plt.draw()
bb = legends[0].legendPatch.get_bbox().inverse_transformed(axes[0].transAxes)

# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='cyan', alpha=0.2)


def extract_man_exp(val):
	#Extract mantissae and exponents for pretty printing!
	mystr="%.2e"%(val)
	mantissa=[]
	exponent=[]
	print mystr
	hit_e=False
	for c in mystr:
		if (c=='e'):
			hit_e=True
		else:
			if (hit_e==False):
				mantissa.append(c)
			else:
				exponent.append(c)
	man_str=''.join(mantissa)
	exp_str=''.join(exponent)
	
	ret=[float(man_str),int(exp_str)]
	return ret
		

def pretty_exp_note(parstring,unitstring,val,f_force=False,last=False):
	if (('e' not in '%g'%(val)) or (f_force==True)):
		ret = parstring+r'$=%.2f$'%(val)+'  ['+unitstring+']'
	else:
		splitval=extract_man_exp(val)
		ret = parstring+r'$=%.2f \times 10^{%d}$'%(splitval[0],splitval[1])+' ['+unitstring+']'
	if (last==False):
		ret+='\n'
	return ret



# place a text box relative to the legend

textlist=[
r'$\phi=\pi/4$  [rad]'+'\n', #Special value for phi
pretty_exp_note(r'$\theta$','rad',theta,f_force=True),
pretty_exp_note(r'$V$','TeV',potential/param.TeV),
pretty_exp_note(r'$\lambda_1$','TeV',l1/param.TeV),
pretty_exp_note(r'$\lambda_2$','TeV',l2/param.TeV,last=True)]
textstr=''.join(textlist)


axes[0].text(0.25, 0.37, textstr, transform=axes[0].transAxes, fontsize=16,verticalalignment='top', bbox=props)
"""

plt.show()


#np.savez("4TeVnodcy_matter"+str(cmdargs[0])+"timeplot"+str(cmdargs[1]),d=distvec,p=pvec)
