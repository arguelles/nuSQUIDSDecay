#! /usr/bin/python

import numpy as np
import sys

numneu=4
fname=sys.argv[1].replace(".dat","")
print fname
txtdata = np.loadtxt(sys.argv[1])

#Assuming three flavors
energy = txtdata[:,0]
nsamples = len(energy)
flux = np.zeros((2,numneu,nsamples))
for ntype in range(0,2):
	for flv in range(0,numneu):
		flux[ntype,flv,:] = txtdata[:,1+numneu*ntype+flv]

np.savez("output/"+fname.replace("mains/","")+".npz", e=energy, flux=flux)
