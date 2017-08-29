#!/usr/bin/env python
#Creates PKL with data for loss figures
import numpy as np
import os
import lfmPostproc as lfmpp
import cPickle as pickle
import glob
import argparse
import h5py

if __name__ == "__main__":
	gId = "Step#1000" #Final step

	parser = argparse.ArgumentParser(description="Generates PKLs from injected wedge data.")
	parser.add_argument('pops',nargs='+',metavar='nval',help="List of pops (n=0,3)")

	args = parser.parse_args()

	RootDir = os.path.expanduser('~') + "/Work/StormPSD/Data/"

	Stubs = ["Trap","Inj0","Inj21","Inj3"]

	nPops = args.pops #List of poopulations to do

	NumP = len(nPops)

	print(nPops)
	#Loop over stubs, and each file within each directory
	for nS in nPops:
		n = np.int(nS)
		fDir = RootDir+Stubs[n]
		h5ps = glob.glob(fDir+"/*.h5part")
		fOut = "%s.pkl"%Stubs[n]
		print("Starting collection %s"%Stubs[n])
		#Create initial holders, form list of arrays

		Xs = []
		Ys = []
		Zs = []
		Xeqs = []
		Yeqs = []
		Teqs = []
		Ks = []
		IDs = []

		
		Nread = 1
		Ntot = len(h5ps)
		#Loop through each file, keep track of number
		Np = 0
		for h5p in h5ps:
			#Open file
			with h5py.File(h5p,'r') as hf:
				#Find particles that were ever lost
				grp = hf.get(gId)
				pid = grp.get("id").value
				inP = grp.get("in").value
				x   = grp.get("x").value
				y   = grp.get("y").value
				z   = grp.get("z").value
				xeq = grp.get("xeq").value
				yeq = grp.get("yeq").value
				Teq = grp.get("Teq").value
				K   = grp.get("keveq").value
			#Done with file, now cut out what we want
			I = (inP<0.5) #Particles that were ost
			Xs.append(x[I])
			Ys.append(y[I])
			Zs.append(z[I])
			Xeqs.append(xeq[I])
			Yeqs.append(yeq[I])
			Teqs.append(Teq[I])
			IDs.append(pid[I])
			Ks.append(K[I])

			print("\t(%d/%d): Found %d particles"%(Nread,Ntot,I.sum()))
			Np = Np + I.sum()
			Nread = Nread + 1
			
		#Save data for this collection
		print("Found %d total particles"%(Np))

		#Remap list of small arrays to one big array
		X   = np.concatenate(Xs)
		Y   = np.concatenate(Ys)
		Z   = np.concatenate(Zs)
		Xeq = np.concatenate(Xeqs)
		Yeq = np.concatenate(Yeqs)
		Teq = np.concatenate(Teqs)
		pID = np.concatenate(IDs)
		K   = np.concatenate(Ks)

		print("Saving collection to %s\n\n"%fOut)

		with open(fOut,"wb") as f:
			pickle.dump(X  ,f)
			pickle.dump(Y  ,f)
			pickle.dump(Z  ,f)
			pickle.dump(Xeq,f)
			pickle.dump(Yeq,f)
			pickle.dump(Teq,f)
			pickle.dump(pID,f)
			pickle.dump(K  ,f)

