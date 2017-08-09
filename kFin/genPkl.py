#Creates PKL with data for initial/final K figures

import numpy as np
import os
import lfmPostproc as lfmpp
import cPickle as pickle
import glob

RootDir = os.path.expanduser('~') + "/Work/StormPSD/Data/"
doInj = True
KCut = 1000 #Cutoff for interesting particles

if (doInj):
	Stubs = ["Inj0","Inj21","Inj3"]
else:
	Stubs = ["Trap","TrapC"]


Ns = len(Stubs)
#Loop over stubs, and each file within each directory

for n in range(Ns):
	fDir = RootDir+Stubs[n]
	h5ps = glob.glob(fDir+"/*.h5part")
	fOut = "%s.pkl"%Stubs[n]
	print("Starting collection %s"%Stubs[n])
	#Create initial holders
	K0 = np.array([])
	Kf = np.array([])
	In = np.array([])
	T0 = np.array([])
	A0 = np.array([])
	Af = np.array([])
	IDs = np.array([])

	for h5p in h5ps:
		print("\tReading %s"%(h5p))
		pid,kev0 = lfmpp.getH5pInit(h5p,"kev")
		pid,kevF = lfmpp.getH5pFin (h5p,"kev")
		pid,inP  = lfmpp.getH5pFin (h5p,"in")
		pid,T0p  = lfmpp.getH5pInit(h5p,"T0p")
		pid,a0p  = lfmpp.getH5pInit(h5p,"alpheq")
		pid,aFp  = lfmpp.getH5pFin (h5p,"alpheq")
		
		#Kill bad energy particles
		I0 = np.isnan(kevF)
		kevF[I0] = 0.0

		#Find remaining particles with K above KCut
		I = (inP>0.5) & (kevF>=KCut)
		#print("\t\tKilling %d particles"%(I0.sum()))
		print("\tFound %d particles"%(I.sum()))		
		K0 = np.append(K0,kev0[I])
		Kf = np.append(Kf,kevF[I])
		T0 = np.append(T0,T0p[I])
		A0 = np.append(A0,a0p[I])
		Af = np.append(Af,aFp[I])
		IDs= np.append(IDs,pid[I])

	#Save data for this collection
	print("Found %d total particles"%(K0.shape))
	print("Saving collection to %s\n\n"%fOut)

	with open(fOut,"wb") as f:
		pickle.dump(K0 ,f)
		pickle.dump(Kf ,f)
		pickle.dump(T0 ,f)
		pickle.dump(A0 ,f)
		pickle.dump(Af ,f)
		pickle.dump(IDs,f)
