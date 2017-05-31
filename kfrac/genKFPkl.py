#Creates PKL with data for initial/final K figures

import numpy as np
import os
import lfmPostproc as lfmpp
import cPickle as pickle
import glob

RootDir = os.path.expanduser('~') + "/Work/StormPSD/Data/"
doInj = True

if (doInj):
	Stubs = ["Inj","Inj1"]
	fOut = "kfInj.pkl"
else:
	Stubs = ["Trap","TrapC"]
	fOut = "kfTrap.pkl"

Ns = len(Stubs)

#Loop over stubs, and each file within each directory
K0 = []
Kf = []
In = []
T0 = []
A0 = []
Af = []

for n in range(Ns):
	fDir = RootDir+Stubs[n]
	h5ps = glob.glob(fDir+"/*.h5part")
	for h5p in h5ps:
		print("\tReading %s"%(h5p))
		pid,kev0 = lfmpp.getH5pInit(h5p,"kev")
		pid,kevF = lfmpp.getH5pFin (h5p,"kev")
		pid,inP  = lfmpp.getH5pFin (h5p,"in")
		pid,T0p  = lfmpp.getH5pInit(h5p,"T0p")
		pid,a0p  = lfmpp.getH5pInit(h5p,"alpheq")
		pid,aFp  = lfmpp.getH5pFin (h5p,"alpheq")

		Npp = kev0.shape[0]
		for p in range(Npp):
			K0.append(kev0[p])
			Kf.append(kevF[p])
			In.append(inP[p])
			T0.append(T0p[p])
			A0.append(a0p[p])
			Af.append(aFp[p])

K0 = np.array(K0)
Kf = np.array(Kf)
In = np.array(In)
T0 = np.array(T0)
A0 = np.array(A0)
Af = np.array(Af)

print("Total Particles = %d\n"%(K0.shape[0]))
#Trim out bad partices
I = ~np.isnan(Kf)
K0 = K0[I]
Kf = Kf[I]
In = In[I]
T0 = T0[I]
A0 = A0[I]
Af = Af[I]

print("After Trimming = %d\n"%(K0.shape[0]))
print("Saving data to %s\n"%(fOut))

with open(fOut,"wb") as f:
	pickle.dump(K0,f)
	pickle.dump(Kf,f)
	pickle.dump(In,f)
	pickle.dump(T0,f)
	pickle.dump(A0,f)
	pickle.dump(Af,f)
