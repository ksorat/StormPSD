#Get time series from a bunch of LFM HDFs
import numpy as np
import glob
import os
import cPickle as pickle
import lfmPreproc as lfmpre

#Wedge info
R = [12,12.5]
Z = [-0.25,0.25]
P = [160,180]
K0 = 10

#HDF directory
lfmDir = os.path.expanduser('~') + "/Work/StormPSD/lfmData"

#Output pickle
fOut = "tsWedge.pkl"

fIns =glob.glob(lfmDir + "/*.hdf")

Nf = len(fIns)
Nf = 10

t = np.zeros(Nf)
Vst = np.zeros(Nf)
kTt = np.zeros(Nf)
Nt = np.zeros(Nf)
Nkt = np.zeros(Nf)

I,dvI = lfmpre.lfmWedge(fIns[0],R=R,P=P,Z=Z)

for i in range(Nf):
	fIn = fIns[i]

	t[i],Vst[i],kTt[i],Nt[i],Nkt[i] = lfmpre.injWedge(fIn,I,dvI,K0=K0)

#Got all data, now sort
Is = np.argsort(t)

t = t[Is]
Vst = Vst[Is]
kTt = kTt[Is]
Nt = Nt[Is]
Nkt = Nkt[Is]

print("Saving data to %s"%fOut)
with open(fOut,"wb") as f:
	pickle.dump(t,f)
	pickle.dump(Vst,f)
	pickle.dump(kTt,f)
	pickle.dump(Nt,f)
	pickle.dump(Nkt,f)
