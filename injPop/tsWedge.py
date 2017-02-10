#Get time series from a bunch of LFM HDFs
import numpy as np
import glob
import os
import lfmPreproc as lfmpre
#Wedge info

R = [12,12.5]
Z = [-0.25,0.25]
P = [160,180]
K0 = 10

#HDF directory
lfmDir = os.path.expanduser('~') + "/Work/StormPSD/lfmData"
fIns =glob.glob(lfmDir + "/*.hdf")

Nf = len(fIns)

t = np.zeros(Nf)
Vst = np.zeros(Nf)
kTt = np.zeros(Nf)
Nt = np.zeros(Nf)
Nkt = np.zeros(Nf)

for i in range(Nf):
	fIn = fIns[i]
	t[i],Vst[i],kTt[i],Nt[i],Nkt[i] = lfmpre.lfmWedge(fIn,R=R,P=P,Z=Z,K0=K0)
