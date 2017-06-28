#Just generate a CSV for the time series data
import numpy as np
import glob
import os
import cPickle as pickle
from sys import exit

#Root dir
rDir = os.path.expanduser('~') + "/Work/StormPSD/"

fStub = ["21","3","0"]
NumW = 3
for nw in range(NumW):
	fPkl = "tsWedge_%s.pkl"%(fStub[nw])
	fTab = "tsWedge_%s.csv"%(fStub[nw])


	print("Loading Wedge TS")
	with open(fPkl,"rb") as f:
		t = pickle.load(f)
		Vst = pickle.load(f)
		kTt = pickle.load(f)
		Nt = pickle.load(f)
		Nkt = pickle.load(f)	
	
	N = t.shape[0]
	dOut = np.zeros((3,N))
	dOut[0,:] = t
	dOut[1,:] = Nt
	dOut[2,:] = kTt
	
	np.savetxt(fTab,dOut.T,delimiter=',')