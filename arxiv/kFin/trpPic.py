import numpy as np
import lfmPostproc as lfmpp
import cPickle as pickle
import os
import datetime
import matplotlib.pyplot as plt
import kCyl as kc
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import lfmViz as lfmv
import lfmPostproc as lfmpp


#Time data
T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"

Base = os.path.expanduser('~') + "/Work/StormPSD/Data/"
Stubs = ["Inj0","Inj21","Inj3"]
NSt = ["0","21","3"]
Leg = ["Midnight","2100","0300"]
Ns = len(Stubs)

Rcut = 8.0
doRCut = True
IDCut = 20000
#IDPer = 20000
NpPlt = 10
NtSk = 1
tMax = 4

#Loop over stubs
K0s = []
KFs = []
T0s = []
Rfs = []
AFs = []
AIs = []

fig = plt.figure(figsize=(8,8))
for n in range(0,3):
	fPkl = "%s.pkl"%Stubs[n]
	print("Reading %s"%(fPkl))
	with open(fPkl,"rb") as f:	
		K0  = pickle.load(f)
		Kf  = pickle.load(f)
		T0  = pickle.load(f)
		A0  = pickle.load(f)
		Af  = pickle.load(f)
		Rf  = pickle.load(f)
		IDs = pickle.load(f)
	if (doRCut):
		Np = len(K0)
		print("\tTotal particles before cut = %d"%Np)

		I = (Rf<=Rcut) & (IDs <= IDCut)
		Np = I.sum()
		print("\tAfter cut at %2.2f, Np = %d"%(Rcut,Np))
		K0  = K0 [I]
		Kf  = Kf [I]
		T0  = T0 [I]
		A0  = A0 [I]
		Af  = Af [I]
		Rf  = Rf [I]
		IDs = IDs[I]

	

	
	for i in range(0,Np):
		IDi = IDs[i]
		h5p = Base + Stubs[n] + "/StormInj_%s.0001.h5part"%(NSt[n])
		# if (IDi>IDPer):
		# 	h5p = Base + Stubs[n] + "/StormInj_%s.0002.h5part"%(NSt[n])
		# 	IDi = IDi-IDPer
		# else:
		# 	h5p = Base + Stubs[n] + "/StormInj_%s.0001.h5part"%(NSt[n])

		t,xeq = lfmpp.getH5pid(h5p,"xeq",IDi)
		t,yeq = lfmpp.getH5pid(h5p,"yeq",IDi)
		t,isIn = lfmpp.getH5pid(h5p,"in",IDi)
		t0 = (isIn>0.5).argmax()-1
		plt.plot(xeq[t0:t0+tMax:NtSk],yeq[t0:t0+tMax:NtSk],linewidth=0.5)
		plt.plot(xeq[t0],yeq[t0],'o')
plt.axis('equal')
plt.show()
