import numpy as np
import lfmPostproc as lfmpp
import cPickle as pickle
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
doInj = True

if (doInj):
	fPkl = "kfInj.pkl"
else:
	fPkl = "kfTrap.pkl"

#Load from pickle
if (os.path.isfile(fPkl)):
	print("Loading data from %s"%(fPkl))
	with open(fPkl,"rb") as f:
		K0 = pickle.load(f)
		Kf = pickle.load(f)
		In = pickle.load(f)
		T0 = pickle.load(f)
		A0 = pickle.load(f)
		Af = pickle.load(f)


Nb0 = 50
NbF = 100
I = (In>0.5)
xH = K0[I]
yH = Kf[I]

vNorm=LogNorm(vmin=1.0,vmax=1.0e+3)

xb = np.logspace(np.log10(10),np.log10(5000),Nb0)
yb = np.logspace(np.log10(10),np.log10(5000),NbF)
#xb = np.linspace(10,1000,Nb0)
#yb = np.linspace(10,10000,NbF)

plt.hist2d(xH,yH,[xb,yb],norm=vNorm)
#plt.hist2d(xH,yH,40,norm=vNorm)
plt.xscale('log')
plt.yscale('log')
plt.colorbar()
plt.show()
