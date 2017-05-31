import numpy as np
import lfmPostproc as lfmpp
import cPickle as pickle
import os

doInj = False

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

