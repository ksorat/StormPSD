#Generates spacecraft spectra

import sys
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
from visit import *
from visit_utils import *
from visit_utils.common import lsearch #lsearch(dir(),"blah")
import pyVisit as pyv
import kCyl as kc
import cPickle as pickle

T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"
kDB = os.path.expanduser('~') + "/Work/StormPSD/Data" + "/Merge/KCyl.xmf"
OrbF = "vaporbit.txt"
oFile = "vaporb.pkl"

Nk = 100 #Number of K samples
Nsk = 10

IBds = [1.0e-4,1.0e+2]
cMap = "magma"

#Get KCyl data
Quiet = True
if (Quiet):
	LaunchNowin()
else:
	Launch()

#Do some defaults
pyv.pvInit()

#Get time data
md0 = GetMetaData(kDB)

dt = md0.times[1] - md0.times[0]
dT0s = md0.times[0]

dT0 = datetime.timedelta(seconds=dT0s)
T0 = dT0 + datetime.datetime.strptime(T0Str,T0Fmt)

tV = np.array(md0.times)
tMin = tV.min()
tMax = tV.max()

#Open spacecraft trajectory data
#Tsc = seconds after T0

Tsc,Xsc,Ysc,Z = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk)

#Prep data
Nsc = Tsc.shape[0]
Isc = np.zeros((Nsc,Nk))


OpenDatabase(kDB)
DefineScalarExpression("L","cylindrical_radius(mesh)")
DefineScalarExpression("K","10^coord(mesh)[2]")
DefineScalarExpression("LogK","log10(K)")
DefineScalarExpression("IScl","2.9979E10*recenter(I)")

#Create Intensity plot to pull from
pyv.lfmPCol(kDB,"IScl",vBds=IBds,Log=True,cMap=cMap)
DrawPlots()

#Get max/min LogK
Query("MinMax","LogK")
kBds = GetQueryOutputValue()
lK = np.linspace(kBds[0],kBds[1],Nk)


#Loop over spacecraft trajecory points
for n in range(Nsc):
	#Get location
	x0 = Xsc[n]
	y0 = Ysc[n]
	R0 = np.sqrt(x0**2.0+y0**2.0)
	if (R0>2.05):
		#Get matching VisIt time slice
		tn = Tsc[n]
		tVsc = np.abs(tV-tn).argmin()
	
		SetTimeSliderState(tVsc)
		#Get data
		lS = (x0,y0,kBds[0])
		lE = (x0,y0,kBds[1])
	
		#Pull out relevant piece
		Lineout(lS,lE,Nk)
		SetActiveWindow(2)
		pInfo = GetPlotInformation()
		Ik = np.array(pInfo['Curve'])[1:-1:2]
	
		Isc[n,:] = Ik
		print(tn)
		#print(Ik)
		DeleteWindow()
		SetActiveWindow(1)
	else:
		Isc[n,:] = 0.0

#Write pickle
with open(oFile, "wb") as f:
	pickle.dump(T0Str,f)
	pickle.dump(T0Fmt,f)
	pickle.dump(Tsc,f)
	pickle.dump(lK,f)
	pickle.dump(Isc,f)
	pickle.dump(Xsc,f)
	pickle.dump(Ysc,f)



