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

T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"
kDB = os.path.expanduser('~') + "/Work/StormPSD/Data" + "/Merge/KCyl.xmf"
Nk = 100 #Number of K samples

IBds = [1.0e-4,1.0e+2]
cMap = "magma"

#Open spacecraft trajectory data
#Tsc = seconds after T0

Tsc = np.linspace(30000,80000,10)
Xsc = 5*np.sin(2*np.pi*Tsc/25000.0)
Ysc = 5*np.cos(2*np.pi*Tsc/25000.0)
Nsc = Tsc.shape[0]
Isc = np.zeros((Nsc,Nk))


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
	#Get matching VisIt time slice
	tV = n+5 #TODO: Fix
	SetTimeSliderState(tV)
	#Get data
	lS = (x0,y0,kBds[0])
	lE = (x0,y0,kBds[1])

	#Pull out relevant piece
	Lineout(lS,lE,Nk)
	SetActiveWindow(2)
	pInfo = GetPlotInformation()
	Ik = np.array(pInfo['Curve'])[1:-1:2]

	Isc[n,:] = Ik
	DeleteWindow()
	SetActiveWindow(1)




