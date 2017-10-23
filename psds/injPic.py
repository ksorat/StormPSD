#Figures for Injection rate

import os
import kCyl as kc
import numpy as np
import cPickle as pickle
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import lfmViz as lfmv


def tWindow(t,Q,dt):
	#Window time series t,Q based on window size dt
	Nt = len(t)
	Qw =  np.zeros(Nt)
	Qw[:] = Q[:]
	J = (Q>0)
	for i in range(Nt):
		t0 = t[i]
		I = (np.abs(t-t0) <= dt)
		IJ = I & J
		if (IJ.sum() > 0):
			Qw[i] = Q[IJ].mean()
		else:
			Qw[i] = 0.0
		
	return Qw

T0Cut = 40000.0
lfmv.ppInit()
doSmoothTS = True
figQ = 300
#Time data
T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"

dR_W = 3 #Wedge radial length [Re]
ReKM = 6.38e+3

tMin = 33600.0
tMax = 189000.0
tsID = ["0","21","3"]
tsLabs = ["Midnight","2100","0300"]
#dt = 3000.0
#dt = 600.0
#dtW = 3000.0
dt = 150.0*2
dtW = dt

NumPop = len(tsID)

aV = []
aT = []
aN = []
aNk = []


for n in range(NumPop):
	
	fPkl = "tsWedge_%s.pkl"%(tsID[n])
	fTab = "tsWedge_%s.csv"%(tsID[n])
	print("Creating CSV for injection %s w/ dt=%f"%(tsID[n],dt))
	with open(fPkl,"rb") as f:
		t = pickle.load(f) #Time [s]
		Vst = pickle.load(f) #Earthward tail velocity [km/s]
		kTt = pickle.load(f) #Thermal energy, kT [keV]
		Nt  = pickle.load(f) #Number density, [#/cm3]
		Nkt = pickle.load(f) #Number density above Kcrit

	N = t.shape[0]
	dOut = np.zeros((3,N))
	dOut[0,:] = t

	wVst = tWindow(t,Vst,dtW)
	wkTt = tWindow(t,kTt,dtW)
	wNt  = tWindow(t,Nt ,dtW)
	wNkt = tWindow(t,Nkt,dtW)

	if (doSmoothTS):
		nScl = (dt*wVst)/(dR_W*ReKM)
		dOut[1,:] = nScl*wNt
		dOut[2,:] = wkTt
		aV.append(wVst)
		aT.append(wkTt)
		aN.append(wNt)
		aNk.append(wNkt)
	else:
		nScl = (dt*Vst)/(dR_W*ReKM)
		dOut[1,:] = nScl*Nt
		dOut[2,:] = kTt
		aV.append(Vst)
		aT.append(kTt)
		aN.append(Nt)
		aNk.append(Nkt)

Tp = kc.Ts2date(t,T0Str)
figSize = (12,8)

plt.figure(figsize=figSize)
for n in range(NumPop):
	plt.plot(Tp,aV[n])
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))			
plt.legend(tsID)
plt.ylabel("Convection Velocity")
#Save and close
plt.savefig("WedgeV.png",dpi=figQ)
plt.close('all')
lfmv.trimFig("WedgeV.png")

plt.figure(figsize=figSize)
for n in range(NumPop):
	plt.plot(Tp,aT[n])
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))			
plt.legend(tsID)
plt.ylabel("Temperature [keV]")
#Save and close
plt.savefig("WedgeT.png",dpi=figQ)
plt.close('all')
lfmv.trimFig("WedgeT.png")

plt.figure(figsize=figSize)
for n in range(NumPop):
	plt.plot(Tp,aN[n])
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))			
plt.legend(tsID)
plt.ylabel("Number Density [1/cm3]")
#Save and close
plt.savefig("WedgeN.png",dpi=figQ)
plt.close('all')
lfmv.trimFig("WedgeN.png")

plt.figure(figsize=figSize)
for n in range(NumPop):
	Vst = aV[n]
	I = (t<=T0Cut)
	Vst[I] = 0.0
	pInj = dt*Vst*aNk[n]/(dR_W*ReKM)
	plt.plot(Tp,pInj,linewidth=0.75)
	T0Str = "2013-03-16T17:10:00Z"
datemin = datetime.datetime.strptime("2013-03-17T04:00:00Z",T0Fmt)
datemax = datetime.datetime.strptime("2013-03-18T06:00:00Z",T0Fmt)

plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))			
plt.gca().set_xlim(datemin,datemax)
plt.legend(tsLabs)
plt.ylabel("Injection Rate")
#Save and close
plt.savefig("WedgeInj.png",dpi=figQ)
plt.close('all')
lfmv.trimFig("WedgeInj.png")



