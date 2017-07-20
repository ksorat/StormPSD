#Make plots showing intensity/pitch angle access for RB due to latitude
import kCyl as kc
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

#Time data
T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"
Labs = ["NULL","Trapped","Injected","Combined"]

tMin = 33600.0
tMax = 189000.0
figQ = 300

#tMax = 195000.0

#RB Opts
rbStrs = ["A","B"]

NumRB = len(rbStrs)

aT = []
aAeC = []
aIscl = []
aL = []
aLCut = []
for nrb in range(NumRB):
	#Create relevant files
	rbStr = rbStrs[nrb]

	OrbF = "vaporbRB" + rbStr.upper() + ".txt"

	#Get projected/non-proj trajectory
	#Tsc = seconds after T0
	#Non-Proj
	Tsc,Xsc,Ysc,Zsc = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk=1,doEq=False)
	Tp,Xp,Yp,Zp = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk=1,doEq=True)
	L = np.sqrt(Xp**2.0+Yp**2.0)
	r = np.sqrt(Xsc**2.0+Ysc**2.0+Zsc**2.0)
	mlat = np.arcsin(Zsc/r)
	cL = np.cos(mlat)
	sL = np.sin(mlat)

	lArg = cL**6.0/np.sqrt(1+3*sL*sL)
	AeC = np.arcsin(np.sqrt(lArg))
	IScl = (AeC-0.5*np.sin(2*AeC))/(0.5*np.pi)
	LCut = (L>=2.0)

	aT.append(Tsc)
	aAeC.append(AeC*180.0/np.pi)
	aIscl.append(IScl)
	aL.append(L)
	aLCut.append(LCut)
	#print("Lat/AeC = %f,%f\n"%(mlat*180.0/np.pi,AeC))

Ia = aLCut[0]
Ib = aLCut[1]

L0 = L.max()
Tp = kc.Ts2date(Tsc,T0Str)
figSize = (12,6)
Leg = ["RB-A","RB-B"]
#Critical pitch angle
plt.figure(figsize=figSize)
plt.plot(Tp,aAeC[0],'r',Tp,aAeC[1],'b')
#plt.plot(Tp,aL[0]/L0,'r--',Tp,aL[1]/L0,'b--')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
plt.ylabel('Maximum Accessible Pitch Angle')
plt.legend(Leg)
plt.savefig("AeC.png",dpi=figQ)

plt.close('all')

#IScl fraction
plt.figure(figsize=figSize)
plt.plot(Tp[Ia],aIscl[0][Ia],'r',Tp[Ib],aIscl[1][Ib],'b')
plt.plot(Tp,aL[0]/L0,'r--',Tp,aL[1]/L0,'b--')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
plt.ylabel('Fraction of total Sin^2 Intensity')
plt.legend(Leg)
plt.savefig("IScl.png",dpi=figQ)
plt.close('all')
