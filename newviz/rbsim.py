#Test of kCyl reader
import kCyl as kc
import os
import numpy as np
import scipy
import scipy.interpolate
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import lfmViz as lfmv
lfmv.ppInit()

doPanel = True
doScl = True


T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"
tMin = 33600.0
tMax = 189000.0
Sig = -0.05

Nk = 35 #Number of K samples
Nsk = 5 #Skip number for trajectory

Stubs = ["KCyl_StormT","KCyl_StormI"]


if (doScl):
	fOut = "RBSim_Scl.png"
	IScl = [50,1]
else:
	fOut = "RBSim_NoScl.png"
	IScl = [1.0,1]

OrbF = "vaporbit.txt"
rbF = "rbspa.cdf"


Npop = len(IScl)

Labs = ["RBSP-A","Trapped","Injected","Trapped+Injected"]
#Now get RBSP data
Trb,Krb,dkrb,Irb = kc.GetRBSP(rbF,T0Str,tMin=tMin,tMax=tMax)

Is = [Irb]
Ts = [Trb]
Ks = [Krb]

#Interpolate from simulation
for n in range(Npop):
	fIn = os.path.expanduser('~') + "/Work/StormPSD/Data" + "/Merge/" + Stubs[n] + ".h5"

	#Interpolate from simulation
	R,P,K,t,I = kc.getCyl(fIn)
	if (Sig>0):
		I = kc.SmoothI(I,sig=Sig)
	Ii = kc.GetInterp(R,P,K,t,I)
	kMin = K.min()
	kMax = K.max()
	Ksc = np.logspace(np.log10(kMin),np.log10(kMax),Nk)
	#Ksc = np.linspace(kMin,kMax,Nk)
	
	#Get trajectory data
	#Tsc = seconds after T0
	
	Tsc,Xsc,Ysc,Z = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk=Nsk)
	Rsc = np.sqrt(Xsc**2.0 + Ysc**2.0)
	Psc = np.arctan2(Ysc,Xsc)
	iP = (Psc<0); Psc[iP] = Psc[iP]+2*np.pi
	
	Nsc = len(Tsc)
	iPts = np.zeros((Nk,4))

	#Loop over spacecraft trajectory, get I(K) at each point
	Isc = np.zeros((Nsc,Nk))
	for i in range(Nsc):
		r = Rsc[i]
		p = Psc[i]
		t = Tsc[i]
		iPts[:,0] = r
		iPts[:,1] = p
		iPts[:,2] = Ksc
		iPts[:,3] = t
		Isc[i,:] = Ii(iPts)

	Isc = IScl[n]*Isc
	Is.append(Isc)
	Ts.append(Tsc)
	Ks.append(Ksc)

Is.append(Is[1]+Is[2])
Ts.append(Tsc)
Ks.append(Ksc)

Lsc = np.sqrt(Xsc**2.0 + Ysc**2.0)
Psc = np.arctan2(Ysc,Xsc)*180.0/np.pi
I = Psc<0; Psc[I] = Psc[I] + 360
MLTsc = np.mod(Psc/15 + 12.0,24)

#Psc = (Psc/15.0)+12
#Do RB/SIM panel figure
if (doPanel):
	figSize = (12,16)
	figQ = 300 #DPI
	figName = "rbsimI.png"
	vMin = 1.0
	vMax = 1.0e+6
	
	cMap = "jet"

	fig = plt.figure(figsize=figSize)
	Np = len(Is)

	#Plots: DST,RBSP,SIM,BLANK,COLOR
	gs = gridspec.GridSpec(1+Np+1,1,height_ratios=[10,10,10,10,1,1])
	vNorm = LogNorm(vmin=vMin,vmax=vMax)


	#Both panels
	#Tp = [kc.Ts2date(Trb,T0Str),kc.Ts2date(Tsc,T0Str)]
	#Kp = [Krb,Ksc]
	#Ip = [Irb,Isc]

	for i in range(Np):
		Ax = fig.add_subplot(gs[i,0])
		Tp = kc.Ts2date(Ts[i],T0Str)

		iPlt = Ax.pcolormesh(Tp,Ks[i],Is[i].T,norm=vNorm,cmap=cMap)
		plt.ylim([50,5.0e+3])
		yStr = "%s\nEnergy [keV]"%Labs[i]
		plt.ylabel(yStr,fontsize="large")
		#
		plt.yscale('log')
		if (i<Np-1):
			plt.setp(Ax.get_xticklabels(),visible=False)
		else:
			#fig.autofmt_xdate()
			Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
			#plt.xlabel("Date",fontsize="large")

	#Do colorbar
	AxC = fig.add_subplot(gs[-1,0])
	cb = mpl.colorbar.ColorbarBase(AxC,cmap=cMap,norm=vNorm,orientation='horizontal')
	cb.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")
	plt.savefig(fOut,dpi=figQ)
	

# plt.pcolormesh(Tsc,Ksc,Isc.T,norm=vNorm,cmap=cMap)
# plt.ylim([50,5.0e+3])
# plt.ylabel("Energy [keV]",fontsize="large")
# plt.xlabel("Date",fontsize="large")
# plt.yscale('log')
# plt.colorbar()