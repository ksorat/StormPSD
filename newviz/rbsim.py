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
import numpy


lfmv.ppInit()

doPanel = True
doLine = True

doScl = False
doA = False
doEq = True

T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"
tMin = 33600.0
tMax = 189000.0
#tMax = 197000.0
Sig = -0.05
#Sig = 0.5

Nsk = 1 #Skip number for trajectory
Nk = 50 #Number of K samples
imeth = "linear"
#imeth = "nearest"

#Nk = 25 #Number of K samples
#Nsk = 10
Stubs = ["KCyl_StormT","KCyl_StormI"]

fOut = "RBSim"


if (doA):
	OrbF = "vaporbRBA.txt"
	rbF = "rbspa.cdf"
	fOut = fOut+"A"
	rbStr = "rbspa"
	Labs = ["RBSP-A","Trapped","Injected","Trapped+Injected"]
else:
	OrbF = "vaporbRBB.txt"
	rbF = "rbspb.cdf"
	fOut = fOut+"B"
	rbStr = "rbspb"
	Labs = ["RBSP-B","Trapped","Injected","Trapped+Injected"]

if (doScl):
	#IScl = [50,1]
	IScl = [12,1]

else:
	fOut = fOut+"_NoScl"
	IScl = [1.0,1]



Npop = len(IScl)


#Now get RBSP data
Trb,Krb,dkrb,Irb = kc.GetRBSP(rbF,T0Str,tMin=tMin,tMax=tMax,rbID=rbStr)

Is = [Irb]
Ts = [Trb]
Ks = [Krb]

#Interpolate from simulation
for n in range(Npop):
	fIn = os.path.expanduser('~') + "/Work/StormPSD/Data" + "/Merge/" + Stubs[n] + ".h5"

	#Interpolate from simulation
	R,P,K,Tkc,I = kc.getCyl(fIn)
	if (Sig>0):
		I = kc.SmoothI(I,sig=Sig)
	Ii = kc.GetInterp(R,P,K,Tkc,I,imeth=imeth)
	kMin = K.min()
	kMax = K.max()
	Ksc = np.logspace(np.log10(kMin),np.log10(kMax),Nk)
	#Ksc = np.linspace(kMin,kMax,Nk)
	
	#Get trajectory data
	#Tsc = seconds after T0
	
	Tsc,Xsc,Ysc,Z = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk=Nsk,doEq=doEq)
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
Np = len(Is)
vMin = 1.0
vMax = 1.0e+6

#Psc = (Psc/15.0)+12
#Do RB/SIM panel figure
if (doPanel):
	figSize = (12,16)
	#figSize = (24,32)
	figQ = 300 #DPI
	figName = "rbsimI.png"
	#vMin = 1.0
	#vMin = 5
	
	cMap = "jet"
	#cMap = "viridis"
	fig = plt.figure(figsize=figSize)

	#Plots: DST,RBSP,SIM,BLANK,COLOR
	gs = gridspec.GridSpec(1+Np+1+1,1,height_ratios=[10,10,10,10,1,5,1])
	vNorm = LogNorm(vmin=vMin,vmax=vMax)


	#Both panels
	#Tp = [kc.Ts2date(Trb,T0Str),kc.Ts2date(Tsc,T0Str)]
	#Kp = [Krb,Ksc]
	#Ip = [Irb,Isc]

	for i in range(Np):
		Ax = fig.add_subplot(gs[i,0])
		Tp = kc.Ts2date(Ts[i],T0Str)
		#Tp = Ts[i]
		iPlt = Ax.pcolormesh(Tp,Ks[i],Is[i].T,norm=vNorm,cmap=cMap)
		plt.ylim([50,5.0e+3])
		yStr = "%s\nEnergy [keV]"%Labs[i]
		plt.ylabel(yStr,fontsize="large")
		#
		plt.yscale('log')
		if (i==0):
			Ax.xaxis.tick_top()
			Ax.xaxis.set_label_position('top')
			Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
		elif (i<Np-1):
			plt.setp(Ax.get_xticklabels(),visible=False)
		else:
			Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))

	#Do colorbar
	AxC = fig.add_subplot(gs[-1,0])
	cb = mpl.colorbar.ColorbarBase(AxC,cmap=cMap,norm=vNorm,orientation='horizontal')
	cb.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")
	
	
	#Do L/MLT
	AxP = fig.add_subplot(gs[-2,0])
	Tp = kc.Ts2date(Tsc,T0Str)
	AxP.plot(Tp,Lsc,'b')
	plt.setp(AxP.get_xticklabels(),visible=False)
	AxP.set_ylabel("L",color='b')
	AxP.tick_params('y',colors='b')

	AxP2 = AxP.twinx()
	AxP2.plot(Tp,MLTsc,'r.')
	AxP2.set_yticks([0,6,12,18,24])
	AxP2.set_ylabel("MLT",color='r')
	AxP2.tick_params('y',colors='r')
	AxP2.set_ylim([0,24])
	#Save
	fOutP = fOut+"_Itk.png"
	plt.savefig(fOutP,dpi=figQ)
	# xT = Ax.get_xticklabels()
	# Nxt = len(xT)
	# for n in range(Nxt):
	# 	xTS = str(xT[n].get_text())
	# 	print(xTS)
	# 	xTS = xTS + "\nB"
	# 	print(xTS)
	# 	xT[n].set_text(unicode(xTS))
	# Ax.set_xticklabels(xT)
	# plt.savefig(fOut,dpi=figQ)
	plt.close('all')
#Show >MeV I
if (doLine):
	K0 = 1000.0
	kR = 200.0
	figSize = (24,8)
	figQ = 300 #DPI
	figName = "rbsimMeV.png"
	plt.close('all')
	vMin = 1.0e-1
	vMax = 1.0e+4
	fig = plt.figure(1,figsize=figSize)
	LW = 1.5
	Ws = 10
	ColI = ['b-','g-','m-','r-']

	for i in range(Np):
		Tp = kc.Ts2date(Ts[i],T0Str)
		#Find critical value of K>1MeV
		kC = (Ks[i]>K0).argmax()
		kC1 = (Ks[i]>=(K0-kR)).argmax()
		kC2 = (Ks[i]<=(K0+kR)).argmin()
		Ic = (Is[i][:,kC1:kC2+1]).mean(axis=1)
		print(i,kC,kC1,kC2)
		if (i==0):
			Ics = Ic
		else:
			Ics = Ic
		#Ic = Is[i][:,kC] #Point value
		# if (i == 0):
		# 	
		# else:
		# 	Ic = (Is[i][:,kC1:kC2]).mean(axis=1)
		plt.semilogy(Tp,Ics,ColI[i],linewidth=LW,label=Labs[i])

	plt.ylim([vMin,vMax])
	plt.ylabel("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")
	plt.legend(loc='lower right')
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
	fOutL = fOut+"_I1.png"
	plt.savefig(fOutL,dpi=figQ)

# plt.pcolormesh(Tsc,Ksc,Isc.T,norm=vNorm,cmap=cMap)
# plt.ylim([50,5.0e+3])
# plt.ylabel("Energy [keV]",fontsize="large")
# plt.xlabel("Date",fontsize="large")
# plt.yscale('log')
# plt.colorbar()