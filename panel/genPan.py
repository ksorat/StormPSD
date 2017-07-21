#Figures for RB comparison
import kCyl as kc
import os
import numpy as np
import scipy
import scipy.interpolate
import scipy.ndimage
from scipy.ndimage.filters import gaussian_filter1d
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import lfmViz as lfmv

lfmv.ppInit()
#Time data
T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"

tMin = 33600.0
tMax = 189000.0

doSmooth = True
Niter = 2

#RB Opts
rbStrs = ["A","B"]

#KC opts
kcStrs = ["KCyl_StormT","KCyl_StormI"]
kcScls = np.pi*4*np.array([1.0,10.0])
LabFS = "large"
TitFS = "large"

Ks = [100,250,500,750,1000]
mSize = 6
mTrk = 2
Ntrk = 10
Nskp = 6

doTracks = True
Ts = np.linspace(35000.0,185000.0,6)
vNorm = LogNorm(vmin=1.0,vmax=1.0e+6)
figSize = (12,12)
figQ = 300
cMap = "jet"
cMap = "viridis"
rbAC = "cyan"
rbBC = "magenta"

#Get started
NumPop = len(kcStrs)
Nk = len(Ks)
Nt = len(Ts)
NumRB = len(rbStrs)

aI = []
#Get KCyls (assuming same grid)
for n in range(NumPop):
	fIn = os.path.expanduser('~') + "/Work/StormPSD/Data" + "/Merge/" + kcStrs[n] + ".h5"
	R,P,K,Tkc,I0 = kc.getCyl(fIn)
	I0 = kcScls[n]*I0

	#Smooth cylinder
	if (doSmooth):
		Irpkt = kc.SmoothKCyl(I0,Niter)
	else:
		Irpkt = I0

	aI.append(Irpkt)

#Total intensity
Irpkt = aI[0]+aI[1]

#Get grid
XX,YY = kc.xy2rp(R,P)

#Get RB trajectories
Xrbs = []
Yrbs = []
for n in range(NumRB):
	rbStr = rbStrs[n]
	OrbF = "vaporbRB" + rbStr.upper() + ".txt"
	Tsc,Xsc,Ysc,Zsc = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk=1,doEq=True)
	Xrbs.append(Xsc)
	Yrbs.append(Ysc)

#Prep for figure
fig = plt.figure(figsize=figSize)
Nkp = 2
HRs = np.ones(Nk+Nkp)
HRs[Nk:] = 0.1

gs = gridspec.GridSpec(Nk+Nkp,Nt,height_ratios=HRs)

Tp = kc.Ts2date(Ts,T0Str)

#Loop over time/energies, find bounds and lin interpolate
for n in range(Nt):
	T0 = Ts[n]
	it1 = (Tkc>=T0).argmax()
	it0 = it1-1
	dt = (T0-Tkc[it0])/(Tkc[it1]-Tkc[it0])

	#Get single KCyl
	IKc = Irpkt[:,:,:,it0] + dt*(Irpkt[:,:,:,it1]-Irpkt[:,:,:,it0])

	#Get RB points at this time
	itRB = np.abs(Tsc-T0).argmin()
	Xrba = Xrbs[0][itRB]
	Yrba = Yrbs[0][itRB]
	Xrbb = Xrbs[1][itRB]
	Yrbb = Yrbs[1][itRB]

	for k in range(Nk):
		K0 = Ks[Nk-k-1]
		if (K0>=1000):
			KLab = "%d MeV"%(np.round(K0/1000))
		else:
			KLab = "%d keV"%(np.round(K0))
		ik1 = (K>=K0).argmax()
		ik0 = ik1-1
		dk = (K0-K[ik0])/(K[ik1]-K[ik0])
		#Get single slice
		Ik = IKc[:,:,ik0] + dt*(IKc[:,:,ik1]-IKc[:,:,ik0])

		Ax = fig.add_subplot(gs[k,n])
		Ax.pcolormesh(XX,YY,Ik,norm=vNorm,cmap=cMap)
		plt.axis('scaled')
		plt.xlim([-12.5,12.5])
		plt.ylim([-15,15])
		lfmv.addEarth2D()
		if (n == 0):
			plt.ylabel(KLab,fontsize=LabFS)
		elif (n == Nt-1):
			plt.ylabel("GSM-Y [Re]",fontsize=LabFS)
			Ax.yaxis.tick_right()
			Ax.yaxis.set_label_position("right")
		else:
			plt.setp(Ax.get_yticklabels(),visible=False)
	
		if (k < Nk-1):
			plt.setp(Ax.get_xticklabels(),visible=False)
		else:
			plt.xlabel('GSM-X [Re]',fontsize=LabFS)
		if (k == 0):
			nT = Tp[n]
			TLab = "%s\n%s"%(str(nT.time()),str(nT.date()))
			plt.title(TLab,fontsize=TitFS)

		#Plot RB points
		Ax.plot(Xrba,Yrba,color=rbAC,marker="o",markersize=mSize)
		Ax.plot(Xrbb,Yrbb,color=rbBC,marker="o",markersize=mSize)
		#Plot RB tracks
		if (doTracks):
			for i in range(Ntrk):
				iRB = itRB-(i+1)*Nskp
				if (iRB<0):
					continue
				else:
					Ax.plot(Xrbs[0][iRB],Yrbs[0][iRB],color=rbAC,marker="o",markersize=mTrk)
					Ax.plot(Xrbs[1][iRB],Yrbs[1][iRB],color=rbBC,marker="o",markersize=mTrk)
AxC = fig.add_subplot(gs[-1,0:Nt/2])
cb = mpl.colorbar.ColorbarBase(AxC,cmap=cMap,norm=vNorm,orientation='horizontal')
cb.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")

plt.savefig("IPans.png",dpi=figQ)
plt.close('all')



