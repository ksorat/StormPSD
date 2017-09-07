#Figures for RB comparison
import kCyl as kc
import pyStorm as pS
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import lfmViz as lfmv

eqDT = 600
eqT0 = 2000.0
vtiDir = os.path.expanduser('~') + "/Work/StormPSD/Data/eqSlc"
def getFld(t,eqStub="eqSlc"):
	tSlc = np.int( (t-eqT0)/eqDT )
	vtiFile = vtiDir + "/" + eqStub + ".%04d.vti"%(tSlc)
	#print("Reading %s"%vtiFile)

	dBz = lfmv.getVTI_SlcSclr(vtiFile).T
	ori,dx,ex = lfmv.getVTI_Eq(vtiFile)
	xi = ori[0] + np.arange(ex[0],ex[1]+1)*dx[0]
	yi = ori[1] + np.arange(ex[2],ex[3]+1)*dx[1]
	#print(xi.shape)
	#print(dBz.shape)

	return xi,yi,dBz
lfmv.ppInit()

#Time data
T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"

tMin = 33600.0
tMax = 189000.0

doSmooth = True
doField = False
Niter = 1

#RB Opts
rbStrs = ["A","B"]

#KC opts
kcStrs = ["KCyl_StormT","KCyl_StormI"]
kcScls = np.array([0.5,2.5])
LabFS = "large"
TitFS = "large"

#Ks = [100,250,500,750,1000]
Ks = [100,250,500,1000,1500]
mSize = 6
mTrk = 2
Ntrk = 10
Nskp = 15
alTrk = 0.5
doTracks = True
DelT = (50+6*60)*60 #Seconds to get to 3/17

#Ts = (60*60)*np.linspace(4,46,8) + DelT
Ts = (60*60)*np.linspace(4,46,8) + DelT

vNorm = LogNorm(vmin=1.0,vmax=1.0e+6)
Vc = np.linspace(-35,35,5)
vcNorm = Normalize(vmin=Vc.min(),vmax=Vc.max())

figSize = (18,12)
figQ = 300
#cMap = "jet"
cMap = "gnuplot2"
cMapC = "RdGy"
cAl = 0.5
cLW = 0.25
#cMap = "viridis"
rbAC = "cyan"
rbBC = "magenta"

#Get started
NumPop = len(kcStrs)
Nk = len(Ks)
Nt = len(Ts)
NumRB = len(rbStrs)


R,P,K,Tkc,Irpkt = pS.TotCyl(doSmooth=True)

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
			KLab = "%s MeV"%(str(K0/1000))
		else:
			KLab = "%s keV"%(str(K0))
		ik1 = (K>=K0).argmax()
		ik0 = ik1-1
		dk = (K0-K[ik0])/(K[ik1]-K[ik0])
		#Get single slice
		Ik = IKc[:,:,ik0] + dt*(IKc[:,:,ik1]-IKc[:,:,ik0])


		#Create plots
		Ax = fig.add_subplot(gs[k,n])
		Ax.pcolormesh(XX,YY,Ik,norm=vNorm,cmap=cMap)
		if (doField):
			#Get dBz data
			xi,yi,dBz = getFld(T0)
			Ax.contour(xi,yi,dBz,Vc,cmap=cMapC,alpha=cAl,linewidth=cLW)
		plt.axis('scaled')
		plt.xlim([-12.5,12.5])
		plt.ylim([-12.5,12.5])
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
					Ax.plot(Xrbs[0][iRB],Yrbs[0][iRB],color=rbAC,marker="o",markersize=mTrk,alpha=alTrk)
					Ax.plot(Xrbs[1][iRB],Yrbs[1][iRB],color=rbBC,marker="o",markersize=mTrk,alpha=alTrk)
if (doField):
	AxCC = fig.add_subplot(gs[-1,Nt/2:])
	cb = mpl.colorbar.ColorbarBase(AxCC,cmap=cMapC,norm=vcNorm,orientation='horizontal')
	cb.set_label("Residual Vertical Field [nT]",fontsize="large")
	AxC = fig.add_subplot(gs[-1,0:Nt/2])

else:
	AxC = fig.add_subplot(gs[-1,Nt/2-Nt/4:Nt/2+Nt/4+1])
cb = mpl.colorbar.ColorbarBase(AxC,cmap=cMap,norm=vNorm,orientation='horizontal')
cb.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")

plt.savefig("IPans.png",dpi=figQ)
plt.close('all')
lfmv.trimFig("IPans.png")



