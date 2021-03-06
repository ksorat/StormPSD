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

#Average over time window
def TWin(Ik,Nw=4):
	Nt = Ik.shape[0]
	IkS = np.zeros(Nt)
	IkS[:] = Ik[:]
	if (Nw>0):
		for i in range(Nt):
			i0 = np.maximum(i-Nw,0)
			i1 = np.minimum(i+Nw,Nt-1)
			IkS[i] = Ik[i0:i1].mean()
	return IkS
#Time data
T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"
Labs = ["NULL","Trapped","Injected","Combined"]

tMin = 33600.0
tMax = 189000.0
#tMax = 195000.0

#RB Opts
rbStrs = ["A","B"]
rbSK = 1 #Skip cadence for RB intensity
rbSKt = 1 #Skip cadence for RB trajectory

#KC opts
doZScl = True #Attenuate for Z
En = 2
kcStrs = ["KCyl_StormT","KCyl_StormI"]
dOm = np.pi*4
#kcScls = np.pi*4*np.array([2.0,2.0])
kcScls = np.array([1/(4*np.pi),2.5])

#Figure opts
doSmoothFig = True
Niter = 1
NiterT = 1
NTWin = 2
#NTWin = 2


doPanelFig = True
doLineFig = True
#KLs = [2000,1000,750,500,250]
KLs = [3000,2000,1000,500,250,50]
vMins = [1.0e-2,1.0e-1,1.0e+0,1.0e+1,1.0e+2,1.0e+3]
vMaxs = [1.0e+2,1.0e+3,1.0e+4,1.0e+5,1.0e+6,1.0e+7]

#KLs = [2500,1000,800,600,200]
figQ = 300 #DPI
cMap = "jet"

NumPop = len(kcStrs)
NumRB = len(rbStrs)
#NumRB = 1

for nrb in range(NumRB):
	aI = []
	aT = []
	aK = []

	#Create relevant files
	rbStr = rbStrs[nrb]

	OrbF = "vaporbRB" + rbStr.upper() + ".txt"
	rbFs = "rbsp" + rbStrs[nrb].lower()
	rbF  = rbFs+".cdf"
	Labs[0] = "RBSP-" + rbStr.upper()

	#Get RB intensity
	print("Reading from %s"%rbStr)
	Trb,Krb,dkrb,Irb = kc.GetRBSP(rbF,T0Str,tMin=tMin,tMax=tMax,rbID=rbFs,rbSK=rbSK)
	#Krb = KrbT[0,:]
	aI.append(Irb)
	aT.append(Trb)
	aK.append(Krb)

	Nec = Krb.shape[0]
	print("Energies = %s"%str(Krb))
	print("RB Cadence = %f"%(Trb[1]-Trb[0]))

	#Get RB trajectory data
	#Tsc = seconds after T0
	#Get full 3D trajectory
	Tsc,Xsc,Ysc,Zsc = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk=rbSKt,doEq=False)

	#Get simulation K-Cyls
	for n in range(NumPop):
		fIn = os.path.expanduser('~') + "/Work/StormPSD/Data" + "/Merge/" + kcStrs[n] + ".h5"
		#fIn = os.path.expanduser('~') + "/Work/StormPSD/Data" + "/MergeWedge/" + kcStrs[n] + ".h5"
		#fIn = os.path.expanduser('~') + "/Work/StormPSD/grab/std/" + kcStrs[n] + ".h5"
		if (n == 0):
			fIn = os.path.expanduser('~') + "/Work/StormPSD/Data" + "/Merge/" + kcStrs[n] + ".h5"
			R,P,K,Tkc,I0 = kc.getCyl(fIn)
		else:
			IStubs = ["_0","_21","_3"]
			fIns = []
			for s in IStubs:
				fIn = os.path.expanduser('~') + "/Work/StormPSD/Data" + "/Merge/" + kcStrs[n] + s + ".h5"
				fIns.append(fIn)
			print(fIns)
			R,P,K,Tkc,I0 = kc.getCyls(fIns)
			
		
		#_,_,_,_,Ntp = kc.getCyl(fIn,fVar="Ntp")

		I0 = kcScls[n]*I0

		#Use different interpolation method
		SimKC = [R,P,K,Tkc,I0]
		rbDat = [Xsc,Ysc,Zsc,Tsc,Krb]

		Is,Isc = kc.InterpSmooth(SimKC,rbDat,Niter=Niter,NiterT=NiterT,doZScl=doZScl,en=En)

		#Save individual contributions
		aI.append(Isc)
		aT.append(Tsc)
		aK.append(Krb)


	#Save combined intensity
	aI.append(aI[1]+aI[2])
	aT.append(aT[1])
	aK.append(aK[1])
	

	#Do figs for each RB
	if (doPanelFig):
		Np = len(aI)
		figSize = (16,16)
		vNorm = LogNorm(vmin=1.0,vmax=1.0e+6)
		figName = "IComp_%s.png"%(rbStr)

		fig = plt.figure(figsize=figSize)
		gs = gridspec.GridSpec(1+Np+1,1,height_ratios=[10,10,10,10,1,1])
		for n in range(Np):
			Ax = fig.add_subplot(gs[n,0])
			Tp = kc.Ts2date(aT[n],T0Str)
			iPlt = Ax.pcolormesh(Tp,aK[n],aI[n].T,norm=vNorm,cmap=cMap)

			#Set axes
			yStr = "%s\nEnergy [keV]"%Labs[n]
			plt.ylabel(yStr,fontsize="large")
			plt.ylim([50,5.0e+3])
			plt.yscale('log')
			if (n==0):
				Ax.xaxis.tick_top()
				Ax.xaxis.set_label_position('top')
				Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
			elif (n<Np-1):
				plt.setp(Ax.get_xticklabels(),visible=False)
			else:
				Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))			
		#Do colorbar
		AxC = fig.add_subplot(gs[-1,0])
		cb = mpl.colorbar.ColorbarBase(AxC,cmap=cMap,norm=vNorm,orientation='horizontal')
		cb.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")

		#Save and close
		plt.savefig(figName,dpi=figQ)		
		plt.close('all')
		lfmv.trimFig(figName)

	if (doLineFig):
		Np = len(KLs)
		Cols = ['r','b','g','m','k']
		#vMins = 1.0*np.ones(Np)
		#vMaxs = (1.0e+5)*np.ones(Np)
		lw1 = 1.5
		lw2 = 0.5

		figSize = (12,12)
		vNorm = LogNorm(vmin=1.0,vmax=1.0e+6)

		Tp = kc.Ts2date(Tsc,T0Str)
		#Create holder for energy interpolant
		Nt = len(Tsc)
		Irbt = np.zeros(Nt)
		Ikct = np.zeros(Nt)
		iPts = np.zeros((Nt,2))
		iPts[:,0] = Tsc

		#Create interpolants
		rbIi = scipy.interpolate.RegularGridInterpolator((Trb,aK[0]),aI[0],method='linear',bounds_error=False)
		kcIi1 = scipy.interpolate.RegularGridInterpolator((Tsc,aK[1]),aI[1],method='linear',bounds_error=False)
		kcIi2 = scipy.interpolate.RegularGridInterpolator((Tsc,aK[2]),aI[2],method='linear',bounds_error=False)
		kcIi = scipy.interpolate.RegularGridInterpolator((Tsc,aK[-1]),aI[-1],method='linear',bounds_error=False)

		figName = "ICompK_%s.png"%(rbStr)

		fig = plt.figure(figsize=figSize)
		gs = gridspec.GridSpec(Np,1)
		for n in range(Np):
			Ax = fig.add_subplot(gs[n,0])
			iPts[:,1] = KLs[n]
			Ik0 = rbIi(iPts)
			Ik1 = kcIi1(iPts)
			Ik2 = kcIi2(iPts)
			Ik3 = kcIi(iPts)
			
			
			Ik1 = TWin(Ik1,Nw=NTWin)
			Ik2 = TWin(Ik2,Nw=NTWin)
			Ik3 = TWin(Ik3,Nw=NTWin)
			
			Ax.semilogy(Tp,Ik0,'k',linewidth=lw1)
			Ax.semilogy(Tp,Ik1,'g',linewidth=lw2)
			Ax.semilogy(Tp,Ik2,'r',linewidth=lw2)
			Ax.semilogy(Tp,Ik3,'m',linewidth=lw1)
			
			plt.ylim([vMins[n],vMaxs[n]])
			#Set axes
			if (n==0):
				Ax.xaxis.tick_top()
				Ax.xaxis.set_label_position('top')
				Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
			elif (n<Np-1):
				plt.setp(Ax.get_xticklabels(),visible=False)
			else:
				Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))			
			#Add labels
			KLab = "%s keV"%(str(KLs[n]))
			Ax.text(0.05,0.75,KLab,transform=Ax.transAxes,fontsize='large')

		plt.suptitle(Labs[0])
		#Save and close
		plt.savefig(figName,dpi=figQ)
		plt.close('all')
		lfmv.trimFig(figName)


if (doSmoothFig):
	#Figure comparing original versus smoothed intensity
	Kc0 = 1000
	#Kc0 = 250.0
	Tc0 = 120000.0
	#ILab = "1MeV Intensity"
	ILab = "%s keV Intensity"%(str(Kc0))
	vNorm = LogNorm(vmin=1.0e-1,vmax=1.0e+4)

	figName = "smPic.png"
	figSize = (8,8)

	pp,rr = np.meshgrid(P,R)
	xx = rr*np.cos(pp)
	yy = rr*np.sin(pp)
	ik0 = np.abs(K-Kc0).argmin()
	it0 = np.abs(Tkc-Tc0).argmin()
	Ixys = [I0[:,:,ik0,it0],Is[:,:,ik0,it0]]
	#Make figure
	fig = plt.figure(figsize=figSize)
	gs = gridspec.GridSpec(3,1,height_ratios=[10,10,1])
	for n in range(2):
		Ax = fig.add_subplot(gs[n,0])
		Ax.pcolormesh(xx,yy,Ixys[n],norm=vNorm,cmap=cMap)
	AxC = fig.add_subplot(gs[-1,0])
	cb = mpl.colorbar.ColorbarBase(AxC,cmap=cMap,norm=vNorm,orientation='horizontal')
	cb.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")
	plt.suptitle(ILab)
	plt.savefig("smoothPic.png",dpi=figQ)
	#plt.show()
	plt.close('all')


