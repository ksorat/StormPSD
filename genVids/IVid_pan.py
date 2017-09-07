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
import cPickle as pickle
from matplotlib.patches import Wedge

vNormP = LogNorm(vmin=1.0,vmax=1.0e+6)
cMapP = "jet"
KBds = [75,4.0e+3]
def pI2D(Ax,T,K,I,Lab="Stupid",doX=False):
	Tp = kc.Ts2date(T,pS.T0Str)
	iPlt = Ax.pcolormesh(Tp,K,I.T,norm=vNormP,cmap=cMapP)
	Ax.set_ylim(KBds)
	Ax.set_yscale('log')
	Ax.yaxis.tick_right()
	Ax.yaxis.set_label_position("right")

	if (not doX):
		plt.setp(Ax.get_xticklabels(),visible=False)
	else:
		Ax.xaxis.tick_top()
		Ax.xaxis.set_label_position("top")
	Ax.set_ylabel("Energy [keV]",fontsize="x-small")
	Ax.set_xlim(dMin,dMax)
	Ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
	Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
	Ax.text(-0.04,0.5,Lab,rotation='vertical',transform=Ax.transAxes,fontsize=12,ha='center')

#K value for main panel
K0 = 1000.0
#K for line plots
Ks = [500,1000,1500]
KLab = ["500 keV","1MeV","1.5MeV"]

dpiQ = 300

#Visual defaults
figSize = (16,7.75)
vNorm = LogNorm(vmin=1.0e-2,vmax=1.0e+3)
cMap = "viridis"

#Field contours
cMapC = "RdGy"
cAl = 0.75
cLW = 1.0
NumC = 13
Vc = np.linspace(-35,35,NumC)
vcNorm = Normalize(vmin=Vc.min(),vmax=Vc.max())

#RB positions/tracks
mSize = 12
alTrk = 0.5
Ntrk = 3
Nskp = 30

#Lines
lwDST = 1.5
lwRB = 1.5
NumT = 500
Nsk = 1


#-------------------
#Get data
#Get total KCyl
R,P,K,Tkc,I0 = pS.TotCyl()
Tsc,Xrb,Yrb,Zrb,Lrb = pS.GetRBs()

#Get RB data
tdst,dst = pS.GetDST()
tdstP = kc.Ts2date(tdst,pS.T0Str)
tI = np.linspace(tdst.min(),tdst.max(),NumT)
IkAs,IkBs = pS.GetRBKt(tI,Ks)
tIp = kc.Ts2date(tI,pS.T0Str)

#Get KCyl->RB trajectories
#2D-I plots
Trbs,Krbs,Irbs = pS.GetRB_I2D()

SimKC = [R,P,K,Tkc,I0]
rbDat = [Tsc,Xrb,Yrb,Zrb,Lrb]
Nk2D = 80
Ksim = np.logspace(1,np.log10(5000),Nk2D)
Tsims,Iksims = pS.GetSim_I2D(SimKC,rbDat,Ksim)


#Prep counters
Nkc = len(Tkc)
k0 = np.abs(K-1000).argmin() #Lazy cut for energy

#Prep for figure
lfmv.ppInit()
plt.close('all')
NRow = 7
NCol = 8
Npx = 4
Npy = 4
HRs = 1.25*np.ones(NRow)
HRs[-1] = 0.2
HRs[-2] = 0.125


#Get grid
XX,YY = kc.xy2rp(R,P)

#Set XLims
dMin = datetime.datetime.strptime(pS.T0Str,kc.T0Fmt) + datetime.timedelta(seconds=pS.tMin)
dMax = datetime.datetime.strptime(pS.T0Str,kc.T0Fmt) + datetime.timedelta(seconds=pS.tMax)

#for n in range(0,Nkc):
nMin = 10
nMax = 1000
nVid = 0

TINY = 1.0e-8
I0[I0<TINY] = TINY

print("Starting video ...")
for n in range(nMin,nMax):
#for n in range(500,501):
	#------------------
	#Setup figure
	fig = plt.figure(figsize=figSize)
	gs = gridspec.GridSpec(NRow,NCol,height_ratios=HRs,left=0.1,right=0.9,top=0.9,bottom=0.1)

	AxCI = fig.add_subplot(gs[-1,0:2])
	AxCF = fig.add_subplot(gs[-1,2:4])

	#Add colorbars for main panel
	cbI = mpl.colorbar.ColorbarBase(AxCI,cmap=cMap,norm=vNorm,orientation='horizontal')
	cbI.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")
	cbF = mpl.colorbar.ColorbarBase(AxCF,cmap=cMapC,norm=vcNorm,orientation='horizontal')
	cbF.set_label("Residual Vertical Field [nT]",fontsize="large")

	AxM  = fig.add_subplot(gs[0:Npx+1,0:Npy])
	AxM.set_aspect('equal')
	AxM.set_xlim(-15.5,12.5)
	AxM.set_ylim(-15,15)

	#Add wedge markers
	wLW = 1.5
	w1 = Wedge((0,0), 12,170,190, width=3.0,fill=False,ec='b',linewidth=wLW)
	w2 = Wedge((0,0), 12,125,145, width=3.0,fill=False,ec='lime',linewidth=wLW)
	w3 = Wedge((0,0), 12,215,235, width=3.0,fill=False,ec='r',linewidth=wLW)
	for w in [w1,w2,w3]:
		AxM.add_artist(w)

	#Add RB labels
	rbFS = 14
	AxM.text(0.75,1.01,'RBSP-B',color=pS.rbBC,transform=AxM.transAxes,fontsize=rbFS)
	AxM.text(0.10,1.01,'RBSP-A',color=pS.rbAC,transform=AxM.transAxes,fontsize=rbFS)


	#Add pcolor panels
	AxRBaa = fig.add_subplot(gs[0,4:])
	AxRBa  = fig.add_subplot(gs[1,4:])
	AxRBba = fig.add_subplot(gs[2,4:])
	AxRBb  = fig.add_subplot(gs[3,4:])

	#Add pcolor colorbar
	AxCI2D = fig.add_subplot(gs[-1,4:])
	cbI = mpl.colorbar.ColorbarBase(AxCI2D,cmap=cMapP,norm=vNormP,orientation='horizontal')
	cbI.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")

	AxDST = fig.add_subplot(gs[4,4:])
	AxNull = fig.add_subplot(gs[6,5])
	AxNull.set_visible(False)

	#------------------
	#Start plotting things
	DateN = datetime.datetime.strptime(pS.T0Str,kc.T0Fmt) + datetime.timedelta(seconds=Tkc[n])
	#print("Writing image at t=%f"%Tkc[n])
	TitS = "1MeV Intensity\n%s"%(str(DateN))
	nRB = np.abs(Tsc-Tkc[n]).argmin()

	#------------------
	#Main panel
	AxM.pcolormesh(XX,YY,I0[:,:,k0,n],norm=vNorm,cmap=cMap)
	AxM.set_title(TitS)
	#Field contours
	Xc,Yc,dBz = pS.getFld(Tkc[n])
	AxM.contour(Xc,Yc,dBz,Vc,cmap=cMapC,alpha=cAl,linewidth=cLW)

	#Plot tracks (do before RB current)
	for i in range(Ntrk):
		iRB = nRB-(i+1)*Nskp
		if (iRB<0):
			continue
		else:
			msTrk = mSize/(i+2.0)
			AxM.plot(Xrb[0][iRB],Yrb[0][iRB],color=pS.rbAC,marker="o",markersize=msTrk,alpha=alTrk)
			AxM.plot(Xrb[1][iRB],Yrb[1][iRB],color=pS.rbBC,marker="o",markersize=msTrk,alpha=alTrk)
	#Plot RB points
	AxM.plot(Xrb[0][nRB],Yrb[0][nRB],color=pS.rbAC,marker="o",markersize=mSize)
	AxM.plot(Xrb[1][nRB],Yrb[1][nRB],color=pS.rbBC,marker="o",markersize=mSize)


	lfmv.addEarth2D(ax=AxM)
	#-----------------------
	#I2D Plots
	pI2D(AxRBaa,Trbs[0],Krbs[0],Irbs[0],Lab='RBSP-A\nActual',doX=True)
	pI2D(AxRBa ,Tsims[0],Ksim,Iksims[0],Lab='RBSP-A\nModel' )
	AxRBaa.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbAC,linewidth=lwRB)
	AxRBa.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbAC,linewidth=lwRB)

	pI2D(AxRBba,Trbs[1],Krbs[1],Irbs[1],Lab='RBSP-B\nActual')
	pI2D(AxRBb ,Tsims[1],Ksim,Iksims[1],Lab='RBSP-B\nModel' )
	AxRBba.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbBC,linewidth=lwRB)
	AxRBb.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbBC,linewidth=lwRB)


	#-----------------------
	#DST plot
	AxDST.plot(tdstP,dst,'k')
	AxDST.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color='k',linewidth=lwDST)
	AxDST.yaxis.tick_right()
	AxDST.yaxis.set_label_position("right")
	AxDST.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
	AxDST.set_xlim(dMin,dMax)
	AxDST.xaxis.set_major_locator(mdates.HourLocator(interval=6))
	AxDST.set_ylabel("DST [nT]")

	fOut = "Data/Vid.%04d.png"%nVid
	plt.savefig(fOut,dpi=dpiQ)
	lfmv.trimFig(fOut)
	plt.close('all')
	nVid = nVid+1
	

