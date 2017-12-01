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

doTest = True
doDark = False

NvSkp = 1 #Video frame skip
nMin = 10
nMax = 1000

vNormP = LogNorm(vmin=1.0,vmax=1.0e+6)
cMapP = "gnuplot2"
cMapP = "gnuplot"

KBds = [75,4.0e+3]
def pI2D(Ax,T,K,I,Lab="Stupid",doX=False,tC='r'):
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
	Ax.text(-0.015,0.6,Lab,color=tC,rotation='vertical',transform=Ax.transAxes,fontsize=16,ha='center')

#K value for main panel
K0 = 1000.0
Rcut = 2.35
dpiQ = 300

#Visual defaults
figSize = (16.0,7.8)
vNorm = LogNorm(vmin=1.0e-2,vmax=1.0e+3)
#vNorm = LogNorm(vmin=1.0e-2,vmax=1.0e+4)
cMap = "viridis"
#cMap = "inferno"

#Field contours
cMapC = "RdGy"
#cMapC = "BrBG"
cAl = 0.75
cLW = 1.0
NumC = 13
dB0 = 35
Vc = np.linspace(-dB0,dB0,NumC)
vcNorm = Normalize(vmin=-dB0,vmax=dB0)

#RB positions/tracks
mSize = 12
alTrk = 0.5
Ntrk = 3
Nskp = 30
lwTRK = 1.5

#Lines
lwDST = 1.5
lwRB = 1.5
NumT = 500
Nsk = 1
cmeLW = 2

#-------------------
#Get data
#Get total KCyl
R,P,K,Tkc,I0 = pS.TotCyl()
Tsc,Xrb,Yrb,Zrb,Lrb = pS.GetRBs()

#Get RB data
tdst,dst = pS.GetDST()
tdstP = kc.Ts2date(tdst,pS.T0Str)
tI = np.linspace(tdst.min(),tdst.max(),NumT)

#Get KCyl->RB trajectories
#2D-I plots
Trbs,Krbs,Irbs = pS.GetRB_I2D()

SimKC = [R,P,K,Tkc,I0]
rbDat = [Tsc,Xrb,Yrb,Zrb,Lrb]
Nk2D = 80
Ksim = np.logspace(1,np.log10(5000),Nk2D)
Tsims,Iksims = pS.GetSim_I2D(SimKC,rbDat,Ksim)

#Get injection rates
nDT = 3
tIR,aJt = pS.wIRate(nDT=nDT)
NumW = len(aJt)
datemin = datetime.datetime.strptime("2013-03-17T04:00:00Z",kc.T0Fmt)
#datemax = datetime.datetime.strptime("2013-03-18T06:00:00Z",kc.T0Fmt)
datemax = datetime.datetime.strptime("2013-03-18T12:00:00Z",kc.T0Fmt)
dCME = datetime.datetime.strptime(pS.CME_T0Str,kc.T0Fmt)


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
rI = (R>Rcut).argmax()

#Set XLims
dMin = datetime.datetime.strptime(pS.T0Str,kc.T0Fmt) + datetime.timedelta(seconds=pS.tMin)
dMax = datetime.datetime.strptime(pS.T0Str,kc.T0Fmt) + datetime.timedelta(seconds=pS.tMax)

nVid = 0

TINY = 1.0e-8
I0[I0<TINY] = TINY
#-------------------
#Flip for dark version
if (doDark):
	plt.style.use('dark_background')
	dstC = 'w'
else:
	dstC = 'k'

print("Starting video ...")
if (doTest):
	nMin = 500
	nMax = 501
	NvSkp = 1
for n in range(nMin,nMax,NvSkp):
	#------------------
	#Setup figure
	fig = plt.figure(figsize=figSize)
	gs = gridspec.GridSpec(NRow,NCol,height_ratios=HRs,left=0.1,right=0.9,top=0.9,bottom=0.1)

	AxCI = fig.add_subplot(gs[-1,0:2])
	AxCF = fig.add_subplot(gs[-1,2:4])

	#Add colorbars for main panel
	cbI = mpl.colorbar.ColorbarBase(AxCI,cmap=cMap,norm=vNorm,orientation='horizontal')
	cbI.set_label("Intensity [cm-2 sr-1 s-1 keV-1]",fontsize="large")
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

	#Tweak pcolors to touch
	x0 = 0.50851064
	y0A = 0.770166
	y0A = 0.785
	#y0B = 0.32212742
	y0B = 0.31
	y0DST = 0.17278123
	y0DST = 0.1715
	dX = 0.4
	dY = 0.15
	dYDST = 0.125

	eps = 0.005/2
	AxRBaa.set_position([x0,y0A,dX,dY])
	AxRBa .set_position([x0,y0A-dY-eps,dX,dY])
	AxRBba.set_position([x0,y0B+dY+eps,dX,dY])
	AxRBb .set_position([x0,y0B,dX,dY])

	#Add pcolor colorbar
	AxCI2D = fig.add_subplot(gs[-1,4:])
	cbI = mpl.colorbar.ColorbarBase(AxCI2D,cmap=cMapP,norm=vNormP,orientation='horizontal')
	cbI.set_label("Intensity [cm-2 sr-1 s-1 keV-1]",fontsize="large")

	AxDST = fig.add_subplot(gs[4,4:])
	AxDST.set_position([x0,y0DST,dX,dYDST])
	AxIR = AxDST.twinx()
	AxIR.set_position([x0,y0DST,dX,dYDST])

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
	Ixy = I0[:,:,k0,n]
	Ixy[0:rI,:] = 0.0
	AxM.pcolormesh(XX,YY,Ixy,norm=vNorm,cmap=cMap)
	AxM.set_title(TitS)
	#Field contours
	Xc,Yc,dBz = pS.getFld(Tkc[n])
	AxM.contour(Xc,Yc,dBz,Vc,cmap=cMapC,alpha=cAl,linewidths=cLW)

	#Plot tracks (do before RB current)
	iRB = max(0,nRB-Ntrk*Nskp)
	AxM.plot(Xrb[0][iRB:nRB],Yrb[0][iRB:nRB],color=pS.rbAC,linewidth=lwTRK)
	AxM.plot(Xrb[1][iRB:nRB],Yrb[1][iRB:nRB],color=pS.rbBC,linewidth=lwTRK)

	#Plot RB points
	AxM.plot(Xrb[0][nRB],Yrb[0][nRB],color=pS.rbAC,marker="o",markersize=mSize,markeredgecolor='k')
	AxM.plot(Xrb[1][nRB],Yrb[1][nRB],color=pS.rbBC,marker="o",markersize=mSize,markeredgecolor='k')


	lfmv.addEarth2D(ax=AxM)
	#-----------------------
	#I2D Plots
	pI2D(AxRBaa,Trbs[0],Krbs[0],Irbs[0],tC=pS.rbAC,Lab='Data',doX=True)
	pI2D(AxRBa ,Tsims[0],Ksim,Iksims[0],tC=pS.rbAC,Lab='Model' )
	AxRBaa.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbAC,linewidth=lwRB)
	AxRBa.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbAC,linewidth=lwRB)

	pI2D(AxRBba,Trbs[1],Krbs[1],Irbs[1],tC=pS.rbBC,Lab='Data')
	pI2D(AxRBb ,Tsims[1],Ksim,Iksims[1],tC=pS.rbBC,Lab='Model' )
	AxRBba.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbBC,linewidth=lwRB)
	AxRBb.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbBC,linewidth=lwRB)

	#Add double labels
	rb1DX = -0.055
	rb1DY = -0.05
	rbFS = 18
	#rbAbox = dict(boxstyle="darrow", ec=pS.rbAC, lw=2)
	AxRBaa.text(rb1DX,rb1DY,"RBSP-A",color=pS.rbAC,rotation='vertical',transform=AxRBaa.transAxes,fontsize=rbFS,va="center",ha='center')
	AxRBba.text(rb1DX,rb1DY,"RBSP-B",color=pS.rbBC,rotation='vertical',transform=AxRBba.transAxes,fontsize=rbFS,va="center",ha='center')


	#-----------------------

	#DST plot
	dstCol = 'darkorange'
	#AxDST.plot(tdstP,dst,dstC)
	AxDST.plot(tdstP,dst,dstCol)
	AxDST.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=dstC,linewidth=lwDST)
	AxDST.yaxis.tick_right()
	AxDST.tick_params('y',colors=dstCol)
	AxDST.yaxis.set_label_position("right")
	AxDST.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
	AxDST.set_xlim(dMin,dMax)
	AxDST.xaxis.set_major_locator(mdates.HourLocator(interval=6))
	AxDST.set_ylabel("Dst [nT]",color=dstCol)

	
	TpIR = kc.Ts2date(tIR,pS.T0Str) 
	for i in range(NumW):
		AxIR.plot(TpIR,aJt[i],color=pS.iCols[i],linewidth=0.5)
	AxIR.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
	AxIR.set_xlim(dMin,dMax)
	AxIR.xaxis.set_major_locator(mdates.HourLocator(interval=6))
	AxIR.set_xlim(dMin,dMax)
	AxIR.set_ylim(0,0.05)
	AxIR.yaxis.tick_left()
	AxIR.yaxis.set_label_position("left")

	AxIR.set_ylabel("Injection Rate")


	#----------------------------
	#Finalize

	fOut = "Data/Vid.%04d.png"%nVid
	plt.savefig(fOut,dpi=dpiQ)
	lfmv.trimFig(fOut,bLenX=20,bLenY=20)
	#lfmv.shaveFig(fOut,bLenX=1,bLenY=0)
	plt.close('all')
	nVid = nVid+1
	

