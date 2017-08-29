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

inPkl = "IVid2d.pkl"

#K value for main panel
K0 = 1000.0
#K for line plots
Ks = [500,1000,1500]
KLab = ["500 keV","1MeV","1.5MeV"]
dpiQ = 100

#Visual defaults
figSize = (16,7.75)
vNorm = LogNorm(vmin=1.0e-2,vmax=1.0e+3)
cMap = "viridis"

#Field contours
cMapC = "RdGy"
cAl = 0.75
cLW = 1.0
NumC = 13
Vc = np.linspace(-30,30,NumC)
vcNorm = Normalize(vmin=Vc.min(),vmax=Vc.max())

#RB positions
mSize = 12

#Line plots
lwDST = 1.5
lwRB = 1.5
NumT = 500

Nsk = 1
MSkl = 4
vNL = (1.0,5.0e+4)
#-------------------
#Get data
#Get total KCyl
R,P,K,Tkc,I0 = pS.TotCyl()
Tsc,Xrb,Yrb,Zrb = pS.GetRBs()

#Get RB data
tdst,dst = pS.GetDST()
tdstP = kc.Ts2date(tdst,pS.T0Str)
tI = np.linspace(tdst.min(),tdst.max(),NumT)
IkAs,IkBs = pS.GetRBKt(tI,Ks)
tIp = kc.Ts2date(tI,pS.T0Str)

#Get KCyl->RB trajectories
SimKC = [R,P,K,Tkc,I0]
rbDat = [Tsc,Xrb,Yrb,Zrb]
tS,sIkAs,sIkBs = pS.GetSimRBKt(SimKC,rbDat,Ks,Nsk)
tSp = kc.Ts2date(tS,pS.T0Str)

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
HRs = np.ones(NRow)
HRs[-1] = 0.25
HRs[-2] = 0.1

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

	AxRBb = fig.add_subplot(gs[2:4,4:])
	AxRBa = fig.add_subplot(gs[0:2,4:])
	AxDST = fig.add_subplot(gs[4:,4:])
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
	#Plot RB points
	AxM.plot(Xrb[0][nRB],Yrb[0][nRB],color=pS.rbAC,marker="o",markersize=mSize)
	AxM.plot(Xrb[1][nRB],Yrb[1][nRB],color=pS.rbBC,marker="o",markersize=mSize)


	lfmv.addEarth2D(ax=AxM)
	#-----------------------
	#RB A K-lines
	AxRBa.semilogy(tIp,IkAs[0],'g',tIp,IkAs[1],'b',tIp,IkAs[2],'r')
	AxRBa.semilogy(tSp,sIkAs[0],'g:',tSp,sIkAs[1],'b:',tSp,sIkAs[2],'r:')
	AxRBa.set_ylim(vNL)
	AxRBa.yaxis.tick_right()
	AxRBa.yaxis.set_label_position("right")
	plt.setp(AxRBa.get_xticklabels(),visible=False)
	AxRBa.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbAC,linewidth=lwRB)
	AxRBa.legend(KLab,bbox_to_anchor=(0.2,1.15))#,loc='upper left')
	AxRBa.set_xlim(dMin,dMax)
	AxRBa.xaxis.set_major_locator(mdates.HourLocator(interval=6))
	AxRBa.set_ylabel("Intensity")
	#-----------------------
	#RB B K-lines
	AxRBb.semilogy(tIp,IkBs[0],'g',tIp,IkBs[1],'b',tIp,IkBs[2],'r')
	AxRBb.semilogy(tSp,sIkBs[0],'g:',tSp,sIkBs[1],'b:',tSp,sIkBs[2],'r:')
	AxRBb.set_ylim(vNL)
	plt.setp(AxRBb.get_xticklabels(),visible=False)
	AxRBb.yaxis.tick_right()
	AxRBb.yaxis.set_label_position("right")
	AxRBb.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbBC,linewidth=lwRB)
	AxRBb.set_xlim(dMin,dMax)
	AxRBb.xaxis.set_major_locator(mdates.HourLocator(interval=6))
	AxRBb.set_ylabel("Intensity")
	#-----------------------
	#DST plot
	AxDST.plot(tdstP,dst,'k')
	AxDST.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color='k',linewidth=lwDST)
	AxDST.yaxis.tick_right()
	AxDST.yaxis.set_label_position("right")
	AxDST.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ'))
	AxDST.set_xlim(dMin,dMax)
	AxDST.xaxis.set_major_locator(mdates.HourLocator(interval=6))
	AxDST.set_ylabel("DST [nT]")

	fOut = "Data/Vid.%04d.png"%nVid
	plt.savefig(fOut,dpi=dpiQ)
	lfmv.trimFig(fOut)
	plt.close('all')
	nVid = nVid+1
	

# if (os.path.isfile(inPkl)):
# 	print("Reading PKL")
# 	with open(inPkl, "rb") as f:
# 		R,P,K,Tkc,I0 = pickle.load(f)
# 		Tsc,Xrb,Yrb,Zrb = pickle.load(f)
# 		tdst,dst = pickle.load(f)
# else:

	# #Save to pkl
	# print("Creating PKL")
	# with open(inPkl,"wb") as f:
	# 	pickle.dump((R,P,K,Tkc,I0),f)
	# 	pickle.dump((Tsc,Xrb,Yrb,Zrb),f)
	# 	pickle.dump((tdst,dst),f)

