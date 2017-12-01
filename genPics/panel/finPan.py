#Figures for main figure panel
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
from matplotlib.patches import Wedge
from matplotlib.patches import ConnectionPatch

lfmv.ppInit()

#Values to plot
Ks = [250,500,750,1000]
Nk = len(Ks)
DelT = (50+6*60)*60 #Seconds to get to 3/17
Nt = 6

Ts = (60*60)*np.array([4,8,12,18,25,30]) + DelT

#----------------
#Image options
figSize = (14,14)
figQ = 300
vNorm = LogNorm(vmin=1.0,vmax=1.0e+6)
cMap = "inferno"
#cMap = "plasma"
cMapC = "RdGy"
#cMapC = "PRGn"
cMapC = "BrBG"
#Contours
doField = True
Nc = 7
dB0 = 27.5
Vc = np.linspace(-dB0,dB0,Nc)
vcNorm = Normalize(vmin=-dB0,vmax=dB0)
#RB tracks
cLW = 0.75
cAl = 0.5
LegFS = "medium"
LabFS = "large"
TitFS = "large"
lwTRK = 1.5
mSize = 8
Ntrk = 3
Nskp = 30
#Line plots
pLW = 1.0
dstLW = 1.5
cmeLW = 2

SMxy = 12.5

#----------------
#Get cylinder data
#Get total KCyl
R,P,K,Tkc,I0 = pS.TotCyl()
Tsc,Xrb,Yrb,Zrb,Lrb = pS.GetRBs()

#Get grid
XX,YY = kc.xy2rp(R,P)

#Cut out inner values
TINY = 1.0e-8
I0[I0<TINY] = TINY

#Get DST
tDST,dst = pS.GetDST()

#Get injection rates
nDT = 3
tIR,aJt = pS.wIRate(nDT=nDT)
NumW = len(aJt)
datemin = datetime.datetime.strptime("2013-03-17T04:00:00Z",kc.T0Fmt)
#datemax = datetime.datetime.strptime("2013-03-18T06:00:00Z",kc.T0Fmt)
datemax = datetime.datetime.strptime("2013-03-18T12:00:00Z",kc.T0Fmt)
dCME = datetime.datetime.strptime(pS.CME_T0Str,kc.T0Fmt)

#----------------


#Prep for figure
fig = plt.figure(figsize=figSize)
Nkp = 2+2
nOff = 2
HRs = np.ones(Nk+Nkp)
HRs[nOff-1] = 0.15
HRs[-1] = 0.15
HRs[-2] = 0.25

gs = gridspec.GridSpec(Nk+Nkp,Nt,height_ratios=HRs,hspace=0.05,wspace=0.05)

#Plot two time panels at top
#pLab = ["DST"]
AxDST = fig.add_subplot(gs[0,:])
TpDST = kc.Ts2date(tDST,pS.T0Str)
AxDST.plot(TpDST,dst,'k',linewidth=dstLW,label="Dst")
AxDST.set_ylabel('Dst [nT]',fontsize=LabFS)
AxDST.legend(fontsize=LegFS,loc='upper right')
AxDST.set_xlim(datemin,datemax)
AxDST.xaxis.tick_top()
AxDST.set_xlabel('Date [UT]',fontsize="x-large")
AxDST.xaxis.set_label_position('top')


#Injection rate time plots
AxIR = AxDST.twinx()
TpIR = kc.Ts2date(tIR,pS.T0Str) 
for i in range(NumW):
	AxIR.plot(TpIR,aJt[i],color=pS.iCols[i],linewidth=pLW,label=pS.iLabs[i])
	#pLab.append(pS.iLabs[i])
AxIR.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
AxIR.axvline(dCME,color=pS.CME_C,linewidth=cmeLW)
AxIR.annotate('CME',color=pS.CME_C,xy=(dCME,0.03),xytext=(20,-75),textcoords='offset pixels',fontsize=LabFS)
AxIR.set_xlim(datemin,datemax)
AxIR.set_ylim(0,0.05)
AxIR.set_ylabel("Injection Rate",fontsize=LabFS)
AxIR.legend(fontsize=LegFS,loc='lower right')
#AxIR.xaxis.tick_top()
#AxIR.set_xlabel('Date [UT]')
#AxIR.xaxis.set_label_position('top')

#Fake axis for matching
AxCon = AxDST.twinx()
mdMin = mdates.date2num(datemin)
mdMax = mdates.date2num(datemax)
AxCon.set_xlim(mdMin,mdMax)
plt.setp(AxCon.get_yticklabels(),visible=False)
AxCon.yaxis.set_tick_params(size=0)
#Convert to date times
Tp = kc.Ts2date(Ts,pS.T0Str)

#Loop over time/energies, find bounds and lin interpolate
for n in range(Nt):

	T0 = Ts[n]
	print("T = %f"%(T0))
	it1 = (Tkc>=T0).argmax()
	it0 = it1-1
	dt = (T0-Tkc[it0])/(Tkc[it1]-Tkc[it0])

	#Get single KCyl
	IKc = I0[:,:,:,it0] + dt*(I0[:,:,:,it1]-I0[:,:,:,it0])
	#Get RB matching point
	nRB = np.abs(Tsc-T0).argmin()

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


		#----------------------
		#Create plots
		Ax = fig.add_subplot(gs[k+nOff,n])
		Ax.pcolormesh(XX,YY,Ik,norm=vNorm,cmap=cMap)
		if (doField):
			#Get dBz data
			xi,yi,dBz = pS.getFld(T0)
			Ax.contour(xi,yi,dBz,Vc,cmap=cMapC,alpha=cAl,linewidths=cLW)
		plt.axis('scaled')
		plt.xlim([-SMxy,SMxy])
		plt.ylim([-SMxy,SMxy])
		lfmv.addEarth2D()
		if (n == 0):
			plt.ylabel(KLab,fontsize=LabFS)
		elif (n == Nt-1):
			plt.ylabel("SM-Y [Re]",fontsize=LabFS)
			Ax.yaxis.tick_right()
			Ax.yaxis.set_label_position("right")
		else:
			plt.setp(Ax.get_yticklabels(),visible=False)
	
		if (k < Nk-1):
			plt.setp(Ax.get_xticklabels(),visible=False)
		else:
			plt.xlabel('SM-X [Re]',fontsize=LabFS)
		if (k == 0):
			#Add connection for time matching
			nT = mdates.date2num(Tp[n])
			xyCon = (nT,0.0)
			xyAx = (0.0,SMxy)

			ConT = ConnectionPatch(xyA=xyAx,xyB=xyCon,coordsA="data",coordsB="data",axesA=Ax,axesB=AxCon,color='k')
			Ax.add_artist(ConT)

		if (k == Nk-1):
			#Plot RB positions for lowest energies
			#Plot tracks
			iRB = max(0,nRB-Ntrk*Nskp)
			Ax.plot(Xrb[0][iRB:nRB],Yrb[0][iRB:nRB],color=pS.rbAC,linewidth=lwTRK)
			Ax.plot(Xrb[1][iRB:nRB],Yrb[1][iRB:nRB],color=pS.rbBC,linewidth=lwTRK)
			#Plot RB points
			Ax.plot(Xrb[0][nRB],Yrb[0][nRB],color=pS.rbAC,marker="o",markersize=mSize,markeredgecolor='k')
			Ax.plot(Xrb[1][nRB],Yrb[1][nRB],color=pS.rbBC,marker="o",markersize=mSize,markeredgecolor='k')

		if ( (k==Nk-1) and (n==0) ):
			#Add wedge markers
			wLW = 1.5
			w1 = Wedge((0,0), 12,170,190, width=3.0,fill=False,ec='b',linewidth=wLW)
			w2 = Wedge((0,0), 12,125,145, width=3.0,fill=False,ec='lime',linewidth=wLW)
			w3 = Wedge((0,0), 12,215,235, width=3.0,fill=False,ec='r',linewidth=wLW)
			for w in [w1,w2,w3]:
				Ax.add_artist(w)

if (doField):
	AxCC = fig.add_subplot(gs[-1,Nt/2:])
	cb = mpl.colorbar.ColorbarBase(AxCC,cmap=cMapC,norm=vcNorm,orientation='horizontal')
	cb.set_label("Residual Vertical Field [nT]",fontsize="large")
	AxC = fig.add_subplot(gs[-1,0:Nt/2])

else:
	AxC = fig.add_subplot(gs[-1,Nt/2-Nt/4:Nt/2+Nt/4+1])
cb = mpl.colorbar.ColorbarBase(AxC,cmap=cMap,norm=vNorm,orientation='horizontal')
#cb.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")
cb.set_label("Intensity [cm$^{-2}$ sr$^{-1}$ s$^{-1}$ keV$^{-1}$]",fontsize="large")


plt.savefig("IPans.png",dpi=figQ)
plt.close('all')
lfmv.trimFig("IPans.png")
