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

lfmv.ppInit()

#Values to plot
Ks = [250,500,750,1000]
Nk = len(Ks)
DelT = (50+6*60)*60 #Seconds to get to 3/17
Nt = 6
#Ts = (60*60)*np.linspace(4,44,Nt) + DelT
Ts = (60*60)*np.array([4,8,12,18,25,30]) + DelT

#----------------
#Image options
figSize = (17.5,12)
figQ = 300
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
cLW = 0.5
cAl = 0.5
LabFS = "large"
TitFS = "large"
lwTRK = 1.5
mSize = 8
Ntrk = 3
Nskp = 30

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


vNorm = LogNorm(vmin=1.0,vmax=1.0e+6)

#Prep for figure
fig = plt.figure(figsize=figSize)
Nkp = 2
HRs = np.ones(Nk+Nkp)
HRs[-1] = 0.1
HRs[-2] = 0.1
gs = gridspec.GridSpec(Nk+Nkp,Nt,height_ratios=HRs)#,hspace=0.1,wspace=0.1)

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
		Ax = fig.add_subplot(gs[k,n])
		Ax.pcolormesh(XX,YY,Ik,norm=vNorm,cmap=cMap)
		if (doField):
			#Get dBz data
			xi,yi,dBz = pS.getFld(T0)
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

		if (k == Nk-1):
			#Plot RB positions for lowest energies
			#Plot tracks
			iRB = max(0,nRB-Ntrk*Nskp)
			Ax.plot(Xrb[0][iRB:nRB],Yrb[0][iRB:nRB],color=pS.rbAC,linewidth=lwTRK)
			Ax.plot(Xrb[1][iRB:nRB],Yrb[1][iRB:nRB],color=pS.rbBC,linewidth=lwTRK)
			#Plot RB points
			Ax.plot(Xrb[0][nRB],Yrb[0][nRB],color=pS.rbAC,marker="o",markersize=mSize)
			Ax.plot(Xrb[1][nRB],Yrb[1][nRB],color=pS.rbBC,marker="o",markersize=mSize)

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
cb.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")

plt.savefig("IPans.png",dpi=figQ)
plt.close('all')
lfmv.trimFig("IPans.png")
