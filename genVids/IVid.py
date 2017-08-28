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

doLoad = True

#K value for main panel
K0 = 1000.0

#Visual defaults
figSize = (14,8)
vNorm = LogNorm(vmin=1.0e-2,vmax=1.0e+3)
cMap = "viridis"

#Field contours
cMapC = "RdGy"
cAl = 1.0
cLW = 1.0
NumC = 9
Vc = np.linspace(-35,35,NumC)
vcNorm = Normalize(vmin=Vc.min(),vmax=Vc.max())

#RB positions
mSize = 12

if (doLoad):
	#Get total KCyl
	R,P,K,Tkc,I0 = pS.TotCyl()
	Tsc,Xrb,Yrb,Zrb = pS.GetRBs()
	tdst,dst = pS.GetDST()
#Prep counters
Nkc = len(Tkc)
k0 = np.abs(K-1000).argmin() #Lazy cut for energy

#Prep for figure
plt.close('all')
lfmv.ppInit()
NRow = 7
NCol = 6
Npx = 4
Npy = 4
HRs = np.ones(NRow)
HRs[-1] = 0.25
HRs[-2] = 0.2

#Get grid
XX,YY = kc.xy2rp(R,P)

fig = plt.figure(figsize=figSize)
gs = gridspec.GridSpec(NRow,NCol,height_ratios=HRs)

AxCI = fig.add_subplot(gs[-1,0:2])
AxCF = fig.add_subplot(gs[-1,2:4])

#Add colorbars for main panel
cbI = mpl.colorbar.ColorbarBase(AxCI,cmap=cMap,norm=vNorm,orientation='horizontal')
cbI.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")
cbF = mpl.colorbar.ColorbarBase(AxCF,cmap=cMapC,norm=vcNorm,orientation='horizontal')
cbF.set_label("Residual Vertical Field [nT]",fontsize="large")

AxM  = fig.add_subplot(gs[:-2,:-2])
AxM.set_aspect('equal')
AxM.set_xlim(-12.5,12.5)
AxM.set_ylim(-12.5,12.5)

AxRBb = fig.add_subplot(gs[2:4,4:])
AxRBa = fig.add_subplot(gs[0:2,4:])
AxDST = fig.add_subplot(gs[4:6,4:])
AxNull = fig.add_subplot(gs[6,5])
AxNull.set_visible(False)

#for n in range(0,Nkc):

for n in range(500,501):

	#print("Writing image at t=%f"%Tkc[n])
	nRB = np.abs(Tsc-Tkc[n]).argmin()

	#------------------
	#Main panel
	AxM.pcolormesh(XX,YY,I0[:,:,k0,n],norm=vNorm,cmap=cMap)
	#Field contours
	Xc,Yc,dBz = pS.getFld(Tkc[n])
	AxM.contour(Xc,Yc,dBz,Vc,cmap=cMapC,alpha=cAl,linewidth=cLW)
	#Plot RB points
	AxM.plot(Xrb[0][nRB],Yrb[0][nRB],color=pS.rbAC,marker="o",markersize=mSize)
	AxM.plot(Xrb[1][nRB],Yrb[1][nRB],color=pS.rbBC,marker="o",markersize=mSize)


	lfmv.addEarth2D(ax=AxM)

	AxDST.plot(tdst,dst)
	#plt.close(1)
	plt.show()
