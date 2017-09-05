#Main comparison versus observational data
#Do 2D I(K,t) and line plots

import kCyl as kc
import pyStorm as pS
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import lfmViz as lfmv

#K-Lines
KLs = [3000,2000,1000,500,250,50]
vMins = [1.0e-2,1.0e-1,1.0e+0,1.0e+1,1.0e+2,1.0e+3]
vMaxs = [1.0e+2,1.0e+3,1.0e+4,1.0e+5,1.0e+6,1.0e+7]

doPanelFig = True
doLineFig = True

#-------------------
#Get data
#Get total KCyl
R,P,K,Tkc,I0 = pS.TotCyl()
#Get RB data
Tsc,Xrb,Yrb,Zrb = pS.GetRBs()
Trbs,Krbs,Irbs = pS.GetRB_I2D()

Labs = ["NULL","Trapped","Injected","Combined"]
rbStrs = ["A","B"]
Nrb = 2


for n in range(Nrb):
	Labs[0] = "RBSP-" + rbStrs[n].upper()
	aT = []; aK = []; aI = []
	#Start w/ RB data
	aT.append(Trbs[n])
	aK.append(Krbs[n])
	aI.append(Irbs[n])


	#-------------------
	#Panel figure
	if (doPanelFig):
		Np = 4
		figSize = (16,16)
		vNorm = LogNorm(vmin=1.0,vmax=1.0e+6)
		cMap = "jet"
		figName = "IComp2D_%s.png"%(rbStrs[n])

		fig = plt.figure(figsize=figSize)
		gs = gridspec.GridSpec(1+Np+1,1,height_ratios=[10,10,10,10,1,1])
		for np in range(1):
			print(aT[np].shape,aK[np].shape,aI[np].shape)
			Ax = fig.add_subplot(gs[np,0])
			Tp = kc.Ts2date(aT[np],pS.T0Str)
			iPlt = Ax.pcolormesh(Tp,aK[np],aI[np].T,norm=vNorm,cmap=cMap)

			#Set axes
			yStr = "%s\nEnergy [keV]"%Labs[np]
			plt.ylabel(yStr,fontsize="large")
			plt.ylim([50,5.0e+3])
			plt.yscale('log')
			if (np==0):
				Ax.xaxis.tick_top()
				Ax.xaxis.set_label_position('top')
				Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
			elif (np<Np-1):
				plt.setp(Ax.get_xticklabels(),visible=False)
			else:
				Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))			

		#Do colorbar
		AxC = fig.add_subplot(gs[-1,0])
		cb = mpl.colorbar.ColorbarBase(AxC,cmap=cMap,norm=vNorm,orientation='horizontal')
		cb.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")

		#Save and close
		plt.savefig(figName,dpi=pS.figQ)
		plt.close('all')
		lfmv.trimFig(figName)