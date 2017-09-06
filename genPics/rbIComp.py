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
KLs = [3000,2000,1000,500,250,75]
vMins = [1.0e-2,1.0e-1,1.0e+0,1.0e+1,1.0e+2,1.0e+3]
vMaxs = [1.0e+2,1.0e+3,1.0e+4,1.0e+5,1.0e+6,1.0e+7]

doPanelFig = True
doLineFig = True
Nk2D = 80

#-------------------
#Get data
#Get 3 KCyls
R,P,K,Tkc,Itot = pS.TotCyl(doSmooth=True)
_,_,_,_,Iinj   = pS.InjCyl(doSmooth=True)
_,_,_,_,Itrp   = pS.TrapCyl(doSmooth=True)

#Get RB data
Tsc,Xrb,Yrb,Zrb,Lrb = pS.GetRBs()
Trbs,Krbs,Irbs = pS.GetRB_I2D()

#Get KCyl->RB trajectories
SimKC_tot = [R,P,K,Tkc,Itot]
SimKC_inj = [R,P,K,Tkc,Iinj]
SimKC_trp = [R,P,K,Tkc,Itrp]
rbDat = [Tsc,Xrb,Yrb,Zrb,Lrb]

K = np.logspace(1,np.log10(5000),Nk2D)
Ts_tot,Ik_tot = pS.GetSim_I2D(SimKC_tot,rbDat,K)
Ts_inj,Ik_inj = pS.GetSim_I2D(SimKC_inj,rbDat,K)
Ts_trp,Ik_trp = pS.GetSim_I2D(SimKC_trp,rbDat,K)
Ts = Ts_tot[0]

#Get line data for RB
rbaIKs,rbbIKs = pS.GetRBKt(Tsc,KLs)
#Get line data from SIM
Tkl,simaIKs,simbIKs = pS.GetSimRBKt(SimKC_tot,rbDat,KLs)

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

	#Add trapped/injected/total
	aT.append(Ts)
	aK.append(K)
	aI.append(Ik_trp[n])

	aT.append(Ts)
	aK.append(K)
	aI.append(Ik_inj[n])

	aT.append(Ts)
	aK.append(K)
	aI.append(Ik_tot[n])

	if (n == 0):
		rbIKs = rbaIKs
		simIKs = simaIKs
	else:
		rbIKs = rbbIKs
		simIKs = simbIKs

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
		for npp in range(Np):
			#print(aT[np].shape,aK[np].shape,aI[np].shape)
			Ax = fig.add_subplot(gs[npp,0])
			Tp = kc.Ts2date(aT[npp],pS.T0Str)
			iPlt = Ax.pcolormesh(Tp,aK[npp],aI[npp].T,norm=vNorm,cmap=cMap)

			#Set axes
			yStr = "%s\nEnergy [keV]"%Labs[npp]
			plt.ylabel(yStr,fontsize="large")
			plt.ylim([75,5.0e+3])
			plt.yscale('log')
			if (npp==0):
				Ax.xaxis.tick_top()
				Ax.xaxis.set_label_position('top')
				Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
			elif (npp<Np-1):
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
	#-------------------
	#Line figure
	if (doLineFig):
		Np = len(KLs)
		figSize = (12,12)
		figName = "ICompKs_%s.png"%(rbStrs[n])
		LWrb = 2.0
		LWsim = 1.0
		#height_ratios=[10,10,10,10,1,1]
		#HRs = 10*np.ones(Np+2)
		#HRs[0] = 1; HRs[-1] = 1

		fig = plt.figure(figsize=figSize)
		gs = gridspec.GridSpec(Np,1,hspace=0.2)#,height_ratios=HRs)
		Trb = kc.Ts2date(Tsc,pS.T0Str)
		Tsim = kc.Ts2date(Tkl,pS.T0Str)

		for npp in range(Np):
			Ax = fig.add_subplot(gs[npp,0])
			Ax.semilogy(Trb,rbIKs[npp],'k',linewidth=LWrb)
			Ax.semilogy(Tsim,simIKs[npp],color='r',linewidth=LWsim)
			Ax.set_ylim([vMins[npp],vMaxs[npp]])

			KLab = "%s keV"%(str(KLs[npp]))
			Ax.text(0.025,0.8,KLab,transform=Ax.transAxes,fontsize='large')

			#Set axes
			if (npp==0):
				Ax.xaxis.tick_top()
				Ax.xaxis.set_label_position('top')
				Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
				Ax.tick_params(axis='x', which='major', pad=15)
			elif (npp<Np-1):
				plt.setp(Ax.get_xticklabels(),visible=False)
			else:
				Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))	
				Ax.tick_params(axis='x', which='major', pad=15)		

		plt.suptitle(Labs[0])

		#Save and close
		plt.savefig(figName,dpi=pS.figQ)
		plt.close('all')
		lfmv.trimFig(figName)


