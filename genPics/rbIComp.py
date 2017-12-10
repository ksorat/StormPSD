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
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LogLocator

lfmv.ppInit()
#K-Lines
KLs = [3000,2000,1000,500,250,75]
vMins = [1.0e-2,1.0e-1,1.0e+0,1.0e+1,1.0e+2,1.0e+3]
vMaxs = [1.0e+2,1.0e+3,1.0e+4,1.0e+5,1.0e+6,1.0e+7]

#print(plt.style.available)
#vNormP = LogNorm(vmin=1.0,vmax=1.0e+6)
vNormP = LogNorm(vmin=1.0,vmax=1.0e+6)
cMapP = "gnuplot2"
#cMapP = "nipy_spectral"
#cMapP = "gist_rainbow"
doTitle = False
doPanelFig = False
doLimPanelFig = True
doLineFig = True
Nk2D = 80

#plt.style.use('ggplot')

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
_  ,simaIKs_trp,simbIKs_trp = pS.GetSimRBKt(SimKC_trp,rbDat,KLs)
_  ,simaIKs_inj,simbIKs_inj = pS.GetSimRBKt(SimKC_inj,rbDat,KLs)

#Labs = ["NULL","Initial","Injected","Combined"]
Labs = ["NULL","Combined","Initial","Injected"]

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
	aI.append(Ik_tot[n])

	aT.append(Ts)
	aK.append(K)
	aI.append(Ik_trp[n])

	aT.append(Ts)
	aK.append(K)
	aI.append(Ik_inj[n])


	if (n == 0):
		rbIKs = rbaIKs
		simIKs = simaIKs
		simIKs_trp = simaIKs_trp
		simIKs_inj = simaIKs_inj
	else:
		rbIKs = rbbIKs
		simIKs = simbIKs
		simIKs_trp = simbIKs_trp
		simIKs_inj = simbIKs_inj
	#-------------------
	#Panel figure
	if (doPanelFig):
		#Change ordering to Data/Comb/Init/Inj

		Np = 4
		figSize = (16,16)
		vNorm = vNormP
		cMap = cMapP
		figName = "IComp2D_%s.png"%(rbStrs[n])

		fig = plt.figure(figsize=figSize)
		gs = gridspec.GridSpec(1+Np+1,1,height_ratios=[10,10,10,10,1,1],hspace=0.05)
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
			Ax.grid(color='silver',axis='y',linewidth=1)#,alpha=0.5)

		#Do colorbar
		AxC = fig.add_subplot(gs[-1,0])
		cb = mpl.colorbar.ColorbarBase(AxC,cmap=cMap,norm=vNorm,orientation='horizontal')
		#cb.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize="large")
		cb.set_label("Intensity [cm$^{-2}$ sr$^{-1}$ s$^{-1}$ keV$^{-1}$]",fontsize="large")
		#'\Large{Intensity}\n\small{[cm$^{-2}$ sr$^{-1}$ s$^{-1}$ keV$^{-1}$]}'
		#Save and close
		plt.savefig(figName,dpi=pS.figQ)
		plt.close('all')
		lfmv.trimFig(figName)
	#-------------------
	#Limited Panel figure
	if (doLimPanelFig):
		doCbar = True
		LimLabs = [Labs[0],"Simulation"]
		Np = 2
		figSize = (16,10)
		vNorm = vNormP
		cMap = cMapP
		figName = "I2D_%s.png"%(rbStrs[n])

		fig = plt.figure(figsize=figSize)
		gs = gridspec.GridSpec(1+Np+1,1,height_ratios=[10,10,1,1])
		npan = 0
		for npp in [0,1]:
			#print(aT[np].shape,aK[np].shape,aI[np].shape)
			Ax = fig.add_subplot(gs[npan,0])
			Tp = kc.Ts2date(aT[npp],pS.T0Str)
			iPlt = Ax.pcolormesh(Tp,aK[npp],aI[npp].T,norm=vNorm,cmap=cMap)

			#Set axes
			yStr = "%s\nEnergy [keV]"%LimLabs[npan]
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
			Ax.tick_params(labelright=True,right=True,which='both')
			npan = npan+1
			
		if (doCbar):
			#Do colorbar
			AxC = fig.add_subplot(gs[-1,0])
			cb = mpl.colorbar.ColorbarBase(AxC,cmap=cMap,norm=vNorm,orientation='horizontal')
			cb.set_label("Intensity [cm-2 sr-1 s-1 keV-1]",fontsize="large")

		#Save and close
		plt.savefig(figName,dpi=pS.figQ)
		plt.close('all')
		lfmv.trimFig(figName)		
	#-------------------
	#Line figure
	if (doLineFig):
		Np = len(KLs)
		figSize = (12,12)
		figName = "ILine_%s.png"%(rbStrs[n])
		LWrb = 2
		LWsim = 1.5
		LWsim_c = 0.75
		Leg = ['Data','Model']
		#height_ratios=[10,10,10,10,1,1]
		#HRs = 10*np.ones(Np+2)
		#HRs[0] = 1; HRs[-1] = 1

		fig = plt.figure(figsize=figSize)
		gs = gridspec.GridSpec(Np,1,hspace=0.15)#,height_ratios=HRs)
		Trb = kc.Ts2date(Tsc,pS.T0Str)
		Tsim = kc.Ts2date(Tkl,pS.T0Str)

		for npp in range(Np):
			Ax = fig.add_subplot(gs[npp,0])
			#Ax.semilogy(Trb,rbIKs[npp],'b',linewidth=LWrb)
			Ax.semilogy(Trb,rbIKs[npp],'b',linewidth=LWsim)
			Ax.semilogy(Tsim,simIKs[npp],color='r',linewidth=LWsim)
			#Ax.semilogy(Tsim,simIKs_trp[npp],color='g',linewidth=LWsim_c)
			#Ax.semilogy(Tsim,simIKs_inj[npp],color='m',linewidth=LWsim_c)
			#Ax.semilogy(Tsim,simIKs[npp],color='r',linewidth=LWsim)

			Ax.set_ylim([vMins[npp],vMaxs[npp]])
			pMin = np.log10(vMins[npp])
			pMax = np.log10(vMaxs[npp])

			majTicks = 10**np.arange(pMin,pMax+1)
			minTicks = []
			#Crazy annoying way of setting minor ticks
			for i in range(4):
				for j in range(2,10):
					minTicks.append(j*majTicks[i])
			minTicks = np.array(minTicks)

			Ax.yaxis.set_ticks(majTicks,minor=False)
			Ax.yaxis.set_ticks(minTicks,minor=True)
			

			K0 = KLs[npp]
			if (K0>=1000):
				KLab = "%s MeV"%(str(K0/1000.0))
			else:
				KLab = "%s keV"%(str(K0))
			Ax.text(0.025,0.8,KLab,transform=Ax.transAxes,fontsize='large')

			#Set axes
			if (npp==0):
				Ax.xaxis.tick_top()
				Ax.xaxis.set_label_position('top')
				Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
				#Ax.tick_params(axis='x', which='major', pad=15)
			elif (npp<Np-1):
				plt.setp(Ax.get_xticklabels(),visible=False)
			else:
				Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))	
				#Ax.tick_params(axis='x', which='major', pad=15)
				
				#Add legend to bottom
				Ax.legend(Leg,loc='upper right',ncol=2)
			Ax.tick_params(labelright=True,right=True,which='both')
			Ax.set_ylabel('\Large{Intensity}\n\small{[cm$^{-2}$ sr$^{-1}$ s$^{-1}$ keV$^{-1}$]}')	
			Ax.set_xlim(Tsim[0],Tsim[-1])
		#SupS = "Intensity Comparison (%s)\n"%(Labs[0]) + r"\textcolor{blue}{Data}/\textcolor{red}{Model}"
		#plt.suptitle("Intensity Comparison (%s)\n"%(Labs[0]) + r'\textcolor{blue}{Data}/\textcolor{red}{Model}')
		
		if (doTitle):
			plt.suptitle("Intensity Comparison (%s)"%(Labs[0]),fontsize="x-large")
		#Save and close
		plt.savefig(figName,dpi=pS.figQ)
		plt.close('all')
		lfmv.trimFig(figName)


