#Generate multi-multipanel figure for kappa/kTScl sweep
import kCyl as kc
import os
import numpy as np
import scipy
import scipy.interpolate
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

#Defaults
#Time data
T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"
Base = "/Work/StormPSD/Data" + "/f0Sweep/"

tMin = 33600.0
tMax = 189000.0
KLs = [2000,1000,500,250]
vMins = [1.0e-1,1.0e+0,1.0e+2,1.0e+3]
vMaxs = [1.0e+3,1.0e+4,1.0e+6,1.0e+7]

figSize = (36,18)
figQ = 300 #DPI
cMap = "jet"

kappas = [3.0,3.5,4.0]
kappaStr = ["A","B","C"]

kTScls = [0.15,0.25,0.4]
kTSclStr = ["a","b","c"]

kcStrs = ["KCyl_StormT","KCyl_StormI"]
kcScls = np.pi*4*np.array([1.0,10.0])

rbStrs = ["A","B"]
rbSK = 2 #Skip cadence for RB intensity
rbSKt = 1 #Skip cadence for RB trajectory
doDip = True #Do dipole projection of trajectory
lw = 2.0

#------------------
Nkap = len(kappas)
Nkt = len(kTScls)
NumPop = len(kcStrs)

#NumRB = len(rbStrs)
NumRB = 1 #Only use A
NKL = len(KLs)

for nrb in range(NumRB):
	#Create relevant files
	rbStr = rbStrs[nrb]

	OrbF = "vaporbRB" + rbStr.upper() + ".txt"
	rbFs = "rbsp" + rbStrs[nrb].lower()
	rbF  = rbFs+".cdf"
	#Labs[0] = "RBSP-" + rbStr.upper()

	#Get RB intensity
	print("Reading from %s"%rbStr)
	Trb,Krb,dkrb,Irb = kc.GetRBSP(rbF,T0Str,tMin=tMin,tMax=tMax,rbID=rbFs,rbSK=rbSK)

	#Get RB trajectory data
	#Tsc = seconds after T0
	Tsc,Xsc,Ysc,Zsc = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk=rbSKt,doEq=False)

	#Create figure
	fig = plt.figure(figsize=figSize)
	gs = gridspec.GridSpec(Nkap,Nkt)
	#Get simulation KCyls, only ALL
	for nk in range(Nkap):
		for nt in range(Nkt):
			aI = []
			aT = []
			aK = []

			#Get simulation K-Cyls
			for n in range(NumPop):
				fIn = os.path.expanduser('~') + Base + kcStrs[n] + "_" + kappaStr[nk] + kTSclStr[nt] + ".h5"
				print("Reading %s"%fIn)

				R,P,K,Tkc,I0 = kc.getCyl(fIn)
				I0 = kcScls[n]*I0
				SimKC = [R,P,K,Tkc,I0]
				rbDat = [Xsc,Ysc,Zsc,Tsc,Krb]
				Is,Isc = kc.InterpSmooth(SimKC,rbDat)

				#Save contributions
				aI.append(Isc)
				aT.append(Tsc)
				aK.append(Krb)

			#Get combined intensity
			Ikc = aI[0]+aI[1]
			gsIJ = gridspec.GridSpecFromSubplotSpec(NKL,1,subplot_spec=gs[nk,nt])
			Tp = kc.Ts2date(Tsc,T0Str)

			#Create holder for energy interpolant
			Nt = len(Tsc)
			iPts = np.zeros((Nt,2))
			iPts[:,0] = Tsc
			#Create interpolants
			rbIi = scipy.interpolate.RegularGridInterpolator((Trb,Krb),Irb,method='linear',bounds_error=False)
			kcIi = scipy.interpolate.RegularGridInterpolator((Tsc,Krb),Ikc,method='linear',bounds_error=False)
			for nl in range(NKL):
				Ax = fig.add_subplot(gsIJ[nl,0])
				#Ax = plt.subplot(gsIJ(nl,0))
				iPts[:,1] = KLs[nl]
				Ik0 = rbIi(iPts)
				Ik1 = kcIi(iPts)
				Ax.semilogy(Tp,Ik0,'k',linewidth=lw)
				Ax.semilogy(Tp,Ik1,'m',linewidth=lw)
				plt.ylim([vMins[nl],vMaxs[nl]])
				#Add labels
				KLab = "%s keV"%(str(KLs[nl]))
				Ax.text(0.05,0.75,KLab,transform=Ax.transAxes,fontsize='medium')
				if (nk==0 and nl==0):
					#Ax.xaxis.tick_top()
					#Ax.xaxis.set_label_position('top')
					#Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
					plt.setp(Ax.get_xticklabels(),visible=False)
				elif (nk<Nkap-1 or nl<NKL-1):
					plt.setp(Ax.get_xticklabels(),visible=False)
				else:
					Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
				if (nl==0):
					TitS = "Kappa = %2.1f\nTe = %2.2fTi"%(kappas[nk],kTScls[nt])
					plt.title(TitS)	
	figName = "Sweep_RB%s.png"%(rbStr)
	plt.savefig(figName,dpi=figQ)
	plt.close('all')
	