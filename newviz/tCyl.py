#Test of kCyl reader
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

doPanel = True

T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"
tMin = 30002.0
tMax = 110000.0

Nk = 40 #Number of K samples
IScl = 1.0/(np.pi)

OrbF = "vaporbit.txt"
rbF = "rbspa.cdf"

Nsk = 1 #Skip number for trajectory
fIn = os.path.expanduser('~') + "/Work/StormPSD/Data" + "/Merge/KCyl.h5"

#Interpolate from simulation
R,P,K,t,I = kc.getCyl(fIn)
#I = kc.SmoothI(I,sig=0.25)
Ii = kc.GetInterp(R,P,K,t,I)
kMin = K.min()
kMax = K.max()
Ksc = np.logspace(np.log10(kMin),np.log10(kMax),Nk)

#Get trajectory data
#Tsc = seconds after T0

Tsc,Xsc,Ysc,Z = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk=Nsk)
Rsc = np.sqrt(Xsc**2.0 + Ysc**2.0)
Psc = np.arctan2(Ysc,Xsc)
iP = (Psc<0); Psc[iP] = Psc[iP]+2*np.pi

Nsc = len(Tsc)

iPts = np.zeros((Nk,4))

#Loop over spacecraft trajectory, get I(K) at each point
Isc = np.zeros((Nsc,Nk))
for i in range(Nsc):
	r = Rsc[i]
	p = Psc[i]
	t = Tsc[i]
	iPts[:,0] = r
	iPts[:,1] = p
	iPts[:,2] = Ksc
	iPts[:,3] = t
	Isc[i,:] = Ii(iPts)
Isc = IScl*Isc

#Now get RBSP data
Trb,Krb,dkrb,Irb = kc.GetRBSP(rbF,T0Str,tMin=tMin,tMax=tMax)

#Get DST data
Tdst,dst = kc.GetRBSP_DST(rbF,T0Str,tMin=tMin,tMax=tMax)

#Do RB/SIM panel figure
if (doPanel):
	figSize = (8,8)
	figQ = 300 #DPI
	figName = "rbsimI.png"
	vMin = 1.0
	vMax = 1.0e+5
	
	cMap = "jet"

	fig = plt.figure(figsize=figSize)
	Np = 2

	#Plots: DST,RBSP,SIM,BLANK,COLOR
	idst = 0; ipan=1
	gs = gridspec.GridSpec(1+Np+2,1,height_ratios=[2,10,10,1,1])
	vNorm = LogNorm(vmin=vMin,vmax=vMax)

	#DST
	Ax = fig.add_subplot(gs[idst,0])
	xdst = kc.Ts2date(Tdst,T0Str)

	Ax.plot(xdst,dst)
	plt.setp(Ax.get_yticklabels(),visible=False)

	#Both panels
	Tp = [kc.Ts2date(Trb,T0Str),kc.Ts2date(Tsc,T0Str)]
	Kp = [Krb,Ksc]
	Ip = [Irb,Isc]
	for i in range(Np):
		Ax = fig.add_subplot(gs[ipan+i,0])
		iPlt = Ax.pcolormesh(Tp[i],Kp[i],Ip[i].T,norm=vNorm,cmap=cMap)
		plt.ylim([50,5.0e+3])
		plt.ylabel("Energy [keV]",fontsize="large")
		#plt.xlabel("Date",fontsize="large")
		plt.yscale('log')

	#Do colorbar
	AxC = fig.add_subplot(gs[-1,0])
	cb = mpl.colorbar.ColorbarBase(AxC,cmap=cMap,norm=vNorm,orientation='horizontal')
	plt.show()

# plt.pcolormesh(Tsc,Ksc,Isc.T,norm=vNorm,cmap=cMap)
# plt.ylim([50,5.0e+3])
# plt.ylabel("Energy [keV]",fontsize="large")
# plt.xlabel("Date",fontsize="large")
# plt.yscale('log')
# plt.colorbar()