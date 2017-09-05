#Specific routines for StormPSD project

import numpy as np
import datetime
import kCyl as kc
import os
import lfmViz as lfmv
import scipy
import scipy.interpolate

T0Str = "2013-03-16T17:10:00Z"
figQ = 300 #DPI

#Globals
Base = os.path.expanduser('~') + "/Work/StormPSD/Data"
VTIDir = Base + "/eqSlc/"
InjWs = [0,21,3]

dtVTI = 150.0
T0VTI = 30000.0

tMin = 33600.0
tMax = 189000.0

#Scale factors for both populations
injScl = 2.5
trapScl = 0.5

#Smoothing defaults
nSm = 1
NTWin = 2

#Defaults for RB data
fOrbA = Base + "/VAP/vaporbRBA.txt"
fOrbB = Base + "/VAP/vaporbRBB.txt"
fRBa = Base + "/VAP/rbspa.cdf"
fRBb = Base + "/VAP/rbspb.cdf"

rbAC = "cyan"
rbBC = "magenta"

#Defaults for weighting info
fRBw = Base + "/VAP/rbWgt.cdf"

def getFld(t):
	tSlc = np.int( (t-T0VTI)/dtVTI )
	vtiFile = VTIDir + "eqSlc.%04d.vti"%(tSlc)

	dBz = lfmv.getVTI_SlcSclr(vtiFile).T
	ori,dx,ex = lfmv.getVTI_Eq(vtiFile)
	xi = ori[0] + np.arange(ex[0],ex[1]+1)*dx[0]
	yi = ori[1] + np.arange(ex[2],ex[3]+1)*dx[1]

	return xi,yi,dBz

def InjCyl():
	fIns = []
	NumW = len(InjWs)
	for i in range(NumW):
		fIn = Base + "/Merge/KCyl_StormI_%d.h5"%(InjWs[i])
		fIns.append(fIn)
	print(fIns)
	R,P,K,Tkc,I0 = kc.getCyls(fIns)
	I0 = injScl*I0
	return R,P,K,Tkc,I0

def TrapCyl():
	fIn = Base + "/Merge/KCyl_StormT.h5"
	R,P,K,Tkc,I0 = kc.getCyl(fIn)
	I0 = trapScl*I0
	return R,P,K,Tkc,I0

def TotCyl(doSmooth=True):
	R,P,K,Tkc,tI = TrapCyl()
	_,_,_,_,iI = InjCyl()
	I0 = tI + iI
	if (doSmooth):
		I0 = kc.SmoothKCyl(R,P,I0,nSm)
	return R,P,K,Tkc,I0

#Get RB positions
def GetRBs():
	Tsc,Xa,Ya,Za = kc.getTraj(fOrbA,T0Str,tMin,tMax,Nsk=1,doEq=True)
	Tsc,Xb,Yb,Zb = kc.getTraj(fOrbB,T0Str,tMin,tMax,Nsk=1,doEq=True)
	Xrb = [Xa,Xb]
	Yrb = [Ya,Yb]
	Zrb = [Za,Zb]

	return Tsc,Xrb,Yrb,Zrb

#Get I(K,t) lines for RB A/B
def GetRBKt(t,Ks):
	#Get RB data
	TrbA,KrbA,_,IrbA = kc.GetRBSP(fRBa,T0Str,tMin,tMax,"rbspa")
	TrbB,KrbB,_,IrbB = kc.GetRBSP(fRBb,T0Str,tMin,tMax,"rbspb")

	#Create interpolants
	Ia = scipy.interpolate.RegularGridInterpolator((TrbA,KrbA),IrbA,method='linear',bounds_error=False,fill_value=0.0)
	Ib = scipy.interpolate.RegularGridInterpolator((TrbB,KrbB),IrbB,method='linear',bounds_error=False,fill_value=0.0)

	IkAs = []
	IkBs = []
	NumK = len(Ks)
	NumT = len(t)
	iPts = np.zeros((NumT,2))
	iPts[:,0] = t
	for i in range(NumK):
		iPts[:,1] = Ks[i]
		IkA = Ia(iPts)
		IkB = Ib(iPts)
		
		IkAs.append(IkA)
		IkBs.append(IkB)
	return IkAs,IkBs

#Get 2D I(K,t) from RB A/B
def GetRB_I2D():
	#Get RB data
	TrbA,KrbA,_,IrbA = kc.GetRBSP(fRBa,T0Str,tMin,tMax,"rbspa")
	TrbB,KrbB,_,IrbB = kc.GetRBSP(fRBb,T0Str,tMin,tMax,"rbspb")
	
	Trbs = [TrbA,TrbB]
	Krbs = [KrbA,KrbB]
	Irbs = [IrbA,IrbB]
	return Trbs,Krbs,Irbs	
#Get simulated I(K,t) lines for RB trajectories
#SimKC = [R,P,K,Tkc,Is]
#rbDat = [Xsc,Ysc,Zsc,Tsc,Ksc]

def GetSimRBKt(SimKC,rbDat,Ks,Nsk=1):
	Tsc,Xrb,Yrb,Zrb = rbDat
	R,P,K,Tkc,Ikc = SimKC

	Ksc = np.array(Ks)
	#Already smoothed if gonna smooth
	Ii = kc.GetInterp(R,P,K,Tkc,Ikc)
	IkA = kc.InterpI_XYZ(Ii,Xrb[0],Yrb[0],Zrb[0],Tsc,Ksc,doScl=True,en=2)
	IkB = kc.InterpI_XYZ(Ii,Xrb[1],Yrb[1],Zrb[1],Tsc,Ksc,doScl=True,en=2)

	sIkAs = []
	sIkBs = []
	NumK = len(Ks)
	for n in range(NumK):
		IknA = TWin(IkA[:,n],Nw=NTWin)
		IknB = TWin(IkB[:,n],Nw=NTWin)
		sIkAs.append(IknA[::Nsk])
		sIkBs.append(IknB[::Nsk])

	Tsc = Tsc[::Nsk]
	return Tsc,sIkAs,sIkBs

#Average over time window
def TWin(Ik,Nw=4):
	Nt = Ik.shape[0]
	IkS = np.zeros(Nt)
	IkS[:] = Ik[:]
	if (Nw>0):
		for i in range(Nt):
			i0 = np.maximum(i-Nw,0)
			i1 = np.minimum(i+Nw,Nt-1)
			IkS[i] = Ik[i0:i1].mean()
	return IkS

#Get DST
def GetDST():
	Ts,dst = kc.GetRBSP_DST(fRBa,T0Str,tMin,tMax,rbID="rbspa")
	return Ts,dst

#Get L3 data for f(L,K) PSD
def GetRBPSD():
	from spacepy import pycdf
	cdf = pycdf.CDF(fRBw)
	print(cdf)
	#Get main data
	Tdt = cdf['rbspa_ect-mageis_l3_time_epoch'][...]
	L = cdf['rbspa_ect-mageis_l3_L'][...]
	K = cdf['rbspa_ect-mageis_l3_FEDU_CORR_T5_channel_energy'][...]
	Ilk = cdf['rbspa_ect-mageis_l3_FEDU_CORR_T5'][...]
	cdf.close()

	#Kill bad energies
	#Ilk[K<=0] = 0.0

	#Silly hard-coded numbers to find orbit
	i0 = 7400
	i1 = 10500

	iMin = L[i0:i1].argmin()
	iMax = L[i0:i1].argmax()

	#New index bounds
	j0 = i0+iMin
	j1 = i0+iMax

	#Cut down to what we want
	L = L[j0:j1]
	K = K[j0:j1,:]
	Ilk = Ilk[j0:j1,:]

	#Find good energy channels
	kMax = K.max(axis=0)
	Ik = (kMax >= 0)

	K = kMax[Ik]
	Ilk = Ilk[:,Ik]

	#More cutting
	I0 = Ilk.max(axis=0)
	Ik = (I0 >= 1.0e-8)

	K = K[Ik]
	Ilk = Ilk[:,Ik]
	return L,K,Ilk