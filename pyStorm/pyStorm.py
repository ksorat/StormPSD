#Specific routines for StormPSD project

import numpy as np
import datetime
import kCyl as kc
import os
import lfmViz as lfmv
T0Str = "2013-03-16T17:10:00Z"

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

#Defaults for RB data
fOrbA = Base + "/VAP/vaporbRBA.txt"
fOrbB = Base + "/VAP/vaporbRBB.txt"
fRBa = Base + "/VAP/rbspa.cdf"
fRBb = Base + "/VAP/rbspb.cdf"

rbAC = "cyan"
rbBC = "magenta"

def getFld(t):
	tSlc = np.int( (t-T0VTI)/dtVTI )
	vtiFile = VTIDir + "eqSlc.%04d.vti"%(tSlc)
	print("Reading %s"%vtiFile)

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

#Get DST
def GetDST():
	Ts,dst = kc.GetRBSP_DST(fRBa,T0Str,tMin,tMax,rbID="rbspa")
	return Ts,dst