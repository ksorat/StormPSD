#Final kappa sensitivity pic

import tpkCyl as kc
#import okCyl as okc
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import kaiViz as kv

Base = "/Work/StormPSD/Data" + "/f0Sweep/"
KStr = "KCyl_Storm"
kcScls = np.pi*4*np.array([1.0,10.0])

def grabKCyl(kStr,ktStr="b"):
	fIn1 = os.path.expanduser('~') + Base + KStr
	fInT = fIn1 + "T_%s%s.h5"%(kStr,ktStr)
	fInI = fIn1 + "I_%s%s.h5"%(kStr,ktStr)
	R,P,K,Tkc,iI0 = kc.getCyl(fInI)
	R,P,K,Tkc,tI0 = kc.getCyl(fInT)
	I0 = kcScls[0]*tI0 + kcScls[1]*iI0
	I0 = kc.SmoothKCyl(R,P,I0,kc.nSm)

	return R,P,K,Tkc,I0

vN = kv.genNorm(1.0,1.0e+6,doLog=True)
cM = "gnuplot2"

kaps = [3.0,3.5,4.0]
kapStr = ["A","B","C"]
ktStr = "b"
Nk2D = 80

#Get RB data
Tsc,Xrb,Yrb,Zrb,Lrb = pS.GetRBs()
Trbs,Krbs,Irbs = pS.GetRB_I2D()
rbDat = [Tsc,Xrb,Yrb,Zrb,Lrb]

Nkap = len(kaps)
K = np.logspace(1,np.log10(5000),Nk2D)


#for i in range(Nkap):
for i in [0]:	
	Rkc,Pkc,Kkc,Tkc,Ikc = grabKCyl(kapStr[i])
	SimKC = [Rkc,Pkc,Kkc,Tkc,Ikc]
	Ts,Iks = kc.GetSim_I2D(SimKC,rbDat,K)

	
