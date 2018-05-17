#Final kappa sensitivity pic

import tpkCyl as kc
#import okCyl as okc
import cPickle as pickle
import pyStorm as pS
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
pData = "kSens.pkl"

def grabKCyl(kStr,ktStr="b"):
	fIn1 = os.path.expanduser('~') + Base + KStr
	fInT = fIn1 + "T_%s%s.h5"%(kStr,ktStr)
	fInI = fIn1 + "I_%s%s.h5"%(kStr,ktStr)
	R,P,K,Tkc,iI0 = kc.getCyl(fInI)
	R,P,K,Tkc,tI0 = kc.getCyl(fInT)
	I0 = kcScls[0]*tI0 + kcScls[1]*iI0
	I0 = kc.SmoothKCyl(R,P,I0,2*pS.nSm)

	return R,P,K,Tkc,I0

vN = kv.genNorm(1.0,1.0e+6,doLog=True)
cM = "gnuplot2"

kaps = [3.0,3.5,4.0]
kapStr = ["A","B","C"]
Nkap = len(kaps)
ktStr = "b"
Nk2D = 80
K = np.logspace(1,np.log10(5000),Nk2D)
rbI = 0 #RB-A

#

if (os.path.isfile(pData)):
	print("Loading data ...")
	with open(pData,"rb") as f:
		Trbs  = pickle.load(f)
		Krbs  = pickle.load(f)
		Irbs  = pickle.load(f)
		Ts    = pickle.load(f)
		Ik2Ds = pickle.load(f)
else:
	print("Generating data ...")
	#Get RB data
	Tsc,Xrb,Yrb,Zrb,Lrb = pS.GetRBs()
	Trbs,Krbs,Irbs = pS.GetRB_I2D()
	rbDat = [Tsc,Xrb,Yrb,Zrb,Lrb]
	Trbs,Krbs,Irbs = pS.GetRB_I2D()
	Ik2Ds = []
	Ik2Ds.append(Irbs[rbI])

	for i in range(Nkap):
		Rkc,Pkc,Kkc,Tkc,Ikc = grabKCyl(kapStr[i])
		SimKC = [Rkc,Pkc,Kkc,Tkc,Ikc]
		Ts,Iks = pS.GetSim_I2D(SimKC,rbDat,K)
		Ts = Ts[rbI]
		Ik = Iks[rbI]
		Ik2Ds.append(Ik)
	print("Writing pickle ...")
	with open(pData, "wb") as f:
		pickle.dump(Trbs ,f)
		pickle.dump(Krbs ,f)
		pickle.dump(Irbs ,f)
		pickle.dump(Ts   ,f)
		pickle.dump(Ik2Ds,f)
		

figSize = (8,16)
fig = plt.figure(figsize=figSize)
HRs = [20,20,20,20,1]

gs = gridspec.GridSpec(5,1,height_ratios=HRs)

Tpc = [Trbs[rbI],Ts,Ts,Ts]
Kpc = [Krbs[rbI],K,K,K]

for i in range(4):
	Ax = fig.add_subplot(gs[i,0])
	Tp = kc.Ts2date(Tpc[i],pS.T0Str)
	iPC = Ax.pcolormesh(Tp,Kpc[i],Ik2Ds[i].T,norm=vN,cmap=cM)
	plt.yscale('log')
	plt.ylim([75,5.0e+3])

plt.show()
