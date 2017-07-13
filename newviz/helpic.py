#Make some helpful pics for geometry

import kCyl as kc
import os
import numpy as np

import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import lfmViz as lfmv

#Time data
T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"
#Labs = ["NULL","Trapped","Injected","Combined"]
figQ = 300
tMin = 33600.0
tMax = 189000.0
#tMax = 195000.0
rbSKt = 1

#RB Opts
rbStrs = ["A","B"]
NumRB = len(rbStrs)

Xs = []
Ys = []
Zs = []
Xps = []
Yps = []
for nrb in range(NumRB):
	#Create relevant files
	rbStr = rbStrs[nrb]

	OrbF = "vaporbRB" + rbStr.upper() + ".txt"
	rbFs = "rbsp" + rbStrs[nrb].lower()
	rbF  = rbFs+".cdf"
	#Labs[0] = "RBSP-" + rbStr.upper()

	#Get RB trajectory data
	#Tsc = seconds after T0
	Tsc,Xsc,Ysc,Z   = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk=rbSKt,doEq=False)
	Tsc,XscP,YscP,Zd = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk=rbSKt,doEq=True)
	Xs.append(Xsc)
	Ys.append(Ysc)
	Zs.append(Z)
	Xps.append(XscP)
	Yps.append(YscP)


Tp = kc.Ts2date(Tsc,T0Str)

#Make time plot
plt.figure(figsize=(12,4))
plt.plot(Tp,Tsc,'r')
plt.ylabel('Seconds')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
plt.savefig("TimeConv.png",dpi=figQ)
plt.close('all')

#Make trajectory plot
Cols = ['r','b','m']
figSize = (12,12)
fig = plt.figure(figsize=figSize)
gs = gridspec.GridSpec(2,1)

for n in range(NumRB):
	Ax = fig.add_subplot(gs[n,0])
	L = np.sqrt(Xps[n]**2.0+Yps[n])
	Ax.plot(Tp,Xs[n],'r',Tp,Ys[n],'b',Tp,Zs[n],'m',Tp,L,'k',Tp,Xps[n],'r--',Tp,Yps[n],'b--')
	TitS = "RB-%s"%(rbStrs[n])
	plt.title(TitS)
	Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
	plt.legend(["X","Y","Z","L"])
plt.savefig("RBTraj.png",dpi=figQ)
plt.close('all')
