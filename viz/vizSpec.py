import sys
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
import cPickle as pickle
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib.patches import Wedge
import matplotlib
import lfmViz as lfmv

lfmv.ppInit()
clr='k'
bck='w'
fontsize=16
labelsize=18
titlesize=18
ann_size=20
linewidth=2.0
matplotlib.rc('text',color=clr)
matplotlib.rc('font',size=fontsize)
matplotlib.rc('axes',facecolor=bck,edgecolor=clr,labelcolor=clr,linewidth=linewidth,titlesize=titlesize)
matplotlib.rc('xtick',color=clr)
matplotlib.rc('ytick',color=clr)
matplotlib.rc('savefig',facecolor=bck,edgecolor=bck)

#Load data
inFile = "vaporb.pkl"

with open(inFile, "rb") as f:
	T0Str = pickle.load(f)
	T0Fmt = pickle.load(f)
	Tsc   = pickle.load(f)
	lK    = pickle.load(f)
	Isc   = pickle.load(f)
	Xsc   = pickle.load(f)
	Ysc   = pickle.load(f)

K = 10.0**lK

T0 = datetime.datetime.strptime(T0Str,T0Fmt)

Nsc = Tsc.shape[0]
scT = []
for i in range(Nsc):
	scdt = T0 + datetime.timedelta(seconds=Tsc[i])
	scT.append(scdt)
scT = np.array(scT)

vNorm = LogNorm(vmin=1.0,vmax=1.0e+4)
cMap = "jet"

Is = np.zeros(Isc.shape)
Is[:,:] = Isc[:,:]
Nx = Isc.shape[0]
N = 4
for i in range(N,Nx-N):
	Is[i,:] = np.mean(Isc[i-N/2:i+N/2,:],axis=0)

figSize = (18,6)
figQ = 300 #DPI
fig = plt.figure(figsize=figSize)

plt.pcolormesh(scT,K,4*Is.T,norm=vNorm,cmap=cMap)
plt.ylim([50,5.0e+3])
plt.ylabel("Energy [keV]",fontsize="large")
plt.xlabel("Date",fontsize="large")
plt.yscale('log')
fig = plt.gcf()
fig.autofmt_xdate()

plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
plt.savefig("simSC.png",dpi=figQ)

#plt.gca().xaxis.set_major_locator(mdates.AutoDateLocator())

#ax = plt.gca()
#ax.fmt_xdata = mdates.DateFormatter('%Y-%m-%d')
#plt.colorbar(orientation="horizontal",aspect=50,label="Intensity (ish)")
