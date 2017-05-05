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
import lfmViz as lfmv

lfmv.ppInit()

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

vNorm = LogNorm(vmin=1.0e-1,vmax=1.0e+6)
cMap = "jet"

plt.pcolormesh(scT,K,Isc.T/(np.pi*4),norm=vNorm,cmap=cMap)
plt.ylim([50,1000])
plt.ylabel("Energy [keV]")
plt.yscale('log')
fig = plt.gcf()
fig.autofmt_xdate()

plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
#plt.gca().xaxis.set_major_locator(mdates.AutoDateLocator())

#ax = plt.gca()
#ax.fmt_xdata = mdates.DateFormatter('%Y-%m-%d')
plt.colorbar(orientation="horizontal",aspect=50,label="Intensity (ish)")
