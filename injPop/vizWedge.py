#Visualize time series pulled from wedge
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import cPickle as pickle
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import scipy.interpolate as interpolate
import lfmPreproc as lfmpre
import datetime
import lfmViz as lfmv
import matplotlib.dates as mdates

def injRate(fPkl):
	#Load data
	with open(fPkl,"rb") as f:
		t = pickle.load(f)
		Vst = pickle.load(f)
		kTt = pickle.load(f)
		Nt = pickle.load(f)
		Nkt = pickle.load(f)	
	#Create probability distribution
	fTot = Vst*Nkt
	f = fTot/fTot.sum()	
	return t,f
lfmv.ppInit()

T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"
tm = 30000
tM = 120000
fPkls = ["tsWedge_0.pkl","tsWedge_21.pkl","tsWedge_3.pkl"]
wLabs = ["Midnight","2100 MLT","0300 MLT"]
Nsk = 10

T0 = datetime.datetime.strptime(T0Str,T0Fmt)


#Make figure
figSize = (20,4)
figQ = 300 #DPI
fig = plt.figure(1,figsize=figSize)

NumW = len(fPkls)

for nw in range(NumW):
	t,f = injRate(fPkls[nw])
	i0 = (t>tm).argmax()
	i1 = (t<tM).argmin()

	tC = t[i0:i1:Nsk]
	fC = f[i0:i1:Nsk]

	#Presentation figure
	tD = []
	Nt = len(tC)
	for i in range(Nt):
		dt = T0 + datetime.timedelta(seconds=tC[i])
		tD.append(dt)
	plt.plot(tD,fC*100,linewidth=0.75)

plt.ylabel('Injection Rate',fontsize='large')
fig.autofmt_xdate()
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
plt.legend(wLabs)
plt.show()
