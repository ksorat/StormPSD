#Makes injection rate pic

import kCyl as kc
import pyStorm as pS
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import lfmViz as lfmv

nDT = 3
figSize = (12,6)
datemin = datetime.datetime.strptime("2013-03-17T04:00:00Z",kc.T0Fmt)
#datemax = datetime.datetime.strptime("2013-03-18T06:00:00Z",kc.T0Fmt)
datemax = datetime.datetime.strptime("2013-03-18T02:00:00Z",kc.T0Fmt)
dCME = datetime.datetime.strptime(pS.CME_T0Str,kc.T0Fmt)

LW = 1
cmeLW = 2

FS = "large"
t,aJt = pS.wIRate(nDT=nDT)
Tp = kc.Ts2date(t,pS.T0Str)

plt.figure(figsize=figSize)
Ax = plt.gca()
NumW = len(aJt)
for i in range(NumW):
	plt.plot(Tp,aJt[i],color=pS.iCols[i],linewidth=LW)
Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))	
#Ax.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbBC,linewidth=lwRB)
Ax.axvline(dCME,color=pS.CME_C,linewidth=cmeLW)
Ax.annotate('CME',color=pS.CME_C,xy=(dCME,0.03),xytext=(20,0),textcoords='offset pixels',fontsize=FS)
Ax.set_xlim(datemin,datemax)
plt.legend(pS.iLabs,fontsize=FS)
plt.ylabel("Injection Rate",fontsize=FS)
plt.savefig("WedgeInj.png",dpi=pS.figQ)
plt.close('all')
lfmv.trimFig("WedgeInj.png")
