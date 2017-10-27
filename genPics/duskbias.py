#Generate comparison of I(L,t,K=K0)
#Shell/time dep. of intensity at given channel

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

#K0 value
LW = 1.5
cmeLW = 2
L0 = 3.0
L1 = 6.0

K0 = 1000
figSize = (12,8)
FS = "large"
Col = ["k","b","g","r"]
Legs = ["Initial","Midnight","2100","0300"]
fOut = "ICon.png"
#-------------------


#----------


t,Ikts = pS.ICons(K0,L0,L1-L0,doMin=True,doSmooth=True)
I0 = Ikts[0][0] #Initial K0 intensity from trapped
Tp = kc.Ts2date(t,pS.T0Str)


plt.figure(figsize=figSize)
Ax = plt.gca()

I0 = 1.0
dCME = datetime.datetime.strptime(pS.CME_T0Str,kc.T0Fmt)
for i in range(4):
	plt.semilogy(Tp,Ikts[i]/I0,color=Col[i],linewidth=LW)

plt.ylim([1.0e+2,2.0e+4])
plt.ylim([5.0e+1,2.5e+4])
plt.ylabel("Intensity [cm$^{-2}$ sr$^{-1}$ s$^{-1}$ keV$^{-1}$]",fontsize="large")
Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))	
#Ax.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbBC,linewidth=lwRB)
Ax.axvline(dCME,color=pS.CME_C,linewidth=cmeLW)
Ax.annotate('CME',color=pS.CME_C,xy=(dCME,0.03),xytext=(20,0),textcoords='offset pixels',fontsize=FS)
#Ax.set_xlim(datemin,datemax)
plt.legend(Legs,loc="upper right",ncol=2,fontsize=FS)

plt.savefig(fOut,dpi=pS.figQ)
plt.close('all')
lfmv.trimFig(fOut)
