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

L0 = 2.1
L1 = 6.0
K0 = 1000.0
K1 = 1100.0

figSize = (12,8)
FS = "large"
Col = ["k","b","g","r"]
Legs = ["Initial","2400","2100","0300"]
Legs = [pS.pLabs[0],pS.iLabs[0],pS.iLabs[1],pS.iLabs[2]]
fOut = "ICons.png"
#-------------------
datemin = datetime.datetime.strptime("2013-03-17T04:00:00Z",kc.T0Fmt)
datemax = datetime.datetime.strptime("2013-03-18T22:00:00Z",kc.T0Fmt)


#----------


t,Ikts = pS.ICons(K0,K1=K1,L=L0,dL=L1-L0,doMin=True,doSmooth=True)
I0 = Ikts[0][0] #Initial K0 intensity from trapped
Tp = kc.Ts2date(t,pS.T0Str)
print("Total I (Initial) I = %f"%(I0))

plt.figure(figsize=figSize)
Ax = plt.gca()

I0 = 1.0
dCME = datetime.datetime.strptime(pS.CME_T0Str,kc.T0Fmt)
IkF = np.zeros(4)
for i in range(4):
	Ik = Ikts[i]

	plt.semilogy(Tp,Ik/I0,color=Col[i],linewidth=LW)
	#print(Ikts[i][-1])
	IkF[i] = Ik[-1]
	print("%s Intensity (Final) = %f"%(Legs[i],IkF[i]))



print("Total I (Final) I = %f"%(IkF.sum()))
print("Enhancement, I/I0 = %f"%(IkF.sum()/Ikts[0][0]))
plt.ylim([1.0e+2,2.0e+4])
plt.ylim([2.0e-1,2.0e+3])
#plt.ylim([1.0e-2,1.0e+4])

plt.ylabel("Intensity [cm$^{-2}$ sr$^{-1}$ s$^{-1}$ keV$^{-1}$]",fontsize="large")
Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
Ax.tick_params(labeltop=False,labelright=True,right=True,which='both')
#Ax.axvline(kc.Date2Num(Tkc[n],pS.T0Str),color=pS.rbBC,linewidth=lwRB)
Ax.axvline(dCME,color=pS.CME_C,linewidth=cmeLW)
Ax.annotate('CME',color=pS.CME_C,xy=(dCME,100),xytext=(10,-10),textcoords='offset pixels',fontsize=FS)
Ax.set_xlim(datemin,datemax)
plt.legend(Legs,loc="lower right",ncol=2,fontsize=FS)

plt.savefig(fOut,dpi=pS.figQ)
plt.close('all')
lfmv.trimFig(fOut)
