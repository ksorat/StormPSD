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
K0 = 1000
figSize = (8,12)
vNorm = LogNorm(vmin=1.0e0,vmax=1.0e+4)
cMap = "viridis"
FS = "large"
#-------------------
#Get data
R,P,K,Tkc,I = pS.TotCyl(doSmooth=True)

k0 = np.abs(K-K0).argmin()
Ip = I.mean(axis=1)
Ilt = Ip[:,k0,:]

#Prep for figure
lfmv.ppInit()
plt.close('all')

#--------------
#L,t
fname = "Ilt.png"
fig = plt.figure(figsize=figSize)
gs = gridspec.GridSpec(3,1,height_ratios=[10,0.1,0.2])
Ax = fig.add_subplot(gs[0,0])
AxC = fig.add_subplot(gs[-1,0])
Ax.set_xlim([2.5,12])

Ax.pcolormesh(R,Tkc,Ilt.T,norm=vNorm,cmap=cMap)
cb = mpl.colorbar.ColorbarBase(AxC,cmap=cMap,norm=vNorm,orientation='horizontal')
cb.set_label("Intensity [cm-2 sr-1 s-1 kev-1]",fontsize=FS)
plt.savefig(fname,dpi=pS.figQ)
plt.close('all')
lfmv.trimFig(fname)

#----------
LW = 1.5
L0 = 4.0
L1 = 8.0
Col = ["k","b","g","r"]
t,Ikts = pS.ICons(K0,L0,L1-L0,doMin=True,doSmooth=True)
I0 = Ikts[0][0] #Initial K0 intensity from trapped
Tp = t


# plt.plot(Tp,Ikts[0]/I0,'k',linewidth=LW)
# plt.plot(Tp,Ikts[1]/I0,'b',linewidth=LW)
# plt.plot(Tp,Ikts[2]/I0,'g',linewidth=LW)
# plt.plot(Tp,Ikts[3]/I0,'r',linewidth=LW)
I0 = 1.0
for i in range(4):
	plt.semilogy(Tp,Ikts[i]/I0,color=Col[i],linewidth=LW)

plt.ylim([10,1.0e+4])
plt.show()
#--------------
#L,t
# fname = "Ixt.png"
# Irp = I[:,:,k0,:]
# Nt = len(Tkc)
# Nr = len(R)
# Nx = 2*Nr+1
# X = np.zeros(Nx)
# Ixt = np.zeros((Nx,Nt))
# P1 = len(P)/2

# for i in range(Nt):
# 	for n in range(Nr):
# 		np = n+Nr+1
# 		X[n] = -R[n]
# 		X[np] = R[n]
# 		Ixt[n,i] = 0.5*(Irp[n,0,i] + Irp[n,-1,i])
# 		Ixt[np,i] = 0.5*(Irp[n,P1,i] + Irp[n,P1-1,i])



# fig = plt.figure(figsize=figSize)
# gs = gridspec.GridSpec(3,1,height_ratios=[10,0.1,0.2])
# Ax = fig.add_subplot(gs[0,0])
# AxC = fig.add_subplot(gs[-1,0])
# Ax.pcolormesh(X,Tkc,Ixt.T,norm=vNorm,cmap=cMap)
# plt.show()

