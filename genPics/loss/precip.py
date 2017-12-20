#Various pics for trapped particle precipitation

import numpy as np
import lfmPostproc as lfmpp
import cPickle as pickle
import os
import datetime
import matplotlib.pyplot as plt
import kCyl as kc
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import lfmViz as lfmv

def xy2phi(X,Y,I):
	P = 180.0*np.arctan2(Y[I],X[I])/np.pi
	I = P<0
	P[I] = P[I] + 360
	return P

#Time data
T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"

fPkl = "Trap.pkl"
lfmv.ppInit()

#Critical values
Rcrit = 2.1
Zcrit = 7.0
TCrit = 35000.0
KCrit0 = 1000.0
KCrit1 = 2500.0

#Load trapped data from PKL
with open(fPkl,"rb") as f:
	X   = pickle.load(f)
	Y   = pickle.load(f)
	Z   = pickle.load(f)
	Xeq = pickle.load(f)
	Yeq = pickle.load(f)
	Teq = pickle.load(f)
	pID = pickle.load(f)
	K   = pickle.load(f)

Ib = ~np.isnan(X) & (K>=KCrit0) & (K<=KCrit1) & (Teq>=TCrit)

X   = X  [Ib] 
Y   = Y  [Ib] 
Z   = Z  [Ib] 
Xeq = Xeq[Ib] 
Yeq = Yeq[Ib] 
Teq = Teq[Ib] 
pID = pID[Ib] 
K   = K  [Ib] 

R = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
Req = np.sqrt(Xeq**2.0 + Yeq**2.0)

Iatm = (R<=Rcrit) & (Teq>=TCrit) & (Req>=Rcrit)
Imp = (np.abs(Z) >= Zcrit) & (Teq>=TCrit)

Tatm = Teq[Iatm]
Tmp  = Teq[Imp]
Patm = xy2phi(Xeq,Yeq,Iatm)
Pmp  = xy2phi(Xeq,Yeq,Imp)


#--------------
#T0 = 46800.0
T0 = 45000.0
TLoss = (Teq-T0)/(60)


figSize = (8,8)
dpiQ = 300

#------------------
#Histogram for time of losses
Tfin = 180
Ntb = 600
Legs = ["Atmosphere","Magnetopause"]
Tpb = np.linspace(0,Tfin,Ntb)
figName = "LossT.png"
fig = plt.figure(num=1,figsize=figSize)
TLs = [TLoss[Iatm],TLoss[Imp]]
plt.hist(TLs,Tpb,alpha=0.75,color=['b','r'],log=True,histtype='stepfilled',cumulative=True)
plt.ylim([1,100000])
plt.xlim([0,Tfin])
plt.xlabel("Time after compression [m]")
plt.ylabel("Cumulative TP's Lost")
plt.title("Losses from Initial Population")
plt.legend(Legs)
plt.savefig(figName,dpi=dpiQ)
lfmv.trimFig(figName)
plt.close('all')

#------------------

# #2D histogram, MP
# figName = "Loss2D_MP.png"
# fig = plt.figure(figsize=figSize)
# vNorm = LogNorm(vmin=1.0,vmax=5.0e+3)
# plt.hist2d(Pmp,Tmp,bins=50,norm=vNorm,cmap="viridis")
# plt.title("Magnetopause Losses")
# plt.xlabel("MLT")
# plt.ylabel("Time")
# plt.colorbar()
# plt.savefig(figName,dpi=dpiQ)
# lfmv.trimFig(figName)
# plt.close('all')


# MS = 4
# plt.plot(Xeq[Iatm],Yeq[Iatm],'b.',markersize=MS)
# plt.plot(Xeq[Imp ],Yeq[Imp ],'r.',markersize=MS)

# If = Iatm & (Req>=2.1)

# T0 = 45000
# Ntb = 1000
# Tpb = np.linspace(0,200,Ntb)
# #plt.hist(TLoss[Iatm],Tpb)

# If = If & (TLoss <= 240)
# tNorm = mpl.colors.Normalize(vmin=0.0,vmax=200.0)
# plt.scatter(Xeq[If],Yeq[If],c=TLoss[If],norm=tNorm,cmap="viridis",alpha=0.5,edgecolors='k')
# plt.axis('equal')
# cbar = plt.colorbar()
# cbar.set_alpha(1)
# cbar.draw_all()
# plt.show()
