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
Zcrit = 7.5
Tcrit = 35000.0
KCrit = 750.0

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

Ib = ~np.isnan(X) & (K>=KCrit)
X   = X  [Ib] 
Y   = Y  [Ib] 
Z   = Z  [Ib] 
Xeq = Xeq[Ib] 
Yeq = Yeq[Ib] 
Teq = Teq[Ib] 
pID = pID[Ib] 
K   = K  [Ib] 

R = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)

Iatm = (R<=Rcrit) & (Teq>=Tcrit)
Imp = np.abs(Z) >= Zcrit

Tatm = Teq[Iatm]
Tmp  = Teq[Imp]
Patm = xy2phi(Xeq,Yeq,Iatm)
Pmp  = xy2phi(Xeq,Yeq,Imp)

figSize = (8,8)
dpiQ = 300
Al = 0.75
Leg = ["Atmosphere","Magnetopause"]
Nt = 25

#Comparison histogram for losses
Tmin = 734944.1
Tmax = 734945.25
TBins = np.linspace(Tmin,Tmax,Nt)
figName = "LossT.png"
fig = plt.figure(figsize=figSize)
mdAtm = kc.Ts2date(np.double(Tatm),T0Str)
mdMP = kc.Ts2date(np.double(Tmp),T0Str)
TLs = [mdates.date2num(mdAtm),mdates.date2num(mdMP)]
plt.hist(TLs,TBins,alpha=Al,log=True)
#plt.hist(mdates.date2num(mdAtm),TBins,alpha=Al,log=True)
#plt.hist(mdates.date2num(mdMP),TBins,alpha=Al,log=True)
Ax = plt.gca()
Ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%MZ\n%m-%d'))
Ax.xaxis.set_major_locator(mdates.HourLocator(np.arange(0,25,4)))

plt.xlim(Tmin,Tmax)
plt.ylim(10,1.0e+5)
#plt.hist(Tmp ,Nt,alpha=Al,log=True)
plt.legend(Leg)
plt.xlabel("Loss Time")
plt.ylabel("Number of Test Particles")
plt.savefig(figName,dpi=dpiQ)
lfmv.trimFig(figName)
plt.close('all')

#2D histogram, ATM
figName = "Loss2D_Atm.png"
fig = plt.figure(figsize=figSize)
vNorm = LogNorm(vmin=1.0,vmax=5.0e+3)
plt.hist2d(Patm,Tatm,bins=50,norm=vNorm,cmap="viridis")
plt.title("Atmospheric Losses")
plt.xlabel("MLT")
plt.ylabel("Time")
plt.colorbar()
plt.savefig(figName,dpi=dpiQ)
lfmv.trimFig(figName)
plt.close('all')

#2D histogram, MP
figName = "Loss2D_MP.png"
fig = plt.figure(figsize=figSize)
vNorm = LogNorm(vmin=1.0,vmax=5.0e+3)
plt.hist2d(Pmp,Tmp,bins=50,norm=vNorm,cmap="viridis")
plt.title("Magnetopause Losses")
plt.xlabel("MLT")
plt.ylabel("Time")
plt.colorbar()
plt.savefig(figName,dpi=dpiQ)
lfmv.trimFig(figName)
plt.close('all')

