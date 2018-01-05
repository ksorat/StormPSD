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
KCrit0 = 750.0
KCrit1 = 5000.0

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

#Get 
P = 180.0*np.arctan2(Yeq,Xeq)/np.pi
I = P<0
P[I] = P[I] + 360.0
#mlt = (P/15.0)+12



#--------------
#T0 = 46800.0
T0 = 45000.0
TLoss = (Teq-T0)/(60)


figSize = (8,8)
dpiQ = 300

#------------------
#Histogram for time of losses (both)
Tfin = 120
Ntb = 600

Legs = ["Magnetopause","Atmosphere"]
Tpb = np.linspace(0,Tfin,Ntb)
figName = "LossT.png"
fig = plt.figure(num=1,figsize=figSize)

TLs = (TLoss[Imp],TLoss[Iatm])
plt.hist(TLs,Tpb,alpha=1.0,log=True,histtype='bar',cumulative=True)
plt.ylim([1,100000])
plt.xlim([0,Tfin])
plt.xlabel("Time after compression [m]")
plt.ylabel("Cumulative TPs")
plt.title("Losses from Initial Population")
plt.legend(Legs,loc='upper left')
plt.savefig(figName,dpi=dpiQ)
lfmv.trimFig(figName)
plt.close('all')

#------------------
#Focused view of precipitation
Tfin = 120
Ntb = 600
Tpb = np.linspace(0,Tfin,Ntb)
figName = "LossTatm.png"
fig = plt.figure(num=1,figsize=figSize)
plt.hist(TLoss[Iatm],Tpb)
plt.xlim([0,Tfin])
plt.xlabel("Time after compression [m]")
plt.ylabel("Instantaneous TPs")
plt.title("Precipitation Losses from Initial Population")
plt.savefig(figName,dpi=dpiQ)
lfmv.trimFig(figName)
plt.close('all')

#------------------
#L/MLT of precipitation
figName = "AtmEQ.png"
fig = plt.figure(num=1,figsize=figSize)

Nr = 30
Np = 48*2
rOut = 6.0*2
vNorm = LogNorm(vmin=1.0,vmax=25)

deg2rad = np.pi/180.0
rE = np.linspace(Rcrit,rOut,Nr+1)
pE = np.linspace(0,2*np.pi,Np+1)
H,_,_ = np.histogram2d(Req[Iatm],deg2rad*P[Iatm],[rE,pE])
Phi,R = np.meshgrid(pE,rE)
xx = R*np.cos(Phi)
yy = R*np.sin(Phi)
plt.pcolormesh(xx,yy,H,norm=vNorm,cmap="viridis")
plt.axis('scaled')
lfmv.addEarth2D()
plt.colorbar()
plt.title("Precipitation Losses")
plt.savefig(figName,dpi=dpiQ)
lfmv.trimFig(figName)
plt.close('all')

#------------------
#L/MLT of mp losses
figName = "mpEQ.png"
fig = plt.figure(num=1,figsize=figSize)

Nr = 60*2
Np = 96*2
rOut = 14.0
vNorm = LogNorm(vmin=1.0,vmax=100)

deg2rad = np.pi/180.0
rE = np.linspace(Rcrit,rOut,Nr+1)
pE = np.linspace(0,2*np.pi,Np+1)
H,_,_ = np.histogram2d(Req[Imp],deg2rad*P[Imp],[rE,pE])
Phi,R = np.meshgrid(pE,rE)
xx = R*np.cos(Phi)
yy = R*np.sin(Phi)
plt.pcolormesh(xx,yy,H,norm=vNorm,cmap="viridis")
plt.axis('scaled')
lfmv.addEarth2D()
plt.colorbar()
plt.title("Magnetopause Losses")
plt.savefig(figName,dpi=dpiQ)
lfmv.trimFig(figName)
plt.close('all')


