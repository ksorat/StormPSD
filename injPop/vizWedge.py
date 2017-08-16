#Figures for Injection rate

import os
import kCyl as kc
import numpy as np
import cPickle as pickle
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import lfmPreproc as lfmpre
import lfmViz as lfmv

def tWindow(t,Q,dt):
        #Window time series t,Q based on window size dt
        Nt = len(t)
        Qw =  np.zeros(Nt)
        Qw[:] = Q[:]
        J = (Q>0)
        for i in range(Nt):
                t0 = t[i]
                I = (np.abs(t-t0) <= dt)
                IJ = I & J
                if (IJ.sum() > 0):
                        Qw[i] = Q[IJ].mean()
                else:
                        Qw[i] = 0.0

        return Qw

lfmv.ppInit()
T0Cut = 40000.0

doSmoothTS = True
figQ = 300
dtW = 150.0

tsID = ["0","21","3"]

NumPop = len(tsID)
T0s = []
Np = 1000000

for n in range(NumPop):

        fPkl = "tsWedge_%s.pkl"%(tsID[n])
        fTab = "tsWedge_%s.csv"%(tsID[n])
      
        with open(fPkl,"rb") as f:
                t = pickle.load(f) #Time [s]
                Vst = pickle.load(f) #Earthward tail velocity [km/s]
                kTt = pickle.load(f) #Thermal energy, kT [keV]
                Nt  = pickle.load(f) #Number density, [#/cm3]
                Nkt = pickle.load(f)

        I = (t<T0Cut)
        Vst[I] = 0.0
        fTot = Vst*Nkt

        if (doSmoothTS):
        	fTotW = tWindow(t,fTot,dtW)
        else:
        	fTotW = fTot

        T0p = lfmpre.genSample(t,fTotW,Np)

        T0s.append(T0p)

Nb = 50
plt.hist(T0s,Nb)
plt.legend(tsID)
plt.xlim([30000,195000])
plt.show()