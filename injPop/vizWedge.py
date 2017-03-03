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

fPkl = "tsWedge.pkl"
doPanel = False
doPDF = True
#Load data
with open(fPkl,"rb") as f:
	t = pickle.load(f)
	Vst = pickle.load(f)
	kTt = pickle.load(f)
	Nt = pickle.load(f)
	Nkt = pickle.load(f)


T = t
figSize = (12,12)
figQ = 300 #DPI
if (doPanel):
	fig = plt.figure(figsize=figSize)
	gs = gridspec.GridSpec(4,1)
	
	Ax = fig.add_subplot(gs[0,0])
	plt.plot(t,Vst*1.0e-3)
	plt.ylabel('Sunward Velocity [km/s]')
	
	Ax = fig.add_subplot(gs[1,0])
	plt.plot(t,kTt*1.0e-3)
	plt.ylabel('Temperature [keV]')
	
	Ax = fig.add_subplot(gs[2,0])
	plt.plot(t,Nt*1.0e-6)
	plt.ylabel('Number Density')
	
	Ax = fig.add_subplot(gs[3,0])
	plt.plot(t,Nkt*1.0e-5)
	plt.ylabel('Number Density\nK>10 [keV]')
	plt.xlabel("Seconds")
	#plt.ylim(1.0e+0,1.0e+6)
	
	plt.suptitle("Nightside Wedge Averages\nSt. Patrick's Day 2013")
	plt.savefig("WedgeAvg.png")
	plt.close('all')

#Create probability distribution using bulk/Nk0
Nr = 10000
fTot = Vst*Nkt
f = fTot/fTot.sum()
# nScl = fTot.sum()
# f = fTot/nScl
# #Use normalized PDF to get CDF 
# F = np.cumsum(f)
# #Interpolate the inverse cdf
# cdf = interpolate.interp1d(t,F)
# icdf = interpolate.interp1d(F,t)

# r = np.random.rand(Nr)
# tR = icdf(r)
tR = lfmpre.genSample(t,fTot,Nr)


#Show both
T0 = t.min()
T1 = t.max()
dt = (T1-T0)/len(t)
if (doPDF):
	fig = plt.figure(figsize=figSize)
	Nb = 100
	NumT,eb,patches = plt.hist(tR,Nb,label="Sampled 10k Particles")
	fScl = NumT.max()/f.max()
	plt.plot(t,f*fScl,'r-',label="MHD")
	plt.xlabel("Time")
	plt.ylabel("High-Energy (K>10 keV) Particle Flux")
	plt.xlim(0,150000)
	plt.legend()
	plt.savefig("tInj.png")

	plt.close('all')