#Create h5 data for pre-existing population

import kCyl as kc
import pyStorm as pS
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import scipy
import scipy.interpolate
import h5py
import lfmViz as lfmv

fOut = "rbInterp.h5"

def SmoothI(I,sig=1.0):
	import scipy
	import scipy.ndimage
	from scipy.ndimage.filters import gaussian_filter

	I = gaussian_filter(I,sig)
	return I

vNorm = LogNorm(vmin=1.0e-1,vmax=1.0e+5)

def genFig(L,K,I,fOut):
	plt.pcolormesh(L,K,I.T,norm=vNorm,cmap="jet")
	plt.xlim([1.5,8])
	plt.ylim([10,5000])
	plt.gca().set_yscale('log')
	cb = plt.colorbar()
	cb.set_label('Intensity\ncm-2 sr-1 s-1 keV-1')
	plt.xlabel('L [Re]')
	plt.ylabel('Energy [keV]')
	plt.savefig(fOut)
	plt.close('all')
	lfmv.trimFig(fOut)

#-----------------
#L2 data
tm = 30000.0
tM = 50000.0
TrbA,KrbA,_,IrbA = kc.GetRBSP(pS.fRBa,pS.T0Str,tm,tM,"rbspa")
Tsc,Xa,Ya,Za = kc.getTraj(pS.fOrbA,pS.T0Str,tm,tM)

L = np.sqrt(Xa**2.0 + Ya**2.0)

i0 = L.argmax()
Tsc = Tsc[0:i0]
L = L[0:i0]

Nt = len(Tsc)
Nk = len(KrbA)

Irb = np.zeros((Nt,Nk))

for i in range(Nt):
	#Find matching time
	j = np.abs(TrbA-Tsc[i]).argmin()
	Irb[i,:] = IrbA[j,:]
genFig(L,KrbA,Irb,"RBPSD_L2.png")

#-----------------
#L3 data

L,K,Ilk = pS.GetRBPSD()
genFig(L,K,Ilk,"RBPSD_L3.png")

#----------------------
#Create interpolation
Sig = 1.1
Icrit = 1
Nl = 128
Nk = 256
IlkS = SmoothI(Ilk,sig=Sig)
LL,KK = np.meshgrid(L,K)
Li = np.linspace(1.5,10,Nl)
Ki = np.logspace(1,np.log10(5000),Nk)
LLi,KKi = np.meshgrid(Li,Ki)

Ii = scipy.interpolate.griddata((LL.flatten(),KK.flatten()),IlkS.T.flatten(),(LLi,KKi),method='linear',fill_value=0)
Ii = Ii.T

#Chop out floor values
Ii[Ii<=Icrit] = 0.0

#Push outward to fill domain
IkMax = Ii.max(axis=1)
l0 = IkMax.argmin()
Ii[l0-1:,:] = Ii[l0-2,:]

#Push down to fill energies
IkMax = Ii.max(axis=0)
k0 = (IkMax>0).argmax()
for i in range(Nl):
	Ii[i,0:k0] = Ii[i,k0]

genFig(Li,Ki,Ii,"RBPSD_Interp.png")

#----------------------
#Scale and convert to I->f
#I units from RB are j = cm-2 sr-1 s-1 keV-1
#Want f = I/p^2, w/ units (keV*s)-3
vc = 2.997e+10 #cm/s
mec2 = 0.511e+3 #keV

fLK = np.zeros((Nl,Nk))
for n in range(Nk):
	gamma = 1 + Ki[n]/mec2
	jScl = ((vc/mec2)**2.0)/(gamma*gamma-1)
	fLK[:,n] = jScl*Ii[:,n]
	#print("K/jScl = %f,%e"%(Ki[n],jScl))

#----------------------
#Write out
with h5py.File(fOut,'w') as hfw:
	hfw.create_dataset("Li" ,data=Li)
	hfw.create_dataset("Ki" ,data=Ki)
	hfw.create_dataset("fLK",data=fLK.T) #Transpose for Fortran

#-----------------
#Analytic approximation

B = -3.5
A = 4.54238496204e+26
en = 2.0
Ian = np.zeros((Nl,Nk))
for i in range(Nl):
	for j in range(Nk):
		L = Li[i]
		K = Ki[j]
		if (K >= 1000.0):
			L0 = 3.5
			Lm = 4.5
		elif (K >= 200.0):
			L0 = 4.75 -1.25*(K/1000.0)
			Lm = 6.37 - 1.88*(K/1000.0)
		else:
			L0 = 2.8 + 7.5*(K/1000.0)
			Lm = 6.0
		gamma = 1 + K/mec2
		fScl = (gamma*gamma-1)/((vc/mec2)**2.0)
		if (L <= L0):
			fIJ = 0.0
		else:
			Lscl = (L-L0)/(Lm-L0)
			eAA = (1 - Lscl**en)
			fExp = np.exp((2.0/en)*eAA)
			fLKscl = Lscl**2.0
			fIJ = fLKscl*fExp
		Ian[i,j] = fScl*A*(K**B)*fIJ

genFig(Li,Ki,Ian,"RBPSD_An.png")