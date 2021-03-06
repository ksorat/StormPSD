#Generate XMLs for PSD calculation
import xml.etree.ElementTree as et
import xml.dom.minidom
import cPickle as pickle
import numpy as np
import os

def tWindow(t,Q,dt,SigQ=0.25):
	import scipy.ndimage
	from scipy.ndimage.filters import gaussian_filter1d
	#First gaussian window to denoise
	#Q = gaussian_filter1d(Q,sigma=SigQ)
	#Window time series t,Q based on window size dt
	Nt = len(t)
	Qw =  np.zeros(Nt)
	for i in range(Nt):
		t0 = t[i]
		I = (np.abs(t-t0) <= dt)
		Qw[i] = Q[I].mean()
	return Qw

doTest = True
doSmoothTS = True

aIDs = ["StormA","StormT","StormI"]

BaseP = "~/Work/StormPSD/Data/"
MaskP = ["Inj0","Inj21","Inj3","Trap","TrapC"]
doTS = [True,True,True,False,False]
tsID = ["0","21","3",None,None]

tsS = "tsWedge_"
dR_W = 3 #Wedge radial length [Re]
ReKM = 6.38e+3

#Uniform parameters
T0 = 33600.0
Tf = 197000.0
dt = 600.0
Rin = 2
Rout = 18
Nth = 16 #Number of threads
doLogR = True
if (doTest):
	dt = 3000.0
	#Tf = T0+10*dt

#Parameters

#Config 2
#kappa = 3.5
#kTScl = 0.25
kappas = [3.0,3.5,4.0]
kappaStr = ["A","B","C"]

kTScls = [0.15,0.25,0.4]
kTSclStr = ["a","b","c"]

Nkap = len(kappas)
Nkt = len(kTScls)

Nr = 60
Np = 48
Nk = 20
kMin = 30.0
kMax = 4100.0
Na = 12


NumPSD = len(aIDs)
NumPop = len(MaskP)

#Start by creating CSV time series files from PKLs
#Cheating by scaling number density by (dt*Ux/dR) to introduce dN particles per step
for n in range(NumPop):
	if (doTS[n]):
		fPkl = "tsWedge_%s.pkl"%(tsID[n])
		fTab = "tsWedge_%s.csv"%(tsID[n])
		print("Creating CSV for injection %s w/ dt=%f"%(tsID[n],dt))
		with open(fPkl,"rb") as f:
			t = pickle.load(f) #Time [s]
			Vst = pickle.load(f) #Earthward tail velocity [km/s]
			kTt = pickle.load(f) #Thermal energy, kT [keV]
			Nt  = pickle.load(f) #Number density, [#/cm3]

		N = t.shape[0]
		dOut = np.zeros((3,N))
		dOut[0,:] = t

		wVst = tWindow(t,Vst,dt)
		wkTt = tWindow(t,kTt,dt)
		wNt  = tWindow(t,Nt ,dt)
		if (doSmoothTS):
			nScl = (dt*wVst)/(dR_W*ReKM)
			dOut[1,:] = nScl*wNt
			dOut[2,:] = wkTt
		else:
			nScl = (dt*Vst)/(dR_W*ReKM)
			dOut[1,:] = nScl*Nt
			dOut[2,:] = kTt
		print("\tMin/Mean/Max nScl = %f,%f,%f"%(nScl.min(),nScl.mean(),nScl.max()))
		#print(nScl)
		np.savetxt(fTab,dOut.T,delimiter=',')

for nk in range(Nkap):
	kappa = kappas[nk]
	for nkt in range(Nkt):
		kTScl = kTScls[nkt]
		IDs = [s + "_" + kappaStr[nk] + kTSclStr[nkt] for s in aIDs]

		for i in range(NumPSD):

			#Create XML
			iDeck = et.Element('Params')

			#Particle/population details
			if (IDs[i] == IDs[0]):
				#All
				doPop = [True,True,True,True,True]
			elif (IDs[i] == IDs[1]):
				#Trapped
				doPop = [False,False,False,True,True]
			else:
				#Injected
				doPop = [True,True,True,False,False]

			pInfo = et.SubElement(iDeck,"particles")
			pInfo.set("species","")
			pInfo.set("equatorial","T")
			pn = 0
			for p in range(NumPop):
				if (doPop[p]):
					pID = "population"+str(pn+1)
					pn = pn + 1
					popInfo = et.SubElement(pInfo,pID)
					popInfo.set("files",BaseP+MaskP[p]+"/*.h5part")
					if (doTS[p]):
						popInfo.set("weighting","kapNTSeries")
						popInfo.set("kappa",str(kappa))
						popInfo.set("kTScl",str(kTScl))
						tsF = tsS + tsID[p] + ".csv"
						popInfo.set("seriesfile",tsF)
					else:
						popInfo.set("weighting","rbtrap")
			#Timing data
			tInfo = et.SubElement(iDeck,"timing")
			tInfo.set("weighting",str(T0))
			cStr = "%s:%s:%s"%(str(T0),str(dt),str(Tf))
			tInfo.set("calculation",cStr)

			#Options
			oInfo = et.SubElement(iDeck,"options")
			oInfo.set("relativistic","T")
			oInfo.set("doFlux","T")

			#Phasespace
			psInfo = et.SubElement(iDeck,"phasespace")
			psR = et.SubElement(psInfo,"r")
			psR.set("N",str(Nr))
			psR.set("min",str(Rin))
			psR.set("max",str(Rout))
			if (doLogR):
				psR.set("log","T")
			psT = et.SubElement(psInfo,"theta")
			psT.set("N","1")
			psP = et.SubElement(psInfo,"phi")
			psP.set("N",str(Np))
			psK = et.SubElement(psInfo,"k")
			psK.set("N",str(Nk))
			psK.set("min",str(kMin))
			psK.set("max",str(kMax))
			psK.set("log","T")
			psA = et.SubElement(psInfo,"alpha")
			psA.set("N",str(Na))
			psS = et.SubElement(psInfo,"psi")
			psS.set("N","1")

			#Output
			ioInfo = et.SubElement(iDeck,"output")
			ioInfo.set("base",IDs[i])
			ioInfo.set("fullevery","0")
			slc = et.SubElement(ioInfo,"slice3D1")
			slc.set("dim1","r")
			slc.set("dim2","phi")
			slc.set("dim3","k")
			slcT = et.SubElement(slc,"theta")
			slcT.set("ignore","T")

			#Write out
			fOut = IDs[i]+".xml"
			#Finished creating XML tree, now write
			xmlStr = xml.dom.minidom.parseString(et.tostring(iDeck)).toprettyxml(indent="    ")
			with open(fOut,"w") as f:
				f.write(xmlStr)


#Generate runners
#Individual PSD
for nk in range(Nkap):
	for nkt in range(Nkt):
		IDs = [s + "_" + kappaStr[nk] + kTSclStr[nkt] for s in aIDs]

		for i in range(NumPSD):
			IDi = IDs[i]
			RunP = "RunPSD_%s.sh"%(IDi)
			with open(RunP,"w") as fID:
				fID.write("#!/bin/bash")
				fID.write("\n\n")
				fID.write("module restore lfmtp\n")
				fID.write("module list\n")
				fID.write("export OMP_NUM_THREADS=%d\n"%Nth)
				ComS = "psd.x %s.xml\n"%IDi
				fID.write(ComS)
				ComS = "mv %s_r_phi_k_Slice3D#1.h5 KCyl_%s.h5\n"%(IDi,IDi)
				fID.write(ComS)
			os.chmod(RunP, 0744)

#Sub all PSDs
RunF = "SubPSDs.sh"
wcS = "12:00"
qS = "regular"
pS = "UJHB0003"
with open(RunF,"w") as fID:
	fID.write("#!/bin/bash")
	fID.write("\n\n")
	for nk in range(Nkap):
		for nkt in range(Nkt):
			IDs = [s + "_" + kappaStr[nk] + kTSclStr[nkt] for s in aIDs]
			for i in range(NumPSD):
				IDi = IDs[i]
				logF = "PSD_%s.log"%(IDi)
				jobS = "PSD%s"%(IDi)
				ComS = "bsub -a poe -P " + pS + " -W " + wcS + " -n 1 -q " + qS + " -J " + jobS + " -e %s -o %s"%(logF,logF)
				ComS = ComS + " \"RunPSD_%s.sh\" "%(IDi) + "\n"
				fID.write(ComS)

os.chmod(RunF, 0744)

