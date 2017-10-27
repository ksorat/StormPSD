#Generate XMLs for PSD calculation
import xml.etree.ElementTree as et
import xml.dom.minidom
import cPickle as pickle
import numpy as np
import os
import pyStorm as pS



doTest = False
doSmoothTS = True

#All,All trapped, injected wedges, all injected
IDs = ["StormA","StormT","StormI_0","StormI_21","StormI_3","StormI"]
doID = [False,True,True,True,True,False]

BaseP = "~/Work/StormPSD/Data/"
MaskP = ["Inj0","Inj21","Inj3","Trap"]
doTS = [True,True,True,False]
tsID = ["0","21","3",None]

tsS = "tsWedge_"
dR_W = 3 #Wedge radial length [Re]
ReKM = 6.38e+3

#Uniform parameters
T0 = 33600.0
Tf = 197000.0
dt = 150.0
Rin = 2.05
Rout = 18
Nth = 36 #Number of threads

doLogR = True
if (doTest):
	dt = 5*dt

dtW = dt
dtFull = 0 #Num of timesteps between full output

#Parameters

#Config 1

kappa = 3.5
kTScl = 0.25
#Nr = 60
#Using Nr=30 for better aspect ratio
Nr = 30
Np = 48
Nk = 20
kMin = 30.0
kMax = 4100.0
Na = 12


NumPSD = len(IDs)
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
			Nkt = pickle.load(f) #High-energy number density

		N = t.shape[0]
		dOut = np.zeros((3,N))
		dOut[0,:] = t

		wVst = pS.twWindow(t,Vst,dtW)
		wkTt = pS.twWindow(t,kTt,dtW)
		wNt  = pS.twWindow(t,Nt ,dtW)
		wNkt = pS.twWindow(t,Nkt,dtW)

		if (doSmoothTS):
			print("Using smoothed TS with dtW = %f"%(dtW))
			nScl = (dt*wVst)/(dR_W*ReKM)
			dOut[1,:] = nScl*wNt
			dOut[2,:] = wkTt
			Ntot = np.sum(dtW*nScl*wNt)
		else:
			nScl = (dt*Vst)/(dR_W*ReKM)
			dOut[1,:] = nScl*Nt
			dOut[2,:] = kTt
			Ntot = np.sum(dtW*nScl*Nt)

		print("\tMin/Mean/Max nScl = %f,%f,%f"%(nScl.min(),nScl.mean(),nScl.max()))
		print("\tTotal particles = %f"%Ntot)
		#print(nScl)
		np.savetxt(fTab,dOut.T,delimiter=',')

for i in range(NumPSD):
	if (not doID[i]):
		print("Skipping %s"%(IDs[i]))
		continue
	print("Creating %s"%(IDs[i]))

	#Create XML
	iDeck = et.Element('Params')

	#Particle/population details
	if (IDs[i] == IDs[0]):
		#All
		doPop = [True,True,True,True]
	elif (IDs[i] == IDs[1]):
		#Trapped
		doPop = [False,False,False,True]
	elif (IDs[i] == "StormI_0"):
		#0 wedge
		doPop = [True,False,False,False]
	elif (IDs[i] == "StormI_21"):
		#21 wedge
		doPop = [False,True,False,False]
	elif (IDs[i] == "StormI_3"):
		#3 wedge
		doPop = [False,False,True,False]
	else:
		#All Injected
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
				#popInfo.set("weighting","rbtrap")
				popInfo.set("weighting","rbinterp")
				popInfo.set("infile","rbInterp.h5")
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
	ioInfo.set("fullevery",str(dtFull))
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
#pS = "P28100045"
pS = "UJHB0010"

#Individual PSD
for i in range(NumPSD):
	if (doID[i]):
		IDi = IDs[i]
		RunP = "RunPSD_%s.pbs"%(IDi)
		with open(RunP,"w") as fID:
			fID.write("#!/bin/bash\n\n")
			fID.write("#PBS -A %s\n"%(pS))
			fID.write("#PBS -N %s\n"%(IDi))
			fID.write("#PBS -j oe\n")
			fID.write("#PBS -q regular\n")
			fID.write("#PBS -l walltime=12:00:00\n")
			fID.write("#PBS -l select=1:ncpus=72:ompthreads=72\n")
			fID.write("\n\n")
			fID.write("source ~/.bashrc\n")
			fID.write("module restore lfmtp\n")
			fID.write("module list\n")
			fID.write("hostname\n")
			fID.write("date\n")
			fID.write("export OMP_NUM_THREADS=%d\n"%Nth)
			fID.write("export KMP_STACKSIZE=128M\n")
			ComS = "omplace -nt $OMP_NUM_THREADS psd.x %s.xml > %s.out\n"%(IDi,IDi)
			#ComS = "psd.x %s.xml > %s.out\n"%(IDi,IDi)
			fID.write(ComS)
			ComS = "mv %s_r_phi_k_Slice3D#1.h5 KCyl_%s.h5\n"%(IDi,IDi)
			fID.write(ComS)
		os.chmod(RunP, 0744)
	
#Sub all PSDs
# RunF = "SubPSDs.sh"
# #wcS = "12:00"
# #qS = "regular"
# wcS = "24:00"
# qS = "geyser"



# with open(RunF,"w") as fID:
# 	fID.write("#!/bin/bash")
# 	fID.write("\n\n")

# 	for i in range(NumPSD):
# 		if (doID[i]):
# 			IDi = IDs[i]
# 			logF = "PSD_%s.log"%(IDi)
# 			jobS = "PSD%s"%(IDi)
# 			ComS = "bsub -a poe -P " + pS + " -W " + wcS + " -n 1 -q " + qS + " -J " + jobS + " -e %s -o %s"%(logF,logF)
# 			ComS = ComS + " \"RunPSD_%s.sh\" "%(IDi) + "\n"
# 			fID.write(ComS)
# os.chmod(RunF, 0744)

