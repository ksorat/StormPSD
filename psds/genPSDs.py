#Generate XMLs for PSD calculation


import xml.etree.ElementTree as et
import xml.dom.minidom

import numpy as np
import os

doTest = True

IDs = ["StormA","StormT","StormI"]

BaseP = "~/Work/StormPSD/Data/"
#MaskP = ["Trap/StormTrap.*.h5part","TrapC/StormTrapC.*.h5part","Inj/StormInj.*.h5part","Inj1/StormInj.*.h5part"]
MaskP = ["Inj0","Inj21","Inj3","Trap","TrapC"]
doTS = [True,True,True,False,False]
tsID = [0,21,3,None,None]

tsS = "tsWedge_"

#Uniform parameters
T0 = 33600.0
Tf = 197000.0
dt = 600.0
Rin = 2
Rout = 14
Nth = 16 #Number of threads
doLogR = True
if (doTest):
	dt = 6000.0
	#Tf = T0+10*dt

#Parameters

#Config 1
kappa = 3.4
kTScl = 0.5
Nr = 18
Np = 24
Nk = 25
kMin = 30.0
kMax = 4000.0
Na = 18

#Config 2
# kappa = 2.6
# kTScl = 0.5
# Rin = 2
# Rout = 16
# Nr = 28
# Np = 24
# Nk = 50
# kMin = 10.0
# kMax = 5000.0
# Na = 10


NumPSD = len(IDs)
NumPop = len(MaskP)

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
				tsF = tsS + str(tsID[p]) + ".csv"
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


#Generate runner
RunF = "RunPSD.sh"
with open(RunF,"w") as fID:
	fID.write("#!/bin/bash")
	fID.write("\n\n")
	fID.write("export OMP_NUM_THREADS=%d\n"%Nth)
	for i in range(NumPSD):
		IDi = IDs[i]
		ComS = "psd.x %s.xml\n"%IDi
		fID.write(ComS)
		ComS = "mv %s_r_phi_k_Slice3D#1.h5 KCyl_%s.h5\n"%(IDi,IDi)
		fID.write(ComS)
