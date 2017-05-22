#Generate XMLs for PSD calculation


import xml.etree.ElementTree as et
import xml.dom.minidom

import numpy as np
import os

IDs = ["StormA","StormT","StormI"]

BaseP = "~/Work/StormPSD/Data/"
MaskP = ["Trap/StormTrap.*.h5part","TrapC/StormTrapC.*.h5part","/Inj/StormInj.*.h5part","/Inj1/StormInj.*.h5part"]
doKap = [False,False,True,True]
tsF = "tsWedge.csv"
kappa = 2.8
kTScl = 0.5

T0 = 33600.0
Tf = 197000.0
dt = 3600.0
Rin = 2
Rout = 20
Nr = 48
Np = 24
Nk = 50
kMin = 10.0
kMax = 5000.0
Na = 20

Nth = 8

NumPSD = len(IDs)
NumPop = len(MaskP)

for i in range(NumPSD):

	#Create XML
	iDeck = et.Element('Params')

	#Particle/population details
	if (IDs[i] == IDs[0]):
		doPop = [True,True,True,True]
	elif (IDs[i] == IDs[1]):
		doPop = [True,True,False,False]
	else:
		doPop = [False,False,True,True]

	pInfo = et.SubElement(iDeck,"particles")
	pInfo.set("species","")
	pInfo.set("equatorial","T")
	for p in range(NumPop):
		if (doPop[p]):
			pID = "population"+str(p+1)
			popInfo = et.SubElement(pInfo,pID)
			popInfo.set("files",BaseP+MaskP[p])
			if (doKap[p]):
				popInfo.set("weighting","kapNTSeries")
				popInfo.set("kappa",str(kappa))
				popInfo.set("kTScl",str(kTScl))
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
