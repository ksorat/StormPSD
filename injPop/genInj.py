#Routine to generate H5 particle blocks for storm run
#Injected population

import numpy as np
import glob
import os
import cPickle as pickle
import lfmPreproc as lfmpre

#Root dir
rDir = os.path.expanduser('~') + "/Work/StormPSD/"
#Output pickle
fPkl = "tsWedge.pkl"

#Wedge info
R = [11.5,12.5]
Z = [-0.25,0.25]
P = [165,195]
Kc = 10 #Cutoff for high-energy

#Run info
T0 = 30000
Tfin = 197000
Ns = 60
Nf = 600
Alpha = [40,140]
K0 = [10,1000]

#PBlock info
Npp = 2500 #Number of particles per block
Nb = 10 #Number of particle blocks
Np = Npp*Nb #Number of particles

#HDF directory/H5 IC data/Input decks/Output
oTag = "StormInj"
lfmDir = rDir + "lfmData"
ipDir = rDir + "Data/ip0/"
inpDir = rDir + "Inps/"
outDir = rDir + "Data/Inj"

#Generate wedge time series if necessary, otherwise load
if (os.path.isfile(fPkl)):
	print("Loading Wedge TS")
	with open(fPkl,"rb") as f:
		t = pickle.load(f)
		Vst = pickle.load(f)
		kTt = pickle.load(f)
		Nt = pickle.load(f)
		Nkt = pickle.load(f)	
else:
	print("Generating Wedge TS")
	fIns =glob.glob(lfmDir + "/*.hdf")

	Nf = len(fIns)

	t = np.zeros(Nf)
	Vst = np.zeros(Nf)
	kTt = np.zeros(Nf)
	Nt = np.zeros(Nf)
	Nkt = np.zeros(Nf)
	#Volume info for wedge
	I,dvI = lfmpre.lfmWedge(fIns[0],R=R,P=P,Z=Z)
	#Loop through files, get wedge data
	for i in range(Nf):
		print("Reading file %d of %d\n"%(i,Nf))
		fIn = fIns[i]
		t[i],Vst[i],kTt[i],Nt[i],Nkt[i] = lfmpre.injWedge(fIn,I,dvI,K0=Kc,T0=T0)
	#Got all data, now sort
	Is = np.argsort(t)
	
	t = t[Is]
	Vst = Vst[Is]
	kTt = kTt[Is]
	Nt = Nt[Is]
	Nkt = Nkt[Is]
	#Save data
	print("Saving data to %s"%fOut)
	with open(fPkl,"wb") as f:
		pickle.dump(t,f)
		pickle.dump(Vst,f)
		pickle.dump(kTt,f)
		pickle.dump(Nt,f)
		pickle.dump(Nkt,f)


#Done wedge data
#----------------------------

#----------------------------
#Generate H5s

#Start by calculating injection times
fTot = Vst*Nkt #High-energy number flux in wedge
T0p = lfmpre.genSample(t,fTot,Np)
T0p = np.maximum(T0p,T0)

#Create particle block
pBlock = lfmpre.ParticleBlock(0,T0,Np)

Rp = lfmpre.genLinR(R[0],R[1],Np)
Pp = (np.pi/180.0)*lfmpre.genLinR(P[0],P[1],Np)

#Leaving time-dep quantities undefined (V,Mu)
Xp = Rp*np.cos(Pp)
Yp = Rp*np.sin(Pp)

pBlock.x = Xp
pBlock.y = Yp
pBlock.z[:] = 0.0
pBlock.kev = lfmpre.genLogR(K0[0],K0[1],Np)
pBlock.alpha = lfmpre.genLinR(Alpha[0],Alpha[1],Np)
pBlock.T0p = T0p
pBlock.xeq = pBlock.x
pBlock.yeq = pBlock.y
pBlock.Teq[:] = T0
pBlock.alpheq = pBlock.alpha
pBlock.keveq = pBlock.kev
pBlock.GC[:] = False #Start off as FP, calculate Mu

#Split into blocks
pList = lfmpre.splitBlock(pBlock,Nb)

#Write out blocks
lfmpre.writeBlocks(ipDir + oTag+".",pList)

#--------------------------
#Generate input decks
for n in range(Nb):
	nId = n+1 #Handle offset
	resFile = ipDir + "%s.%04d.h5part"%(oTag,nId)
	ideckFile = "%s%s.%04d.xml"%(inpDir,oTag,nId)
	iDeck = lfmpre.genDeck(spc="e",tagStr=oTag,outDir=outDir)
	iDeck = lfmpre.itDeck(iDeck,T0=T0,Tf=Tfin,dt=1.0,dtS=60,dtF=600,iMeth="dynamic")
	iDeck = lfmpre.resDeck(iDeck,resFile,outTag="StormInj",doApp=False,Ts=T0)
	iDeck = lfmpre.streamDeck(iDeck)
	lfmpre.writeDeck(iDeck,fOut=ideckFile)
