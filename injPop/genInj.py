#Routine to generate H5 particle blocks for storm run
#Injected population

import numpy as np
import glob
import os
import cPickle as pickle
import lfmPreproc as lfmpre
from sys import exit
#import matplotlib.pyplot as plt

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

#ipID = 0 #Which injected population: (0,3,21)
#Nswp = Sweep number, generates block 1+(Ns-1)*Nb,Ns*Nb
#Nswp = 1
ipIDs = [0,3,21]
Nswps = [5,6]


#Universal constants
#-------------------
#Root dir
rDir = os.path.expanduser('~') + "/Work/StormPSD/"
doWedgeOnly = False #Only create wedge TS
Kc = 10 #Cutoff for high-energy
T0Cut = 40000.0 #Cutoff for earliest injection time

#PBlock info / Sweep has blocks, blocks have particles

Npp = 20000 #Number of particles per block
Nb = 100 #Number of particle blocks

#Run info
T0 = 30000
Tfin = 197000
Ns = 150
Nf = 150
dtW = 150.0
doSmooth = True

#Alpha = [40,140]
Alpha = [1,179]
R = [9.0,12.0]
Z = [-0.25,0.25]
K0 = [10,2000]


for ipID in ipIDs:
	for Nswp in Nswps:
		print(ipID)
		print(Nswp)

		#Output pickle
		fPkl = "tsWedge_%d.pkl"%(ipID)

		if (ipID == 0):
			#Standard wedge
			#Wedge info
			P = [170,190]
		elif (ipID == 1):
			#Fat wedge w/ same injection data as ip0
			fPkl = "tsWedge_0.pkl"
			#Wedge info
			R = [10.0,11.0]
			P = [135,225]
		elif (ipID == 21):
			P = [125,145]
		elif (ipID == 3):
			P = [215,235]

		Np = Npp*Nb #Number of particles

		#Some derived quantities
		rSeed = 31337+Nswp
		np.random.seed(rSeed)
		id0 = 1+(Nswp-1)*Np #First particle ID
		pb0 = 1+(Nswp-1)*Nb #First block ID

		print("Generating injected popuation")
		print("\tWedge %d"%(ipID))
		print("\tParticles %d to %d"%(id0,id0+Np-1))
		print("\tBlocks %d to %d"%(pb0,pb0+Nb-1))
		print("\tUsing random # seed %d"%rSeed)

		#HDF directory/H5 IC data/Input decks/Output
		oTag = "StormInj_%d"%ipID
		lfmDir = rDir + "lfmData"
		ipDir = rDir + "Data/ip%d/"%ipID
		inpDir = rDir + "Inps/"
		outDir = rDir + "Data/Inj%d/"%ipID

		#Create directories if necessary
		if not os.path.exists(ipDir):
			os.makedirs(ipDir)
		if not os.path.exists(outDir):
			os.makedirs(outDir)

		#Generate wedge time series if necessary, otherwise load
		if (os.path.isfile(fPkl)):
			print("Loading Wedge TS %s"%(fPkl))
			with open(fPkl,"rb") as f:
				t = pickle.load(f)
				Vst = pickle.load(f)
				kTt = pickle.load(f)
				Nt = pickle.load(f)
				Nkt = pickle.load(f)	

		else:
			print("Generating Wedge TS")
			fIns =glob.glob(lfmDir + "/*.hdf")

			t,Vst,kTt,Nt,Nkt = lfmpre.tsWedge(fIns,R,P,Z)

			#Save data
			print("Saving data to %s"%fPkl)
			with open(fPkl,"wb") as f:
				pickle.dump(t,f)
				pickle.dump(Vst,f)
				pickle.dump(kTt,f)
				pickle.dump(Nt,f)
				pickle.dump(Nkt,f)

		if (doWedgeOnly):
			exit()
		#Done wedge data
		#----------------------------

		#----------------------------
		#Generate H5s

		#Start by calculating injection times
		I = (t<T0Cut)
		Vst[I] = 0.0
		fTotTmp = Vst*Nkt #High-energy number flux in wedge
		if (doSmooth):
			fTot = tWindow(t,fTotTmp,dtW)
		else:
			fTot = fTotTmp

		T0p = lfmpre.genSample(t,fTot,Np)
		T0p = np.maximum(T0p,T0)

		#Create particle block
		pBlock = lfmpre.ParticleBlock(0,T0,Np,id0)

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
		lfmpre.writeBlocks(ipDir + oTag+".",pList,pb0=pb0)

		#--------------------------
		#Generate input decks
		for n in range(Nb):
			nId = n+pb0 #Handle offset
			resFile = ipDir + "%s.%04d.h5part"%(oTag,nId)
			ideckFile = "%s%s.%04d.xml"%(inpDir,oTag,nId)
			iDeck = lfmpre.genDeck(spc="e",tagStr=oTag,outDir=outDir)
			iDeck = lfmpre.itDeck(iDeck,T0=T0,Tf=Tfin,dt=1.0,dtS=Ns,dtF=Nf,iMeth="dynamic",kMin="1.0")
			iDeck = lfmpre.resDeck(iDeck,resFile,outTag=oTag,doApp=False,Ts=T0)
			iDeck = lfmpre.streamDeck(iDeck)
			lfmpre.writeDeck(iDeck,fOut=ideckFile)
