import sys
import os
import numpy as np
import datetime
from visit import *
from visit_utils import *
from visit_utils.common import lsearch #lsearch(dir(),"blah")
import pyVisit as pyv

titS = "St. Patrick's Storm"
outVid ="StPattys.mp4"

#Field
fBds = [1,250]
fMap = "viridis"
#Particles
pBds = [100,2000]
pMap = "Reds"
pSz = 5

Base = os.path.expanduser('~') + "/Work/StormPSD/Data"
pFile = "Storm.Min3D.h5part"
Quiet = True


dbP = Base + "/H5p/" + pFile
dbF = Base + "/eqSlc/eqSlc.*.vti database"
dbs = [dbF, dbP]

if (Quiet):
	LaunchNowin()
else:
	Launch()

DefineScalarExpression("ev","max(1000*kev,1.0e-8)")
DefineScalarExpression("logev","log10(ev)")

DefineScalarExpression("isIn","in")

#Do some defaults
pyv.lfmExprs()
pyv.pvInit()

md0 = GetMetaData(dbs[0])
mdH5p = GetMetaData(dbs[1])

dt = md0.times[1] - md0.times[0]
T0 = md0.times[0] - md0.times[0] #Reset to zero


OpenDatabase(dbs[0])
OpenDatabase(dbs[1])

#Create database correlation
CreateDatabaseCorrelation("P2Fld",dbs,0)
	
#Create fields/particle plots
ActivateDatabase(dbs[0])
pyv.lfmPCol(dbs[0],"Bmag",vBds=fBds,pcOpac=0.7,Inv=False,Log=True,cMap=fMap)
pyv.chopInner2D()

ActivateDatabase(dbs[1])
pyv.lfmPScat(dbs[1],v4="kev",vBds=pBds,cMap=pMap,Log=False,Inv=False,pSize=pSz)
pyv.onlyIn()


#Gussy things up
tit = pyv.genTit( titS=titS)
plXs = [0.03,0.03]
plYs = [0.9,0.4]
plTits = ["Field Strength [nT]", "Particle Energy [keV]"]
pyv.cleanLegends(plXs,plYs,plTits)
pyv.setAtts()

#Let's see what we got
DrawPlots()

# #Do time loop
pyv.doTimeLoop(T0=T0,dt=dt,Save=True,tLabPos=(0.3,0.05),Trim=True)

pyv.makeVid(Clean=True,outVid=outVid,tScl=1)
DeleteAllPlots()
CloseDatabase(dbs[0])
CloseDatabase(dbs[1])
pyv.killAnnotations()
os.system("mkdir tmpVid")

