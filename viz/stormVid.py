import sys
import os
import numpy as np
import datetime
from visit import *
from visit_utils import *
from visit_utils.common import lsearch #lsearch(dir(),"blah")
import pyVisit as pyv

titS = "St. Patrick's Storm 2013"
outVid ="StPattys.mp4"

T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"

vidScl = 4 #>1 to slow down
#Field
fBds = [1,250]
fMap = "viridis"

#Particles (T)rap, (I)injected
#Using files Base/(Trap/Inj)/Storm(Trap/Inj).h5part

pMapT = "Reds"
pMapI = "Blues"

pSzI = 8
pSzT = 4
#pBds = [4,6]
pVar = "kev"
pLab = "Particle Energy [keV]"
pBds = [10,1000]

Base = os.path.expanduser('~') + "/Work/StormPSD/Data"

Quiet = True

dbPI = Base + "/H5p/StormInj.Min3D.h5part"
dbPT = Base + "/H5p/StormTrap.Min3D.h5part"
dbF = Base + "/eqSlc/eqSlc.*.vti database"
dbs = [dbF,dbPT,dbPI]

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

#Get time info
#dT0s = 2000 #Start time [s] from T0
dT0s = md0.times[0]
dT0 = datetime.timedelta(seconds=dT0s)
T0 = dT0 + datetime.datetime.strptime(T0Str,T0Fmt)


OpenDatabase(dbs[0])
OpenDatabase(dbs[1])
OpenDatabase(dbs[2])

#Create database correlation
CreateDatabaseCorrelation("P2Fld",dbs,0)
	
#Create fields/particle plots
ActivateDatabase(dbs[0])
pyv.lfmPCol(dbs[0],"Bmag",vBds=fBds,pcOpac=0.7,Inv=False,Log=True,cMap=fMap)
pyv.chopInner2D()

ActivateDatabase(dbs[1])
pyv.lfmPScat(dbs[1],v4=pVar,vBds=pBds,cMap=pMapT,Log=False,Inv=False,pSize=pSzT)
pyv.onlyIn()

ActivateDatabase(dbs[2])
pyv.lfmPScat(dbs[2],v4=pVar,vBds=pBds,cMap=pMapI,Log=False,Inv=False,pSize=pSzI)
pyv.onlyIn()


#Gussy things up
tit = pyv.genTit( titS=titS)
plXs = [0.03,0.05,0.01]
plYs = [0.9,0.4,0.4]

plTits = ["Field Strength [nT]", "",pLab]

pyv.cleanLegends(plXs,plYs,plTits)
pyv.setAtts()

AnLab = 'Plot0002'
GetAnnotationObject(AnLab).drawTitle = 0
GetAnnotationObject(AnLab).drawLabels = 0

#Let's see what we got
DrawPlots()

#Do time loop
pyv.doTimeLoop(T0=T0,dt=dt,Save=True,tLabPos=(0.3,0.025),Trim=True)
pyv.makeVid(Clean=True,outVid=outVid,tScl=vidScl)
DeleteAllPlots()
CloseDatabase(dbs[0])
CloseDatabase(dbs[1])
pyv.killAnnotations()
os.system("mkdir tmpVid")

