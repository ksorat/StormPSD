import sys
import os
import numpy as np
import datetime
from visit import *
from visit_utils import *
from visit_utils.common import lsearch #lsearch(dir(),"blah")
import pyVisit as pyv

Quiet = True

titS = "St. Patrick's Storm 2013"
outVid ="StPattys.mp4"

T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"

vidScl = 2 #>1 to slow down
#Field
fBds = [0.1,250]
fMap = "viridis"

Base = os.path.expanduser('~') + "/Work/StormPSD/Data"

dbSlc = Base+"/eqSlc/eqSlc.*.vti database"
dbLn  = Base+"/blines/blines.*.vtp database"
dbPI = Base + "/H5p/StormInj.Min3D.h5part"
dbPT = Base + "/H5p/StormTrap.Min3D.h5part"
#dbPT = Base + "/Trap/StormTrap.0001.h5part"
#dbPI = Base + "/Inj/StormInj.0001.h5part"

dbs = [dbSlc,dbLn,dbPI,dbPT]

pMapT = "Cool"
pMapI = "Cool"
pSzI = 4
pSzT = 4

#pBds = [4,6]
pVar = "kev"
pLab = "Particle Energy [keV]"
pBdsI = [-5,0]
pBdsT = [10000,50000]

if (Quiet):
	LaunchNowin()
else:
	Launch()

#Do some defaults
pyv.lfmExprs()
pyv.pvInit()
DeleteExpression("Bmag")

DefineScalarExpression("slcMag","sqrt(Bx*Bx+By*By+Bz*Bz)")
DefineScalarExpression("pZero","kev*0.0")
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
OpenDatabase(dbs[3])

#Create database correlation
CreateDatabaseCorrelation("SLP",dbs,0)

#Plot field lines
ActivateDatabase(dbs[1])
AddPlot("Pseudocolor","Bmag")
pcOp = GetPlotOptions()
pcOp.legendFlag = 0
SetPlotOptions(pcOp)


#Plot equatorial slice
ActivateDatabase(dbs[0])
pyv.lfmPCol(dbs[0],"slcMag",vBds=fBds,pcOpac=1.0,Inv=False,Log=True,cMap=fMap,Legend=False)
pyv.chopInner2D()
pyv.to3D(opNum=1)

#Plot particles
#Injected
ActivateDatabase(dbs[2])
pyv.lfmPScat(dbs[2],v3="pZero",v4=pVar,vBds=pBdsI,cMap=pMapI,Log=False,Inv=False,pSize=pSzI,Legend=False)
pyv.onlyIn()

#Trapped
ActivateDatabase(dbs[3])
pyv.lfmPScat(dbs[3],v3="pZero",v4=pVar,vBds=pBdsT,cMap=pMapT,Log=False,Inv=False,pSize=pSzT,Legend=False)
pyv.onlyIn()

#Cleanup
v3d = GetView3D()
v3d.viewNormal = (0.883195, -0.0484622, 0.466495)
v3d.focus = (-1.04995, -2.09808e-05, 0.0199747)
v3d.viewUp = (-0.467268, -0.00537966, 0.884099)
v3d.viewAngle = 30
v3d.parallelScale = 24.5542
v3d.nearPlane = -49.1085
v3d.farPlane = 49.1085
v3d.imagePan = (0, 0)
v3d.imageZoom = 1.21
v3d.perspective = 1
v3d.eyeAngle = 2
v3d.centerOfRotationSet = 0
v3d.centerOfRotation = (-1.04995, -2.09808e-05, 0.0199747)
v3d.axis3DScaleFlag = 0
v3d.axis3DScales = (1, 1, 1)
v3d.shear = (0, 0, 1)
v3d.windowValid = 1
SetView3D(v3d)

#Cleanup

pyv.killAnnotations()

pyv.setAtts()


DrawPlots()
pyv.doTimeLoop(T0=T0,dt=dt,Save=True,tLabPos=(0.45,0.25),Trim=True)
pyv.makeVid(Clean=True,outVid=outVid,tScl=vidScl)
DeleteAllPlots()

