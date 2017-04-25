import sys
import os
import numpy as np
import datetime
from visit import *
from visit_utils import *
from visit_utils.common import lsearch #lsearch(dir(),"blah")
import pyVisit as pyv

def AddMol(db,cMap="Cool",vBds=[0,1],rScl=1.0e-4,doLeg=True):
	ActivateDatabase(db)
	AddPlot("Molecule","kev")
	mOp = GetPlotOptions()
	mOp.drawBondsAs = 0
	mOp.atomSphereQuality = 2
	mOp.continuousColorTable = cMap
	mOp.minFlag = 1
	mOp.maxFlag = 1
	mOp.scalarMin = vBds[0]
	mOp.scalarMax = vBds[1]
	mOp.scaleRadiusBy = 3
	mOp.radiusVariable = "kevRad"
	mOp.radiusScaleFactor = rScl
	if (doLeg):
		mOp.legendFlag = 1
	else:
		mOp.legendFlag = 0

	SetPlotOptions(mOp)

	#Move from x,y,z->xeq,yeq,0
	AddOperator("Displace",0)
	dOp = GetOperatorOptions(0)
	dOp.variable = "dR"
	SetOperatorOptions(dOp)

def AddLine(x0=(0,0,0),x1=(1,1,1),lw=2):
	axL = CreateAnnotationObject("Line3D")
	axL.point1 = tuple(x0)
	axL.point2 = tuple(x1)
	axL.width = lw

def SetWin(W,H):
	#Set view
	ResizeWindow(1,W,H)

	# v3d = GetView3D()
	# v3d.viewNormal = (0.525379, -0.801214, 0.286414)
	# v3d.focus = (-1.04995, -2.09808e-05, 0.0199747)
	# v3d.viewUp = (-0.183865, 0.221756, 0.957611)
	# v3d.viewAngle = 30
	# v3d.parallelScale = 24.5542
	# v3d.nearPlane = -49.1085
	# v3d.farPlane = 49.1085
	# v3d.imagePan = (0, 0)
	# v3d.imageZoom = 2.14359
	# v3d.perspective = 1
	# v3d.eyeAngle = 2
	# v3d.centerOfRotationSet = 0
	# v3d.centerOfRotation = (-1.04995, -2.09808e-05, 0.0199747)
	# v3d.axis3DScaleFlag = 0
	# v3d.axis3DScales = (1, 1, 1)
	# v3d.shear = (0, 0, 1)
	# v3d.windowValid = 1
	# SetView3D(v3d)

	v3d = GetView3D()
	v3d.viewNormal = (0.431009, -0.865563, 0.255014)
	v3d.focus = (-1.04995, -2.09808e-05, 0.0199747)
	v3d.viewUp = (-0.117092, 0.226573, 0.96693)
	v3d.viewAngle = 30
	v3d.parallelScale = 24.5542
	v3d.nearPlane = -49.1085
	v3d.farPlane = 49.1085
	v3d.imagePan = (0.0144521, 0.0167345)
	v3d.imageZoom = 3.13843
	v3d.perspective = 1
	v3d.eyeAngle = 2
	v3d.centerOfRotationSet = 0
	v3d.centerOfRotation = (-1.04995, -2.09808e-05, 0.0199747)
	v3d.axis3DScaleFlag = 0
	v3d.axis3DScales = (1, 1, 1)
	v3d.shear = (0, 0, 1)
	v3d.windowValid = 1
	SetView3D(v3d)

H = 900
W = 1800

Quiet = True
Prod  = True
rIn = 2.0
rOpac = 150 #[0,255]

Nsk = 1

doScat = False #Scatter/molecule
doTwoP = False



titS = "St. Patrick's Storm 2013"
outVid ="fldP.mp4"

T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"

vidScl = 2 #>1 to slow down
#Field
#Total
# fBds = [0.1,250]
# fMap = "viridis"
# fVar = "sclMag"
# fInv = False
# fLog = True
# fLab = "Field Strength [nT]"

#Residual
fBds = [-25,25]
fMap = "RdGy"
fVar = "dBz"
fInv = True
fLog = False
fLab = "dBz [nT]"

Base = os.path.expanduser('~') + "/Work/StormPSD/Data"

dbSlc = Base+"/eqSlc/eqSlc.*.vti database"
dbLn  = Base+"/blines/blines.*.vtp database"

dbPI = Base + "/H5p/StormInj.Min3D.h5part"
dbPT = Base + "/H5p/StormTrap.Min3D.h5part"

dbs = [dbSlc,dbLn,dbPI,dbPT]


#pBds = [4,6]
pVar = "kev"
pLab = "Particle Energy [keV]"
kMax = 1000
kMaxR = 1500
kMinR = 250

#rScl = 5.0e-5
rScl = 7.5e-5
pBdsI = [0,kMax]
pBdsT = [0,kMax]
pMapT = "Reds"
#pMapI = "Cool"
pMapT = "Winter"
pSzI = 4
pSzT = 4

if (doTwoP):
	doLegP = True
else:
	doLegP = False
	pMapT = pMapI

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

DefineScalarExpression("fL","Bmag/Bmag - 0.5")
DefineVectorExpression("dR","{-x+xeq,-y+yeq,-z}")
DefineScalarExpression("kevRad","max(min(%f,kev),%f)"%(kMaxR,kMinR))

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
AddPlot("Pseudocolor","fL")
pcOp = GetPlotOptions()
pcOp.legendFlag = 0
pcOp.colorTableName = "xray"
pcOp.minFlag = 1
pcOp.maxFlag = 1
pcOp.min = 0
pcOp.max = 1
# pcOp.lineType = 1
# pcOp.tubeRadiusBBox = 0.005
# pcOp.opacityType = 2
# pcOp.opacity = 0.5
SetPlotOptions(pcOp)
AddOperator("Tube")
tOp = GetOperatorOptions(0)
tOp.radiusFractionBBox = 0.005
SetOperatorOptions(tOp)


#Plot equatorial slice
ActivateDatabase(dbs[0])
pyv.lfmPCol(dbs[0],fVar,vBds=fBds,pcOpac=1.0,Inv=fInv,Log=fLog,cMap=fMap,Legend=True)
#pyv.chopInner2D()
pyv.to3D(opNum=0)

#Plot particles
if (doScat):
	#Injected
	ActivateDatabase(dbs[2])
	pyv.lfmPScat(dbs[2],v3="pZero",v4=pVar,vBds=pBdsI,cMap=pMapI,Log=False,Inv=False,pSize=pSzI,Legend=True)
	pyv.onlyIn()
	
	#Trapped
	ActivateDatabase(dbs[3])
	pyv.lfmPScat(dbs[3],v3="pZero",v4=pVar,vBds=pBdsT,cMap=pMapT,Log=False,Inv=False,pSize=pSzT,Legend=doLegP)
	pyv.onlyIn()
else:
	AddMol(dbs[2],cMap=pMapI,vBds=pBdsI,rScl=rScl,doLeg=True)
	pyv.onlyIn(opNum=1)

	AddMol(dbs[3],cMap=pMapT,vBds=pBdsT,rScl=rScl,doLeg=doLegP)
	pyv.onlyIn(opNum=1)

	

#Cleanup
# plXs = [0.03,0.03,0.05,0.01]
# plYs = [0.9,0.9,0.25,0.25]

tH = 0.025
lH = 0.015
stH = 0.03
if (doTwoP):
	plXs = [0.75,0.75,0.05,0.02]
	plYs = [0.1,0.2,0.9,0.9]
else:
	plXs = [0.75,0.75,0.035,0.02]
	plYs = [0.1,0.2,0.9,0.9]

plTits = ["",fLab, "", "Particle Energy [keV]"]
pyv.cleanLegends(plXs,plYs,plTits,plHt=lH)

#Clean field strength


AnLab = 'Plot0001'
GetAnnotationObject(AnLab).orientation = 3
GetAnnotationObject(AnLab).yScale = 0.75
GetAnnotationObject(AnLab).xScale = 0.75
GetAnnotationObject(AnLab).fontHeight = tH

#Clean first particle energy
AnLab = 'Plot0003'
GetAnnotationObject(AnLab).drawTitle = 0
GetAnnotationObject(AnLab).drawLabels = 0
GetAnnotationObject(AnLab).xScale = 0.75
GetAnnotationObject(AnLab).fontHeight = tH
AnLab = 'Plot0002'
GetAnnotationObject(AnLab).xScale = 0.75
GetAnnotationObject(AnLab).fontHeight = tH

pyv.setAtts(pHeight=H,pWidth=W)
pyv.SetBaseColors()

SetWin(W,H)

#----------------------------------
#Add Earth
ActivateDatabase(dbs[0])
AddPlot("Contour","RadAll")
cOp = GetPlotOptions()
cOp.legendFlag = 0
cOp.contourMethod = 1
cOp.contourValue = (rIn)
cOp.colorType = 0
cOp.singleColor = (0, 204, 255, rOpac) #Last number opacity (out of 255) 
SetPlotOptions(cOp)
AddOperator("Revolve",0)

#----------------------------------
#Add axis lines
LW = 2
# AddLine(x1=(12.5,0,0),lw=LW)
# AddLine(x1=(0,20,0),lw=LW)
# AddLine(x1=(0,0,7.5),lw=LW)
# pyv.genTit(titS="Sunward",Pos=(0.75,0.4),height=lH)
# pyv.genTit(titS="Dusk",Pos=(0.7,0.675),height=lH)
AnOp = GetAnnotationAttributes()
AnOp.axes3D.triadFlag = 1
SetAnnotationAttributes(AnOp)

#----------------------------------
#Set toggles
ToggleLockViewMode()
ToggleMaintainViewMode()

DrawPlots()

SetTimeSliderState(100)
SaveWindow()
sys.exit()

if (Prod):
	pyv.doTimeLoop(Ninit=1,T0=T0,dt=dt,Ns=Nsk,Save=True,tLabPos=(0.2,0.8),tH=stH,Trim=False)
	pyv.makeVid(Clean=True,outVid=outVid,tScl=vidScl)
	DeleteAllPlots()
else:
	print("Here")
	#OpenGUI()
	SetWin(W,H)

