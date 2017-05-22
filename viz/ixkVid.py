import sys
import os
import numpy as np
import datetime
from visit import *
from visit_utils import *
from visit_utils.common import lsearch #lsearch(dir(),"blah")
import pyVisit as pyv


#Config
#----------
pHeight = 1200
pWidth = 1600

Quiet = True 
doAvg = True

outVid = "ixk.mp4"
T0Str = "2013-03-16T17:10:00Z"

vidScl = 4 #>1 to slow down
IBds = [1.0e-2,1.0e+5]
Nk = 50
Nx = 100

Base = os.path.expanduser('~') + "/Work/StormPSD/Data"
kDB = Base + "/Merge/KCyl.xmf"
cMap = "inferno"
xLab = "X [Re]"
yLab = "Log(K) [keV]"

#------------
#Make Video
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"

if (Quiet):
	LaunchNowin()
else:
	Launch()


#Do some defaults
pyv.pvInit()

#Get time data
md0 = GetMetaData(kDB)

dt = md0.times[1] - md0.times[0]
dT0s = md0.times[0]

dT0 = datetime.timedelta(seconds=dT0s)
T0 = dT0 + datetime.datetime.strptime(T0Str,T0Fmt)

OpenDatabase(kDB)
DefineScalarExpression("L","cylindrical_radius(mesh)")
DefineScalarExpression("X","coord(mesh)[0]")
DefineScalarExpression("Y","coord(mesh)[1]")
DefineScalarExpression("K","10^coord(mesh)[2]")
DefineScalarExpression("Phi","cylindrical_theta(mesh)")
DefineScalarExpression("LogK","log10(K)")
DefineScalarExpression("IScl","f")

#DefineScalarExpression("IScl","2.9979E10*recenter(I)")

pyv.lfmPCol(kDB,"IScl",vBds=IBds,Log=True,cMap=cMap)

if (doAvg):
	#Cut out meridional slab
	AddOperator("Threshold")
	tOp = GetOperatorOptions(0)
	tOp.listedVarNames = ("default","Y")
	tOp.lowerBounds = (-1e+37,-0.5)
	tOp.upperBounds = (1e+37 , 0.5)
	SetOperatorOptions(tOp)
	
	#Bin into X
	AddOperator("DataBinning")
	dbOps = GetOperatorOptions(1)
	dbOps.numDimensions = 1 #2D
	dbOps.dim1Var = "X"
	dbOps.dim1NumBins = Nx
	dbOps.dim1SpecifyRange = 1
	dbOps.dim1MinRange = -15.0
	dbOps.dim1MaxRange = 15.0
	
	dbOps.dim2Var = "LogK"
	dbOps.dim2SpecifyRange = 1
	dbOps.dim2MinRange = 1.0
	dbOps.dim2MaxRange = 3.75
	
	dbOps.dim2NumBins = Nk
	dbOps.reductionOperator = 0 #Sum
	dbOps.varForReduction = "IScl"
	SetOperatorOptions(dbOps)
else:
	AddOperator("Slice")
	sOp = GetOperatorOptions(0)
	sOp.axisType = 1	
	SetOperatorOptions(sOp)
	pyv.SetWin2D((-15,15,1.0,3.6))
ToggleFullFrameMode()

#OpenGUI()

#Gussy things up
plXs = [0.03]
plYs = [0.4]

plTits = ["Intensity\ns-1 cm-2 keV-1"]

pyv.cleanLegends(plXs,plYs,plTits)
pyv.setAtts(pHeight=pHeight,pWidth=pWidth,xLab=xLab,yLab=yLab)

anAt = AnnotationAttributes()
anAt.axes2D.xAxis.grid = 1
anAt.axes2D.yAxis.grid = 1
anAt.foregroundColor = (0, 204, 255, 255)
anAt.userInfoFlag = 0
anAt.databaseInfoFlag = 0
anAt.timeInfoFlag = 0
anAt.axes2D.xAxis.title.userTitle = 1
anAt.axes2D.yAxis.title.userTitle = 1
anAt.axes2D.xAxis.title.title = xLab
anAt.axes2D.yAxis.title.title = yLab
SetAnnotationAttributes(anAt)

#Let's see what we got
DrawPlots()

#Do time loop
pyv.doTimeLoop(Ninit=1,T0=T0,dt=dt,Save=True,tLabPos=(0.3,0.025),Trim=True)
pyv.makeVid(Clean=True,outVid=outVid,tScl=vidScl)

DeleteAllPlots()
CloseDatabase(kDB)

