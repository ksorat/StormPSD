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
Quiet = False
outVid = "ilk.mp4"
T0Str = "2013-03-16T17:10:00Z"

vidScl = 4 #>1 to slow down
IBds = [1.0e-4,1.0e+2]
Nk = 50
Nl = 40

Base = os.path.expanduser('~') + "/Work/StormPSD/Data"
kDB = Base + "/Merge/KCyl.xmf"
cMap = "magma"
xLab = "L [Re]"
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
DefineScalarExpression("K","10^coord(mesh)[2]")
DefineScalarExpression("LogK","log10(K)")
DefineScalarExpression("IScl","2.9979E10*recenter(I)")

#Get time evolution of # of test particles
pyv.lfmPCol(kDB,"Ntp")
DrawPlots()
QueryOverTime("Variable Sum")
SetActiveWindow(2)
pinfo = GetPlotInformation()
C = np.array(pinfo['Curve'])
TPt = C[1:-1:2] #Number of TPs
DeleteWindow()

SetActiveWindow(1)
DeleteAllPlots()

pyv.lfmPCol(kDB,"IScl",vBds=IBds,Log=True,cMap=cMap)
AddOperator("DataBinning")
dbOps = GetOperatorOptions(0)
dbOps.numDimensions = 1 #2D
dbOps.dim1Var = "L"
dbOps.dim1NumBins = Nl
dbOps.dim1SpecifyRange = 1
dbOps.dim1MinRange = 2.0
dbOps.dim1MaxRange = 18.0

dbOps.dim2Var = "LogK"
dbOps.dim2SpecifyRange = 1
dbOps.dim2MinRange = 1.0
dbOps.dim2MaxRange = 3.6

dbOps.dim2NumBins = Nk
dbOps.reductionOperator = 0 #Sum
dbOps.varForReduction = "IScl"
SetOperatorOptions(dbOps)

ToggleFullFrameMode()

#Gussy things up
#tit = pyv.genTit( titS=titS)
plXs = [0.03]
plYs = [0.4]

plTits = ["Intensity"]

pyv.cleanLegends(plXs,plYs,plTits)
pyv.setAtts(xLab=xLab,yLab=yLab)

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
pyv.doTimeLoop(T0=T0,dt=dt,Save=True,tLabPos=(0.3,0.025),Trim=True)
pyv.makeVid(Clean=True,outVid=outVid,tScl=vidScl)

DeleteAllPlots()
CloseDatabase(kDB)

