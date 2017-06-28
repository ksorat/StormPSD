#2D video showing positions of RB-A/B and intensity

import sys
import os
import numpy as np
import datetime
from visit import *
from visit_utils import *
from visit_utils.common import lsearch #lsearch(dir(),"blah")
import pyVisit as pyv

Quiet = False
doVid = False
doField = True

outVid = "irb2d.mp4"
Base = os.path.expanduser('~') + "/Work/StormPSD/Data/"
fPSD = Base + "Merge/KCyl_StormA.xmf"
fEq  = Base + "psdEq/eqSlc.*.vti database"

kSlc = 1000
cMapI = "viridis"
Ibds = [1.0e-1,5.0e+3]

cMapF = "RdGy"
fbds = [-35,35]
Nc = 17

if (Quiet):
	LaunchNowin()
else:
	Launch()

pyv.lfmExprs()

#Open databases
OpenDatabase(fPSD)
OpenDatabase(fEq)

md0 = GetMetaData(fPSD)
dt = md0.times[1] - md0.times[0]
T0 = md0.times[0]# - md0.times[0] #Reset to zero

#Create correlation
CreateDatabaseCorrelation("ifCor",[fPSD,fEq],0)

#Create intensity pcolor
pyv.lfmPCol(fPSD,"f",cMap=cMapI,vBds=Ibds,Log=True)
AddOperator("Slice")
sOp = GetOperatorOptions(0)
sOp.originType = 1
sOp.originIntercept = np.log10(kSlc)
SetOperatorOptions(sOp)
pyv.addThreshold("f",Ibds[0],1.0e+12,opNum=1)

#Create field contours
if (doField):
	pyv.plotContour(fEq,"dBz",cmap=cMapF,vBds=fbds,Nlevels=Nc,lineWidth=1)
	cOp = GetPlotOptions()

#Set background and such
pyv.setAtts()
anAt = AnnotationAttributes()
anAt.backgroundColor = (0,0,0,0)
anAt.foregroundColor = (0, 204, 255, 255)
#anAt.axes2D.xAxis.grid = 1
#anAt.axes2D.yAxis.grid = 1

anAt.userInfoFlag = 0
anAt.databaseInfoFlag = 0
anAt.timeInfoFlag = 0
anAt.axes2D.xAxis.title.userTitle = 1
anAt.axes2D.yAxis.title.userTitle = 1
anAt.axes2D.xAxis.title.title = "X [Re]"
anAt.axes2D.yAxis.title.title = "Y [Re]"
SetAnnotationAttributes(anAt)

#tit = pyv.genTit( titS=titS,Pos=(0.35,0.955),height=0.02)
#pyv.cleanLegends(plXs,plYs,plTits)

pyv.SetWin2D([-15,12.5,-20,20])

DrawPlots()

if (doVid):
	print("Blah")
else:
	SetTimeSliderState(50)