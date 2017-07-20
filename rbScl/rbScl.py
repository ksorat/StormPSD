#Make plots showing intensity/pitch angle access for RB due to latitude
import kCyl as kc
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

#Time data
T0Str = "2013-03-16T17:10:00Z"
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"
Labs = ["NULL","Trapped","Injected","Combined"]

tMin = 33600.0
tMax = 189000.0
#tMax = 195000.0

#RB Opts
rbStrs = ["A","B"]

NumRB = len(rbStrs)

for nrb in range(NumRB):
	#Create relevant files
	rbStr = rbStrs[nrb]

	OrbF = "vaporbRB" + rbStr.upper() + ".txt"

	#Get projected/non-proj trajectory
	#Tsc = seconds after T0
	#Non-Proj
	Tsc,Xsc,Ysc,Zsc = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk=1,doEq=False)
	Tp,Xp,Yp,Zp = kc.getTraj(OrbF,T0Str,tMin,tMax,Nsk=1,doEq=True)
