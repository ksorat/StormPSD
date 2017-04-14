#Various routines to deal with K-Cylinders from PSDs
import numpy as np
import datetime
Re = 6.38e+3 #Earth radius [km]

#Get trajectory from orbit file, return time in seconds after T0 (datetime)
def getTraj(oFile,T0S):
	T0Fmt = "%Y-%m-%dT%H:%M:%SZ"
	#Expecting format: Year,Month,Day,Hour,Minute,Second, SMX [KM], SMY [KM], SMZ [KM]
	VA = np.loadtxt(oFile,skiprows=1).T
	Nt = VA.shape[1]

	T = np.zeros(Nt)
	X = VA[6,:]/Re
	Y = VA[7,:]/Re
	Z = VA[8,:]/Re

	T0 = datetime.datetime.strptime(T0S,T0Fmt)
	for i in range(Nt):
		tSlc = VA[0:6,i]
		tSC = "%d-%d-%dT%d:%d:%dZ"%(tuple(tSlc))
		Ti = datetime.datetime.strptime(tSC,T0Fmt)
		dt = Ti-T0
		T[i] = dt.total_seconds()