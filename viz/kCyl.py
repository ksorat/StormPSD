#Various routines to deal with K-Cylinders from PSDs
import numpy as np
import datetime
Re = 6.38e+3 #Earth radius [km]


#Expecting format: Year,Month,Day,Hour,Minute,Second, SMX [KM], SMY [KM], SMZ [KM]
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"

#Get trajectory from orbit file, return time in seconds after T0 (datetime)
#Chop out times outside of tMin,tMax range
def getTraj(oFile,T0S,tMin=None,tMax=None,Nsk=1,doEq=False):
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

	if (doEq):
		#Project to equatorial plane along dipole line
		for i in range(Nt):
			x = X[i]
			y = Y[i]
			z = Z[i]
			r = np.sqrt(x**2.0+y**2.0+z**2.0)
			lamb = np.arcsin(z/r)
			phi = np.arctan2(y,x)
			req = r/(np.cos(lamb)**2.0)
			xeq = req*np.cos(phi)
			yeq = req*np.sin(phi)
			#print("xy / xyeq = %f,%f / %f,%f\n"%(x,y,xeq,yeq))
			#print("d(xy) = %f"%(np.sqrt( (x-xeq)**2.0 + (y-yeq)**2.0)))
			X[i] = xeq
			Y[i] = yeq
			Z[i] = 0.0

	if (tMin is not None):
		I = (T>=tMin) & (T<=tMax) 
		T = T[I]
		X = X[I]
		Y = Y[I]
		Z = Z[I]
		
	T = T[0:-1:Nsk]
	X = X[0:-1:Nsk]
	Y = Y[0:-1:Nsk]
	Z = Z[0:-1:Nsk]

	return T,X,Y,Z

#Reads in data cylinder
def getCyl(fIn,fVar="f"):
	import h5py
	
	with h5py.File(fIn,'r') as hf:
		#Start by getting number of steps
		KeyL = list(hf.keys())
		LStp = [s for s in KeyL if 'Step' in s]
		Vars = [s for s in KeyL if 'Step' not in s]
		Nt = len(LStp)
		#print(Vars)

		#Get dimension sizes
		R  = np.array(hf.get("Cr").value)
		P  = np.array(hf.get("Cphi").value)
		K  = np.array(hf.get("Ck").value)
		
		# Ri = np.array(hf.get("Ir").value)
		# Pi = np.array(hf.get("Iphi").value)
		# Ki = np.array(hf.get("Ik").value)

		P = P*np.pi/180.0

		Nr = R.size
		Np = P.size
		Nk = K.size
		t = np.zeros(Nt)
		print("Reading KCyl from %s of size (R,p,K,t) = (%d,%d,%d,%d)"%(fIn,Nr,Np,Nk,Nt))
		I = np.zeros((Nr,Np,Nk,Nt))
		for n in range(0,Nt):
			gId = "Step#%d"%(n)
			grp = hf.get(gId)
			I[:,:,:,n] = grp.get(fVar).value.T
			t[n] = grp.attrs.get("time")
		return R,P,K,t,I

#Smooth intensity data
def SmoothI(I,sig=1.0):
	import scipy
	import scipy.ndimage
	from scipy.ndimage.filters import gaussian_filter

	I = gaussian_filter(I,sig)
	return I
#Create interpolator for K-Cylinder
def GetInterp(R,P,K,t,I,imeth="linear"):
	import scipy
	import scipy.interpolate
	Irpkt = scipy.interpolate.RegularGridInterpolator((R,P,K,t),I,method=imeth,bounds_error=False,fill_value=0)

	return Irpkt

#Get Intensity from RBSP CDF
def GetRBSP(fIn,T0S,tMin=None,tMax=None,rbID="rbspa"):
	from spacepy import pycdf
	cdf = pycdf.CDF(fIn)
	print(cdf)
	#Get main data
	T    = cdf[rbID+'_ect-mageis_l2_ele_time_epoch'][...]
	kRB  = cdf[rbID+'_ect-mageis_l2_ele_FESA_channel_energy'][...]
	dkRB = cdf[rbID+'_ect-mageis_l2_ele_FESA_channel_width'][...]
	Itk  = cdf[rbID+'_ect-mageis_l2_ele_FESA'][...]

	cdf.close()
	#Pull data
	kRB = kRB[0,:]
	dkRB = dkRB[0,:]

	#Constrain time domain
	I,Ts = CutTime(T,T0S,tMin=tMin,tMax=tMax)
	Itk = Itk[I,:]

	return Ts,kRB,dkRB,Itk

def GetRBSP_DST(fIn,T0S,tMin=None,tMax=None,rbID="rbspa"):
	from spacepy import pycdf
	cdf = pycdf.CDF(fIn)

	#Get main data
	T   = cdf['DST_time_epoch'][...]
	dst = cdf['DST_DST'][...]

	#Constrain time domain
	I,Ts = CutTime(T,T0S,tMin=tMin,tMax=tMax)
	dst = dst[I]

	return Ts,dst

#Index of inside time domain
def CutTime(T,T0S,tMin=None,tMax=None):
	T0 = datetime.datetime.strptime(T0S,T0Fmt)
	Nt = len(T)
	Ts = np.zeros(Nt)
	for i in range(Nt):
		dt = T[i]-T0
		Ts[i] = dt.total_seconds()
	I = (Ts>=Ts.min())

	if (tMin is not None):
		I = I & (Ts>=tMin)
	if (tMax is not None):
		I = I & (Ts<=tMax)
	Ts = Ts[I]

	return I,Ts

#Convert from s after T0 to datetimes
def Ts2date(Ts,T0S):
	T0 = datetime.datetime.strptime(T0S,T0Fmt)
	Td = []
	Nt = len(Ts)
	for i in range(Nt):
		tdi = T0+datetime.timedelta(seconds=Ts[i])
		Td.append(tdi)
	Td = np.array(Td)
	return Td

#Given I(t,K) and K0 return total intensity above K0
# def ICum(K,K0,I):
# 	kC = (K>K0).argmax()
# 	sK = np.sqrt(K)
# 	Nk = K.shape[0]
# 	Nt = I.shape[0]
# 	Ic = np.zeros(Nt)
# 	for n in range(Nt):
# 		It = I[n,:]
# 		Iscl = sK*It
# 		Ic = 