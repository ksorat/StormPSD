#Various routines to deal with K-Cylinders from PSDs
import numpy as np
import datetime
Re = 6.38e+3 #Earth radius [km]
Rmin = 1.9 #Minimum worthwhile radius

#Expecting format: Year,Month,Day,Hour,Minute,Second, SMX [KM], SMY [KM], SMZ [KM]
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"

#Smoothing parameters for center, cross, diagonals
a0 = np.exp(0)
#aC = np.exp(-0.5)
#aD = np.exp(-1.0)
aC = np.exp(-1.0)
aD = np.exp(-2.0)

aScl2D = a0+4*aC+4*aD
aScl1D = a0+2*aC

#Given sin^n(alpha) dep. on intensity calculate fraction based on accessible Alpha
def getIScl(Ac,en=2.0):
	Na = 360
	A = np.linspace(0,0.5*np.pi,Na)
	da = A[1]-A[0]
	Ia = np.sin(A)**en
	Ic = np.zeros(Ia.shape)
	Nt = len(Ac)
	I0 = Ia.sum()
	
	It = np.zeros(Nt)
	for n in range(Nt):
		Ic[:] = Ia[:]
		Icut = (A>Ac[n])
		Ic[Icut] = 0.0
		It[n] = Ic.sum()/I0
	

	return It



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

#Reads in and sums multiple data cylinders
#Assume same domain
def getCyls(fIns,fVar="f"):
	N = len(fIns)
	R,P,K,t,I = getCyl(fIns[0],fVar)
	for i in range(1,N):
		Rp,Pp,Kp,tp,Ip = getCyl(fIns[i],fVar)
		I = I + Ip
	return R,P,K,t,I	

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

#Create interpolator for K-Cylinder
def GetInterp(R,P,K,t,I,imeth="linear",fillVal=0.0):
	import scipy
	import scipy.interpolate
	Irpkt = scipy.interpolate.RegularGridInterpolator((R,P,K,t),I,method=imeth,bounds_error=False,fill_value=fillVal)
	
	return Irpkt


#Interpolate intensities at trajectories
def InterpI(Ii,Xsc,Ysc,Tsc,K):
	Rsc = np.sqrt(Xsc**2.0+Ysc**2.0)
	Psc = np.arctan2(Ysc,Xsc)
	iP = (Psc<0); Psc[iP] = Psc[iP]+2*np.pi

	Nsc = len(Tsc)
	Nk = len(K)
	#Create arrays
	iPts = np.zeros((Nk,4))
	Isc = np.zeros((Nsc,Nk))
	for i in range(Nsc):
		#Evaluate for all energies at this point
		iPts[:,0] = Rsc[i]
		iPts[:,1] = Psc[i]
		iPts[:,2] = K
		iPts[:,3] = Tsc[i]
		Isc[i,:] = Ii(iPts)
		# if (Rsc[i]<=Rmin):
		# 	Isc[i,:] = 0.0
		# else:
		# 	Isc[i,:] = Ii(iPts)

	return Isc

#Interpolate intensities at dipole-projection and attenuate model intensity for latitude
def InterpI_XYZ(Ii,X,Y,Z,Tsc,K,doScl=True,en=2.0):
	Nsc = len(Tsc)
	Nk = len(K)

	#Calculate dipole projections from 3D positions
	R = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
	Lam = np.arcsin(Z/R)
	Psc = np.arctan2(Y,X)
	iP = (Psc<0); Psc[iP] = Psc[iP]+2*np.pi

	Req = R/(np.cos(Lam)**2.0)
	Xeq = Req*np.cos(Psc)
	Yeq = Req*np.sin(Psc)

	#Calculate attention factor from latitude
	cL = np.cos(Lam)
	sL = np.sin(Lam)
	lArg = cL**6.0/np.sqrt(1+3*sL*sL)
	AeC = np.arcsin(np.sqrt(lArg)) #Critical equatorial alpha for this latitude
	#I0 = (AeC-0.5*np.sin(2*AeC))/(0.5*np.pi) #Attenuation factor
	I0 = getIScl(AeC,en=en)

	#Create arrays
	iPts = np.zeros((Nk,4))
	Isc = np.zeros((Nsc,Nk))
	for i in range(Nsc):
		#Get 3D position, project to plane along dipole
		#Evaluate for all energies at this point
		iPts[:,0] = Req[i]
		iPts[:,1] = Psc[i]
		iPts[:,2] = K
		iPts[:,3] = Tsc[i]
		if (Req[i]<=Rmin):
			Isc[i,:] = 0.0
		else:
			Isc[i,:] = Ii(iPts)
		if (doScl):
			Isc[i,:] = Isc[i,:]*I0[i]
		#print("Scaling by %f"%(I0[i]))
	return Isc

#Interpolate intensities from KCyl onto RB trajectory
#SimKC = [R,P,K,Tkc,Is]
#rbDat = [Xsc,Ysc,Zsc,Tsc,Ksc]

def InterpSmooth(SimKC,rbDat,Niter=1,NiterT=1,doZScl=True,en=2.0):
	Xsc,Ysc,Zsc,Tsc,Ksc = rbDat
	R,P,K,Tkc,Ikc = SimKC

	
	IkcS = np.zeros(Ikc.shape)
	IkcS[:] = Ikc[:]
	#IkcS = Ikc
	#Now apply smoothing window
	Nr = len(R)
	Np = len(P)
	Nt = len(Tsc)
	Nk = len(Ksc)

	print("Interpolating from KCyl onto T,K grid of size (%d,%d)\n"%(Nt,Nk))
	print("\tdtRB = %f"%(Tsc[1]-Tsc[0]))
	IkcS = SmoothKCyl(R,P,IkcS,Niter=Niter)

	Ii = GetInterp(R,P,K,Tkc,IkcS)
	Isc = InterpI_XYZ(Ii,Xsc,Ysc,Zsc,Tsc,Ksc,doScl=doZScl,en=en)


	IscS = np.zeros(Isc.shape)
	IscS[:] = Isc[:]
	for n in range(NiterT):
		IscS = SmoothIterT(IscS)

	return IkcS,IscS



def SmoothKCyl(R,P,Ikc,Niter=1):
	IkcS = np.zeros(Ikc.shape)
	IkcS[:] = Ikc[:]
	Nr = len(R)
	Np = len(P)
	xx = np.zeros((Nr,Np))
	yy = np.zeros((Nr,Np))
	for i in range(Nr):
		for j in range(Np):
			xx[i,j] = R[i]*np.cos(P[j])
			yy[i,j] = R[i]*np.sin(P[j])

	for n in range(Niter):
		IkcS = SmoothIter(IkcS,xx,yy)
	return IkcS

#Does one iteration of smoothing on Ikc
def SmoothIter(Ikc,xx,yy,doT=True):
	Nr = Ikc.shape[0]
	Np = Ikc.shape[1]
	Nk = Ikc.shape[2]
	Nt = Ikc.shape[3]

	IkcS = np.zeros(Ikc.shape)
	IkcS[:] = Ikc[:]


	for i in range(1,Nr-1):
		for j in range(Np):
			iM = i-1
			iP = i+1
			jM = j-1
			jP = j+1
			if (j==0):
				jM = Np-1
			if (j==Np-1):
				jP = 0

			# a00,aCj,aCip,aCim,aDp,aDm = getWeights(xx,yy,i,j,iM,iP,jM,jP)
			# aScl = a0 + 2*aCj + aCip + aCim + 2*aDp + 2*aDm

			# IkcS[i,j,:,:] = (a00*Ikc[i,j,:,:] + 
			# 				 aCj*(Ikc[i,jP,:,:] + Ikc[i,jM,:,:]) +
			# 				 aCip*Ikc[iP,j,:,:] + aCim*Ikc[iM,j,:,:] +
			# 				 aDp*(Ikc[iP,jP,:,:] + Ikc[iP,jM,:,:]) +
			# 				 aDm*(Ikc[iM,jP,:,:] + Ikc[iM,jM,:,:])
			# 				)/aScl

			IkcS[i,j,:,:] = (a0*Ikc[i,j,:,:] + 
							 aC*(Ikc[iP,j,:,:] + Ikc[iM,j,:,:] + Ikc[i,jP,:,:] + Ikc[i,jM,:,:]) +
							 aD*(Ikc[iP,jP,:,:] + Ikc[iP,jM,:,:] + Ikc[iM,jP,:,:] + Ikc[iM,jM,:,:])
							 )/aScl2D

	Ikc[:] = IkcS[:]
	for n in range(1,Nt-1):
		IkcS[:,:,:,n] = (a0*Ikc[:,:,:,n]+aC*Ikc[:,:,:,n-1]+aC*Ikc[:,:,:,n+1])/aScl1D	
	return IkcS

#Weight function
#Order: a0,aCj,aCip,aCim,aDp,aDm
def getWeights(xx,yy,i,j,iM,iP,jM,jP):
	x0 = xx[i,j]
	y0 = yy[i,j]
	Is = [i,i,iP,iM,iP,iM]
	Js = [j,jP,j,j,jP,jP]
	#Use main diagonal for L
	L = np.sqrt( (xx[iM,jM]-xx[iP,jP])**2.0 + (yy[iM,jM]-yy[iP,jP])**2.0)
	L = np.sqrt( (xx[i,jP]-xx[i,jM])**2.0 + (yy[i,jP]-yy[i,jM])**2.0)
	A = np.zeros(6)
	for n in range(6):
		xp = xx[Is[n],Js[n]]
		yp = yy[Is[n],Js[n]]
		t = np.sqrt( (xp-x0)**2.0 + (yp-y0)**2.0 )/L
		if (t<=1):
			Wt = (1-t**3.0)**3.0
			#Wt = 0.75*(1-t**2.0)
		else:
			Wt = 0.0
		A[n] = Wt
	return A[0],A[1],A[2],A[3],A[4],A[5]

#Does one time iteration of smoothing on IscS
def SmoothIterT(Isc):
	Nt = Isc.shape[0]
	IscS = np.zeros(Isc.shape)
	IscS[:] = Isc[:]
	for n in range(1,Nt-1):
		IscS[n,:] = (a0*Isc[n,:] + aC*Isc[n+1,:] + aC*Isc[n-1,:])/aScl1D
	return IscS
#Get Intensity from RBSP CDF
def GetRBSP(fIn,T0S,tMin=None,tMax=None,rbID="rbspa",rbSK=1,CutR=True):
	from spacepy import pycdf
	cdf = pycdf.CDF(fIn)
	print(cdf)

	#Get main data
	T    = cdf[rbID+'_ect-mageis_l2_ele_time_epoch'][...]
	kRB  = cdf[rbID+'_ect-mageis_l2_ele_FESA_channel_energy'][...]
	dkRB = cdf[rbID+'_ect-mageis_l2_ele_FESA_channel_width'][...]
	Itk  = cdf[rbID+'_ect-mageis_l2_ele_FESA'][...]
	L    = cdf[rbID+'_ect-mageis_l2_ele_L'][...]
	cdf.close()

	#Pull out empty channels, assuming at bottom
	Nt = kRB.shape[0]
	dkT = np.zeros(Nt)
	for i in range(Nt):
		dkT[i] = np.abs((kRB[i,:] - kRB[0,:])).max()
	MdkT = dkT.max()
	if (MdkT>1.0e-8):
		i0 = dkT.argmax()
		print("Max dkT = %f"%(dkT.max()))
		print("i0 = %d"%(i0))
		print("T = %s"%(T[i0]))
		print("K0 = %s"%(kRB[0,:]))
		print("Kd = %s"%(kRB[i0,:]))
	
	kMax = kRB.max(axis=0)
	k0 = (kMax>0).argmax()
	kRB  = kRB [:,k0:]
	dkRB = dkRB[:,k0:]
	Itk  = Itk [:,k0:]

	Itk[kRB<=0] = 0.0
	#print("RB Energies = %s"%str(kRB[0,:]))
	#print("RB Widths = %s"%str(dkRB[0,:]))
	#Pull data

	kRB = kRB[0,:]
	dkRB = dkRB[0,:]


	#Constrain time domain
	if (CutR):
		I,Ts = CutTime(T,T0S,tMin=tMin,tMax=tMax)
		Itk = Itk[I,:]
		L = L[I]
		Itk[L<=Rmin,:] = 0.0

	#Pull out duplicate times
	Nt = Ts.shape[0]
	dt = np.zeros(Nt)
	dt[0:-1] = Ts[1:]-Ts[0:-1]
	I = (dt>=1.0e-8)
	Ts = Ts[I]
	Itk = Itk[I,:]
	
	#Ts = Ts[L>=Rin]

	#Enforce skipping if necessary
	if (rbSK>1):
		Ts = Ts[0::rbSK]
		Itk = Itk[0::rbSK,:]

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

#Get pcolor bounds from R,Phi
def xy2rp(R,P):
	Nr = len(R)
	Np = len(P)

	Pi = np.linspace(0,2*np.pi,Np+1)
	Ri = np.zeros(Nr+1)
	for n in range(Nr-1):
		dr = 0.5*(R[n+1]-R[n])
		Ri[n] = R[n] - dr

	Ri[Nr-1] = R[Nr-1] + dr
	R0 = np.round(Ri[0])
	R1 = np.round(Ri.max())

	Ri = np.logspace(np.log10(R0),np.log10(R1),Nr+1)

	XX = np.zeros((Nr+1,Np+1))
	YY = np.zeros((Nr+1,Np+1))

	for n in range(Nr+1):
		for m in range(Np+1):
			XX[n,m] = Ri[n]*np.cos(Pi[m])
			YY[n,m] = Ri[n]*np.sin(Pi[m])
	return XX,YY	

def ResampleCyl(Ikc,Ntp,Ncut=4):
	Nr = Ikc.shape[0]
	Np = Ikc.shape[1]
	Nk = Ikc.shape[2]
	Nt = Ikc.shape[3]

	IkcS = np.zeros(Ikc.shape)
	IkcS[:] = Ikc[:]

	for t in range(Nt):
		for k in range(Nk):
			for i in range(1,Nr-1):
				for j in range(Np):
					n0 = Ntp[i,j,k,t]
					if (n0 >= Ncut):
						continue
					else:
						#Replace with TP-weighted average
						iM = i-1
						iP = i+1
						jM = j-1
						jP = j+1
						if (j==0):
							jM = Np-1
						if (j==Np-1):
							jP = 0
						NScl = Ntp[iP,j ,k,t] + Ntp[iM,j ,k,t] +Ntp[i ,jP,k,t] + Ntp[i ,jM,k,t]
						
						if (NScl>0):
							IkcS[i,j,k,t] = Ntp[iP,j ,k,t]*Ikc[iP,j ,k,t] + Ntp[iM,j ,k,t]*Ikc[iM,j ,k,t] + Ntp[i ,jP,k,t]*Ikc[i ,jP,k,t] + Ntp[i ,jM,k,t]*Ikc[i ,jM,k,t]
							IkcS[i,j,k,t] = IkcS[i,j,k,t]/NScl
	return IkcS
