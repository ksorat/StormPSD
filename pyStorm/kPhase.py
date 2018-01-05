#Various routines to deal with full phase spaces
import numpy as np
import datetime

Re = 6.38e+3 #Earth radius [km]
Rmin = 1.9 #Minimum worthwhile radius

#Expecting format: Year,Month,Day,Hour,Minute,Second, SMX [KM], SMY [KM], SMZ [KM]
T0Fmt = "%Y-%m-%dT%H:%M:%SZ"


def getPSD(fIn,fVar="f"):
	import h5py
	
	with h5py.File(fIn,'r') as hf:
		#Start by getting number of steps
		KeyL = list(hf.keys())
		LStp = [s for s in KeyL if 'Step' in s]
		Vars = [s for s in KeyL if 'Step' not in s]
		Nt = len(LStp)
		#print(Vars)
		R  = np.array(hf.get("Cr").value)
		P  = np.array(hf.get("Cphi").value)
		K  = np.array(hf.get("Ck").value)
		A  = np.array(hf.get("Calpha").value)

		dG = hf.get("dG").value.T
		dG = dG.squeeze()

		Nr = R.size
		Np = P.size
		Nk = K.size
		Na = A.size

		V = np.zeros((Nr,Np,Nk,Na,Nt))
		t = np.zeros(Nt)
		
		StpS = [int(s.replace('Step#','')) for s in LStp]
		StpS.sort()
		for n in range(0,Nt):
			
			gId = "Step#%d"%(StpS[n])
			grp = hf.get(gId)
			Vt = grp.get(fVar).value.T
			V[:,:,:,:,n] = Vt.squeeze()
			t[n] = grp.attrs.get("time")

	return R,P,K,A,t,dG,V