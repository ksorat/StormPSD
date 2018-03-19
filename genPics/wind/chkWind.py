import numpy as np
import glob
import lfmInterp as lfm
import lfmPreproc as lfmpre
from pyhdf.SD import SD, SDC
import os

lfmDir = os.path.expanduser('~') + "/Work/lfmHDFs/StPatty"
fIns = glob.glob(lfmDir+"/*17T*.hdf")
print(lfmDir)
#print(fIns)
Nf = len(fIns)

t = np.zeros(Nf)
N = np.zeros(Nf)

for i in range(Nf):
	fIn = fIns[i]
	hdffile = SD(fIn)
	print("Reading %s"%(fIn))
	t[i] = hdffile.attributes().get('time')
	D3 = lfm.getHDFScl(hdffile,"rho")
	print("\t T = %f"%(t[i]))
	#print("fIn = %s / T = %f"%(fIn,t[i]))