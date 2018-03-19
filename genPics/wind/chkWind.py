import numpy as np
import glob
import lfmPreproc as lfmpre
from pyhdf.SD import SD, SDC

lfmDir = "os.path.expanduser('~')" + "/Work/lfmHDFs"
fIns = glob.glob(lfmDir+"/*17T*.hdf")

Nf = len(fIns)

t = np.zeros(Nf)
N = np.zeros(Nf)

for i in range(Nf):
	fIn = fIns[i]
	hdffile = SD(fIn)
	t[i] = hdffile.attributes().get('time')
	D3 = lfm.getHDFScl(hdffile,"rho")
	print("fIn = %s / T = %f"%(fIn,t[i]))