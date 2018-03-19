import numpy as np
import glob
import lfmInterp as lfm
import lfmPreproc as lfmpre
from pyhdf.SD import SD, SDC
import os

#lfmDir = os.path.expanduser('~') + "/Work/lfmHDFs/StPatty"
lfmDir = "/glade/p/hao/wiltbemj/StPatty13/LR60-Quad-15s-AEH"
fIns = glob.glob(lfmDir+"/*mhd_2013-03-17T0[56]-?[028]-*.hdf")
print(lfmDir)
#print(fIns)
Nf = len(fIns)

t = np.zeros(Nf)
N = np.zeros(Nf)

for i in range(Nf):
	fIn = fIns[i]
	print("Reading %s"%(fIn))
	hdffile = SD(fIn,SDC.READ)
	#print(hdffile.datasets())
	t[i] = hdffile.attributes().get('time')
	D3 = lfm.getHDFScl(hdffile,"rho")
	print("\t T = %f"%(t[i]))
	print(D3.shape)
	#print("fIn = %s / T = %f"%(fIn,t[i]))
	hdffile.end()