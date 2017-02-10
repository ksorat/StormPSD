#Get time series from a bunch of LFM HDFs
import numpy as np
import glob
import os

#HDF directory
lfmDir = os.path.expanduser('~') + "/Work/StormPSD/lfmData"
fIns =glob.glob(lfmDir + "/*.hdf")
