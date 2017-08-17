#!/usr/bin/env python
import argparse
import numpy as np
import math
import lfmViz as lfmv
import lfmPostproc as lfmpp

#Various constants
Qe = 4.8032e-10 #Electron charge [CGS]
Me = 9.1094e-28 #Electron mass [g]
Re = 6.38e+8 #Earth radius [cm]
clightCMS = 2.9979e+10 #Speed of light [cm/s]

EBscl = (Qe*Re/Me)/(clightCMS**2.0)

def sclMag(Bx,By,Bz):
	#Incoming data is in dim-less units (EB from LFMTP)
	#Scale to nT
	eb2cgs = 1/EBscl #Convert EB to CGS
	G2nT = 10**5.0 #Convert Gauss to nT
	scl = (eb2cgs*G2nT)
	return scl*Bx,scl*By,scl*Bz

def sclElec(Ex,Ey,Ez):
	#Incoming data is in dim-less units (EB from LFMTP)
	#Scale to mili-Volts over meters [mV/m]

	eb2cgs = 1/EBscl #Convert EB to CGS
	G2T = 10**(-4.0)
	V2mV = 10**3.0
	clight = 2.9979e+8 #Speed of light, [m/s]
	#Convert E -> "Gauss" then Tesla, xC_light to V/m then mV/m
	scl = eb2cgs*G2T*clight*V2mV
	return scl*Ex,scl*Ey,scl*Ez

if __name__ == "__main__":
	#Set defaults
	Stub = "eqSlc"

	fldStr = 'B'
	labS = ["dBx","dBy","dBz","Ex","Ey","Ez"]

	parser = argparse.ArgumentParser(description="Reads EB table and particle output to construct matching cadence equatorial slices for visualization.")
	parser.add_argument('tab',nargs=1,type=str,metavar='ebtab.txt',help="EB table to read VTI info from")
	parser.add_argument('dat',nargs=1,type=str,metavar='H.100K.h5part',help="Particle output data to match cadence to")
	parser.add_argument('--meridional', nargs=1, type=float, metavar=90, help='Output slice at given meridional angle (0 or 90) rather than equatorial slice')
	parser.add_argument('--output', nargs=1, type=str, metavar='slice', help='Base name for output files')
	parser.add_argument('--thick', nargs=1, type=int, metavar=9, help='Thickness of slice (# of cells)')
	parser.add_argument('--ns', nargs=1, type=int, metavar=1, help='Start index')
	parser.add_argument('--ne', nargs=1, type=int, metavar=1, help='End index')


	#Finished getting arguments, parse and move on
	args = parser.parse_args()
	ebFile = args.tab[0]
	datFile = args.dat[0]

	meridional = args.meridional
	if meridional is not None:
		meridional = meridional[0]

	if args.output is not None:
		Stub = args.output[0]

	thick = 1
	if args.thick is not None:
		thick = args.thick[0]

	#Start by reading data
	Tv,vtis = lfmv.readTab(ebFile) #VTI's and their associated time

	Tps, pIds = lfmpp.getH5p(datFile,"id") #Get particle times
	Nt = len(Tps)


	ns = args.ns[0]
	ne = args.ne[0]

	Nstart = ns
	Nend = max(Nt,ne)

	print("Total number of slices = %d"%Nt)
	print("Writing between %d and %d"%(Nstart,Nend))
	#Get grid data from first VTI file in ebtab
	o3d,dx3d,ex3d = lfmv.getVTI_Grid(vtis[0])

	Ndim = len(dx3d)
	if (Ndim == 3):
		Nz = np.abs( ex3d[5]-ex3d[4] )
		if (Nz == 0):
			Ndim = 2

	print("Using %d dimensional data"%(Ndim))
	#for i in range(Nt):
	for i in range(Nstart,Nend):
		ti = Tps[i]
		if (math.isnan(ti)):
			print("Exiting on nan")
			exit()
		#Construct output file name
		vtiFileOut = Stub + ".%04d"%(i) + ".vti"
		print("Writing %s @ t = %3.2f ..."%(vtiFileOut,ti))
		if (Ndim == 3):
			dBx,dBy,dBz = lfmv.interpTime(Tv,vtis,ti,'B',thick=thick,meridional=meridional)
			Ex,Ey,Ez    = lfmv.interpTime(Tv,vtis,ti,'E',thick=thick,meridional=meridional)

			#Scale
			#B [EB] -> B [nT]
			#E [EB] -> E [mV/m]
			dBx,dBy,dBz = sclMag(dBx,dBy,dBz)
			Ex,Ey,Ez    = sclElec(Ex,Ey,Ez)
		else:
			#Ndim == 2
			dBx = lfmv.interpTimeSclr(Tv,vtis,ti,'dBx')
			dBy = lfmv.interpTimeSclr(Tv,vtis,ti,'dBy')
			dBz = lfmv.interpTimeSclr(Tv,vtis,ti,'dBz')
			Ex  = lfmv.interpTimeSclr(Tv,vtis,ti,'Ex')
			Ey  = lfmv.interpTimeSclr(Tv,vtis,ti,'Ey')
			Ez  = lfmv.interpTimeSclr(Tv,vtis,ti,'Ez')

		Qs = [dBx,dBy,dBz,Ex,Ey,Ez]
		lfmv.writeVTI_Eq(Qs,o3d,dx3d,ex3d,vtiFileOut,labS,t=ti,thick=thick,meridional=meridional)


