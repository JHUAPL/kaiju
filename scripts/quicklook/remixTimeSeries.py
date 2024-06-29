#!/usr/bin/env python

################ first figure out the time ################

#standard python
import sys
import os
import argparse
import pickle
from argparse import RawTextHelpFormatter

# Numpy and matplotlib
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#Kaipy and related
import h5py
import kaipy.kaiH5 as kaiH5
import kaipy.remix.remix as remix
import kaipy.kaiViz as kv
from kaipy.kaiTools import MJD2UT

MainS = """Creates a summary plot of global ionospheric paramters including
CPCP, HP, FAC.  Creates a PNG image file and python pickle file with the data.
"""

ftag = 'msphere'
ptag = '.'
parser = argparse.ArgumentParser(description=MainS,
	formatter_class=RawTextHelpFormatter)
parser.add_argument('-id',type=str,metavar='runid',default=ftag,
	help='RunID of data (default: %(default)s)')
parser.add_argument('-path',type=str,metavar='path',default=ptag,
	help='Path to directory containing REMIX files (default: %(default)s)')

args = parser.parse_args()
#Open the file and read the time information
remixFile = os.path.join(args.path,args.id+'.mix.h5')
assert(os.path.isfile(remixFile))
nsteps,sIds=kaiH5.cntSteps(remixFile)
mjd = kaiH5.getTs(remixFile,sIds,aID='MJD')
utall = MJD2UT(mjd)

#Setup the dictionary data structures to handle the timeseries information
hemispheres = ('NORTH','SOUTH')
cpcp = {'NORTH':[],'SOUTH':[]}
hp = {'NORTH':[],'SOUTH':[]}
ipfac = {'NORTH':[],'SOUTH':[]}
ijh = {'NORTH':[],'SOUTH':[]}

#Read the data and calculate the integrated quantities
ri = 6500.0e3
for hemi in hemispheres:
	ion = remix.remix(remixFile,sIds.min())
	ion.init_vars(hemi)
	x = ion.ion['X']
	y = ion.ion['Y']
	areaMixGrid = ion.calcFaceAreas(x,y)*ri*ri
	for Id in sorted(sIds):
		ion = remix.remix(remixFile,Id)
		ion.init_vars(hemi)
		# Cross polar cap potential
		cpcp[hemi].append(np.max(ion.variables['potential']['data'])-
						 np.min(ion.variables['potential']['data']))
		# Integrated Postive FAC
		fac = ion.variables['current']['data']
		fac[fac < 0] = 0.0
		pfac = areaMixGrid*fac[:,:]
		ipfac[hemi].append(pfac.sum()/1.0e12)
		# Hemispheric Power
		flux = ion.variables['flux']['data']
		energy = ion.variables['energy']['data']
		hpcalc = areaMixGrid*energy[:,:]*flux[:,:]
		# Convert from keV/cm^2 to mW/m^2 to GW
		hp[hemi].append(hpcalc.sum()*1.6e-21)
		#jh = ion.variables['joule']['data']
		jh = ion.joule()*1.e3
		jhcalc = areaMixGrid*jh[:,:]
		ijh[hemi].append(jhcalc.sum()/1.0e12) #mW to GW


cpcp['units']='kV'
cpcp['name']=r'$\Phi$'
hp['units']='GW'
hp['name']='HP'
ipfac['units']='MA'
ipfac['name']='FAC'
ijh['units']='GW'
ijh['name']='JH'


#Plot the figure
figsize = (10,10)
fig = plt.figure(figsize=figsize)
gs = fig.add_gridspec(4,1)
Ax1 = fig.add_subplot(gs[0,0])
Ax2 = fig.add_subplot(gs[1,0],sharex=Ax1)
Ax3 = fig.add_subplot(gs[2,0],sharex=Ax1)
Ax4 = fig.add_subplot(gs[3,0],sharex=Ax1)
Ax1.plot(utall[1:],cpcp['NORTH'][1:])
Ax1.plot(utall[1:],cpcp['SOUTH'][1:])
kv.SetAxLabs(Ax1,None,cpcp['name']+' ['+cpcp['units']+']')
Ax2.plot(utall[1:],ipfac['NORTH'][1:])
Ax2.plot(utall[1:],ipfac['SOUTH'][1:])
kv.SetAxLabs(Ax2,None,ipfac['name']+' ['+ipfac['units']+']',doLeft=False)
Ax3.plot(utall[1:], hp['NORTH'][1:])
Ax3.plot(utall[1:], hp['SOUTH'][1:])
kv.SetAxLabs(Ax3,None,hp['name']+' ['+hp['units']+']',doLeft=True)
Ax4.plot(utall[1:],ijh['NORTH'][1:])
Ax4.plot(utall[1:],ijh['SOUTH'][1:])
kv.SetAxLabs(Ax4,"UT",ijh['name']+' ['+ijh['units']+']',doLeft=False)
kv.SetAxDate(Ax4)
Ax1.legend(hemispheres,loc='best')
Ax1.set_title(remixFile)
plt.subplots_adjust(hspace=0)
fn = os.path.join(ptag,'remixTimeSeries.png')
kv.savePic(fn)

#Save the results to python pickle file
fn = os.path.join(ptag,'remixTimeSeries.pkl')
fh = open(fn,'wb')
pickle.dump([cpcp,ipfac,hp,ijh,utall],fh)
fh.close()

#As well as a HDF5 file
with h5py.File(os.path.join(ptag,'remixTimeSeries.h5'),'w') as f:
	dset = f.create_dataset('MJD',data=mjd)
	for hemi in hemispheres:
		dset = f.create_dataset('cpcp'+hemi,data=cpcp[hemi])
		dset.attrs['units'] = cpcp['units']
		dset.attrs['name'] = cpcp['name']
		dset = f.create_dataset('hp'+hemi,data=hp[hemi])
		dset.attrs['units'] = hp['units']
		dset.attrs['name'] = hp['name']
		dset = f.create_dataset('ipfac'+hemi,data=ipfac[hemi])
		dset.attrs['units'] = ipfac['units']
		dset.attrs['name'] = ipfac['name']
		dset = f.create_dataset('ijh'+hemi,data=ijh[hemi])
		dset.attrs['units'] = ijh['units']
		dset.attrs['name'] = ijh['name']





