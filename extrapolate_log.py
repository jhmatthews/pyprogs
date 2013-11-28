########################################################################
#
#	University of Southampton, James Matthews, 130625
#
#	"extrapolate_log.py"
#		
#	This code is intended to extrapolate topbase to higher energies
#
#########################################################################

'''
The basic procedure for this code is as follows:
1. list the filenames you want to create
2. work out a fit to the topbase date
3. extrapolate this fit to a higher energy

usage:
	python extrapolate_log.py [-e max_energy] [-f function mode]
returns:
	topbase files in form topbase_h1_phot_extrap.py
'''
	
## import modules
import os, sys
import numpy as np
import pylab
from constants import *
import extrapolate_sub as sub
#from scipy.optimize import curve_fit
#import scipy.optimize as optimization


## we read the functional form and the maximum energy from the command line
## they are set to defaults if not read
func, NEW_MAX_ENERGY= sub.choose_mode(sys.argv)



## open the standard73 atomic data file in order to read topbase filenames
atomic_filename='data/standard73'

## read standard73 into array
atomic_data_files=np.loadtxt(atomic_filename, dtype='string', comments='#', delimiter=None, converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0)


## write topbase filenames to an array
print '''Your filenames are:
	Old       |   	    New
------------------------------------------------'''
topbase_data_files, new_data_files=[],[]
for i in range(len(atomic_data_files)):
	oldfile=atomic_data_files[i]

	## is it a topbase photoionization record?
	if 'topbase' in oldfile and 'phot' in oldfile:

		## if so, create new filename
		newfile=atomic_data_files[i][14:-3]+'_extrap.py'
		topbase_data_files.append(oldfile)
		new_data_files.append(newfile)
		print oldfile[14:],' | ',newfile




## we now have 2 arrays of old and new filenames
## let's do the fit

n_files=len(topbase_data_files)

## quick check
if n_files!=len(new_data_files): 
	print 'Error: differing array lengths'
	sys.exit(0)


## define our new maximum energy in eV
NEW_MAX_ENERGY=100000.0	#use 100KeV for the moment


for i_file in range(n_files):

	## do each filename in turn
	filename_read=topbase_data_files[i_file]
	file_write=open(new_data_files[i_file], 'w')
	sys.stderr.write('Reading '+filename_read+'...')

	## first read the summary records into a set of arrays
	Z,ion,islp,l,E0,num_records = sum_records = np.loadtxt( filename_read, 
			dtype={'names': ('Z', 'ion', 'islp', 'l', 'E0', 'np'), 
			'formats': ('i4', 'i4', 'i4', 'i4', 'float', 'i4')}, 
                        comments='PhotTop ', delimiter=None, converters=None, 
                        skiprows=0, usecols=(1,2,3,4,5,6), unpack=True, ndmin=0)
	
	## then read the actual cross sections
	energy, XS = np.loadtxt(filename_read, dtype='float', 
                        comments='PhotTopS', delimiter=None, converters=None, 
                        skiprows=0, usecols=(1,2), unpack=True, ndmin=0)


	sys.stderr.write('done.\n')
	sys.stderr.write('Extrapolating data...')


	## Now do the fit

	nline=0	## counter for records

	energy_new=[]
	XS_new=[]
	nentries=[]
	
	## we now iterate 
	
	for i in range(len(Z)):

		n_p=int(num_records[i])
		nmore=50	## the number of entries we extrapolate by

		## create temporary arrays based on individual record
		temp_energy = energy[nline:nline+n_p]
		temp_XS = XS[nline:nline+n_p]
		temp_energy_comp = energy[nline:nline+n_p]
		temp_XS_comp  = XS[nline:nline+n_p]
		
		
		Emin, Emax = temp_energy[0], temp_energy[n_p-1]
		XSmax=temp_XS[n_p-1]
		#coefficient = XSmax / (Emax ** -3) 
		if NEW_MAX_ENERGY < Emax: 
			NEW_MAX_ENERGY=Emax 
			nmore=0
		dE = (NEW_MAX_ENERGY - Emax) / nmore
		dlogE = np.log10( (NEW_MAX_ENERGY - Emax) ) / nmore 


		ldx = np.log10 ( temp_energy[n_p-1] ) - np.log10 (temp_energy[n_p-2] )
                ldy = np.log10 ( temp_XS[n_p-1] ) - np.log10 (temp_XS[n_p-2] )
		grad = ldy / ldx

		logXSmax = np.log10( XSmax ) 
		logEmax = np.log10( Emax ) 
		
		nentries.append(nmore+n_p)
		## now extrapolate using gradient in log log space at last 2 points in topbase entry
		for i_energy in range(1, nmore+1):

			E = Emax + 10.0**( i_energy * dlogE) 

			logXSfit =  logXSmax + ( grad * (np.log10(E) - logEmax) )
			XSfit= 10.0 ** logXSfit
			
			## append the values to temporary arrays for this topbase record
			temp_XS = np.append ( temp_XS, XSfit )
			temp_energy = np.append ( temp_energy, E )


		## plot up to check fits
		if i_file==3 and Z[i] == 8 and ion[i] == 5:
			pylab.loglog(temp_energy_comp, temp_XS_comp, 'k--')
			pylab.loglog(temp_energy, temp_XS, c='g')
			pylab.xlabel('Energy eV')
			pylab.ylabel('XS cm^-2')
			print 'summary', Z[i],ion[i],islp[i],l[i],E0[i]

		## put this record in the arrays of all records for this file
		energy_new.append(temp_energy)
		XS_new.append(temp_XS)

		nline = nline + n_p

	sys.stderr.write('done.\n')
	pylab.xlim(10,10000)
	pylab.ylim(1e-29, 1e-15)
	pylab.show()
	




	## Finally, write the files- no numpy modules for flexibility

	sys.stderr.write('Writing to file...')

	for i in range(len(Z)):

		## write the summary records
		file_write.write('PhotTopS  %1d  %1d %3d    %1d     %.6f  %2d\n' %
                                    ( Z[i], ion[i], islp[i], l[i], E0[i], nentries[i] ))

		n_p=int(num_records[i])

		## write the actual XSs
		for j in range(0,nentries[i]):
			
			file_write.write('PhotTop     %.6f %.3e\n' % (energy_new[i][j], XS_new[i][j]) )


	file_write.close()

	sys.stderr.write('done.\n')


## now everything is done!
		

	
		
		






