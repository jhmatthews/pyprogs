########################################################################
#
#	University of Southampton, James Matthews, 130625
#
#	"extrapolate_sub.py"
#		
#	subroutines for extrapolate.py. individually described
#
#########################################################################

## import modules
from math import sqrt, fabs
from scipy.optimize import curve_fit
import scipy.optimize as optimization
import numpy as np

############################################################################
## DOFIT FUNCTION


def dofit(fmode, x, y, nmore, MAX_EXTRAP, n_p):


	'''fitting function for XSs
	
	extrapoltes data in numpy arrays (x, y)
	using functional form of func to maximum 
	of x=MAX_EXTRAP, using nmore extra points
	n_p is the entry in the array we start at
	specifically used for topbase photoionization
	XSs extrapolation
	''' 

	xmin, xmax = x[0], x[n_p-1]

	ymax = y[n_p-1]

	dx = (MAX_EXTRAP - xmax) / nmore

	popt, pcov = optimization.curve_fit(func, x, y, x0=None, sigma=None)

	for i_x in range(1, nmore):

			dx = np.log10( ( MAX_EXTRAP - Emax) ) / nmore 
			X_val = xmax + 10.0**( i_x * dx) 
			temp_energy=np.append ( x, X_val)
			
			if fmode==1:	
				yfit =  func( E, popt[0], popt[1])
			if fmode==2:
				yfit =  func( E, popt[0], popt[1], popt[2])
			if fmode==3:
				yfit =  func( E, popt[0], popt[1], popt[2], popt[3])
			if fmode==4:
				yfit =  func( E, popt[0], popt[1], popt[2], popt[3], popt[4])
			if fmode==5:
				yfit =  func( E, popt[0], popt[1])

			yfit =  func( E, popt)

			y = np.append( y, yfit )

	return x, y



############################################################################
## CHOOSEMODE FUNCTION



def choose_mode(arg_array):

	Plot = False
	fmode=1
	MAX_ENERGY=10000

	if len(arg_array)>1:
		for i in range( len(arg_array) ):
			if arg_array[i] == '-e':
				MAX_ENERGY = float ( arg_array[ i + 1 ] )
			if arg_array[i] == '-f':
				fmode= int ( arg_array[ i + 1 ] )
			if arg_array[i] == 'plot':
				Plot = True
	else:
		fmode=1
		MAX_ENERGY=10000

	return fmode, MAX_ENERGY, Plot


############################################################################
## MINIMISATION FUNCTION

def gmin(f,a,c,tol=3.0e-8):
    """Golden section minimisation.

    Find minimum of f(x) known to have a single minimum
    between a and c, using golden section search. Returns
    minimum position, xmin, and minimum value, f(xmin).
    Suggest set tol to 1.0e-4 for single precision, 3.0e-8
    for double precision."""

    if (c<=a):
        raise ValueError, 'gmin(f,a,c) requires c>a'

    w=(3.0-sqrt(5.0))/2.0 # golden ratio

    b=a+w*(c-a) # set initial abscissas in interval

    bp=c-w*(c-a)

    while(fabs(c-a)>tol*(c+a)/2.0): # loop to refine
        if(f(b)<f(bp)):            # bracketing of minimum
            c=bp
            bp=b
            b=a+w*(c-a)
        else:
            a=b
            b=bp
            bp=c-w*(c-a)

    # finished looping, return best answer
    fb=f(b)
    fbp=f(bp)

    if (fb<fbp): return b, fb
    else: return bp, fbp


############################################################################

def read_topbase(filename):
	'''
	read in XS info from Topbase XS data in Python format

	:INPUT:
		filename 		string
						atomic data filename e.g. topbase_h1_phot.py
	:OUTPUT:
		top 			topbase class instance
						topbase class instance containing information 
						for this filename
	'''

	# read in summary records
	Z,ion,islp,E0,l,num_records = sum_records = np.loadtxt( filename, 
			dtype={'names': ('Z', 'ion', 'islp', 'l', 'E0', 'np'), 
			'formats': ('i4', 'i4', 'i4', 'i4', 'float', 'i4')}, 
                        comments='PhotTop ', delimiter=None, converters=None, 
                        skiprows=0, usecols=(1,2,3,4,5,6), unpack=True, ndmin=0)
	
	# then read the actual cross sections
	energy, XS = np.loadtxt(filename, dtype='float', 
                        comments='PhotTopS', delimiter=None, converters=None, 
                        skiprows=0, usecols=(1,2), unpack=True, ndmin=0)

	# create blank array to store topbase class instances
	top = np.ndarray( len(Z),dtype=np.object)

	nline = 0
	
	for i in range(len(top)):

		n_p = int(num_records[i])
		nmax = nline + n_p
		
		top[i] = topbase_class (Z[i], int(ion[i]), int(islp[i]), int(l[i]), E0[i], n_p, energy[nline:nmax], XS[nline:nmax])

		nline = nmax
		
	# top is an array of class instances like the topbase_ptr in PYTHONRT
	return top


def write_topbase(top, filename):

	'''write array of topbase class instances to file'''

	import numpy as np

	file_write = open( filename, 'w')

	for i in range(len(top)):

		## write the summary records
		file_write.write('PhotTopS  %1d  %1d %3d    %1d     %.6f  %2d\n' %
                                    ( top[i].Z, top[i].ion, top[i].islp, top[i].l, top[i].E0, top[i].np ))

		n_p = top[i].np

		## write the actual XSs
		for j in range( n_p ):
			
			file_write.write('PhotTop     %.6f %.3e\n' % ( top[i].energy[j], top[i].XS[j]) )


	file_write.close()
		
	return 0


class topbase_class:
	'''
	This is a class for topbase photoionization data

	Z 		atomic number
	ion 	ion stage
	islp 	islp number (2s+1, l, parity)
	E0 		threshold energy eV 
	l 		level
	np 		number of entries
	energy 	energies 
	XS 		cross sections
	'''	

	def __init__(self, nz, ne, islp_init, E0_init, linit, np_init, energies, cross_sections):
		self.Z = nz
		self.ion = ne
		self.islp = islp_init
		self.l = linit
		self.E0 = E0_init 
		self.np = np_init
		self.energy = energies
		self.XS = cross_sections







def get_odd_XS():
	'''
	get_odd_XS is a real kluge, it just returns an array of XS identified by eye
	as being odd. These XS are then fudged with a nu**-3 extrapolation.

	:OUTPUT:

		XS_array 		topbase class instance
						array of XS to fudge in topbase_class format

	'''



	# Manually identified list of unusual / anomalous XSections
	XSlist = ['6 2 220', '7 3  220', '7 4 121', '7 6 311', \
	'8 3 341', '8 4 211', '8 4 421', \
	'6 2 231', '7 3 221', '7 4 321', \
	'8 3 141', '8 3 350', '8 4 220', \
	'8 5 321', '6 5 111', '7 3 231', \
	'7 5 200', '8 3 150', '8 3 351', \
	'8 4 231', '8 6 200', '6 5 311', \
	'7 3 411', '7 6 111', '8 3 151', \
	'8 4 200', '8 4 411', '8 7 111']

	import numpy as np

	# number of odd XSections
	n_odd = len(XSlist)

	# create empty object array to store topbase class instances
	XS_array = np.ndarray( n_odd,dtype=np.object)

	Z, ion, islp = [],[],[]

	# cycle over each XS and place in class instance
	for i in range(n_odd):


		# split string
		data = XSlist[i].split()

		
		# get information
		Z.append (int( data[0]))
		ion.append( int(data[1]))
		islp.append(int(data[2]))

	
		# get topbase class instance
	XS_array = topbase_class (np.array(Z), np.array(ion), np.array(islp), 0, 0, 0, [], [])


	# all done, return array
	return XS_array



def check_odd_XS(XS_array, top_record):
	'''
	Routine which checks if the top_record class in question
	is in the flagged list

	:INPUT:
		XS_array 		topbase class instance
						lists of unusual XSections

		top_record 		topbase class instance
						record to compare

	:OUTPUT:
		totalmatch		Bool
						1/True if in list, 0/False if not
	'''

	import numpy as np

	# create boolean arrays for Z, islp and ion
	Zmatch = ( XS_array.Z == top_record.Z )
	ionmatch = ( XS_array.ion == top_record.ion )
	islpmatch = ( XS_array.islp == top_record.islp )

	# create total boolean array (note * means 'and')
	totalmatch = np.sum(ionmatch * islpmatch * Zmatch)


	# error check
	if totalmatch != 1 and totalmatch!=0:
		print "ERROR: Does not equal 1 or zero, exiting"
		print totalmatch
		sys.exit()

	return totalmatch






