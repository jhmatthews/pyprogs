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
import classes as cls

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

	if len(arg_array)>1:
		for i in range( len(arg_array) ):
			if arg_array[i] == '-e':
				MAX_ENERGY = float ( arg_array[ i + 1 ] )
			if arg_array[i] == '-f':
				fmode= int ( arg_array[ i + 1 ] )
	else:
		fmode=1
		MAX_ENERGY=10000

	return fmode, MAX_ENERGY


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


def get_odd_XS():
	'''
	get_odd_XS is a real kluge, it just returns an array of XS identified by eye
	as being odd. These XS are then fudged with a nu**-3 extrapolation.

	Returns:

		XS_array 		object array
						array of XS to fudge in topbase_class format

	'''



	# Manually identified list of unusual / anomalous XSections
	XSlist = ["6 2 231", "7 4 121", "7 3 411", \
	          "7 4 321", "8 3 150", "8 3 141", \
	          "8 3 151", "8 3 341", "8 3 351", \
	          "8 3 350", "8 4 211", "8 4 220", \
	          "8 4 221", "8 5 321", "8 5 121"]

	# number of odd XSections
	n_odd = len(XSlist)

	# create empty object array to store topbase class instances
	XS_array = np.ndarray( n_odd,dtype=np.object)

	# cycle over each XS and place in class instance
	for i in range(n_odd):


		# split string
		data = XSlist[i].split()

		
		# get information
		Z = data[0]
		ion = data[0]
		islp = data[0]

		
		# get topbase class instance
		XS_array[i] = cls.topbase_class (Z, ion, islp, 0, 0, 0, []], []])


	# all done, return array
	return XS_array












