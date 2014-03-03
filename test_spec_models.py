'''
test_spec_models contains routines for analysing the in cell spectrum 

Using a wind_save file which has been run in ionization mode 7 
(this information will even be saved in other ionization modes)
then you can look at what 

just run the following 
>> ipython
>> run test_spec_models.py 
>> run_py_wind("77", "filename")
>> xbands = get_xbands(filename, NPLASMA, NBANDS)
>> plot_cell(xbands, cell_nos)
'''


import csv, sys, os, array, warnings, subprocess
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pywind_sub as ps
import time
from constants import H_OVER_K



def run_py_wind (vers, fname):
	'''
	run version vers of py_wind on file fname.wind_save
	'''
	isys = os.system('py_wind'+vers+' '+fname+' < band_cmds > tempfile &')
	time.sleep(3)
	return isys


#ix,iz,x,z,lx,lz,xmax,zmax=ps.get_wind_geom(fname+'.ioncH1.dat')


#nplasmas = ps.pywind_read(fname+'.pnum.dat',ix,iz)


def read_pywind_smart(filename, return_inwind=False):
	'''
	read a py_wind file using np array reshaping and manipulation
	'''

	# first, simply load the filename 
	d = np.loadtxt(filename, comments="#", dtype = "float", unpack = True)

	# our indicies are already stored in the file- we will reshape them in a sec
	zindices = d[-1]
	xindices = d[-2]

	# we get the grid size by finding the maximum in the indicies list 99 => 100 size grid
	zshape = int(np.max(zindices) + 1)
	xshape = int(np.max(zindices) + 1)


	# reshape our indices arrays
	xindices = xindices.reshape(xshape, zshape)
	zindices = zindices.reshape(xshape, zshape)

	# now reshape our x,z and value arrays
	x = d[0].reshape(xshape, zshape)
	z = d[1].reshape(xshape, zshape)

	values = d[2].reshape(xshape, zshape)

	# these are the values of inwind PYTHON spits out
	inwind = d[3].reshape(xshape, zshape)

	# create an inwind boolean to use to create mask
	inwind_bool = (inwind >= 0)
	mask = (inwind < 0)

	# finally we have our mask, so create the masked array
	masked_values = np.ma.masked_where ( mask, values )

	#print xshape, zshape, masked_values.shape

	#return the transpose for contour plots.
	if return_inwind:
		return x, z, masked_values.T, inwind_bool.T
	else:
		return x, z, masked_values.T


def read_pywind_logspace(filename):

	'''
	convert to logspace coordinates
	'''

	x,z,s = read_pywind_smart(filename)

	x = np.log10(x)
	z = np.log10(z)

	return x,z,s




class xband:
	'''
	a class for storing banded estimators from a cell
	'''
	def __init__(self, nxtot, xj, ave_freq, specmod, alpha, pl_log_w, exp_w, exp_t, fmin, fmax, J_array):
		self.ave_freq = ave_freq
		self.nxtot = nxtot
		self.xj = xj
		self.specmod = specmod
		self.alpha = alpha
		self.pl_log_w = pl_log_w
		self.exp_w = exp_w
		self.exp_t = exp_t
		self.fmin = fmin
		self.fmax = fmax
		self.J = J_array




def f_exp(T, W, nu):

	'''
	exponential spectral model
	'''

	return W * np.exp(- H_OVER_K * nu / T)

def log_pl(alpha, logW, nu):
	'''
	Log PL spectral model
	'''

	return logW + alpha*np.log10(nu)

def pl(alpha, logW, nu):
	'''
	PL spectral model
	'''
	return (10.0**logW) * (nu ** alpha)







def get_xbands(filename, NPLASMA, NBANDS):

	'''
	array of class instances with shape of grid read in
	index as 
	'''

	x, y, nplasmas = read_pywind_smart( "%s.pnum.dat" % (fname) )

	xshape = len(x)
	yshape = len(y)

	ix, iy = np.indices ((xshape, yshape))

	xbands = np.ndarray ( (NBANDS), dtype = np.object )

	x1, y1, pnum = read_pywind_smart( "%s.pnum.dat" % (fname) )

	print np.max(pnum)

	for i in range(NBANDS):

		x1, y1, _nxtot = read_pywind_smart( "%s.nxtot00%i.dat" % (fname, i) )

		x2, y2, _xj = read_pywind_smart( "%s.xj00%i.dat" % (fname, i) )

		x3, y3, _ave_freq, inwind = read_pywind_smart( "%s.xave_freq00%i.dat" % (fname, i), return_inwind=True )

		x1, y1, _specmod = read_pywind_smart( "%s.spec_mod_type00%i.dat" % (fname, i) )

		x1, y1, _exp_t = read_pywind_smart( "%s.exp_temp00%i.dat" % (fname, i) )

		x1, y1, _exp_w = read_pywind_smart( "%s.exp_w00%i.dat" % (fname, i) )

		x1, y1, _pl_log_w = read_pywind_smart( "%s.pl_log_w00%i.dat" % (fname, i) )

		x1, y1, _pl_alpha = read_pywind_smart( "%s.pl_alpha00%i.dat" % (fname, i) )

		x1, y1, _fmax = read_pywind_smart( "%s.fmax_mod00%i.dat" % (fname, i) )

		x1, y1, _fmin = read_pywind_smart( "%s.fmin_mod00%i.dat" % (fname, i) )




		nxtot = _nxtot[inwind]
		ave_freq = _ave_freq[inwind]
		xj = _xj[inwind]
		specmod = _specmod[inwind]
		exp_t = _exp_t[inwind]
		exp_w = _exp_w [inwind]
		pl_log_w = _pl_log_w[inwind]
		pl_alpha = _pl_alpha[inwind]
		fmax = _fmax[inwind]
		fmin = _fmin[inwind]

		#print nxtot.shape, _nxtot.shape

		p = pnum[inwind]

		#for j in range(len(p)):
		#	print p[j], j, nxtot[j]

		#print nxtot



		xbands[i] = xband ( nxtot, ave_freq, xj, specmod, pl_alpha, pl_log_w, exp_w, exp_t, fmin, fmax, np.zeros(len(nxtot)))

		print xbands[i].nxtot.shape

	return xbands





def cell_band(xbands, nplasma, log_f = True, fmin = 1e14, fmax = 1e20):

	'''
	get the in cell spectrum for a given cell
	returns continuous array
	'''

	nbands = len(xbands)

	fmin = xbands[0].fmin[nplasma]
	fmax = xbands[nbands-1].fmax[nplasma]

	log_fmin = np.log10(fmin)
	log_fmax = np.log10(fmax)

	n_freqs = 10000

	freqs =  np.logspace(log_fmin, log_fmax, num = n_freqs)


	J_nu = np.zeros(n_freqs)

	for i in range(nbands):

		specmod = xbands[i].specmod[nplasma]

		fmin = xbands[i].fmin[nplasma]
		fmax = xbands[i].fmax[nplasma]

		freqs_low = (freqs < xbands[i].fmax[nplasma])
		freqs_high = (freqs > xbands[i].fmin[nplasma])

		todo = freqs_low * freqs_high
		

		if specmod == 1:
			fnu = todo * 10.0**log_pl(xbands[i].alpha[nplasma], xbands[i].pl_log_w[nplasma], freqs) 
			print 'parameters:', xbands[i].alpha[nplasma], xbands[i].pl_log_w[nplasma]

		elif specmod == 2:
			fnu = todo * f_exp(xbands[i].exp_t[nplasma], xbands[i].exp_w[nplasma], freqs)
			print 'parameters:', xbands[i].exp_t[nplasma], xbands[i].exp_w[nplasma]

		else: 
			fnu = 0.0* freqs

		# convert nan to zero
		fnu[np.isnan(fnu)] = 0.0

		print "\n\nBAND %i" % (i+1)
		print "%8.4e to %8.4e" % (xbands[i].fmin[nplasma], xbands[i].fmax[nplasma])
		print "spectral model %i" % specmod
		print fnu, xbands[i].nxtot[nplasma]

		J_nu += fnu

	return freqs, J_nu




def plot_cell(xbands, nplasma, lstyle = "-"):

	'''
	plot the spectrum in nplasma (tuple)
	'''

	import pylab

	colors = ["c", "m", "y", "k", "r", "b", "g"]

	for n in range(len(nplasma)):
		f, j = cell_band(xbands, nplasma[n])

		pylab.loglog(f,j,  lstyle)

	pylab.xlabel("log Frequency")
	pylab.ylabel("log Jnu")

	pylab.show()










