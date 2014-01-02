#!/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''
	chianti_sub.py

Subroutines for create_he_chianti.py
'''
import numpy as np
import sys

ANGSTROM  = 1.e-8
C  = 2.997925e10
RYDBERG = 13.60569253		# rydberg in eV
HEV = 4.13620e-15			# Planck's constant in eV 



####################################################
'''
CHIANTI ROUTINES
'''
####################################################



def read_chianti_data ( level_filename="h_1.elvlc", radiative_filename="h_1.wgfa"):

	'''
	read in a chianti atomic database and place in chianti class instance
	Default filenames are for hydrogen
	
	:INPUT
		level_filename		string
							filename for level information
		radiative_filename	string
							filename for radiative information
							
	:OUTPUT
		level	object array
				array of chianti_level class instances
		rad		object array
				array of chianti_rad class instances
	'''
	
	# read in data, place in array
	level_array_read = np.loadtxt (level_filename, comments = "#", dtype = "string")
	rad_array_read = np.loadtxt (radiative_filename, comments = "#", dtype = "string")

	# create blank arrays of type object to store class instances
	level = np.ndarray( len(level_array_read),dtype=np.object)
	rad = np.ndarray( len(rad_array_read),dtype=np.object)
	
	# create level class instances
	for i in range(len(level_array_read)):
	
		index = int (level_array_read[i][0])
		config = int (level_array_read[i][1])
		notation = str (level_array_read[i][2])
		spin = int (level_array_read[i][3]) 
		l = int (level_array_read[i][4])
		l_symbol = str (level_array_read[i][5])
		j = float (level_array_read[i][6])
		multiplicity = int (level_array_read[i][7])
		E_obs = float (level_array_read[i][8])
		E_obs2 = float (level_array_read[i][9])
		E_th = float (level_array_read[i][10])
		E_th2 = float (level_array_read[i][11])
		n = int (notation[0])
		
		level[i] = chianti_level (index,  config, notation, spin, l, l_symbol, 
		                              j, multiplicity, E_obs, E_obs2, E_th, E_th2, n)
	
	# create wgfa class instances
	for i in range(len(rad_array_read)):
		ll = int (rad_array_read[i][0])
		lu = int (rad_array_read[i][1])
		wave = ANGSTROM * float (rad_array_read[i][2])
		if wave!=0:
			freq = C / ( wave ) 
		else:
			freq = 0
		osc = float (rad_array_read[i][3])
		A = float (rad_array_read[i][4])
		note_low = str(rad_array_read[i][5])
		note_up = str(rad_array_read[i][8])
		J_low = float(rad_array_read[i][6][2:])
		J_up = float(rad_array_read[i][9][2:])
		
		rad[i] = chianti_rad(ll, lu, wave, freq, osc, A, note_low, note_up, J_low, J_up)
		
	return level, rad



def chianti_to_lev (level_info, E_thres, Z, ion):
	'''
	create a level class instance for writing to file 
	from chianti data 

	:INPUT:
		level_info 			array like 
							array of chianti level info class instances

		E_thres 			float 
							ionization energy of ion 

	:OUTPUT:
		levels 				array like 
							array of level class instances 


	'''

	nlevels = len(level_info)

	if Z == ion:
		cont  = 1
	else:
		cont = 0

	levels = np.ndarray (nlevels, dtype = np.object)


	for n in range(0, nlevels + cont):

		if n < nlevels:						# non continuum level
			lvl = n+1

			energy_ev = level_info[n].E_obs2 * RYDBERG

			ionpot = energy_ev - E_thres 

			g = level_info[n].multiplicity

			nstring = level_info[n].notation


		elif n == nlevels:				# then this is the continuum level
			
			ion = ion + 1

			lvl = 1

			energy_ev = E_thres

			ionpot = 0.0

			g = 1

			nstring = "cnt"

		# set rad_rates	
		if n == 0 or n == nlevels + cont:
			rad_rate = 1.00e+21
		else: 
			rad_rate = 1.00e-9


		levels[n] = level (Z, ion, lvl, ionpot, energy_ev, g, rad_rate, "()", nstring)

	return levels 




def chianti_to_line (line_info, chianti_levels,  Z, ion):
	'''
	create a line class instance for writing to file 
	from chianti data 

	:INPUT:
		line_info 			array like 
							array of chianti line info class instances

		E_thres 			float 
							ionization energy of ion 

	:OUTPUT:
		lines 				array like 
							array of line class instances 


	'''

	nlines = len(line_info)

	lines = np.ndarray (nlines, dtype = np.object)

	for i in range(nlines):

		ll, il  = where_split_level(chianti_levels, line_info[i].note_low, line_info[i].J_low)

		lu, iu  = where_split_level(chianti_levels, line_info[i].note_up, line_info[i].J_up)

		Eu = chianti_levels[iu].E_obs2 * RYDBERG
		El = chianti_levels[il].E_obs2 * RYDBERG

		print iu, il, lu, ll 

		gl = 2*chianti_levels[il].J + 1
		gu = 2*chianti_levels[iu].J + 1

		lines[i] = line (Z, ion, line_info[i].wave/ANGSTROM, line_info[i].freq, line_info[i].osc, gl, gu, El, Eu, ll, lu)

	return lines 


def where_split_level(chianti_levels, notation, J):

	nlev = None 
	ifind = None 

	for i in range(len(chianti_levels)):

		if chianti_levels[i].J == J and chianti_levels[i].notation == notation:

			ifind = i
			nlev = i+1

	if nlev == None:
		print "No level found"

	return nlev, ifind 



####################################################
'''
READ AND WRITE PYTHON ROUTINES 
'''
####################################################




def write_line_file(lin, filename):
	'''
	write information from array of line class instances to a filename of format

	LinMacro    1   1   1215.33907 0.41620  2  8     0.00000  10.19883     1     2

	:INPUT:
		line 		array like
					array of line class instances

		filename 	string 
					name of file to write, e.g. he20_levels.py

	'''

	header = '''# data from chianti
# z = element, ion= ionstage, f = osc. str., gl(gu) = stat. we. lower(upper) level
# el(eu) = energy lower(upper) level (eV), ll(lu) = lvl index lower(upper) level
#        z ion       lambda      f         gl  gu    el          eu        ll   lu
#
'''

	out = open(filename, "w")

	out.write(header)

	for i in range(len(lin)):

		l = lin[i]

		out.write("LinMacro   %i  %i   %11.5f  %.5f  %i   %i   %.5f   %.5f  %i  %i\n" %
			       (l.z, l.ion, l.wavelength, l.osc, l.g_l, l.g_u, l.e_l, l.e_u, l.ll, l.lu))

	out.close()

	return (0)

def read_line_info(filename):

	'''read in line info from Python's atomic data e.g. h20_lines'''

	line_array_read = np.loadtxt(filename, comments='#', unpack = True, dtype = 'string')
	
	line_array_read = np.transpose(line_array_read)

	lines = np.ndarray( len(line_array_read),dtype=np.object)
	
	for i in range(len(lines)):
		z = float (line_array_read[i][1])
		ion = float (line_array_read[i][2])
		wave = ANGSTROM * float (line_array_read[i][3])
		freq = C / ( wave ) 
		osc = float (line_array_read[i][4])
		gl = int (line_array_read[i][5])
		gu = int (line_array_read[i][6])
		el = float(line_array_read[i][7])
		eu = float(line_array_read[i][8])
		ll =  int (line_array_read[i][9])
		lu =  int (line_array_read[i][10])
		
		lines[i] = line(z, ion, wave, freq, osc, gl, gu, el, eu, ll, lu)
		
	# line is an array of class instances like the line_ptr in PYTHONRT
	return lines


def write_level_file(lev, filename, append = False):

	'''
	write information from array of level class instances to a filename of format

	LevMacro   2  2   1  -78.73318   0.00000   2   1.00e+21  ()  1s

	:INPUT:
		lev 		array like
					array of level class instances

		filename 	string 
					name of file to write, e.g. he20_levels.py

	'''
	
	header = '''# data from chianti
#         z ion lvl ion_pot   ex_energy   g  rad_rate  
#
'''
	if append:
		out = open(filename, "a")
	else:
		out = open(filename, "w")	

	out.write(header)

	iprev = lev[0].ion

	for i in range(len(lev)):

		l = lev[i]

		if l.ion != iprev:
			out.write("#\n#\n")

		out.write("LevMacro   %i  %i  %3i  %.5f   %.5f   %3i   %8.2e  %s  %s\n" %
			       (l.z, l.ion, l.lvl, l.ionpot, l.E, l.g, l.rad_rate, l.bracks, l.nnstring))

		iprev = l.ion

	out.close()

	return (0)



def read_level_info(filename):

	'''read in level info from Python's atomic data e.g. h20_level'''

	level_array_read = np.loadtxt(filename, comments='#', unpack = True, dtype = 'string')
	
	level_array_read = np.transpose(level_array_read)

	lev = np.ndarray( len (level_array_read), dtype=np.object)
	
	for i in range(len(lev)):
		z = int (level_array_read[i][1])
		ion = int (level_array_read[i][2])
		lvl = int (level_array_read[i][3])
		ionpot = float (level_array_read[i][4]) 
		E = float (level_array_read[i][5])
		g = int (level_array_read[i][6])
		rad_rate = float (level_array_read[i][7])
		brack =  str (level_array_read[i][8])
		nnstring =  str  (level_array_read[i][9])
		
		lev[i] = level(z, ion, lvl, ionpot, E, g, rad_rate, brack, nnstring)
		
	# level is an array of class instances like the line_ptr in PYTHONRT
	return lev


def convert_h_to_hydrogenic_class(lev, Z, ion_threshold):

	'''
	Applies a Z^2 correction to convert hydrogen energy levels 
	into Z atomic number hydrogenic ion levels 

	:INPUT:
		lev 			array like
						array of level class instances

		Z 				int
						atomic number 

		ion_threshold	float 
						ionization energy to continuum
	'''

	for i in range(len(lev)):

		lev[i].ionpot *= Z*Z

		lev[i].E = ion_threshold + lev[i].ionpot

		lev[i].Z = Z 

		lev[i].ion = Z 

	return (lev)






# line class: analogous to line ptr in python. contains freq, oscillator strength, 
class line:
	'''This is a class analogous to line ptr in python'''
	def __init__(self, _z, _ion, _wavelength, _freq, _osc, _g_l, _g_u, _ll, _lu):
		self.z = _z
		self.ion = _ion
		self.wavelength = _wavelength
		self.freq = _freq
		self.osc = _osc
		self.g_l = _g_l
		self.g_u = _g_u
		self.ll = _ll
		self.lu = _lu



####################################################
'''
NIST ROUTINES 
'''
####################################################


def read_nist(filename, z, ion):

	nist_read = open (filename, "r")

	nist_levels = []

	lvl = 1

	for line in nist_read:

		data = line.split()

		if len(data)>0:

			if data[0] != "#":

				print data

				if len(data) == 5:

					note = data[0] 

					j = int (data[2])

					E = float (data[3])

					g = 2*j + 1

					note_prev = note



				elif len(data) == 3:

					note = note_prev

					j = int (data[0])

					E = float (data[1])

					g = 2*j + 1


				nist_levels.append ( nist (z, ion, lvl, note, j, E, g) )

				lvl += 1


	nist_read.close()			
	nist_levels = np.array(nist_levels)


	return nist_levels



def nist_to_lev(nistclass, threshold_energy):

	nlevels = len(nistclass)

	lev = np.ndarray( nlevels,dtype=np.object)

	for i in range(nlevels):

		nnstring = nistclass[i].notation
		bracks = "()"

		E = nistclass[i].E

		ionpot = E - threshold_energy 

		rad_rate = 1.00e-09

		lev[i] = level (nistclass[i].z, nistclass[i].ion, nistclass[i].lvl, ionpot, E, nistclass[i].g, rad_rate, bracks, nnstring)


	return lev



def compress_levels (nistclass):

	'''compresses levels in NIST data

	'''

	VERY_BIG = 1e50

	FINE_ENERGY = 1e-10


	nlevels_original = len(nistclass)

	fine_energy = 0.1

	note_prev, conf_prefix_old  = nistclass[0].notation, nistclass[0].notation[:-1]

	gsum, Esum = nistclass[0].g, nistclass[0].E* nistclass[0].g

	nistclass_new = []		# this is where we store the class instances

	lvl = 1

	Eold = -VERY_BIG


 
	# too compress the levels, we can choose to compress them by
	# a few criteria
	# 1 - if they have the same n,l but different J
	# 2 - if they have the same n, but different J and l
	# 3 - i

	for i in range(nlevels_original):

		conf_prefix = nistclass[i].notation[:-1]

		delta_E = nistclass[i].E - Eold 

		print conf_prefix 

		close = ( delta_E < FINE_ENERGY and conf_prefix == conf_prefix_old )

		if nistclass[i].notation == note_prev or close == True :

			gsum += nistclass[i].g 

			Esum += nistclass[i].E * nistclass[i].g

			z = nistclass[i].z

			ion = nistclass[i].z

			print nistclass[i].notation, i, Esum, gsum 


		else:

 
			Enew = Esum / gsum 


			nistclass_new.append ( nist (nistclass[i-1].z, nistclass[i-1].ion, lvl, nistclass[i-1].notation, "#", Enew, gsum) )

			gsum, Esum = nistclass[i].g, nistclass[i].E * nistclass[i].g

			lvl += 1



		note_prev =  nistclass[i].notation

		conf_prefix_old = conf_prefix

		Eold = nistclass[i].E


	nistclass_new = np.array (nistclass_new)

	print lvl 

	return nistclass_new

def read_topbase_levels(filename):

	'''read in level info from original topbase levels and 
	place in class instance array'''

	level_array_read = np.loadtxt(filename, comments='#', unpack = True, dtype = 'string')
	
	level_array_read = np.transpose(level_array_read)

	lev = np.ndarray( len (level_array_read), dtype=np.object)
	
	for i in range(len(lev)):

		ne = int(level_array_read[i][2])
		z = int (level_array_read[i][1])
		ion = int (level_array_read[i][0])

		lvl = int (level_array_read[i][3])
		
		E = float (level_array_read[i][-4])
		g = int (level_array_read[i][-2])
		rad_rate = float (level_array_read[i][-1])

		nnstring =  str  (level_array_read[i][-5])
		
		lev[i] = level(z, ion, lvl, 0.0, E, g, rad_rate, "", nnstring)
		
	# level is an array of class instances like the line_ptr in PYTHONRT
	return lev


def read_topbase_xs(filename, suffix="Mac"):

	'''this needs to be done manually, annoyingly, as format of
	topbase has no PhotTop or PhotMac string'''

	f = open(filename, "r")

	i = -1

	top = []

	for line in f:

		data = line.split()

		if len(data) == 2:
			top_dummy.energy[j] = float(data[0]) * RYDBERG
			top_dummy.XS[j] = 1e-18 * float(data[1])

			if j > top_dummy.np:
				print "error, j > top[i].np"

			j += 1

		elif len(data) == 7:
			i += 1 
			j = 0

			if i>0: top.append(top_dummy)


			Z = int(data[1])
			ion = int(data[2])
			islp = int(data[3])
			l = int(data[4])
			E = -1 * float(data[5]) * RYDBERG
			num = int(data[6])

			top_dummy = topbase_class (Z, ion, islp, l, E, num, np.ndarray( num, dtype=np.object), np.ndarray( num, dtype=np.object))

		else: 
			print "Didn't understand data, exiting"
			print line
			sys.exit()

	top = np.array(top)

	return top




def write_topbase_xs(top, filename, suffix="Mac"):

	'''write array of topbase class instances to file'''


	file_write = open( filename, 'w')

	for i in range(len(top)):

		## write the summary records
		file_write.write('Phot%sS  %1d  %1d %3d    %1d     %.6f  %2d\n' %
                                    ( suffix, top[i].Z, top[i].ion, top[i].islp, top[i].l, top[i].E0, top[i].np ))

		n_p = top[i].np

		## write the actual XSs
		for j in range( n_p ):
			
			file_write.write('Phot%s     %.6f %.3e\n' % (suffix,  top[i].energy[j], top[i].XS[j]) )


	file_write.close()
		
	return 0






####################################################
'''
CLASSES
'''
####################################################




class nist:
	'''This is a class to store nist levels '''
	def __init__(self, _z, _ion, _lvl, _notation, _J, _E, _g):
		self.z = _z
		self.ion = _ion
		self.lvl = _lvl
		self.notation = _notation
		self.J = _J
		self.E = _E
		self.g = _g



class nist_line:
	'''this is a class to store nist_line_data, unfinished'''
	def __init__(self, _z, _ion, _lvl, _notation, _J, _E, _g):
		self.z = _z
		self.ion = _ion
		self.lvl = _lvl
		self.notation_upper = _notation
		self.J = _J
		self.E = _E
		self.g = _g
		self.z = _z
		self.ion = _ion
		self.wavelength = _wavelength
		self.freq = _freq
		self.osc = _osc
		self.g_l = _g_l
		self.g_u = _g_u
		self.ll = _ll
		self.lu = _lu



# line class: analogous to line ptr in python. contains freq, oscillator strength, 
class line:
	'''This is a class analogous to line ptr in python'''
	def __init__(self, _z, _ion, _wavelength, _freq, _osc, _g_l, _g_u, _e_l, _e_u, _ll, _lu):
		self.z = _z
		self.ion = _ion
		self.wavelength = _wavelength
		self.freq = _freq
		self.osc = _osc
		self.g_l = _g_l
		self.g_u = _g_u
		self.e_l = _e_l 
		self.e_u = _e_u 
		self.ll = _ll
		self.lu = _lu
		
		
		
# line class: analogous to line ptr in python. contains freq, oscillator strength, 
class level:
	'''Stores information from a Python levels file'''
	def __init__(self, _z, _ion, _lvl, _ionpot, _E, _g, _rad_rate, _bracks, _nnstring):
		self.z = _z
		self.ion = _ion
		self.lvl = _lvl
		self.ionpot = _ionpot
		self.E = _E
		self.g = _g
		self.rad_rate = _rad_rate
		self.bracks = _bracks
		self.nnstring = _nnstring
		
class chianti_level:
	'''
	Stores chianti level data from file format .elvlc
	See http://www.chiantidatabase.org/cug.pdf Page 8 for description
	'''
	def __init__(self, _index,  _config, _notation, _spin, _l, _l_symbol, _j, _multiplicity, _E_obs, _E_obs2, _E_th, _E_th2, _n):
		self.index = _index
		self.config = _config
		self.notation = _notation
		self.spin = _spin
		self.l = _l
		self.l_symbol = _l_symbol 
		self.J = _j
		self.multiplicity = _multiplicity
		self.E_obs = _E_obs
		self.E_obs2 = _E_obs2
		self.E_th = _E_th
		self.E_th2 = _E_th2
		self.n = _n
		
class chianti_rad:
	'''
	Stores radiative chianti information from wgfa file.
	Contains the wavelengths, gf and A values of the transitions and the indices initial and
	final level corresponding to the indices of the levels as given in the .elvlc file
	See http://www.chiantidatabase.org/cug.pdf Page 9 for description
	'''
	def __init__(self, _ll, _lu, _wave, _freq, _osc, _A, _note_low, _note_up, _J_low, _J_up):
		self.ll = _ll
		self.lu = _lu
		self.wave = _wave
		self.freq = _freq
		self.osc = _osc
		self.A = _A
		self.note_low = _note_low
		self.note_up = _note_up
		self.J_low = _J_low
		self.J_up = _J_up



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

	def __init__(self, nz, ne, islp_init, linit, E0_init, np_init, energies, cross_sections):
		self.Z = nz
		self.ion = ne
		self.islp = islp_init
		self.l = linit
		self.E0 = E0_init 
		self.np = np_init
		self.energy = energies
		self.XS = cross_sections