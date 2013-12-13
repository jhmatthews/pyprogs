'''
	University of Southampton -- JM -- 29 November 2913

				py_macro_data.py

Synopsis:
	This program reads in a chianti data file and converts it to the correct form

Usage:
	python py_macro_data.py [chianti_prefix] [python_prefix]

	chianti_prefix 		e.g  h_1
	python prefix 		e.g. h20

Arguments:
	needs the following chianti filenames

	x_x.elvlc
	x_x.wgfa

Returns:
	creates the following python filenames

	x_x_levels.py
	x_x_lines.py
'''

import numpy as np 
import os, sys 



def read_nist_data ( level_filename="he_1_levels.nist", line_filename="he_1_lines.nist"):

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
	level_array_read = np.loadtxt (level_filename, comments = "%", dtype = "string")
	rad_array_read = np.loadtxt (line_filename, comments = "%", dtype = "string")

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
		
		level[i] = cls.chianti_level (index,  config, notation, spin, l, l_symbol, 
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
		
		rad[i] = cls.chianti_rad(ll, lu, wave, freq, osc, A, note_low, note_up, J_low, J_up)
		
	return level, rad


# get prefixes for input and output filenames
chianti_prefix = sys.argv[1]
python_prefix = sys.argv[2]

# get chianti filenames
chianti_l_file = chianti_prefix + ".clvlc"
chianti_r_file = chianti_prefix + ".wgfa"


# read chianti data from file
level_info, rad_info = read_chianti_data ( level_filename = chianti_l_file , radiative_filename = chianti_r_file )



'''
Format for PYTHON level filename:

#         z ion lvl ion_pot   ex_energy   g  rad_rate
LevMacro  1  1  1 -13.59843    0.000000   2  1.00e+21 () n=1

Format for PYTHON lines filename:

# z = element, ion= ionstage, f = osc. str., gl(gu) = stat. we. lower(upper) level
# el(eu) = energy lower(upper) level (eV), ll(lu) = lvl index lower(upper) level
#        z ion       lambda      f         gl  gu    el          eu        ll   lu
LinMacro    1   1          1215.33907             0.41620     2     8             0.00000            10.19883     1     2
LinMacro    1   1          1025.44253             0.07910     2    18             0.00000            12.08750     1     3
LinMacro    1   1           972.27104             0.02899     2    32             0.00000            12.74854     1     4

Need to convert from Chianti format to PYTHON format
'''





















