'''
Contains classes I use as convenient reference points for atomic data
'''




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
		
		
class level:
	'''Stores information from a Python levels file'''
	def __init__(self, _z, _ion, _lvl, _ionpot, _E, _g, _rad_rate, _nnstring, _islp):
		self.z = _z
		self.ion = _ion
		self.lvl = _lvl
		self.ionpot = _ionpot
		self.E = _E
		self.g = _g
		self.rad_rate = _rad_rate
		self.nnstring = _nnstring
		self.islp = _islp

# line class: analogous to line ptr in python. contains freq, oscillator strength, 
class toplevel:
	'''Stores information from a Python levels file'''
	def __init__(self, _z, _ion, _lvl, _ionpot, _E, _g, _rad_rate, _nnstring, _islp, _ilv):
		self.z = _z
		self.ion = _ion
		self.lvl = _lvl
		self.ionpot = _ionpot
		self.E = _E
		self.g = _g
		self.rad_rate = _rad_rate
		self.nnstring = _nnstring
		self.islp = _islp
		self.ilv = _ilv


class topline:
	'''This is a class analogous to line ptr in python'''
	def __init__(self, _z, _ion, _wavelength, _freq, _gf, _osc, _g_l, _g_u, _e_l, _e_u, _ll, _lu, _conf_u, _conf_l, _ulv, _llv, _islpl, _islpu):
		self.z = _z
		self.ion = _ion
		self.wavelength = _wavelength
		self.freq = _freq
		self.gf = _gf
		self.osc = _osc
		self.g_l = _g_l
		self.g_u = _g_u
		self.e_l = _e_l 
		self.e_u = _e_u 
		self.ll = _ll
		self.lu = _lu
		self.conf_u = _conf_u
		self.conf_l = _conf_l
		self.ulv = _ulv
		self.llv = _llv
		self.islpl = _islpl
		self.islpu = _islpu



		
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



class top_photo_mac:
	'''
	a topbase class for writing

	Z 		atomic number
	ion 	ion stage
	ll		lower level
	E0 		threshold energy eV 
	lu 		upper level
	np 		number of entries
	energy 	energies 
	XS 		cross sections
	'''
	def __init__(self,  _z, _ion, _ll, _lu, _E0, n_p, _energy, _XS):
		self.Z = _z
		self.ion = _ion
		self.ll = _ll
		self.lu = _lu
		self.E0 = _E0
		self.np = n_p
		self.energy = _energy
		self.XS = _XS



class top_photo:
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