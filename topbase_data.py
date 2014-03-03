#!/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''
topbase_data

Routines for dealing with topbase data and often converting into 
format for the radiative transfer code Python
'''


import os, sys
import numpy as np
import atomic_classes as cls
import chianti_sub as sub 

ANGSTROM  = 1.e-8
C  = 2.997925e10
RYDBERG = 13.60569253		# rydberg in eV
HEV = 4.13620e-15			# Planck's constant in eV 









class line():

	''' 
	A load of functions for manipulating and writing lines data from topbase
	'''

	def read_topbase_lines(self, filename):


		'''read in line info from original topbase levels and 
		place in class instance array'''

		f = open(filename, "r")

		linelist = []
		
		for line in f:

			data = line.split()

			if data[0] != "#":

				z = int (data[1])
				ne = int(data[2])

				ion = z + 1 - ne

				gf = data[-4]		# oscillator strength weighted by lower level

				islp = int(data[3])
				jslp = int(data[4])

				ilv = int(data[5])
				jlv = int(data[6])

				gi = float(data[-2])
				gj = float(data[-1])

				wave = float(data[-3])
				freq = C / (ANGSTROM * wave)

				# read string configurations ... a bit complicates
				if len(data[7]) == 2:
					iconf = data[7]+" "+data[8]
					i_conf = 9
				elif len(data[7]) == 3:
					iconf = data[7]
					i_conf = 8
				else:
					iconf = data[7][0:2]+" "+data[7][2:]
					i_conf = 8

				if len(data[i_conf]) == 2:
					jconf = data[i_conf]+" "+data[i_conf+1]
				elif len(data[i_conf]) == 3:
					jconf = data[i_conf]
				else:
					jconf = data[i_conf][0:2]+" "+data[i_conf][2:]


				if gf[0] == "-" or float(gf) < 0:

					gf = - float(gf)

					gu = gj
					gl = gi

					ulv = jlv
					llv = ilv

					conf_u = jconf
					conf_l = iconf

					islpl = islp
					islpu = jslp

				else:

					gf = float(gf)

					gu = gi
					gl = gj

					ulv = ilv
					llv = jlv

					conf_u = iconf
					conf_l = jconf

					islpl = jslp
					islpu = islp


				osc = gf / gl 			# non weighted oscillator strength, gf/glower

				#print conf_u, "->", conf_l


				line_i = cls.topline (z, ion, wave, freq, gf, osc, gl, gu, None, None, None, None, conf_u, conf_l, ulv, llv, islpl, islpu)
				
				linelist.append( line_i )

		linelist = np.array(linelist)

		return linelist


	def match_lines(self, linelist, levellist):

		'''
		matches lines with level list , returns line rather than topline atomic_classes
		'''

		new_line = []

		for i in range(len(linelist)):

			ll = None
			lu = None

			z = linelist[i].z
			ion = linelist[i].ion 
			osc = linelist[i].osc
			gl = linelist[i].g_l
			gu = linelist[i].g_u
			wave = linelist[i].wavelength
			freq = linelist[i].freq
			islp = linelist[i].islpu
			jslp = linelist[i].islpl

			for j in range(len(levellist)):

				if levellist[j].z == z and levellist[j].ion == ion:

					if levellist[j].nnstring == linelist[i].conf_u and islp == levellist[j].islp:

					 	lu = levellist[j].lvl
					 	eu = levellist[j].E

					elif levellist[j].nnstring == linelist[i].conf_l and jslp == levellist[j].islp:

					 	ll = levellist[j].lvl
					 	el = levellist[j].E

			if eu < el:
				print "Error, el>eu, %s to %s not found in levels" % (linelist[i].conf_u, linelist[i].conf_l)
			if ll == None or lu == None:
				print "Error, line %i, %s to %s not found in levels" % (i, linelist[i].conf_u, linelist[i].conf_l)

			new_line.append( cls.line(z, ion, wave, freq, osc, gl, gu, el, eu, ll, lu) )


		return np.array(new_line)







class level():

	''' 
	A load of functions for manipulating and writing level data from topbase
	'''


	def read_topbase_levels(self, filename):

		'''read in level info from original topbase levels and 
		place in class instance array'''

		f = open(filename, "r")

		lev = []
		
		for line in f:

			data = line.split()

			if data[0] != "#":

				ne = int(data[2])
				z = int (data[1])
				ion = z + 1 - ne

				lvl = int (data[0])

				ilv = int (data[4])

				islp = int (data[3])
				
				E = float (data[-2]) * RYDBERG
				ionpot = float (data[-3]) * RYDBERG
				g = float (data[-1])
				#rad_rate = float (level_array_read[i][-1])

				if len(data) == 10:
					nnstring = data[5]+" "+data[6]
				elif len(data) == 9:
					if data[5] == "1s2": 
						nnstring= data[5]
					else:
						nnstring = data[5][0:2]+" "+data[5][2:]

				
				lev.append( cls.toplevel(z, ion, lvl, ionpot, E, g, 0.0, nnstring, islp, ilv) )

		lev = np.array(lev)
			
		# level is an array of class instances like the line_ptr in PYTHONRT
		return lev


	def modify_upper_energies(self, toplev, threshold, ion = 2):

		top = []

		for i in range(len(toplev)):
			if toplev[i].ion == ion:
				toplev[i].E = toplev[i].E + threshold 

			top.append(toplev[i])

		return np.array(toplev)






	def sort_toplev(self, class_array, attr_string="E"):

		'''sorts array of class instances into correct order (by lower level)'''

		import operator

		# sort by lower level as default, provided as keyword arg
		sorted_class_array = sorted(class_array, key=operator.attrgetter(attr_string))

		# return numpy array of class instances
		return np.array(sorted_class_array)


#	def split_by_ion(self, toplev, z):
#
#		topnew = []
#
#		ilast = 0
#
#		for i in range(len(toplev)):











class photo():

	''' 
	A load of functions for manipulating and writing photoionzioation data from topbase
	'''

	def read_original_topbase_xs(self, filename, suffix="Mac"):

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
				ion = Z - int(data[2]) + 1 
				islp = int(data[3])
				l = int(data[4])
				E = -1 * float(data[5]) * RYDBERG
				num = int(data[6])

				top_dummy = cls.top_photo (Z, ion, islp, l, E, num, np.ndarray( num, dtype=np.object), np.ndarray( num, dtype=np.object))

			elif data[0]!="#": 
				print "Didn't understand data, exiting"
				print line
				sys.exit()

		top = np.array(top)

		return top


	def read_topbase(self, filename, mode="Top"):
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
		Z,ion,islp,l, E0, num_records = sum_records = np.loadtxt( filename, 
				dtype={'names': ('Z', 'ion', 'islp', 'l', 'E0', 'np'), 
				'formats': ('i4', 'i4', 'i4', 'i4', 'float', 'i4')}, 
	                        comments='Phot%s ' % mode , delimiter=None, converters=None, 
	                        skiprows=0, usecols=(1,2,3,4,5,6), unpack=True, ndmin=0)
		
		# then read the actual cross sections
		energy, XS = np.loadtxt(filename, dtype='float', 
	                        comments='Phot%sS' % mode, delimiter=None, converters=None, 
	                        skiprows=0, usecols=(1,2), unpack=True, ndmin=0)

		#ion = int(Z) - int(ne) +1

		# create blank array to store topbase class instances
		top = np.ndarray( len(Z),dtype=np.object)

		sr = np.transpose(sum_records)

		nline = 0
		
		for i in range(len(top)):

			n_p = int(num_records[i])
			nmax = nline + n_p

			#ion = int(Z[i]) - int(ne[i]) + 1

			
			top[i] = cls.top_photo(Z[i], int(ion[i]), int(islp[i]), int(l[i]), E0[i], n_p, energy[nline:nmax], XS[nline:nmax])

			#if ion[i] == 2 and int(islp[i]) / 100 != 2:
			#	print "Error, non doublet in He II! ion %i islp[i] %i" % (ion, int(islp[i]))
			#	print sr[i]
			#	print ne[i], ion[i] 

			nline = nmax
			
		# top is an array of class instances like the topbase_ptr in PYTHONRT
		return top






	def prepare_topbase_xs(self, top, levels, top_levels, suffix="Mac", levmax = 1e50):

		'''write array of topbase class instances to file in macro atom format, with lower and upper levels specified'''



		lu = 1 # we assume ionization to ground state of +ion


		print len(top), len(levels), len(top_levels)

		count = 0
		topnew = []

		for i in range(len(top)):


			#print top_levels[i].z, top_levels[i].nnstring



			for j in range(len(levels)):


				#print top_levels[i].nnstring, "//" ,levels[j].nnstring

				if top_levels[i].nnstring == levels[j].nnstring and top_levels[i].z == levels[j].z and top_levels[i].ion == levels[j].ion:


					if levels[j].lvl<=levmax:

						ll = levels[j].lvl 
						ion = levels[j].ion

						#print "match", ll, top_levels[i].z, levels[j].nnstring

						## write the summary records
						n_p = top[i].np
	 

						if n_p>100: 
							n_p = 100

						if n_p>0:
							topnew_class = cls.top_photo_mac ( top[i].Z, top[i].ion, ll, lu, top[i].E0, n_p, top[i].energy, top[i].XS)
							topnew.append ( topnew_class)

							count +=1

		return np.array(topnew)

	def prepare_topbase_xs2(self, top, toplevels, suffix="Mac", levmax = 1e50):

		'''write array of topbase class instances to file in macro atom format, with lower and upper levels specified'''



		lu = 1 # we assume ionization to ground state of +ion



		count = 0
		topnew = []
		cc = 0

		print len(top), len(toplevels)

		for i in range(len(top)):



			for j in range(len(toplevels)):



				if toplevels[j].z == top[i].Z and toplevels[j].ion == top[i].ion and toplevels[j].ilv == top[i].l and toplevels[j].islp== top[i].islp:


					if toplevels[j].lvl<=levmax:

						ll = toplevels[j].lvl 
						ion = toplevels[j].ion

						#print "match", ll, top_levels[i].z, levels[j].nnstring

						## write the summary records
						n_p = top[i].np
	 

						#if n_p>100: 
						#	n_p = 100

						if n_p>0:
							topnew_class = cls.top_photo_mac ( top[i].Z, top[i].ion, ll, lu, top[i].E0, n_p, top[i].energy, top[i].XS)
							topnew.append ( topnew_class)

							count +=1


		return np.array(topnew)






	def read_top_macro(self, filename, mode="Mac"):
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
		Z,ion,ll,lu, E0, num_records = sum_records = np.loadtxt( filename, 
				dtype={'names': ('Z', 'ion', 'islp', 'l', 'E0', 'np'), 
				'formats': ('i4', 'i4', 'i4', 'i4', 'float', 'i4')}, 
	                        comments='Phot%s ' % mode , delimiter=None, converters=None, 
	                        skiprows=0, usecols=(1,2,3,4,5,6), unpack=True, ndmin=0)
		
		# then read the actual cross sections
		energy, XS = np.loadtxt(filename, dtype='float', 
	                        comments='Phot%sS' % mode, delimiter=None, converters=None, 
	                        skiprows=0, usecols=(1,2), unpack=True, ndmin=0)

		# create blank array to store topbase class instances
		top = np.ndarray( len(Z),dtype=np.object)

		nline = 0
		
		for i in range(len(top)):

			n_p = int(num_records[i])
			nmax = nline + n_p
			
			top[i] = cls.top_photo_mac (Z[i], int(ion[i]), int(ll[i]), int(lu[i]), E0[i], n_p, energy[nline:nmax], XS[nline:nmax])

			nline = nmax
			
		# top is an array of class instances like the topbase_ptr in PYTHONRT
		return top



	def write_top_macro(self, topnew, filename, suffix = "Mac", append = False, levmax = 100, z=None, ion = None):

		if append:
			write_string = 'a'
		else:
			write_string = 'w'
		
		file_write = open( filename, write_string)

		topnew = sub.sort_class (topnew)

		#print len(topnew), len(topnew2), count 

		for i in range(len(topnew)):


			zbool = (z == topnew[i].Z or z == None)
			ibool = (ion == topnew[i].ion or ion == None)

			if zbool and ibool:

				if topnew[i].ll <=levmax and topnew[i].lu <=levmax:

					file_write.write('Phot%sS  %1d  %1d %3d    %1d     %.6f  %2d\n' %
				                                    ( suffix, topnew[i].Z, topnew[i].ion, topnew[i].ll, topnew[i].lu, topnew[i].E0, topnew[i].np ))


					## write the actual XSs
					for j in range( topnew[i].np ):
							
						file_write.write('Phot%s     %.6f %.3e\n' % (suffix,  topnew[i].energy[j], topnew[i].XS[j]) )



		file_write.close()
			
		return 0









'''
Z = 2
ion_threshold = 24.310389
dE = 0.001
LEVMAX = 100
ion_threshold = 78.733180




phot = photo()
lev = level() 
lin = line()


topphot = phot.read_topbase("topbase_he_xs.py")					# read topbase Xsections - this should be smoothed with py_top_phot

toplev = lev.read_topbase_levels("topbase_he_levels")			# read topbase_levels




toplev = lev.modify_upper_energies (toplev, 24.310389, ion = 2)		# modify He II energies to be higher than He I by threshold

toplev = lev.sort_toplev(toplev)

ilv = 0
ilast = toplev[0].ion

for i in range(len(toplev)):

	ilv +=1

	toplev[i].lvl = i + 1
	#topphot[i].lvl = i + 1

	if toplev[i].ion != ilast:
		ilv = 0

	ilast = toplev[i].ion


for i in range(len(topphot)):
	if topphot[i].ion == 1: print topphot[i].ion, i


top_macro = phot.prepare_topbase_xs2(topphot, toplev, levmax = LEVMAX)		# bit weird, matching with itself


linelist = lin.read_topbase_lines("topbase_he_lines")

newlines = lin.match_lines (linelist, toplev)



phot.write_top_macro(top_macro, "he_top_phot_macro.py", z = 2, ion = 1)
#sub.write_level_file (toplev, "he_top_levels.py", z = 2, ion = 1)
#sub.write_line_file(newlines, "he_top_lines.py", z = 2, ion = 1)



#sub.write_line_file(newlines, "he_top_lines.py", z = 2, ion = 2, append = True)


'''










