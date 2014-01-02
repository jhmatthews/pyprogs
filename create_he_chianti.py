#!/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''
	create_he_chianti.py

This is a routine to turn Chianti data 

Chianti files needed:

he_1.wgfa
he_1.elvlc
he_2.wgfa
he_2.elvlc

'''
import chianti_sub as sub
import os, sys
import numpy as np

RYDBERG = 13.60569253	# rydberg in eV
HEV = 4.13620e-15	# Planck's constant in eV 


Z = 2
ion_threshold = 24.587387512


# first create level file

level_info, line_info = sub.read_chianti_data(level_filename="he_1.elvlc", radiative_filename="he_1.wgfa")

lev_he1 = sub.chianti_to_lev (level_info, ion_threshold, Z, 1)

sub.write_level_file (lev_he1, "test")



lev = sub.read_level_info ("data/atomic_macro/h10_levels.py")

lev_he2 = sub.convert_h_to_hydrogenic_class (lev, Z, 78.733180)

sub.write_level_file (lev_he2, "test", append = True)




line = sub.chianti_to_line (line_info, level_info, Z, 1)
sub.write_line_file (line, "test_lines")

lines = sub.read_line_info ("data/atomic_macro/h10_lines.py")




topphot = sub.read_topbase_xs("topbase_he_xs")

sub.write_topbase_xs(topphot, "he_top_phot_macro.py")








#le = calc_levels("he", 2, [25], [1], [24.587387512])

#print le 

#sub.write_level_file(le, "test")

'''
Z = 2
ion_threshold = 24.587387512
		
#nist = sub.read_nist("he_1.nist", 2, 1)
#lev = sub.nist_to_chianti (nist, 24.587387512)
#sub.write_level_file (lev, "alllevels")

nist = sub.read_nist("he_1.nist", 2, 1)

nist = sub.compress_levels (nist)

lev = sub.nist_to_lev (nist, ion_threshold)

sub.write_level_file (lev, "test")



Z = 2
ion_threshold = 24.587387512

lev = sub.read_level_info ("data/atomic_macro/h10_levels.py")

lev = sub.convert_h_to_hydrogenic_class(lev, Z, 78.733180)


sub.write_level_file (lev, "test", append = True)


sub.read_chianti_data("he_1.elvlc")'''






