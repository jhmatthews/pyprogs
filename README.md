pyprogs
=======

Scripts used for Python [the RT code]. Often, naturally, confusingly, written in Python [the language]


### py_top_phot:

C script which was traditionally used to smooth over cross sections. 

'''
#Compilation: 
cd py_top_phot
make

#Usage:
The program uses rdpar to get parameters

  rdstr ("Root.for.all.files", root);      Rootname for all of the input files.  The program assumes
					   that the input data file will have the extension .dat
						and therefore the output file will have the extension
						.py
  rdint ("No.of.output.points", &nout);    number of output x-sections for each input x section
  rdflo ("Fractional.smoothing", &delta);  The fractional size, e.g. 0.5 for the boxcar smooth
'''
