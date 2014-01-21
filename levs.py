'''
code to plot departure coefficients for a few cells of choice
from a wind save file
'''

import csv, sys, os, array, warnings, subprocess
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pywind_sub as ps




nlevels = 80


# write out command file for pywind
cmdfile = open("levcmds", "w")

if sys.argv[-1] == "a":
	cmdfile.write("1\nB\nn\nt\nE\nM\n1\n0\n")
else:
	cmdfile.write("1\nB\nn\nt\nE\nM\n1\n1\n")

for i in range(0,nlevels):
	cmdfile.write(str(i) + "\n")

cmdfile.write("q\n")

cmdfile.close()




# run pywind
filename = sys.argv[1]
ncells = len(sys.argv[2:]) - 1

if sys.argv[-1] != "x":
	print "RUNNING PYWIND..."
	os.system("py_wind %s < levcmds > tempout" % filename)
	if sys.argv[-1] != "a":
		ncells += 1





# array of cells we want
nplasma = np.zeros (ncells)
for i in range(ncells):
	nplasma[i] = int(sys.argv[2+i]) 






# read in densities and plasma cell numbers
nplasmas = np.loadtxt("%s.pnum.dat" % (filename) , comments = "#", usecols = (2, 3), unpack = True)[0]
ne = np.loadtxt("%s.ne.dat" % (filename) , comments = "#", usecols = (2, 3), unpack = True)[0]
te = np.loadtxt("%s.te.dat" % (filename) , comments = "#", usecols = (2, 3), unpack = True)[0]
conv = np.loadtxt("%s.conv_whole.dat" % (filename) , comments = "#", usecols = (2, 3), unpack = True)[0]






# create ne arrays
ne_arr=[]
te_arr = []
conv_arr = []
for j in range(len(nplasma)):
	n = nplasma[j]
	i_find = np.where( nplasmas == n)
	ne_arr.append(ne[i_find[0][0]])
	te_arr.append(te[i_find[0][0]])
	conv_arr.append(conv[i_find[0][0]])






# create zeros arrays for each cell
dep_plot = []
for j in range(len(nplasma)):
	dep_plot.append (np.zeros(nlevels))
dep_plot = np.array(dep_plot)








for i in range(nlevels):

	# read in de coefficients
	dep_coefs = np.loadtxt("%s.dep_coef_%i.dat" % (filename, i) , comments = "#", usecols = (2, 3), unpack = True)
	dep_coefs = dep_coefs[0]



	for j in range(len(nplasma)):

		n = nplasma[j]

		i_find = np.where( nplasmas == n)

		if len(i_find) > 1 or len(i_find[0])>1:
			print "error, multiple matches..."


		#print dep_coefs[i_find[0][0]]
		d = dep_coefs[i_find[0][0]]

		#if d < 100:
		dep_plot[j][i] = d


n_range = np.arange(1, nlevels + 1)


fig = plt.figure(figsize=(8.3,11.7),dpi=80)
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

h = 0
he1 = 11
he2 = 64 
all_lev = 75

for j in range(len(nplasma)):
	ax1.plot(n_range[h:he1], dep_plot[j][h:he1], label = "n: %i ne: %4.2e te: %4.2e c: %i" % (nplasma[j], ne_arr[j], te_arr[j], conv_arr[j]) )
	ax2.plot(n_range[he1:he2] - he1, dep_plot[j][he1:he2], label = "n: %i ne: %4.2e te: %4.2e c: %i" % (nplasma[j], ne_arr[j], te_arr[j], conv_arr[j]) )
	ax3.plot(n_range[he2:all_lev] - he2, dep_plot[j][he2:all_lev], label = "n: %i ne: %4.2e te: %4.2e c: %i" % (nplasma[j], ne_arr[j], te_arr[j], conv_arr[j]) )

ax1.set_yscale("log")
ax2.set_yscale("log")
ax3.set_yscale("log")

ax1.set_title("H1")
ax2.set_title("He1")
ax3.set_title("He2")

if sys.argv[-1] == "a":
	ax2.set_ylabel("Lev Pop")
else:
	ax2.set_ylabel("Dep Coef")
ax3.set_xlabel ("nlevel_macro")

ax1.legend(loc = 2)

if sys.argv[-1] == "a":
	savename = "levpop_%i.png" % nplasma[j]
else:
	savename = "depcoef_%i.png" % nplasma[j]


plt.savefig(savename)
os.system("open -a preview %s" % savename )

#os.system("rm -f levcmds")
#os.system("rm -f tempout")


