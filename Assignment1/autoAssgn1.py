#!/usr/bin/python

import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np

UCellSize = np.arange(4,11)
nAtoms = 4*np.power(UCellSize, 3)

# Edit input file
fnames = ["md", "lmd"]
for j in range(2):
	fname = fnames[j]
	infile = fname + '.in'
	cfile = fname + '.c'
	exfile = fname
	fp = open(infile, 'r+')
	strOld = fp.read()
	fp.close()

	time = []
	for i in range(4,11):
		fp  = open(infile, 'w+')
		newSize = str(i) + ' '
		strNew =  3*newSize + '\n'+ strOld[9:]
		fp.write(strNew)
		fp.close()
	
		# Run MD program
		subprocess.call(['gcc', '-O', '-o', exfile, cfile, '-lm'])
		comm = './' + exfile + '<' + infile + '>out.txt'
		os.system('./md<md.in>out.txt')
		time.append(float(open('out.txt','r').read()))

	print(time)

# Generate convergence plots

#plt.loglog(nAtoms, time, 'r-o')
#plt.show()
