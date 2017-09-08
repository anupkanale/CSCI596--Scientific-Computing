#!/usr/bin/python

import subprocess

# Edit input file
fp = open('md.in', 'r+')
strOld = fp.read()
fp.close()

fp  = open('md.in', 'w+')
for i in range(4,10)
	newSize = str(i) + ' '
	strNew =  3*newSize + '\n' + strOld[9:]
	
	fp.write(strNew)
	
	# Run MD program
	subprocess.call('gcc -O -o md md.c -lm')
	subprocess.call('./md < md.in > out.txt')
	
	
	
# Write data into text files


# Generate convergence plots
