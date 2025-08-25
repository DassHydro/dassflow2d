from . import metrics
from . import minimization
from . import read_results
from . import results


import numpy as np
import os

# return an array of size n of real numbers in the Fortran format
def real_array(n):
	return np.empty(dtype='double', shape=(n,),order='F')

# return an array of size n of integer numbers in the Fortran format
def int_array(n):
	return np.empty(dtype='int32', shape=(n,),order='F')
	
# return a 2d array of size (n,m) of real numbers in the Fortran format
def real_array2d(n,m):
	return np.empty(dtype='double', shape=(n,m),order='F')

# return a 2d array of size (n,m) of integer numbers in the Fortran format
def int_array2d(n,m):
	return np.empty(dtype='int32', shape=(n,m),order='F')

# read input from file and save in param
def read_input_py(filename, param):
	#put input file in a more compact format
	os.system("sed -e '/^!/d;/^$/d;s/!.*//g;s/[ \t]//g' %s > input_python.post"%(filename))
	fname = open("input_python.post",'r')
	data = fname.readlines()[1:-1]
	fname.close()
	
	indata = []
	
	for d in data:
		# remove the = symbol in d
		d1 = d.split('=')[1]
		# remove the \n symbol in the end of d1
		d1 = d1[:-1]
		# remove the , symbol in the end of d1 if it's present
		if d1[-1]==',':
			d1 = d1[:-1]
		
		try:
			# first try to convert d1 into an integer
			indata.append(int(d1))
		except ValueError:
			# try to convert d1 into float
			try:
				indata.append(float(d1))
			except ValueError:
				# check if d1 is in the Fortran exponential form of real number
				if 'd-' in d1:
					#replace d- by e-
					d1=d1.replace('d','e')
					indata.append(float(d1))
				else:
					#d1 is str for some names (mesh, scheme name ...), remove the simple quotes
					d1 = d1[1:-1]
					indata.append(d1)
