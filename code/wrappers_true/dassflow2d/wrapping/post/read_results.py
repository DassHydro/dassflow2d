import dassflow2d as df2d
import numpy as np
from os import listdir, getcwd
from os.path import isfile,join
	
def get_result_files():
	current_path = getcwd()+"/res"
	ti = df2d.m_common.get_tc0()
	te = df2d.m_common.get_ts()
	file_list = []
	time_list = []
	for f in listdir(current_path):
		if isfile(join(current_path,f)) and f[len(f)-4:]==".dat": #select .dat files from the current working directory
			
			
			file_list.append(f)
			
			# get the time instance from the file name
			time = f[:len(f)-4].split('_')[1]
			try :
				time = float(time)
			except ValueError:
				if time == "initial":
					time = ti # write the correct time instance using get functions from dassflow ! 
				elif time == "final":
					time = te
			
			# save the time instance
			time_list.append(time)
	return file_list, time_list

# read data from .dat file
def load_dof(f): 
	if not f[len(f)-4:] == ".dat": #check extension of the given file
		raise ValueError('The given file isnot a .dat')
	
	data = {'i':[],'x':[],'y':[],'bathy':[],'h':[],'zs':[],'manning':[],'u':[],'v':[]}
	keys = list(data.keys())
	filename = getcwd()+ "/res/" + f
	with open(filename, 'r') as txt: 
		for line in txt.readlines():
			if line[1] == "#":
				continue
			line = line.split()
			for k in range(len(line)):
				data_temp = line[k]
				if k == 0:
					data[keys[k]].append(int(data_temp))
				else:
					data[keys[k]].append(float(data_temp))
	return data
		
