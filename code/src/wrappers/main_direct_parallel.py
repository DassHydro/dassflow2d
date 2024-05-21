#=======================================================#
# Source librairies
#=======================================================#
import dassflow2d as df2d
import numpy as np
import sys
import os
from mpi4py import MPI

df2d.wrapping.m_mpi.init_mpi()

# store main path
dassflow_dir= os.path.abspath(os.path.join(__file__ ,"../../.."))
code_dir =  f"{dassflow_dir}/code"
bin_dir = f"{code_dir}/bin_A"

##########
# Initialise bin
##########

os.chdir(code_dir)

os.system("make cleanres")			   # removes all in bin_dir/res directory
os.system("make cleanmsh")			   # removes all in bin_dir/msh directory
os.system("make cleanmin")			   # removes all in bin_dir/min directory

if os.path.isfile(f"rm {bin_dir}/restart.bin"):
	os.system(f"rm {bin_dir}/restart.bin")   # removes all in bin_dir/msh directory

os.chdir(bin_dir)

# ------------ Define Default values------------------

df2d.wrapping.read_input(f"{bin_dir}/input.txt")

model = df2d.dassflowmodel(bin_dir =  bin_dir, hdf5_path = f"{bin_dir}/res/simu.hdf5" , run_type = "direct", clean = True, custom_config = None)

model.init_mesh()

model.kernel.dof  = df2d.wrapping.m_model.unk(model.kernel.mesh)
model.kernel.dof0 = model.kernel.dof     

df2d.wrapping.call_model.init_fortran(model.kernel)
df2d.wrapping.call_model.run(model.kernel, arg = "direct")
df2d.wrapping.call_model.clean_model(model.kernel)  
