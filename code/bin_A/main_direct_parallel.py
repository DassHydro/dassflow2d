#=======================================================#
# Source librairies
#=======================================================#
import dassflow2d as df2d
import numpy as np
import sys
import os


run_type ="direct" # if runtype = rundirect runmin, classic treatment, if runtype = runminpython

df2d.wrapping.m_mpi.init_mpi()
rank = df2d.wrapping.m_mpi.get_proc() # get the rank and number of processors
nproc = df2d.wrapping.m_mpi.get_np()
mpi = [rank, nproc]
print("mpi=", mpi)

# store main path
dassflow_dir="/home/leo/DISTANT/dassflow2d"
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
