#=======================================================#
# Source librairies
#=======================================================#
import dassflow2d as df2d
import numpy as np
import sys
import os



run_type ="min" # if runtype = rundirect runmin, classic treatment, if runtype = runminpython


df2d.wrapping.m_mpi.init_mpi()
rank = df2d.wrapping.m_mpi.get_proc() # get the rank and number of processors
nproc = df2d.wrapping.m_mpi.get_np()
mpi = [rank, nproc]
print("mpi=", mpi)
# store main path
bin_dir        = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A"
hdf5_path      = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/res/simu.hdf5"
        # run type
clean = True
print("clean ok")
        # ------------ Define Default values------------------
config = df2d.core.config.Config()
print("config ok")

df2d.wrapping.read_input(f"{bin_dir}/input.txt") 
print("read input ok")
config.get() 
print("before intimodel ")
kernel=df2d.wrapping.call_model.Model()
print("after init model ")
os.chdir(bin_dir)
kernel.mesh = df2d.wrapping.m_mesh.msh()
df2d.wrapping.call_model.init_solver(kernel)
meshing = df2d.core.meshing.Meshing(mesh_fortran = kernel.mesh)
print("mmeshing ok")



os.chdir(bin_dir)
kernel.dof = df2d.wrapping.m_model.unk(meshing.mesh_fortran)
kernel.dof0 = df2d.wrapping.m_model.unk(meshing.mesh_fortran)        
print("dof0 ok")

df2d.wrapping.call_model.init_fortran(kernel)
df2d.wrapping.call_model.run(self = kernel, arg = run_type)
df2d.wrapping.call_model.clean_model(kernel)  
