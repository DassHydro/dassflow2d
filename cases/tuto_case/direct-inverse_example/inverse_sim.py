import dassflow2d as df2d
import matplotlib.pyplot as plt
import shutil, os
import numpy as np
import csv
from mpi4py import MPI

# TODO : adapt this into ? model.meshing object + mesh_extent_box from a real dassflow mesh/setup (with wrapped object)
class mesh_box:
    def __init__(self, xmin, ymax, ncol, nrow, dx, dy, x_cell, y_cell):
        self.xmin = xmin
        self.ymax = ymax
        self.ncol = ncol
        self.nrow = nrow
        self.dx = dx
        self.dy = dy
        self.x_cell = x_cell
        self.y_cell = y_cell


##############################################################################################################

code_dir =  os.getcwd() #os.path.join(dassflow_dir,"code")
bin_dir = os.path.join(code_dir,"bin_A")

##############################################################################################################

##########
# Initialise bin
##########

os.chdir(code_dir)

os.system("make cleanres")			   # removes all in bin_dir/res directory
os.system("make cleanmsh")			   # removes all in bin_dir/msh directory
os.system("make cleanmin")			   # removes all in bin_dir/min directory

#if os.path.isfile(f"rm {bin_dir}/restart.bin"):
#	os.system(f"rm {bin_dir}/restart.bin")   # removes all in bin_dir/msh directory

os.chdir(bin_dir)

##########
# MPI
##########

df2d.wrapping.m_mpi.init_mpi()
rank = df2d.wrapping.m_mpi.get_proc() # get the rank and number of processors
nproc = df2d.wrapping.m_mpi.get_np()
mpi = [rank, nproc]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

##########
# Model
##########
input_params={ "mesh_name":'mesh.geo',
              "ts":"1000",
			   "use_obs":'1', #<- set to 1 
			   "use_UVobs":'0',
			   "use_Zobs":'1', #<- we use water surface elevation observations
               
			   "w_obs":'1',
               "w_vtk":'2',
               "w_gnuplot":'1',
               "w_tecplot":'0',
      		   "c_hydrograph":'0',
               "c_manning":'1', #<- we infer the (homogeneously distributed) friction
               "dtw":'100',
               "dta":"150",
               
               "bc_infil":"0",
			   "bc_rain":"0",}

df2d.wrapping.read_input(os.path.join(bin_dir,"input.txt"))
my_model = df2d.dassflowmodel(bin_dir =  bin_dir, hdf5_path = os.path.join(bin_dir,"res","simu.hdf5") , run_type = "direct", clean = True, custom_config=input_params)

#Re initialize mesh from new updated geo file
my_model.init_mesh()
#my_model.meshing.plot()
# Input parameters for the generation of observations

Config = df2d.core.config.Config()
Config.set(custom_config = input_params)

##########
# Friction
##########

#Create Python class by calling wrapped initialise routines
my_model.kernel.my_friction =  df2d.wrapping.m_model.friction_data(my_model.kernel.mesh)
#Allocate and get initial values from Fortran
my_model.kernel.my_friction.nland = 1
df2d.wrapping.call_model.init_friction(my_model.kernel)

#Provide values, on top of initial ones from Fortran initialization routine, in Python structure

my_model.kernel.my_friction.manning[:] = 0.1 #<- value changed from target model
my_model.kernel.my_friction.manning_beta[:] = 0
my_model.kernel.my_friction.land[:] = 1 #<- homogeneous friction

##########
# Hydraulic states
##########

my_model.kernel.dof  = df2d.wrapping.m_model.unk(my_model.kernel.mesh)
my_model.kernel.dof0 = my_model.kernel.dof
#my_model.kernel.dof0.h[:] = 0  #Disregarded if ic.bin is provided

##############################################################################################################

##########################################
# Run
##########################################

df2d.wrapping.call_model.init_fortran(my_model.kernel)

df2d.wrapping.call_model.run(my_model.kernel, arg = "min")

df2d.wrapping.call_model.clean_model(my_model.kernel)

