import dassflow2d as df2d
import matplotlib.pyplot as plt
import shutil, os, sys
import numpy as np
import csv
import random
from mpi4py import MPI


path_to_SP = os.path.abspath(os.path.join(os.path.dirname(__file__),'../..'))

path_to_Outils = os.path.join(path_to_SP,'Outils')

sys.path.append(path_to_SP)

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

from Outils.out import *
from Outils.mesh_generator import *

##########
# Initialise bin
##########

os.chdir(code_dir)

os.system("make cleanres")			   # removes all in bin_dir/res directory
os.system("make cleanmsh")			   # removes all in bin_dir/msh directory
os.system("make cleanmin")			   # removes all in bin_dir/min directory

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
# Mesh
##########

L = 1000
dx = 10
type = 'channel'

mesh_name = "channel.geo"

ts = 1400
dtw = ts*0.2

use_porosity = 1

input_params={ "mesh_type" : 'basic', 
              "mesh_name": mesh_name,
              "lx" : L,
              "ly" : 9.0,
              "nx" : 101,
              "ny" : 4,

                "ts": ts ,
                "use_obs":'0',
                "use_UVobs":'0',
                "use_Zobs":'0',
                "use_porosity": use_porosity,
                
                "w_obs":'1',
                "w_vtk":'1', 
                "w_gnuplot":'1',
                "w_tecplot":'0',
                "adapt_dt":'1',

                "dt":'0.1',

                "dtw": dtw ,
                "dta":"100",
                
                "bc_infil":"0",
                "bc_rain":"0",}

df2d.wrapping.read_input(os.path.join(bin_dir,"input.txt"))
my_model = df2d.dassflowmodel(bin_dir =  bin_dir, hdf5_path = os.path.join(bin_dir,"res","simu.hdf5") , run_type = "direct", clean = True, custom_config=input_params)
# Reinitialize mesh from new updated geo file
my_model.init_mesh()
# my_model.meshing.plot() 

# Input parameters for the generation of observations
Config = df2d.core.config.Config()
Config.set(custom_config = input_params)
nc = my_model.kernel.mesh.nc    # number of cells 
nland = nc

##########
# Friction
##########

#Create Python class by calling wrapped initialise routines
my_model.kernel.my_friction =  df2d.wrapping.m_model.friction_data(my_model.kernel.mesh)
#Allocate and get initial values from Fortran
my_model.kernel.my_friction.nland = nland
df2d.wrapping.call_model.init_friction(my_model.kernel)

#Provide values, on top of initial ones from Fortran initialization routine, in Python structure
my_model.kernel.my_friction.manning[:] = 0.0
my_model.kernel.my_friction.manning_beta[:] = 0
my_model.kernel.my_friction.land[:] = 1

##########
# Porosity
##########

#Create Python class by calling wrapped initialise routines
my_model.kernel.my_porosity = df2d.wrapping.m_model.porosity_data(my_model.kernel.mesh)

#Allocate and get initial values from Fortran
my_model.kernel.my_porosity.nland = nland
df2d.wrapping.call_model.init_porosity(my_model.kernel)

#Provide values, on top of initial ones from Fortran initialization routine, in Python structure
for i in range(nland) :
    my_model.kernel.my_porosity.phi[i] = 0.5
    my_model.kernel.my_porosity.hbanks[i] = 18.0 
    my_model.kernel.my_porosity.gamma[i] = 2.0 
my_model.kernel.my_porosity.land[:] = np.arange(1, nland + 1)

##########
# Hydraulic states
##########

my_model.kernel.dof  = df2d.wrapping.m_model.unk(my_model.kernel.mesh)
my_model.kernel.dof0 = my_model.kernel.dof

for i in range(nc) :
    my_model.kernel.dof0.h[i] = 1.0 + random.randint(0,1)*0.1 

my_model.kernel.dof0.u[:] = 0.0
my_model.kernel.dof0.v[:] = 0.0


##########################################
# Run
##########################################

df2d.wrapping.call_model.init_fortran(my_model.kernel)
df2d.wrapping.call_model.run(my_model.kernel, arg = "direct")