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

from Outils.out import *
from Outils.mesh_generator import *

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
# Mesh
##########


#mesh_name = gen_mesh(type,L,dx)
mesh_name = "channel.geo"


##########
# Model
##########

ts = 3000
dtw = ts*0.2

# REMARK : In this test-case, you must set nland equal to the number of cells in your mesh

use_porosity = 1

input_params={ "mesh_name": mesh_name ,
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

#Re initialize mesh from new updated geo file
my_model.init_mesh()
#my_model.meshing.plot()
# Input parameters for the generation of observations

Config = df2d.core.config.Config()
Config.set(custom_config = input_params)

##########
# Friction
##########

nc = my_model.kernel.mesh.nc

nland = nc

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

nc = my_model.kernel.mesh.nc

#Create Python class by calling wrapped initialise routines

my_model.kernel.my_porosity = df2d.wrapping.m_model.porosity_data(my_model.kernel.mesh)

#Allocate and get initial values from Fortran

my_model.kernel.my_porosity.nland = my_model.kernel.my_friction.nland
df2d.wrapping.call_model.init_porosity(my_model.kernel)

#Provide values, on top of initial ones from Fortran initialization routine, in Python structure

for i in range(nland) :

    my_model.kernel.my_porosity.phi[i] = 0.5

my_model.kernel.my_porosity.land[:] = range(1,nland+1)



##########
# Hydraulic states
##########

my_model.kernel.dof  = df2d.wrapping.m_model.unk(my_model.kernel.mesh)
my_model.kernel.dof0 = my_model.kernel.dof

for i in range(nc) :
    my_model.kernel.dof0.h[i] = 1. + random.randint(0,1)*0.1

my_model.kernel.dof0.u[:] = 0.
my_model.kernel.dof0.v[:] = 0.

##############################################################################################################

##########################################
# Run
##########################################

#my_model.meshing.plot_dev("cell")

df2d.wrapping.call_model.init_fortran(my_model.kernel)

df2d.wrapping.call_model.run(my_model.kernel, arg = "direct")

if (os.path.isdir("./obs")):
    shutil.rmtree('./obs')
#shutil.copytree("./res/obs", "./obs")



##############################################################################################################

########################
# Outputs from python
########################

# Meaning of graphe = [ 'h' , 'u' , 'v' , 'qx' , 'qy' ]

graphe = [1,0,0,0,0]

# Remark : reference not usable here

display_ref = 0

h0 = my_model.kernel.dof0.h[:nc]
u0 = my_model.kernel.dof0.u[:nc]
v0 = my_model.kernel.dof0.v[:nc]
h = my_model.kernel.dof.h[:nc]
u = my_model.kernel.dof.u[:nc]
v = my_model.kernel.dof.v[:nc]

plot_dat(graphe,display_ref,mesh_name,ts,h,u,v,h0,u0,v0)

# To save pictures : save = 1

save = 0

# Put : 'initial' or 'final'

time = 'final'

#plot_vtk(code_dir,'h',time,ts,save)
#plot_vtk(code_dir,'u',time,ts,save)
#plot_vtk(code_dir,'v',time,ts,save)
plot_vtk(code_dir,'zs','initial',ts,save)
plot_vtk(code_dir,'zs',time,ts,save)
#plot_vtk(code_dir,'porosity',time,ts,save)


df2d.wrapping.call_model.clean_model(my_model.kernel)

# Remove mesh_file

delete_mesh(mesh_name)
