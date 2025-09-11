import dassflow2d as df2d
import matplotlib.pyplot as plt
import shutil, os, sys
import numpy as np
import csv
from mpi4py import MPI
import json


def main():

	# Retrieve and parse additional arguments from sys.argv
	num_procs = int(sys.argv[1]) if len(sys.argv) > 1 else 1
	batch_params_str = sys.argv[2] if len(sys.argv) > 2 else '{}'
	batch_params = json.loads(batch_params_str)
	print(batch_params)
	call_dassflow(batch_params)
     
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

def call_dassflow(batch_params):
    ##############################################################################################################

    code_dir =  os.getcwd()
    bin_dir = batch_params["bin_folder"] 
    ##############################################################################################################

    ##########
    # Initialise bin
    ##########

    os.chdir(code_dir)

    #os.system("make cleanres")			   # removes all in bin_dir/res directory
    #os.system("make cleanmsh")			   # removes all in bin_dir/msh directory
    #os.system("make cleanmin")			   # removes all in bin_dir/min directory

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
    # Input Parameters
    ##########

    mesh = batch_params["mesh_name"]
    ts = batch_params["ts"]
    dtw = batch_params["dtw"]

    input_params={ "mesh_name": mesh,
                "ts": ts,
                "use_obs":'0',
                "use_UVobs":'0',
                "use_Zobs":'1',
                "use_porosity":'1',
                
                "w_obs":'1',
                "w_vtk":'0',
                "w_gnuplot":'1',
                "w_tecplot":'0',

                "adapt_dt":'1',
                'g':'9.81',
                'dt':'0.05',

                "dtw": dtw,
                "dta":"600",
                
                "heps":"0.0001",
                "eps_min":"0.00001"}

    df2d.wrapping.read_input(os.path.join(bin_dir,"input.txt"))

    ##########
    # Model
    ##########
    print("PARTIE MODEL")

    my_model = df2d.dassflowmodel(bin_dir =  bin_dir, hdf5_path = os.path.join(bin_dir,"res","simu.hdf5") , run_type = "direct", clean = False, custom_config=input_params)

    my_model.init_mesh()

    Config = df2d.core.config.Config()
    Config.set(custom_config = input_params)

    ##########
    # Friction
    ##########

    my_model.kernel.my_friction =  df2d.wrapping.m_model.friction_data(my_model.kernel.mesh)
    my_model.kernel.my_friction.nland = 1

    df2d.wrapping.call_model.init_friction(my_model.kernel)
    
    my_model.kernel.my_friction.land[:] = 1
    my_model.kernel.my_friction.manning[:] = batch_params["friction_prior"]
    my_model.kernel.my_friction.manning_beta[:] = 0


    ##########
    # Porosity
    ##########

    my_model.kernel.my_porosity =  df2d.wrapping.m_model.porosity_data(my_model.kernel.mesh)
    #porosity_values = np.loadtxt(os.path.join(bin_dir,batch_params["porosity_values_file_name"]), skiprows = 1, delimiter = ';', usecols=(1))
    #porosity_lands =  np.loadtxt(os.path.join(bin_dir,batch_params["porosity_land_file_name"]), skiprows = 1, delimiter = ';', usecols=(1))
    #porosity_gamma =  np.loadtxt(os.path.join(bin_dir,batch_params["porosity_gamma_file_name"]), skiprows = 1, delimiter = ';', usecols=(1))
    #porosity_hbanks =  np.loadtxt(os.path.join(bin_dir,batch_params["porosity_hbanks_file_name"]), skiprows = 1, delimiter = ';', usecols=(1))

    #my_model.kernel.my_porosity.nland = len(porosity_gamma)
    my_model.kernel.my_porosity.nland = 1
    df2d.wrapping.call_model.init_porosity(my_model.kernel)

    my_model.kernel.my_porosity.land[:] = 1 # porosity_lands
    #my_model.kernel.my_porosity.phi[:] = porosity_values
    my_model.kernel.my_porosity.gamma[:] = 2 # porosity_gamma
    my_model.kernel.my_porosity.hbanks[:] = 10 # porosity_hbanks
    

    ##########
    # Hydraulic states
    ##########

    my_model.kernel.dof  = df2d.wrapping.m_model.unk(my_model.kernel.mesh)
    my_model.kernel.dof0 = my_model.kernel.dof

    my_model.kernel.dof0.h[:] = 0.1
    my_model.kernel.dof0.u[:] = 0.
    my_model.kernel.dof0.v[:] = 0.

    ##############################################################################################################

    ##########################################
    # Run
    ##########################################

    #my_model.meshing.plot_dev("cell")

    df2d.wrapping.call_model.init_fortran(my_model.kernel)

    df2d.wrapping.call_model.run(my_model.kernel, arg = "direct")

    #if (os.path.isdir("./obs")):
        #shutil.rmtree('./obs')
    #shutil.copytree("./res/obs", "./obs")

    ##############################################################################################################

    df2d.wrapping.call_model.clean_model(my_model.kernel)



main()
