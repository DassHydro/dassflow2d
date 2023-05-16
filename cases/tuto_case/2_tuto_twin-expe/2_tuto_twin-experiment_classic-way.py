##################################################################
##################################################################
# PERFORM TWIN EXPERIMENT (DIRECT RUN + INFERENCE RUN)
#
# Method "old way" using command system and pre-existing files
##################################################################
##################################################################


#=======================================================#
#  define librairies librairies
#=======================================================#

import dassflow2d as df2d              # dassflow2d : main package
import os                        # os, for shell command execution (mainly for file manipulation)
import numpy as np
import matplotlib.pyplot as plt
# for os.chdir() command, to open current birectory as bin_directory, defined in parameters, just below

from subprocess import call # terminal comand 

import h5py # source hdf5 files


# method that copy pre-existing files to bin dir, used in the script
# source is the path of the source file and target is the path where to copy the file
# it directly use the shell
def cp_inference_file(source, target):
    source = f'inference_files/{source}'
    target = f'{target}'
    call(['rm', target]) 
    call(['cp', source, target]) # Linux
              
    
# copy all reference files (for direct run)
def cp_reference_file():

    cp_inference_file(source = 'hydrograph_target.txt' ,  
                      target = 'hydrograph.txt')
    
    cp_inference_file(source = 'mesh_target.txt' ,  
                      target = 'automaticaly_generated_mesh.txt')
    
    cp_inference_file(source = 'land_uses_target.txt' ,  
                      target = 'land_uses.txt')
    
    cp_inference_file(source = 'input_direct.txt' ,  
                      target = 'input.txt')

# cp_dir copy files from source directory to target directory.
# note that the directory is copied (not only the files within)
def cp_dir(source, target):
  call(['cp', '-a', source, target]) # Linux

# save hdf5 (results) files out of /res directory
def cp_hdf5_file(source, target):
    source = f'res/{source}'
    target = f'hdf5/{target}'
#    call(['mkdir', "hdf5"])
#    call(['rm', target]) 
    call(['cp', source, target]) # Linux


# define paths
dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"
bin_dir = f"{dassflow_dir}/code/bin_A"


# delete all files in your simulation directory before starting
os.system(f"rm -r {dassflow_dir}/code/bin_A/*")
# Copy recursively the files provided in DassFlow case repository into your own simulation directory **code/bin_A/**.
os.system(f"cp -r {dassflow_dir}//cases/tuto_case/2_tuto_twin-expe/bin_A/* {dassflow_dir}/code/bin_A")
os.chdir( f"{dassflow_dir}/code/")
os.system(f"make cleanres cleanmin ")

#=======================================================#
# set up true configuration
#=======================================================#

os.chdir(bin_dir)
cp_reference_file()

#=======================================================#
# direct run of the  model
#=======================================================#A


# initialise fortran instance, and python corrponding data
direct_model = df2d.dassflowmodel(bin_dir = bin_dir, hdf5_path = f"{dassflow_dir}/code/bin_A/res/simu.hdf5", run_type = "direct", clean = True)
# then intialise meshing
direct_model.init_all()
# define initial conditions
direct_model.kernel.dof0.h[:] = 1
direct_model.kernel.dof0.u[:] = 0
direct_model.kernel.dof0.v[:] = 0
direct_model.kernel.dof = direct_model.kernel.dof0

direct_model.run()



direct_model.save_all() # save simulation results in hdf5 files

cp_hdf5_file(source = "simu.hdf5", target = "true.hdf5") #save hdf5 (results) files out of /res directory
cp_dir('./res/obs', '.') # copy observation files


#----------- plot
#direct_model.plot_var(what = "bathy", when = "initial", title_plot = "bahtymetry") 
#direct_model.plot_var(what = "h", when = "initial", title_plot = "h")
#direct_model.plot_var(what = "u", when = "initial", title_plot = "u")
#direct_model.plot_var(what = "v", when = "initial", title_plot = "v")
#direct_model.plot_var(what = "manning_alpha", when = "initial", title_plot = "Manning-Strickler coefficient")


a = h5py.File(name = f"{bin_dir}/res/simu.hdf5")



for i in range(a["output/spatial_var/bathy"].shape[1]):
    b=a["output/spatial_var/bathy"][:,i] 
    zs = a["output/spatial_var/h"][:,i] +  a["output/spatial_var/bathy"][:,0] 
    
    plt.scatter(x = np.arange(len(zs )),
                y = zs, 
                linewidths = 1 )
    
    plt.scatter(x = np.arange(len(b )),
                y = b, 
                c="red", 
                linewidths = 1 )
    plt.show()





#bc = df2d.m_model.get_bc()
df2d.wrapping.call_model.clean_model(direct_model.kernel)         # deallocate correctly (necessary action)


# clean eventual hdf5 file open
import gc
for obj in gc.get_objects():   # Browse through ALL objects
    if isinstance(obj, h5py.File):   # Just HDF5 files
        try:
            obj.close()
        except:
            pass # Was already closed

###########################################################
#===========================================================
# RUN INFERENCE
#===========================================================
###########################################################


#----------------------#
#  Define Parameters
#----------------------#

os.chdir( f"{dassflow_dir}/code/")
os.system(f"make cleanres cleanmin ")



os.chdir(bin_dir)
os.system(f"rm restart.bin ")


# /!\ warning :: bathymetry not inferable ? --> lilian = optim not find optimum
print("CHOOSE INFERENCE TYPE (1 hydrograph, 2 land_use, 3 = bathy)")
inference_type = input("Enter 1,2 or 3 \n")

if inference_type == "1":
    # -- infer hydrograph --
    cp_reference_file()
    cp_inference_file(source = 'hydrograph_prior.txt' ,
                      target = 'hydrograph.txt')

    cp_inference_file(source = 'input_hydrograph.txt' ,
                      target = 'input.txt')
elif inference_type == "2":
    # -- same for manning -- #
    cp_reference_file()
    cp_inference_file(source = 'land_uses_prior.txt' ,
                      target = 'land_uses.txt')

    cp_inference_file(source = 'input_land_uses.txt' ,
                      target = 'input.txt')

elif inference_type == "3":
    # -- infer hydrograph --#
    cp_reference_file()
    cp_inference_file(source = 'mesh_prior.txt' ,
                          target = 'automaticaly_generated_mesh.txt.txt')

    cp_inference_file(source = 'input_bathy.txt' ,
                          target = 'input.txt')




#=======================================================#
# Inference
#=======================================================#


my_model = df2d.dassflowmodel(bin_dir = bin_dir, hdf5_path = f"{dassflow_dir}/code/bin_A/res/simu.hdf5", run_type = "min")

# then intialise meshing
my_model.init_all()
# define initial conditions
my_model.kernel.dof0.h[:] = 1
my_model.kernel.dof0.u[:] = 0
my_model.kernel.dof0.v[:] = 0
my_model.kernel.dof = my_model.kernel.dof0
my_model.run() # only inference is performed

# peroform simulation with infered parameted, to perform the last lines of the script.

## save simulation results in hdf5 files
#my_model.save_res()
#
#cp_hdf5_file(source = "simu.hdf5", target = "inference_manning.hdf5")
#
#inf = h5py.File(name =  f"{bin_dir}/hdf5/inference_manning.hdf5", mode = "r")
#
#
##=======================================================#
## LOAD HDF5 files
##=======================================================#
#
## get grid object
#my_model.build_grid()
#grid = my_model.grid
#grid.plot(show_edges =True, cpos = "xy")
#
#true = h5py.File(name = f"{bin_dir}/hdf5/true.hdf5", mode = "r")
#
#plot_title="null"
#if inference_type == "1":
#    # -- infer hydrograph --
#    true =  [20,20]
#    inf = np.loadtxt("min/hydrograph_001.006")[:,1]
#    prior = np.loadtxt("min/hydrograph_001.000")[:,1]
#    plot_title="Inference of hydrograph parameter: Q(m3/s)"
#elif inference_type == "2":
#    # -- same for manning -- #
#    true =  true["spatial_var/manning_alpha"][:,0]
#    inf = np.loadtxt("min/manning.052")[:,1]
#    prior = np.loadtxt("min/manning.000")[:,1]
#    plot_title="inference of manning parameter: n"
#
#elif inference_type == "3":
#    true =  true["spatial_var/bathy"][:,0]
#    inf = np.loadtxt("min/bathy.039")[:,1]
#    prior = np.loadtxt("min/bathy.000")[:,1]
#    plot_title="inference of bathymetry parameter (m)"
#
#
#
#plt.scatter( x = np.arange(len(true)), y = true, label = "target")
#plt.scatter( x = np.arange(len(inf)), y = inf, label = "infered")
#plt.scatter( x = np.arange(len(prior)), y =prior, label ="prior")
#plt.title(plot_title)
#plt.legend()

df2d.wrapping.call_model.clean_model(my_model.kernel)         # deallocate correctly (necessary action)
