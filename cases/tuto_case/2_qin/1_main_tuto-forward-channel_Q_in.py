####################################################################################################################
####################################################################################################################
# PERFORM A DIRECT SIMULATION WITH  DASSFLOW2D AND GENERATE OBSERVED RESULTS
# Q in
#
# In addition, compared to the "lake at rest" test case, here we generate observed data to perform a twin experiment
# ----> w_obs=1  in input.txt and obs.txt file is provided
####################################################################################################################
####################################################################################################################

#=======================================================#
# Source librairies
#=======================================================#

import dassflow2d as df2d
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpi4py import MPI
import h5py                # for save file
import pandas as pd        # for read tables

#=======================================================#
# copy existing case files
#=======================================================#

df2d.wrapping.m_mpi.init_mpi()

dassflow_dir = os.path.abspath(os.path.join(__file__ ,"../../../.."))

case_dir = os.path.join(f"{dassflow_dir}","cases/tuto_case/2_qin")
run_dir = os.path.join(f"{dassflow_dir}","code/bin_A")

print("Dassflow directory is understood as: ", dassflow_dir)
print("Case is copied from: ", case_dir)
print("Running directory is: ", run_dir)

# delete all files in your simulation directory before starting
os.system(f"rm -r {dassflow_dir}/code/bin_A/*")
# Copy recursively the files provided in DassFlow case repository into your own simulation directory **code/bin_A/**.
os.system(f"cp -r {dassflow_dir}/cases/tuto_case/2_qin/bin_A/* {dassflow_dir}/code/bin_A")

os.chdir( f"{dassflow_dir}/code/")
os.system("make cleanres cleanmin")

#=======================================================#
# initialise + run + save results
#=======================================================#

df2d.wrapping.read_input(f" {dassflow_dir}/code/bin_A/input.txt")

my_model = df2d.dassflowmodel(bin_dir = run_dir,hdf5_path = f"{run_dir}/res/simu.hdf5", run_type = "direct",clean=True) # initialise fortran/python instance

my_model.init_all() # allocate and initialise many fortran variables

my_model.run() # run model

my_model.save_all() # save simulation results in hdf5 files

#Acces stored data (vectors)
# bathy = my_model.outputs.result.bathy#[0.0][["bathy"]]
# my_scalar = bathy
# h = my_model.outputs.result.h
# labels = dict(xlabel='X [m]', ylabel='Y [m]', zlabel='')
# plotter = my_model.meshing.plot(my_scalar,
#                                      title_scale_bar ="Zb [m] ", 
#                                      title_plot = "Bathymetry elevation on 2D mesh grid", 
#                                      xlabel = labels["xlabel"],
#                                      ylabel = labels["ylabel"],
#                                      zlabel = labels["zlabel"])

# plotter = my_model.meshing.plot(my_scalar=h)


#Read hdf5 save files (tables)
my_hdf5_file =  h5py.File(f"{run_dir}/res/simu.hdf5", "r")

print(list(my_hdf5_file.keys()))
print(list(my_hdf5_file["input"].keys()))
print(list(my_hdf5_file["output"].keys()))
print(list(my_hdf5_file["input"]["meshing"].keys()))
print(list(my_hdf5_file["input"]["boundary"].keys()))
print(list(my_hdf5_file["output"]["result"].keys()))

bathy = my_hdf5_file["output"]["result"]["bathy"]

allx =[]
allz = []
for i in range(my_model.meshing.mesh_fortran.nc):
    x = my_model.meshing.mesh_fortran.cell[i].grav.x
    y = my_model.meshing.mesh_fortran.cell[i].grav.y

    if(y==50.0):
        allx.append(x)
        allz.append(bathy[i])

plt.plot(allx, allz, 'k')
plt.title("Longitudinal bathymetry")
plt.show()


#Plot water depth at each saved time with Meshing.mesh_pyvista
for i in range(0,my_hdf5_file["output"]["result"]["h"].shape[1]):
    tmp = my_hdf5_file["output"]["result"]["h"][:,i]

    my_model.meshing.mesh_pyvista.plot(scalars = tmp, show_edges=True, cpos= "xy", notebook =False)


#Plot v at a given time through Meshing.plot
v = my_hdf5_file["output"]["result"]["v"][:,2]
plotter = my_model.meshing.plot(my_scalar = v,
                                     title_scale_bar ="", 
                                     title_plot = f"Initial", 
                                     xlabel = "X [m]", 
                                     ylabel = "Y [m]") # for a local run remove notebook option or set notebook=False 
#plotter.show() # remove jupyter_backend if needed

#=======================================================#
# Save results as observed data
#=======================================================#
os.system("rm ../../code/bin_A/obs/*")
os.system("cp ../../code/bin_A/res/obs/* bin_A/obs/")
