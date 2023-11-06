##########################################################
##########################################################
# PERFORM A DIRECT SIMULATION WITH  DASSFLOW2D
# LAKE AT REST
#
# Introduction to basic commands to run and visualise simulation reuslt
##########################################################
##########################################################


#=======================================================#
# Source librairies
#=======================================================#
import dassflow2d as df2d
import numpy as np
import os
import importlib
from mpi4py import MPI
#=======================================================#
# copy case file
#=======================================================#

dassflow_dir = "/home/pagarambois/Documents/dassflow2d"
# or absolute path : dassflow_dir = "/home/pagarambois/Documents/Distant/dassflow2d-wrap/"

print(f"Printing case files \n from {dassflow_dir}/cases/tuto_case/0_lake-at-rest/bin_A/* \n  to  {dassflow_dir}/code/bin_A ")

# delete all files in your simulation directory before starting
os.system(f"rm -r {dassflow_dir}/code/bin_A/*")
# Copy recursively the files provided in DassFlow case repository into your own simulation directory **code/bin_A/**.
os.system(f"cp -r {dassflow_dir}/cases/tuto_case/1_lake-at-rest/bin_A/* {dassflow_dir}/code/bin_A")

os.chdir( f"{dassflow_dir}/code/")
os.system("make cleanres cleanmin")

#=======================================================#
# initialise + run +save results
#=======================================================#
df2d.wrapping.read_input(f" {dassflow_dir}/code/bin_A/input.txt")
my_model = df2d.dassflowmodel(bin_dir =  f"{dassflow_dir}/code/bin_A", hdf5_path = f"{dassflow_dir}/code/bin_A/res/simu.hdf5" , run_type = "direct", clean = True) # initialise fortran/python instance
# initialise all fortran kernel values and source them into dassflowmodel object
#my_model.init_all()
my_model.init_mesh()
my_model.kernel.dof  = df2d.wrapping.m_model.unk(my_model.kernel.mesh)
my_model.kernel.dof0 = my_model.kernel.dof
my_model.kernel.dof0.h[:] = 1
# run fortran kernel
my_model.init_fortran()
my_model.run()

my_model.save_all() # save simulation results in hdf5 files

#=======================================================#
# Post-processing
#=======================================================#

# Plot of the 2D bathymetry (input parameter of the 2D shallow water model) with package plot function

plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                            what = "bathy", 
                                            title_plot = "Bathymetry elevation",
                                            notebook = False )# for a local run remove notebook option or set notebook=False 
                                
# plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed

# Plot the friction parameter field
plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                            what = "manning_alpha", 
                                            title_scale_bar ="n [m-1/3.s] ", 
                                            title_plot = "Friction parameter (Manning coefficient)", 
                                            notebook = False )# for a local run remove notebook option or set notebook=False 
                                
# plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed

# Plot intial flow conditions 

plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                    what = "h", 
                                    when = 0,
                                    title_scale_bar ="h [m] ", 
                                    title_plot = "Initial water depth", 
                                    notebook=False) # for a local run remove notebook option or set notebook=False 
                        
# plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed

plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                    what = "zs", 
                                    when = 0,
                                    title_scale_bar ="zs [m] ", 
                                    title_plot = "Initial water surface elevation", 
                                    notebook=False) # for a local run remove notebook option or set notebook=False 
                        
# plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed

plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                    what = "u", 
                                    when = 0,
                                    title_scale_bar ="u [m/s] ", 
                                    title_plot = "Initial velocity u along x", 
                                    notebook=False) # for a local run remove notebook option or set notebook=False 
                        
# plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed


plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                    what = "v", 
                                    when = 0,
                                    title_scale_bar ="v [m/s] ", 
                                    title_plot = "Initial velocity v along y", 
                                    notebook=False) # for a local run remove notebook option or set notebook=False 
                        
# plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed

# Plot flow depth at a given time
plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                            what = "h", 
                                            when = 3,
                                            title_scale_bar ="h [m] ", 
                                            title_plot = f"Water depth at time = {my_model.outputs.result.all_time[3]}  s ", 
                                            notebook=False) # for a local run remove notebook option or set notebook=False 
                                
# plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed

# Simulation time steps at which variables have been written 
print(my_model.outputs.result.all_time)

print("previous plot for t = ", my_model.outputs.result.all_time[3])

# Compute velocity magnitude

u = my_model.outputs.result.u
v = my_model.outputs.result.v 
norm_vel = np.sqrt(u**2+v**2)

# Print the shape of the output velocity fields
print("shape of velocity array is : ", np.shape(norm_vel))
print("Maximum velocity magnitude at ecah time step is : ", np.amax(norm_vel,axis=0))


# Plot velocity magnitude at final time step
plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista, 
                                            my_scalar = norm_vel[:,-1],
                                            title_scale_bar ="norm(u,v) [m/s] ",
                                            title_plot = f"Velocity magnitude at final time",
                                            notebook=False) # for a local run remove notebook option or set notebook=False 
                                
# plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed      

# Plot water surface elevation at final time step

plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                            what = "zs",
                                            when = -1,
                                            title_scale_bar ="zs [m] ", 
                                            title_plot = f"Water surface elevation at final time", 
                                            notebook=False) # for a local run remove notebook option or set notebook=False 
                                
# plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed    
