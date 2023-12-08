#!/usr/bin/env python
# coding: utf-8

# In[1]:


#####################################################################
#####################################################################
# PERFORM A DIRECT SIMULATION WITH  DASSFLOW2D
# LAKE AT REST
#
# Introduction to basic commands of run and visualisation of results
#####################################################################
#####################################################################


#=======================================================#
# Source librairies
#=======================================================#

import dassflow2d as df2d
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#=======================================================#
# copy of case files
#=======================================================#

os.chdir('../../')
dassflow_dir = os.getcwd() # DassFlow directory (you can also impose your absolute path)

#dassflow_dir="/home/pagarambois/Documents/Distant/dassflow2d"
os.chdir(dassflow_dir)
print("DassFlow directory is: ", dassflow_dir)

# Define directory where case is run
# (its name 'bin_A' is imposed in  {dassflow_dir}/code/makefile.inc : CASEDIR='bin_A')
run_dir = f"{dassflow_dir}/code/bin_A/"

# Define directory containing case data
case_data_dir = f"{dassflow_dir}/cases/tuto_case/1_lake-at-rest/bin_A/"

# Clean run directory
os.system(f"rm -r {run_dir}*")

# Copy case data to runing directory
os.system(f"cp -r {case_data_dir}* {run_dir}") # Copy of case files from existing case to bin_A

# Move to code directory and clean bin directory
os.chdir( f"{dassflow_dir}/code/")
os.system("make cleanres cleanmin") # Clean forward run and minimization results


# In[2]:


#=======================================================#
# Initialisation
#=======================================================#

# input file reading (simulation settings)
df2d.wrapping.read_input(f"{run_dir}/input.txt")

# Creation of dassflowmodel object using case data:
df2d.wrapping.m_mpi.init_mpi() #set the number of processes to 1
my_model = df2d.dassflowmodel(bin_dir =  run_dir, hdf5_path = f"{dassflow_dir}/code/bin_A/res/simu.hdf5" , run_type = "direct", clean = True, custom_config = None)

my_model.config.get()

# Initializion of the Fortran kernel (dassflow Python library is obtained by wrapping Fortran source code)
#initialise all fortran kernel values and source them into dassflowmodel object

my_model.init_all()

my_model.kernel.dof0.h[:] = 1
my_model.kernel.dof.h[:] = 1


# In[3]:


#=======================================================#
# Run Fortran kernel and save results
#=======================================================#

my_model.run()
my_model.save_all()


# In[4]:


#=======================================================#
# Vizualize parameters and results
#=======================================================#


# In[5]:


# Plot of the 2D bathymetry (input parameter of the 2D shallow water model) with package plot function

#plotter = my_model.plot_var(my_model.meshing.mesh_pyvista,
#                                             what = "bathy",
#                                             title_plot = "Bathymetry elevation")# for a local run remove notebook option or set notebook=False

#plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed


# In[6]:


# Plot the friction parameter field
plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                             what = "manning_alpha",
                                             title_scale_bar ="n [m-1/3.s] ",
                                             title_plot = "Friction parameter (Manning coefficient)",
                                             notebook = True )# for a local run remove notebook option or set notebook=False

plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed


# In[7]:


# Plot intial flow conditions

plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                             what = "h",
                                             when = 0,
                                             title_scale_bar ="h [m] ",
                                             title_plot = "Initial water depth",
                                             notebook=True) # for a local run remove notebook option or set notebook=False

plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed

plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                             what = "zs",
                                             when = 0,
                                             title_scale_bar ="zs [m] ",
                                             title_plot = "Initial water surface elevation",
                                             notebook=True) # for a local run remove notebook option or set notebook=False

plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed

plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                             what = "u",
                                             when = 0,
                                             title_scale_bar ="u [m/s] ",
                                             title_plot = "Initial velocity u along x",
                                             notebook=True) # for a local run remove notebook option or set notebook=False

plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed


plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                             what = "v",
                                             when = 0,
                                             title_scale_bar ="v [m/s] ",
                                             title_plot = "Initial velocity v along y",
                                             notebook=True) # for a local run remove notebook option or set notebook=False

plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed


# In[8]:


# Plot flow depth at a given time
plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                             what = "h",
                                             when = 3,
                                             title_scale_bar ="h [m] ",
                                             title_plot = f"Water depth at time = {my_model.outputs.result.all_time[3]}  s ",
                                             notebook=True) # for a local run remove notebook option or set notebook=False

plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed

# Simulation time steps at which variables have been written
print(my_model.outputs.result.all_time)

print("previous plot for t = ", my_model.outputs.result.all_time[3])


# In[9]:


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
                                             notebook=True) # for a local run remove notebook option or set notebook=False

plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed

# Plot water surface elevation at final time step

plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                             what = "zs",
                                             when = -1,
                                             title_scale_bar ="zs [m] ",
                                             title_plot = f"Water surface elevation at final time",
                                             notebook=True) # for a local run remove notebook option or set notebook=False

plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed


# In[10]:


print(my_model.config)

print("The numerical scheme that has been used to solve the 2D shallow water equations is:")
print("Temporal scheme is: ", my_model.config["temp_scheme"])
print("Spatial scheme is: ", my_model.config["spatial_scheme"])


# In[11]:


#clean model
df2d.wrapping.call_model.clean_model(my_model.kernel)

