#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 16:43:57 2022

@author: livillenave
"""
# PERFORM FULL VISUALISATION

#=======================================================#
# Source librairies
#=======================================================#
import dassflow2d as df2d
import numpy as np
import os
import importlib
import pyvista as pv

#=======================================================#
# copy case file
#=======================================================#

dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"
# or absolute path : dassflow_dir = "/home/pagarambois/Documents/Distant/dassflow2d-wrap/"

print(f"Printing case files \n from {dassflow_dir}/cases/tuto_casee/2_tuto_twin-expe/bin_A/* \n  to  {dassflow_dir}/code/bin_A ")

# delete all files in your simulation directory before starting
#os.system(f"rm -r {dassflow_dir}/code/bin_A/*")
# Copy recursively the files provided in DassFlow case repository into your own simulation directory **code/bin_A/**.
#os.system(f"cp -r {dassflow_dir}/cases/tuto_case/2_tuto_twin-expe/bin_A/* {dassflow_dir}/code/bin_A")
os.chdir( f"{dassflow_dir}/code/")
os.system("make cleanres cleanmin")

#=======================================================#
# initialise + run 
#=======================================================#

my_model = df2d.dassflowmodel(bin_dir =  f"{dassflow_dir}/code/bin_A" , run_type = "direct", clean =True)  # initialise fortran/python instance
# then intialise meshing
my_model.init_mesh()
# initialise dof structure
my_model.init_dof()	
# initialise remaining structures
my_model.init_fortran()
my_model.run()

my_model.outputs.save_res(name_hdf5_file="test")

my_model.boundary.plot(what = "values")
my_model.boundary.plot(what = "meshing")
#=======================================================#
# Visualise configuration 
#=======================================================#
my_model.config.get()
print(my_model.config)

#=======================================================#
# Visualise meshing 
#=======================================================#
pv.plot(my_model.grid, show_bounds=True, show_edges=True, cpos = "xy")

my_model.plot_meshing("node")
my_model.plot_meshing("cell")
my_model.plot_meshing("edge")


#=======================================================#
# Visualise bc 
#=======================================================#
my_model.plot_bc_mesh(title_plot = "test", save_plot =True)
my_model.plot_bc_values(save_plot =True)

#=======================================================#
# Visualise initial conditions
#=======================================================#
#save result in hdf5 is a necessary step to access variables
my_model.save_res(filename = "simu")


my_model.plot_var(path_hdf5_file=f"{my_model.bin_dir}/res/simu.hdf5",
                  what = "h",
                  when = "initial") # 0 instead of "initial" also work
my_model.plot_var(path_hdf5_file=f"{my_model.bin_dir}/res/simu.hdf5",
                  what = "zs",
                  when = "initial") # 0 instead of "initial" also work

my_model.plot_var(path_hdf5_file=f"{my_model.bin_dir}/res/simu.hdf5",
                  what = "u",
                  when = "initial") # 0 instead of "initial" also work

my_model.plot_var(path_hdf5_file=f"{my_model.bin_dir}/res/simu.hdf5",
                  what = "v",
                  when = "initial") # 0 instead of "initial" also work

#=======================================================#
# Visualise parameters
#=======================================================#
#hdf5 is a necessary to access variables
my_model.plot_var(path_hdf5_file=f"{my_model.bin_dir}/res/simu.hdf5",
                  what = "bathy",
                  when = "initial") # 0 instead of "initial" also work

my_model.plot_var(path_hdf5_file=f"{my_model.bin_dir}/res/simu.hdf5",
                  what = "manning_alpha",
                  when = "initial") # 0 instead of "initial" also work




my_model.save_res() # save simulation results in hdf5 files
my_model.build_grid()  # build a pyvista.unstructuredgrid object, which is used for plots



bc_ref = df2d.wrapping.call_model.boundaries_copy(python_interface.boundary.bc)

bc_ref.hyd[0].q[:] = 333.33333333333
python_interface.boundary.bc.hyd[0].q[:]  = bc_ref.hyd[0].q[:]

