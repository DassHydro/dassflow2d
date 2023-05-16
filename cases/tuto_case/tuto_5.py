#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 14:42:24 2022

@author: livillenave
"""

import dassflow2d as df2d
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt

#----paths ----

# main directory
dassflow_dir="/home/livillenave/Documents/distant/dassflow2d-wrap/"
# code directory, where compilation happens
code_dir =  f"{dassflow_dir}/code/"
# bin directory, where simulation happens
bin_dir = f"{code_dir}/bin_A/"


#---------------------open fortran interfacing ------------------------------
# open interface with fortran kernel
model = df2d.dassflowmodel(dassflow_dir=dassflow_dir, bin_dir=bin_dir )
# then intialise meshing
model.init_mesh()

# the perform run, using the code, fully wrapped within python calls:
model.init_dof()
model.init_fortran()

#model.boundary.get_values(model)
model.boundary.bc
# The boundary conditions are made accessible
counts = model.boundary.get_metadata()

# this correspond to the id of the edge corresponding to the bc
model.boundary.corresp["discharg1"][1] 

# the values of bc are made accessible here
# for first hydrograph
print(model.boundary.bc.hyd[0].t)
print(model.boundary.bc.hyd[0].q)


model.boundary.bc.hyd[0].q = 100

bc_ref = df2d.wrapping.call_model.boundaries_copy(model.boundary.bc)
bc_ref.hyd[0].q = 800
# for second hydrograph, etc...
#print(model.boundary.bc.hyd[1].t)
#print(model.boundary.bc.hyd[1].q)


# plot bc using the mesh and the boundary correspondance 
model.boundary.plot(what = "meshing", save_plot = False)
model.boundary.plot(what = "values", save_plot = False)

model.run()
model.Output
model.outputs.save_res(hdf5_path="/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/res/simu.hdf5")
# you also can visualise values easily