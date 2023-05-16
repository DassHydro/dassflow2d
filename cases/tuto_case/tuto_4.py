#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 11:23:25 2022

@author: livillenave
"""



# in this tutorial, we show how to acess mesh (NO: and update boundary condictions)
# the same logic as for the first 3 tutorials is applied to end up to the 
#provided framework
#

import dassflow2d as df2d
import pyvista as pv
import numpy as np


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

# now, you stored the mesh object
# as a pyvista unstructured grid object,
# you can now go to pyvista documentation to take full advantage of the object
# some plot methods are already implemented:
# but the most basic  functionalities are also showed below
    # ----------- basic plots --------------------
    # basic plots
model.meshing.mesh_pyvista.plot()
model.meshing.mesh_pyvista.plot(show_edges=True, cpos= "xy", notebook =False)


    # basic plots using plotter instances
sargs = dict(vertical = True, title = "title scale bar", position_y=0.3, position_x = 0.8  )
plotter = pv.Plotter(notebook = False)
# mesh plot
plotter.add_mesh(model.meshing.mesh_pyvista, scalars = np.random.rand(model.kernel.mesh.nc), scalar_bar_args=sargs,
		 show_edges = True)
plotter.view_xy()
plotter.show_bounds()
plotter.add_title("title plot", font_size=18, color=None, font=None, shadow=False)
plotter.show()


# ----------- Existing plots --------------------

model.meshing.plot()

model.meshing.plot_dev(what = "cell")
model.meshing.plot_dev(what = "node")
model.meshing.plot_dev(what = "edge")



# the perform run, using the code, fully wrapped within python calls:
model.init_dof()
model.init_fortran()

df2d.wrapping.call_model.run(model.kernel, arg = "direct")
df2d.wrapping.call_model.clean_model(model.kernel)

