#==============================================================#
#==============================================================#
# show various mesh method (for ploting and return variables)
#==============================================================#
#==============================================================#

import dassflow2d as df2d
import h5py
import os




#=======================================================#
#  EXAMPLE
#=======================================================#
input_param= { "ts":2000 } 
bin_dir =  "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A"
        
import os
os.chdir(bin_dir)
# initialise fortran instance, and python corrponding data
my_model = df2d.DassFlowModel(bin_dir = bin_dir, 
                              run_type = "direct") # run_type can be min or direct (grad ?)

# see df2d.DassFlowModel.build_grid()
pyvista_arrays =  my_model.build_grid(what = "arrays") # [cells, cell_type, points])
pyvista_grid =  my_model.build_grid(what = "pyvista") # [cells, cell_type, points])

# see help(df2d.DassFlowModel.plot_meshing)
edge = my_model.plot_meshing("edge")
cell = my_model.plot_meshing("cell")
node = my_model.plot_meshing("node")


my_model.update_fortran()
my_model.run()

my_model.save_res()

df2d.call_model.clean_model(my_model.model)

