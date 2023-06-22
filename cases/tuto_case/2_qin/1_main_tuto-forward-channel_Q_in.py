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
#=======================================================#
# copy existing case files
#=======================================================#

#os.chdir('/home/pagarambois/Documents/Distant/dasshydro/cases/tuto_case/2_qin')
#dassflow_dir = "/home/pagarambois/Documents/Distant/dassflow2d-wrap"

dassflow_dir = os.getcwd() # os.path.abspath(os.curdir)
if dassflow_dir.split('/')[-1] != "dasshydro":
   os.chdir("../../..")
   dassflow_dir = os.getcwd()

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

my_model = df2d.dassflowmodel(bin_dir = run_dir,hdf5_path = f"{run_dir}/res/simu.hdf5", run_type = "direct",clean=True) # initialise fortran/python instance
my_model.init_all() # allocate and initialise many fortran variables
my_model.run() # run model
my_model.save_all() # save simulation results in hdf5 files

bathy = my_model.outputs.all_res[0.0][["bathy"]]
my_scalar = bathy
h = my_model.outputs.all_res[0.0][["h"]]
labels = dict(xlabel='X [m]', ylabel='Y [m]', zlabel='')
plotter = my_model.meshing.plot(my_scalar,
                                     title_scale_bar ="Zb [m] ", 
                                     title_plot = "bathymetry elevation on 2D mesh grid", 
                                     axis_labels = labels)

plotter = my_model.meshing.plot(my_scalar=h)

allx =[]
allz = []
for i in range(my_model.meshing.mesh_fortran.nc):
    x =my_model.meshing.mesh_fortran.cell[i].grav.x
    y = my_model.meshing.mesh_fortran.cell[i].grav.y
    print(y)
    if(y==50.0):
        allx.append(x)
        allz.append(my_model.outputs.all_res[0.0][["bathy"]].iloc[i-1])

time_out = my_model.outputs.all_times[2]
v = my_model.outputs.all_res[time_out][["v"]]
        



#plot en passant par autre fonction... via struct de vars...
# my_model.

for key, value in my_model.outputs.all_res.items():
    tmp = value["h"]

    # my_model.meshing.plot()

    my_model.meshing.mesh_pyvista.plot(scalars = tmp, show_edges=True, cpos= "xy", notebook =False)

time = 0.0
h0 = my_model.outputs.all_res[time][["h"]]
u0 = my_model.outputs.all_res[time][["u"]]
v0 = my_model.outputs.all_res[time][["v"]]

plotter = my_model.meshing.plot(my_scalar = U,
                                     title_scale_bar ="", 
                                     title_plot = f"Initial", 
                                     xlabel = "X [m]", 
                                     ylabel = "Y [m]") # for a local run remove notebook option or set notebook=False 
plotter.show() # remove jupyter_backend if needed

#=======================================================#
# Save results as observed data
#=======================================================#
os.system("rm ../../code/bin_A/obs/*")
os.system("cp ../../code/bin_A/res/obs/* bin_A/obs/")
