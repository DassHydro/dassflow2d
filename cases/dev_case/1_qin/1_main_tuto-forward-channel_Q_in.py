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

#=======================================================#
# copy existing case files
#=======================================================#

dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"


print(f"Printing case files \n from {dassflow_dir}/cases/tuto_case/0_lake-at-rest/bin_A/* \n  to  {dassflow_dir}/code/bin_A ")

# delete all files in your simulation directory before starting
os.system(f"rm -r {dassflow_dir}/code/bin_A/*")
# Copy recursively the files provided in DassFlow case repository into your own simulation directory **code/bin_A/**.
os.system(f"cp -r {dassflow_dir}/cases/tuto_case/1_qin/bin_A/* {dassflow_dir}/code/bin_A")
os.chdir( f"{dassflow_dir}/code/")
os.system("make cleanres cleanmin")
#=======================================================#
# initialise + run +save results
#=======================================================#

my_model = df2d.DassFlowModel(bin_dir =  f"{dassflow_dir}/code/bin_A" , run_type = "direct") # initialise fortran/python instance
my_model.update_fortran() # allocate and initialise many fortran variables
my_model.run() # run model
my_model.save_res() # save simulation results in hdf5 files
my_model.build_grid()  # build a pyvista.unstructuredgrid object, which is used for plots




#=======================================================#
# Save results as observed data
#=======================================================#
os.system("rm ../../code/bin_A/obs/*")
os.system("cp ../../code/bin_A/res/obs/* bin_A/obs/")
