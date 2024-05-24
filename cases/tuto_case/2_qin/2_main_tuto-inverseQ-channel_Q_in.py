####################################################################################################################
####################################################################################################################
# PERFORM AN INFERENCE WITH  DASSFLOW2D
# Q in
#
# In addition, compared to the "lake at rest" test case, here we generate observed data to perform a twin experiment
# ----> w_obs=1, use_obs = 1  in input.txt and obs.txt file is provided
####################################################################################################################
####################################################################################################################


#=======================================================#
# Source librairies
#=======================================================#

import dassflow2d as df2d
import numpy as np
import os
from mpi4py import MPI

#=======================================================#
# copy existing case files
#=======================================================#

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

#=======================================================#
# Direct simulation and save results
#=======================================================#
# initialise fortran instance, and python corrponding data
df2d.wrapping.m_mpi.init_mpi()

df2d.wrapping.read_input(f" {dassflow_dir}/code/bin_A/input.txt")

my_model = df2d.dassflowmodel(bin_dir = run_dir,hdf5_path = f"{run_dir}/res/simu.hdf5", run_type = "direct",clean=True)

my_model.init_all() # allocate and initialise many fortran variables

# run model
my_model.run()
# save simulation results in hdf5 files
my_model.save_all()
#
#my_model.build_grid() # necessary for plots # builds callable objects
#my_model.grid


#=======================================================#
# Prepare twin experiment for hydrograph inference from water levels observations
#=======================================================#

os.system(f"rm {dassflow_dir}/code/bin_A/obs/*")
os.system(f"cp {dassflow_dir}/code/bin_A/res/obs/* bin_A/obs/")

os.system(f"rm {dassflow_dir}/code/bin_A/hydrograph.txt")                                    # delete the "true" hydrograph used in a run above to generate water level observations
os.system(f"cp {dassflow_dir}/code/bin_A/hydrograph_first_guess.txt  bin_A/hydrograph.txt")  # define first guess on hydrograph for inference

print("Observation files and first guess hydrograph copied for twin experiment")
wait = input("Press Enter to continue.")


#=======================================================#
# Inference of Q_in from water levels observations
#=======================================================#

df2d.wrapping.read_input(f" {dassflow_dir}/code/bin_A/input_inverse.txt")

my_model_inferQ_in = df2d.dassflowmodel(bin_dir = run_dir, hdf5_path = f"{run_dir}/res/simu.hdf5", run_type = "min",clean=True)
my_model_inferQ_in.run()

#=======================================================#
# Post-processing
#=======================================================#


my_model.plot_var(what = "bathy", when = "initial", title_plot = "Bathymetry")
my_model.plot_var(what = "h", when = "initial", title_plot = "INITIAL h")
my_model.plot_var(what = "zs", when = "initial", title_plot = "INITIAL zs")
# etc ...
my_model.plot_var(what = "h", when = 0, title_plot = "INITIAL h")
my_model.plot_var(what = "h", when = 1, title_plot = "h at second time step")
