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

#=======================================================#
# copy existing case files
#=======================================================#
os.system("rm -r ../../code/bin_A/*")
os.system("cp -r ./2_simple_channel_mcdonald/* ../../code/bin_A")

#=======================================================#
# Direct simulation and save results
#=======================================================#
# initialise fortran instance, and python corrponding data
my_model = df2d.DassFlowModel(bin_dir = "../../code/bin_A", arg = "direct")
# run model
my_model.run()
# save simulation results in hdf5 files
my_model.save_res()

my_model.build_grid() # necessary for plots # builds callable objects
my_model.grid


#=======================================================#
# Prepare twin experiment for hydrograph inference from water levels observations
#=======================================================#

os.system("rm ../../code/bin_A/obs/*")
os.system("cp ../../code/bin_A/res/obs/* bin_A/obs/")

os.system("rm ../../code/bin_A/hydrograph.txt")                                    # delete the "true" hydrograph used in a run above to generate water level observations
os.system("cp ../../code/bin_A/hydrograph_first_guess.txt  bin_A/hydrograph.txt")  # define first guess on hydrograph for inference

print("Observation files and first guess hydrograph copied for twin experiment")
wait = input("Press Enter to continue.")



#=======================================================#
# Inference of Q_in from water levels obserationza
#=======================================================#
my_model_inferQ_in = df2d.DassFlowModel(bin_dir = "../../code/bin_A",
                                run_type = "min")



#=======================================================#
# Post-processing
#=======================================================#



my_model.plot_var(what = "bathy", when = "initial", title_plot = "bahtymetry")
my_model.plot_var(what = "h", when = "initial", title_plot = "INITIAL h")
my_model.plot_var(what = "zs", when = "initial", title_plot = "INITIAL zs")
# etc ...
my_model.plot_var(what = "h", when = 0, title_plot = "INITIAL h")
my_model.plot_var(what = "h", when = 1, title_plot = "h at second time step")
