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

#=======================================================#
# copy case file
#=======================================================#

dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"
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

my_model = df2d.dassflowmodel(bin_dir =  f"{dassflow_dir}/code/bin_A", hdf5_path = f"{dassflow_dir}/code/bin_A/res/simu.hdf5" , run_type = "direct", clean = True) # initialise fortran/python instance
# initialise all fortran kernel values and source them into dassflowmodel object
my_model.init_all()

my_model.kernel.dof.h[:] = my_model.kernel.dof0.h[:]=1
# run fortran kernel
my_model.run()

my_model.save_all() # save simulation results in hdf5 files



#=======================================================#
# Post-processing
#=======================================================#

		
# for indication about plot_var method:
# help(my_model.plot_var)
# ~ print("Would you like to plot some model outputs")
# ~ args = input("Press Y or N to continue.") # ajouter exit si pas Y ou N

# ~ if args == "Y" or args == "y":
    # ~ my_model.plot_var(what = "bathy", when = "initial", title_plot = "bahtymetry", save_plot=True, filename = "./res/bathy")
    # ~ my_model.plot_var(what = "h", when = "initial", title_plot = "INITIAL h", save_plot=True, filename = "./res/h_0")
    # ~ my_model.plot_var(what = "zs", when = "initial", title_plot = "INITIAL zs", save_plot=True, filename = "./res/zs_0")
    # ~ my_model.plot_var(what = "u", when = "initial", title_plot = "INITIAL u", save_plot=True, filename = "./res/u_0")
    # ~ my_model.plot_var(what = "v", when = "initial", title_plot = "INITIAL v", save_plot=True, filename = "./res/v_0")
    
    # ~ # etc ...
    # ~ my_model.plot_var(what = "h", when = 0, title_plot = "INITIAL h")
    # ~ my_model.plot_var(what = "h", when = 1, title_plot = "h at second time step", save_plot=True,filename = "./res/h_fin")
    
    # ~ # result at the end of the simulation:
    # ~ my_model.plot_var(what = "vel", when = "final", title_plot = "norm(u,v) at final time step", save_plot=True,filename = "./res/velocity_fin")
    # ~ my_model.plot_var(what = "zs", when = "final", title_plot = "Zs(m) at final time step", save_plot=True,filename = "./res/Zs_fin")
    

# ~ print("Would you like to plot temporal evolution of the free surface")
# ~ args = input("Press Y or N to continue.") # ajouter exit si pas Y ou N

# ~ if args == "Y" or args == "y":
    
    # ~ for i in range(11):
        # ~ my_model.plot_var(what = "zs", when = i, title_plot = f" h  at wriite timestep= {i}") # water eight
        


# ~ print("Would you like to plot temporal evolution of the velocitSy")
# ~ args = input("Press Y or N to continue.") # ajouter exit si pas Y ou N

# ~ if args == "Y" or args == "y":
    
    # ~ for i in range(11):
        # ~ my_model.plot_var(what = "vel", when = i, title_plot = f" norm(u,v) at at wriite timestep= {i}")
        

