#=======================================================#
# Source librairies
#=======================================================#
import dassflow2d as df2d
import numpy as np
import sys
import os



run_type = "direct"
bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A"




#=======================================================#
#  Force configuration for direct, and prepare inference configuration
#=======================================================#
config_direct = {"use_obs":0, "w_obs":1  }
config_inf    = {"use_obs":1, "w_obs":1, "c_hydrograph":1  } 

#=======================================================#
#  run direct fortran model
#=======================================================#
os.chdir(bin_dir)
os.system("rm hydrograph.txt")
os.system("cp hydrograph_true.txt hydrograph.txt")

model = df2d.dassflowmodel(bin_dir = os.getcwd(), 
                           hdf5_path=os.getcwd() + "/res/simu.hdf5", 
                           run_type = "direct", 
                           clean = True)

model.config.set(config_dictionary_to_force = config_direct)
model.config.get()

model.init_all()

model.kernel.dof.h[:]=1
dof_init = df2d.wrapping.call_model.dof_copy(model.kernel.dof)
model.kernel.dof0 = df2d.wrapping.call_model.dof_copy(dof_init)


model.run()
model.save_all()

model.obs = df2d.core.obs.Obs(bin_dir = model.bin_dir,
                              ts = model.config["ts"])
model.obs.sim_to_true()

df2d.wrapping.call_model.clean_model(model.kernel)

os.chdir("../.")
os.system(f"cp -R {bin_dir} {bin_dir}_true")

#=======================================================#
#  run inf fortran model
#=======================================================#
os.chdir(bin_dir)

os.system("rm hydrograph.txt")
os.system("cp hydrograph_prior.txt hydrograph.txt")

model = df2d.dassflowmodel(bin_dir = os.getcwd(), 
                           hdf5_path=os.getcwd() + "/res/simu.hdf5", 
                           run_type = "min", 
                           clean = True)

model.config.set(config_dictionary_to_force = config_inf)
model.config.get()

model.init_all()

model.kernel.dof  = df2d.wrapping.call_model.dof_copy(dof_init)
model.kernel.dof0 = df2d.wrapping.call_model.dof_copy(dof_init)

model.run()
model.save_all()

model.obs = df2d.core.obs.Obs(bin_dir = model.bin_dir,
                              ts = model.config["ts"])
model.obs.sim_to_true()

