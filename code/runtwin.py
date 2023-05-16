#=======================================================#
# Source librairies
#=======================================================#
import dassflow2d as df2d
import numpy as np
import sys
import os


infer_hydrogram = 1      # sys.argv[1]
infer_bathy = 0          #sys.argv[2]
infer_manning_alpha = 0  #sys.argv[3]
bin_dir= os.getcwd()

print(infer_hydrogram)
print(infer_bathy)
print(infer_manning_alpha)
#=======================================================#
#  Force configuration for direct, and prepare inference configuration
#=======================================================#
config_direct = {"use_obs":0, "w_obs":1, "mesh_name":"mesh.geo"}
config_inf    = {"use_obs":1, "w_obs":1}

#if infer_hydrogram == 1:
#	config_inf["c_hydrograph"]=infer_hydrogram
#if infer_bathy == 1:
#	config_inf["c_bathy"]=infer_bathy
#if infer_manning_alpha == 1:
#	config_inf["c_manning"]=infer_manning_alpha
 

#=======================================================#
#  run direct fortran model
#=======================================================#
os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A")

#if infer_hydrogram == 1:
#	os.system("rm hydrograph.txt")
#	os.system("cp hydrograph_true.txt hydrograph.txt")
#if infer_bathy == 1:
#	os.system(f"rm {model.config['mesh_name']}")
#	os.system(f"cp {model.config['mesh_name'][:-4]}_true.geo {model.config['mesh_name']}")
#if infer_manning_alpha == 1:
#	os.system("rm land_uses.txt")
#	os.system("cp land_uses_true.txt land_uses.txt")

model = df2d.dassflowmodel(bin_dir = os.getcwd(), 
                           hdf5_path=os.getcwd() + "/res/simu.hdf5", 
                           run_type = "direct", 
                           clean = False)

model.config.set(config_dictionary_to_force = config_direct)
model.config.get()
model.init_all()
model.boundary.plot(what ="meshing")
model.meshing.plot_dev(what = "edge")

model.kernel.dof.h[:]=1
dof_init = df2d.wrapping.call_model.dof_copy(model.kernel.dof)
model.kernel.dof0 = df2d.wrapping.call_model.dof_copy(dof_init)
model.run()

model.output = df2d.core.output.Output(bin_dir= os.getcwd(), ts = model.config["ts"],  
                                       boundary_metadata= model.boundary.get_metadata() )
model
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

#if infer_hydrogram == 1:
#	os.system("rm hydrograph.txt")
#	os.system("cp hydrograph_prior.txt hydrograph.txt")
#if infer_bathy == 1:
#	os.system(f"rm {model.config['mesh_name']}")
#	os.system(f"cp {model.config['mesh_name'][:-4]}_prior.geo {model.config['mesh_name']}")
#if infer_manning_alpha == 1:
#	os.system("rm land_uses.txt")
#	os.system("cp land_uses_prior.txt land_uses.txt")

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

