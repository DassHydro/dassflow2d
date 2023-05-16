#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 10:27:19 2023

@author: livillenave
"""

import smash
import numpy as np
import matplotlib.pyplot as plt
import os
import dassflow2d as df2d

# write hydrograph.txt
# qsim smash : smash.Output.qsim object
# write_dir: path where to write hydrograph.txt file
# dt: smash timestep to transpose dates to timerange in seconds from zero
def write_hydrograph(qsim_smash, write_dir, dt ):
    time = np.arange(start = 0, stop = dt * np.shape(qsim_smash)[1], step = dt) 
    with open(f'{write_dir}/hydrograph.txt', 'w') as f:
        f.write("#comment\n")
        f.write("#comment\n")
        f.write("#comment\n")
        f.write(str(np.shape(qsim_smash)[0])+"\n")
        for my_id in range(np.shape(qsim_smash)[0]):
            f.write("#comment\n")
            f.write("#comment\n")
            f.write("#comment\n")
            f.write(f"{np.shape(qsim_smash)[1]}\n")
            for i in range(len(time)):
                f.write(f"{time[i]} {qsim_smash[my_id, i]}\n")
                
        f.write("#comment\n")
                
                
            
def positive_or_negative():
    if np.random.random() < 0.5:
        return 1
    else:
        return -1

def positive_or_negative_ndarrray(shape_ndarray_to_mimic):
    res= np.ndarray(shape =shape_ndarray_to_mimic)
    for i in range(shape_ndarray_to_mimic[0]):
        for j in range(shape_ndarray_to_mimic[1]):
            a = positive_or_negative()
            res[i,j]= a
    return(res)
                
  

def ifelse(condition, return_true, return_false):
    if condition:
        return(return_true)
    else:
        return(return_false)              
                
# =============================== #
# SMASH PARAMETER
# =============================== #

#------------------------#
# config                 #
#------------------------#
nb_timestep=24,                          # number of timestep
setup = {
        "dt": 900,
        "start_time": "2020-01-01 00:00",
        "end_time": "2020-01-01 06:00",
    }
# mesh definition
dx = 1_000
(nrow, ncol) = (10, 10)
mesh = {
"dx": dx,
"nrow": nrow,
"ncol": ncol,
"ng": 1,
"nac": nrow * ncol,
"area": nrow * ncol * (dx ** 2),
"gauge_pos": np.array([9, 9], dtype=np.int32),
"code": np.array(["Practice_case"])}

#------------------------#
#  FLOW DIRECTION AND DRAINED AREA    #
#------------------------#
mesh["flwdir"] = np.array(
    [    
    [4, 5, 5, 5, 5, 5, 5, 5, 5, 5],    
    [3, 4, 5, 5, 5, 5, 5, 5, 5, 5],    
    [3, 3, 4, 5, 5, 5, 5, 5, 5, 5],    
    [3, 3, 3, 4, 5, 5, 5, 5, 5, 5],    
    [3, 3, 3, 3, 4, 5, 5, 5, 5, 5],    
    [3, 3, 3, 3, 3, 4, 5, 5, 5, 5],
    [3, 3, 3, 3, 3, 3, 4, 5, 5, 5],    
    [3, 3, 3, 3, 3, 3, 3, 4, 5, 5],    
    [3, 3, 3, 3, 3, 3, 3, 3, 4, 5],    
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 4],    
    ], dtype=np.int32)

mesh["flwacc"] = np.array(
    [      [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],    
           [1, 4, 2, 2, 2, 2, 2, 2, 2, 2],    
           [1, 2, 9, 3, 3, 3, 3, 3, 3, 3],    
           [1, 2, 3, 16, 4, 4, 4, 4, 4, 4],    
           [1, 2, 3, 4, 25, 5, 5, 5, 5, 5],    
           [1, 2, 3, 4, 5, 36, 6, 6, 6, 6],    
           [1, 2, 3, 4, 5, 6, 49, 7, 7, 7],    
           [1, 2, 3, 4, 5, 6, 7, 64, 8, 8],    
           [1, 2, 3, 4, 5, 6, 7, 8, 81, 9],    
           [1, 2, 3, 4, 5, 6, 7, 8, 9, 100]], dtype=np.int32 )


plt.imshow(mesh["flwacc"])
plt.show();plt.close()
plt.imshow(mesh["flwdir"])
plt.show();plt.close()
#------------------------#
#  path    #
#------------------------#
ind_path = np.unravel_index(np.argsort(mesh["flwacc"], axis=None),    
     mesh["flwacc"].shape)
mesh["path"] = np.zeros(shape=(2, mesh["flwacc"].size),    
    dtype=np.int32)
mesh["path"][0, :] = ind_path[0]    
mesh["path"][1, :] = ind_path[1]  

#------------------------#
# rainfall generation    #
#------------------------#
prcp = np.zeros(shape=nb_timestep, dtype=np.float32)
tri = np.linspace(0, 12, 10)
prcp[0:10] = tri
prcp[9:19] = np.flip(tri)

plt.plot(prcp);plt.show();plt.close()
# --------------------------------------- #
# Potential evapotranspiration generation #
# --------------------------------------- #
pet=0

#----------------------------------------#
#       Spatialised parameters
#       (gr-a structure)
#----------------------------------------#

# --- INITIALISE ---#
parameters = dict()
# The maximum capacity of the transfer storage.    [mm]
parameters["cft"] = np.zeros(shape = mesh["flwacc"].shape)
# The maximum capacity of the production storage. [mm]
parameters["cp"] = np.zeros(shape = mesh["flwacc"].shape)
# The non-conservative exchange parameter.        [mm/dt]
parameters["exc"] = np.zeros(shape = mesh["flwacc"].shape)
# The linear routing parameter.                   [minuts]
parameters["lr"] = np.zeros(shape = mesh["flwacc"].shape)


#--- SET VALUES ---#
# >>> SET CFT
parameters["cft"][:,:] = 10
parameters["cft"][(0,1,2,3,4,5,6,7,8,9),(0,1,2,3,4,5,6,7,8,9)] = 1000
parameters["cft"][(6,7,8,9),(6,7,8,9)] = 0.001

# >>> SET CP
parameters["cp"][:,:] = 200
parameters["cp"][(0,1,2,3,4,5,6),(0,1,2,3,4,5,6)] = 300
parameters["cp"][(7,8,9),(7,8,9)] = 250

# set exc to zero (default value)
parameters["exc"][:,:] = 0

# set linear routing to default value
parameters["lr"][:,:] = 5

plt.imshow(parameters["cft"]);plt.show();plt.close()
plt.imshow(parameters["cp"]);plt.show();plt.close()
plt.imshow(parameters["exc"]);plt.show();plt.close()
plt.imshow(parameters["lr"]);plt.show();plt.close()



def copy_dict(init_dict):
    new_dict = dict()
    for key, values in init_dict.items():
        new_dict[key] = values.copy()
    return(new_dict)


copy_dict(init_dict = parameters)

#----------------------------------------#
# gather created data into a unique dictionary
#----------------------------------------#
data_smash_true = dict()
data_smash_true["model1"] = {"setup":setup.copy(), "mesh":mesh.copy(), "rain":prcp.copy(), "pet":pet, "param":copy_dict(init_dict = parameters)} 
data_smash_true["model2"] = {"setup":setup.copy(), "mesh":mesh.copy(), "rain":prcp.copy(), "pet":pet, "param":copy_dict(init_dict = parameters)}
data_smash_true["model2"]["param"]["cp"] = data_smash_true["model1"]["param"]["cp"]- 50 

# --- param --- #
prior_cp_noise = 100
        
data_smash_prior = dict()
data_smash_prior["model1"] = {"setup":setup.copy(), "mesh":mesh.copy(), "rain":prcp.copy(), "pet":pet, "param":copy_dict(init_dict = parameters)} 
data_smash_prior["model2"] = {"setup":setup.copy(), "mesh":mesh.copy(), "rain":prcp.copy(), "pet":pet, "param":copy_dict(init_dict = parameters)}
data_smash_prior["model1"]["param"]["cp"] = data_smash_true["model1"]["param"]["cp"] + (np.random.rand(data_smash_true["model1"]["param"]["cp"].shape[0], 
                data_smash_true["model1"]["param"]["cp"].shape[1])  * prior_cp_noise  * positive_or_negative_ndarrray(data_smash_true["model1"]["param"]["cp"].shape ) )
data_smash_prior["model2"]["param"]["cp"] = data_smash_true["model2"]["param"]["cp"] + (np.random.rand(data_smash_true["model1"]["param"]["cp"].shape[0], 
                data_smash_true["model2"]["param"]["cp"].shape[1])  * prior_cp_noise * positive_or_negative_ndarrray(data_smash_true["model1"]["param"]["cp"].shape) )

#----------------------------------------#
# run true smash
#----------------------------------------#

# CATCHMENT 1
smash_model1 = smash.Model(data_smash_true["model1"]["setup"], data_smash_true["model1"]["mesh"])
smash_model1.input_data.prcp = np.broadcast_to(data_smash_true["model1"]["rain"], smash_model1.input_data.prcp.shape)
smash_model1.input_data.pet = data_smash_true["model1"]["pet"]
smash_model1.parameters.cft = data_smash_true["model1"]["param"]["cft"]
smash_model1.parameters.cp = data_smash_true["model1"]["param"]["cp"]
smash_model1.parameters.exc = data_smash_true["model1"]["param"]["exc"]
smash_model1.parameters.lr = data_smash_true["model1"]["param"]["lr"]
smash_model1.run(inplace=True)
# add to dictionary true qsim
data_smash_true["model1"]["qsim"] = smash_model1.output.qsim.copy()

shape_qsim_smash_model1 = smash_model1.output.qsim.shape
del smash_model1


# CATCHMENT 2
smash_model2 = smash.Model(data_smash_true["model2"]["setup"], data_smash_true["model2"]["mesh"])
smash_model2.input_data.prcp = np.broadcast_to(data_smash_true["model2"]["rain"], 
                                              smash_model2.input_data.prcp.shape)
smash_model2.input_data.pet = data_smash_true["model2"]["pet"]
smash_model2.parameters.cft = data_smash_true["model2"]["param"]["cft"]
smash_model2.parameters.cp = data_smash_true["model2"]["param"]["cp"]
smash_model2.parameters.exc = data_smash_true["model2"]["param"]["exc"]
smash_model2.parameters.lr = data_smash_true["model2"]["param"]["lr"]
smash_model2.run(inplace=True)
# add to dictionary true qsim
data_smash_true["model2"]["qsim"] = smash_model2.output.qsim.copy()
del smash_model2


# ------------------------------------- #
# run true Dassflow based on generated hydrograms before
# ------------------------------------- #

# >>> define paths
dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"
bin_dir = f"{dassflow_dir}/code/bin_A/"
#ref_case_dir = 0

#>>> import 1D2D case from tuto case
os.chdir(f"{dassflow_dir}/code/")
os.system(f"rm -r {bin_dir}/*")    
os.system(f"cp -r {dassflow_dir}/cases/tuto_case/6_1D2D_less_stations/* {bin_dir}/")  

# >>> traduce smash information
config_direct_hy = dict() 
config_direct_hy["ts"] = (data_smash_true["model1"]["qsim"].shape[1]-1)* (data_smash_true["model1"]["setup"]["dt"])
config_direct_hy["dtw"] = config_direct_hy["dtp"]  = config_direct_hy["ts"]/100
config_direct_hy["w_vtk"] = 0


# write the two modeled hydrograph
qsim_smash = np.ndarray(shape = (2, data_smash_true["model1"]["qsim"].shape[1]))
qsim_smash[0,:] = data_smash_true["model1"]["qsim"] + 10
qsim_smash[1,:] = data_smash_true["model2"]["qsim"] + 10
write_hydrograph(qsim_smash=qsim_smash, write_dir = bin_dir, dt = data_smash_true["model1"]["setup"]["dt"])

# >>> Run model
os.chdir(f"{bin_dir}")
dassflow_model_true = df2d.dassflowmodel(bin_dir =  f"{bin_dir}", 
                                    hdf5_path= f"{bin_dir}/res/simu.hdf5" , 
                                    run_type = "direct",
                                    clean = True) # TO DEBUG: DO NOT WORK WITH CLEAN + TRUE

for config_key in config_direct_hy.keys():
    if config_key in dassflow_model_true.config.keys():
        dassflow_model_true.config[config_key]=config_direct_hy[config_key]
dassflow_model_true.config.set()
dassflow_model_true.init_all()

dassflow_model_true.kernel.dof0.h[:] = dassflow_model_true.kernel.dof.h[:] = 1

h_init = dassflow_model_true.kernel.dof.h.copy()
u_init = dassflow_model_true.kernel.dof.u.copy()
v_init = dassflow_model_true.kernel.dof.v.copy()


dassflow_model_true.run()


# Plot simulation
#time_keys = dassflow_model_true.outputs.all_times
#for id_time in time_keys:
#    my_scalar = np.asanyarray(dassflow_model_true.outputs.all_res[id_time]["h"])
#    dassflow_model_true.meshing.mesh_pyvista.plot(cpos = "xy", show_edges = False, scalars = my_scalar)

# done 
import os
os.chdir(bin_dir)
os.system(f"cp -r ./res/obs/* ./obs/")

df2d.wrapping.call_model.clean_model(dassflow_model_true.kernel)
#============================================================#
# RUN SMASH PRIOR
#============================================================#


# CATCHMENT 1
smash_model1 = smash.Model(data_smash_prior["model1"]["setup"], data_smash_prior["model1"]["mesh"])
smash_model1.input_data.prcp = np.broadcast_to(data_smash_prior["model1"]["rain"], smash_model1.input_data.prcp.shape)
smash_model1.input_data.pet = data_smash_prior["model1"]["pet"]
smash_model1.parameters.cft = data_smash_prior["model1"]["param"]["cft"]
smash_model1.parameters.cp = data_smash_prior["model1"]["param"]["cp"]
smash_model1.parameters.exc = data_smash_prior["model1"]["param"]["exc"]
smash_model1.parameters.lr = data_smash_prior["model1"]["param"]["lr"]
smash_model1.run(inplace=True)
# add to dictionary true qsim
data_smash_prior["model1"]["qsim"] = smash_model1.output.qsim.copy()
del smash_model1


# CATCHMENT 2
smash_model2 = smash.Model(data_smash_prior["model2"]["setup"], 
                           data_smash_prior["model2"]["mesh"])
smash_model2.input_data.prcp = np.broadcast_to(data_smash_prior["model2"]["rain"], 
                                              smash_model2.input_data.prcp.shape)
smash_model2.input_data.pet = data_smash_prior["model2"]["pet"]
smash_model2.parameters.cft = data_smash_prior["model2"]["param"]["cft"]
smash_model2.parameters.cp = data_smash_prior["model2"]["param"]["cp"]
smash_model2.parameters.exc = data_smash_prior["model2"]["param"]["exc"]
smash_model2.parameters.lr = data_smash_prior["model2"]["param"]["lr"]
smash_model2.run(inplace=True)
# add to dictionary true qsim
data_smash_prior["model2"]["qsim"] = smash_model2.output.qsim.copy()
del smash_model2


        
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('True and Prior simulation, SMASH model')
ax1.set_title("TRUE")
ax1.plot( data_smash_true["model1"]["qsim"][0], c="red", label = "catchemnt 1")
ax1.plot( data_smash_true["model2"]["qsim"][0], c="blue", label = "catchemnt 2")
ax1.set(xlabel='id timestep', ylabel='Q (m3/s)')
ax2.set_title("PRIOR")
ax2.plot( data_smash_prior["model1"]["qsim"][0], c="red", label = "catchemnt 1")
ax2.plot( data_smash_prior["model2"]["qsim"][0], c="blue", label = "catchemnt 2"  )
ax2.set(xlabel='id timestep')
ax2.legend()
# Set yaxis for both plot
my_ylim1 = ax1.get_ylim() ; my_ydist1=my_ylim1[1]-my_ylim1[0]
my_ylim2 = ax2.get_ylim() ; my_ydist2=my_ylim2[1]-my_ylim2[0]
axref = ifelse(condition = my_ydist1>my_ydist2, return_true= ax1, return_false = ax2)
plt.setp(ax1, ylim=axref.get_ylim())
plt.setp(ax2, ylim=axref.get_ylim())


#===============================#
# Prior run dassflow
#==============================#

# write the two modeled hydrograph
qsim_smash = np.ndarray(shape = (2, data_smash_prior["model1"]["qsim"].shape[1]))
qsim_smash[0,:] = data_smash_prior["model1"]["qsim"] + 10
qsim_smash[1,:] = data_smash_prior["model2"]["qsim"] + 10
write_hydrograph(qsim_smash=qsim_smash, write_dir = bin_dir, dt = data_smash_prior["model1"]["setup"]["dt"])

# >>> Run model
os.chdir(f"{bin_dir}")
dassflow_model_prior = df2d.dassflowmodel(bin_dir =  f"{bin_dir}", 
                                    hdf5_path= f"{bin_dir}/res/simu.hdf5" , 
                                    run_type = "direct",
                                    clean = True) # TO DEBUG: DO NOT WORK WITH CLEAN + TRUE

for config_key in config_direct_hy.keys():
    if config_key in dassflow_model_prior.config.keys():
        dassflow_model_prior.config[config_key]=config_direct_hy[config_key]
dassflow_model_prior.config.set()
dassflow_model_prior.init_all()


dassflow_model_prior.kernel.dof.h[:] = dassflow_model_prior.kernel.dof0.h[:] = h_init[:]
dassflow_model_prior.kernel.dof.u[:] = dassflow_model_prior.kernel.dof0.u[:] = u_init[:]
dassflow_model_prior.kernel.dof.v[:] = dassflow_model_prior.kernel.dof0.v[:] = v_init[:]

dassflow_model_prior.run()


# Plot simulation
#time_keys = dassflow_model_prior.outputs.all_times
#for id_time in time_keys:
#    my_scalar = np.asanyarray(dassflow_model_prior.outputs.all_res[id_time]["h"])
#    dassflow_model_prior.meshing.mesh_pyvista.plot(cpos = "xy", show_edges = False, scalars = my_scalar)


df2d.wrapping.call_model.clean_model(dassflow_model_prior.kernel)
#
#time_keys_prior = dassflow_model_prior.outputs.all_times
#time_keys_true = dassflow_model_true.outputs.all_times
#
#all_times_res = np.ndarray(shape =  (dassflow_model_prior.meshing.mesh_fortran.nc, len(time_keys_prior)) )
#for id_time in range(len(time_keys_prior)):
#    my_scalar = np.asanyarray(dassflow_model_prior.outputs.all_res[time_keys_prior[id_time]]["h"]) -     np.asanyarray(dassflow_model_true.outputs.all_res[time_keys_true[id_time]]["h"])
#    all_times_res[:, id_time] = my_scalar
#   # dassflow_model_prior.meshing.mesh_pyvista.plot(cpos = "xy", show_edges = False, scalars = my_scalar)
#
#tmp = np.linalg.norm(all_times_res, ord=2, axis=None, keepdims=False)
#tmp = np.linalg.norm(all_times_res, ord=2, axis=1, keepdims=False)


df2d.wrapping.call_model.clean_model(dassflow_model_true.kernel)
#==============================#
# INFERENCE run dassflow
#==============================#
config_inf_hy = config_direct_hy
config_inf_hy["use_obs"]=1
config_inf_hy["use_c_hydrograph"]=1


# write the two modeled hydrograph
qsim_smash = np.ndarray(shape = (2, data_smash_prior["model1"]["qsim"].shape[1]))
qsim_smash[0,:] = data_smash_prior["model1"]["qsim"] + 10
qsim_smash[1,:] = data_smash_prior["model2"]["qsim"] + 10 

write_hydrograph(qsim_smash=qsim_smash, write_dir = bin_dir, dt = data_smash_prior["model1"]["setup"]["dt"])

# >>> Run model
os.chdir(f"{bin_dir}")
dassflow_model_prior = df2d.dassflowmodel(bin_dir =  f"{bin_dir}", 
                                    hdf5_path= f"{bin_dir}/res/simu.hdf5" , 
                                    run_type = "min",
                                    clean = True) # TO DEBUG: DO NOT WORK WITH CLEAN + TRUE

for config_key in config_direct_hy.keys():
    if config_key in dassflow_model_prior.config.keys():
        dassflow_model_prior.config[config_key]=config_direct_hy[config_key]
dassflow_model_prior.config.set()
dassflow_model_prior.init_all()

dassflow_model_prior.kernel.dof.h[:] = dassflow_model_prior.kernel.dof0.h[:] = h_init[:]
dassflow_model_prior.kernel.dof.u[:] = dassflow_model_prior.kernel.dof0.u[:] = u_init[:]
dassflow_model_prior.kernel.dof.v[:] = dassflow_model_prior.kernel.dof0.v[:] = v_init[:]

dassflow_model_prior.run()

dassflow_model_prior.min = df2d.core.min.Min(bin_dir = dassflow_model_prior.bin_dir) # to updatage package (done but to push)
dassflow_model_prior.min.source_all()

plt.close()
dassflow_model_prior.min.plot()


id_ite_max=int(max(dassflow_model_prior.min.all_ite))

qinf1 = dassflow_model_prior.min.param["hydrograph"][id_ite_max][1][:,1] -10
qinf2 = dassflow_model_prior.min.param["hydrograph"][id_ite_max][2][:,1] -10

        
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.suptitle('True and Prior simulation, SMASH model')
ax1.set_title("TRUE")
ax1.plot( data_smash_true["model1"]["qsim"][0], c="red", label = "catchemnt 1")
ax1.plot( data_smash_true["model2"]["qsim"][0], c="blue", label = "catchemnt 2")
ax1.set(xlabel='id timestep', ylabel='Q (m3/s)')

ax2.set_title("PRIOR")
ax2.plot( data_smash_prior["model1"]["qsim"][0], c="red", label = "catchemnt 1")
ax2.plot( data_smash_prior["model2"]["qsim"][0], c="blue", label = "catchemnt 2"  )
ax2.set(xlabel='id timestep')
ax2.legend()

# DASSFLOW INFERENCE
ax3.set_title("INFERED DASSFLOW")
ax3.plot( qinf1 , c="red", label = "catchemnt 1")
ax3.plot( qinf2, c="blue", label = "catchemnt 2"  )
ax3.set(xlabel='id timestep')

# Set yaxis for both plot
my_ylim1 = ax1.get_ylim() ; my_ydist1=my_ylim1[1]-my_ylim1[0]
my_ylim2 = ax2.get_ylim() ; my_ydist2=my_ylim2[1]-my_ylim2[0]
my_ylim3 = ax3.get_ylim() ; my_ydist3=my_ylim3[1]-my_ylim3[0]
axref = ifelse(condition = my_ydist1>my_ydist2, return_true= ax1, return_false = ax2)
plt.setp(ax1, ylim=axref.get_ylim())
plt.setp(ax2, ylim=axref.get_ylim())
plt.setp(ax3, ylim=axref.get_ylim())




        
fig, (ax1, ax2 ) = plt.subplots(1, 2)
fig.set_figwidth = 20
fig.set_figheight = 20
fig.suptitle('True and Prior simulation, SMASH model')
ax1.set_title("Catchment 1")
ax1.plot( data_smash_true["model1"]["qsim"][0], c="red", marker ='.', label = "True")
ax1.plot(  data_smash_prior["model1"]["qsim"][0], c="green", label = "Prior")
ax1.plot(  qinf1, c="blue", label = "Inf")
ax1.set(xlabel='id timestep', ylabel='Q (m3/s)')

ax2.set_title("Catchment 2")
ax2.plot(  data_smash_true["model2"]["qsim"][0], marker ='.',  c="red", label = "True")
ax2.plot( data_smash_prior["model2"]["qsim"][0], c="green", label = "Prior"  )
ax2.plot( qinf2, c="blue", label = "Inf"  )
ax2.set(xlabel='id timestep')
ax2.legend()

# Set yaxis for both plot
my_ylim1 = ax1.get_ylim() ; my_ydist1=my_ylim1[1]-my_ylim1[0]
my_ylim2 = ax2.get_ylim() ; my_ydist2=my_ylim2[1]-my_ylim2[0]
axref = ifelse(condition = my_ydist1>my_ydist2, return_true= ax1, return_false = ax2)
plt.setp(ax1, ylim=axref.get_ylim())
plt.setp(ax2, ylim=axref.get_ylim())
plt.show()
plt.close()


plt.plot( data_smash_true["model1"]["qsim"][0] + data_smash_true["model2"]["qsim"][0], c="blue", label = "TRUE")
plt.plot( data_smash_prior["model1"]["qsim"][0] +  data_smash_prior["model2"]["qsim"][0], c="green", label = "PRIOR")
plt.plot( qinf1+qinf2, c="red", label = "INFERED DASSFLOW")
plt.legend()
plt.show()
plt.close()


#==============================#
# INFERENCE SMASH using observation created by dassflow
#==============================#

import copy
data_smash_inf = copy.deepcopy(data_smash_prior)
epsmin=1*10**-3
qinf1[qinf1<epsmin] = epsmin
qinf2[qinf2<epsmin] = epsmin


qobs1 = np.ndarray(shape = shape_qsim_smash_model1)
qobs2 = np.ndarray(shape = shape_qsim_smash_model1)

qobs1[:] = qinf1.copy() 
qobs2[:] = qinf2.copy()
#----------------
# >>> Run
#--------------

# CATCHMENT 1
smash_model1 = smash.Model(data_smash_inf["model1"]["setup"], data_smash_inf["model1"]["mesh"])

smash_model1.input_data.qobs = qobs1.copy()
smash_model1.input_data.prcp = np.broadcast_to(data_smash_inf["model1"]["rain"], smash_model1.input_data.prcp.shape)
smash_model1.input_data.pet = data_smash_inf["model1"]["pet"]
smash_model1.parameters.cft = data_smash_inf["model1"]["param"]["cft"].copy()
smash_model1.parameters.cp = data_smash_inf["model1"]["param"]["cp"].copy()
smash_model1.parameters.exc = data_smash_inf["model1"]["param"]["exc"].copy()
smash_model1.parameters.lr = data_smash_inf["model1"]["param"]["lr"].copy()

smash_model1.optimize(control_vector="cp", inplace=True);
smash_model1.optimize(control_vector="cp",mapping= "distributed", inplace=True);
smash_model1.run(inplace=True)
# add to dictionary true qsim
data_smash_inf["model1"]["qsim"] = smash_model1.output.qsim.copy()
data_smash_inf["model1"]["qobs"] = smash_model1.input_data.qobs.copy()
data_smash_inf["model1"]["param"]["cp"] = smash_model1.parameters.cp.copy()
del smash_model1


# CATCHMENT 2
smash_model2 = smash.Model(data_smash_inf["model2"]["setup"], 
                           data_smash_inf["model2"]["mesh"])
smash_model2.input_data.qobs = qobs2.copy()
smash_model2.input_data.prcp = np.broadcast_to(data_smash_inf["model2"]["rain"].copy(), 
                                              smash_model2.input_data.prcp.shape)
smash_model2.input_data.pet = data_smash_inf["model2"]["pet"]
smash_model2.parameters.cft = data_smash_inf["model2"]["param"]["cft"].copy()
smash_model2.parameters.cp = data_smash_inf["model2"]["param"]["cp"].copy()
smash_model2.parameters.exc = data_smash_inf["model2"]["param"]["exc"].copy()
smash_model2.parameters.lr = data_smash_inf["model2"]["param"]["lr"].copy()
smash_model2.optimize(control_vector="cp", inplace=True);
smash_model2.optimize(control_vector="cp",mapping= "distributed", inplace=True);
smash_model2.run(inplace=True)
# add to dictionary true qsim
data_smash_inf["model2"]["qsim"] = smash_model2.output.qsim.copy()
data_smash_inf["model2"]["qobs"] = smash_model2.input_data.qobs.copy()
data_smash_inf["model2"]["param"]["cp"] = smash_model2.parameters.cp.copy()
del smash_model2




        
fig, (ax1, ax2 ) = plt.subplots(1, 2)
fig.set_figwidth = 20
fig.set_figheight = 20
fig.suptitle('True and Prior simulation, SMASH model')
ax1.set_title("Catchment 1")
ax1.plot( data_smash_true["model1"]["qsim"][0], c="red", marker ='.', label = "True")
ax1.plot(  data_smash_prior["model1"]["qsim"][0], c="green", label = "Prior")
ax1.plot(  data_smash_inf["model1"]["qobs"][0], c="brown", marker ='x', label = "Observation DASSFLOW produced")
ax1.plot(  data_smash_inf["model1"]["qsim"][0], c="blue", label = "Infered smash")
ax1.set(xlabel='id timestep', ylabel='Q (m3/s)')

ax2.set_title("Catchment 2")
ax2.plot(  data_smash_true["model2"]["qsim"][0], marker ='.',  c="red", label = "True")
ax2.plot( data_smash_prior["model2"]["qsim"][0], c="green", label = "Prior"  )
ax2.plot(  data_smash_inf["model2"]["qobs"][0], marker ='x', c="brown", label = "Observation DASSFLOW produced")
ax2.plot(   data_smash_inf["model2"]["qsim"][0], c="blue", label = "Infered smash"  )
ax2.set(xlabel='id timestep')
ax2.legend()

# Set yaxis for both plot
my_ylim1 = ax1.get_ylim() ; my_ydist1=my_ylim1[1]-my_ylim1[0]
my_ylim2 = ax2.get_ylim() ; my_ydist2=my_ylim2[1]-my_ylim2[0]
axref = ifelse(condition = my_ydist1>my_ydist2, return_true= ax1, return_false = ax2)
plt.setp(ax1, ylim=axref.get_ylim())
plt.setp(ax2, ylim=axref.get_ylim())
plt.show()
plt.close()









        
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.suptitle('True and Prior simulation, SMASH model')
ax1.set_title("TRUE")
pos1 = ax1.imshow(data_smash_true["model1"]["param"]["cp"])
fig.colorbar(pos1, ax=ax1)

ax2.set_title("PRIOR")
pos2 = ax2.imshow(data_smash_prior["model1"]["param"]["cp"])
fig.colorbar(pos2, ax=ax2)

ax3.set_title("INFERED DASSFLOW")
pos3 = ax3.imshow(data_smash_inf["model1"]["param"]["cp"])
fig.colorbar(pos3, ax=ax3)



        
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.suptitle('True and Prior simulation, SMASH model')
ax1.set_title("TRUE")
ax1.imshow(data_smash_true["model2"]["param"]["cp"])
ax2.set_title("PRIOR")
ax2.imshow(data_smash_prior["model2"]["param"]["cp"])
# DASSFLOW INFERENCE
ax3.set_title("INFERED DASSFLOW")
ax3.imshow(data_smash_inf["model2"]["param"]["cp"])
ax3.set(xlabel='id timestep')

# Set yaxis for both plot
my_ylim1 = ax1.get_ylim() ; my_ydist1=my_ylim1[1]-my_ylim1[0]
my_ylim2 = ax2.get_ylim() ; my_ydist2=my_ylim2[1]-my_ylim2[0]
my_ylim3 = ax3.get_ylim() ; my_ydist3=my_ylim3[1]-my_ylim3[0]
axref = ifelse(condition = my_ydist1>my_ydist2, return_true= ax1, return_false = ax2)
plt.setp(ax1, ylim=axref.get_ylim())
plt.setp(ax2, ylim=axref.get_ylim())
plt.setp(ax3, ylim=axref.get_ylim())



