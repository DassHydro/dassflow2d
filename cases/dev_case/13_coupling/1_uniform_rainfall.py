#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 10:20:46 2023

@author: livillenave
"""

def print_xybc(model):
    corresp =  model.boundary.corresp
    
    for elem in corresp["discharg1"].values():
        id_cell = elem["id_cell"]
        print(id_cell)
        
        for my_id in id_cell:
            gravity = model.meshing.mesh_fortran.cell[my_id-1].grav
            print(gravity.x, gravity.y)
            

        


import matplotlib as mpl

mpl.rcParams['figure.figsize'] = [12.0, 12.0]
mpl.rcParams['figure.dpi'] = 80
mpl.rcParams['savefig.dpi'] = 100

mpl.rcParams['font.size'] = 16
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'large'
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['lines.linewidth'] = 2.5
mpl.rcParams['lines.markersize'] = 10

#================================================#
#  SCRIPT PATH                             #
#================================================#
dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"
case_dir = dassflow_dir +"/cases/dev_case/13_coupling"
bin_dir = dassflow_dir + "/bin_A"
#================================================#
#  params                          #
#================================================#
nb_timestep_smash = 48
dt_smash = 3600

dx_dassflow = 10 # in m
ly = 10  # width of channel in y direction
ny = 2   # number of nodes in y direction
         # smash resolution is hardcoded 
#================================================#
#  source librairies                             #
#================================================#

import numpy as np
import matplotlib.pyplot as plt
import dassflow2d as df2d
import os
import copy

# import developed librairies
os.chdir(dassflow_dir + "/Tools/1_pre-treatment/2_shapefile_to_gmsh/organised_version/libs")
#from GIS_libs import *
#from gmsh_libs import *
from smash_coupling_libs import *
#from dassflow_mesh_manipulation_libs import *


runfile(f'{case_dir}/libs/coupling_1.py')
runfile(f'{case_dir}/libs/libs.py')

        
        
#================================================#
#  Generate SMASH CASE                           #
#================================================#
data = build_case(nb_timestep=nb_timestep_smash,               
                  start_time="2020-01-01 00:00", 
                  end_time="2020-01-03 00:00")

data["mesh"] = data["mesh"].copy()
data["mesh"]["flwacc"] = data["mesh"]["drained_area"]
data["mesh"], inflows = get_hydrological_coupling(mesh = data["mesh"],
                                          treshold_drained_area=36)
plt.imshow(data["mesh"]["categorical"])
data["mesh"] =  add_xy_coord(data["mesh"])

plt.imshow(data["mesh"]["x_coord"])
plt.colorbar()
plt.title(" x coordinate cell")
plt.show()
plt.imshow(data["mesh"]["y_coord"])
plt.colorbar()
plt.title(" y coordinate cell")
plt.show()


#=============================================================================#
# perform SMASH simulation
#  a) run model_1000_3bc
#  b) extract inflows
#  c) extract net rain on river network
#  d) extract  streamflow at the gauge (outlet of the catchment)
#=============================================================================#
data["rain"] = data["rain"] + 25
data["rain"][:] = 30
#---------------------------------------------------------------------------------------------------#
#        a) run model
#---------------------------------------------------------------------------------------------------#

# Dassflow has been built manually
data["setup"]["save_net_prcp_domain"] = True
data["setup"]["save_qsim_domain"] = True
smash_model = smash.Model(data["setup"], data["mesh"])

smash_model.input_data.prcp = np.broadcast_to(data["rain"], smash_model.input_data.prcp.shape)
smash_model.input_data.pet = data["pet"]
smash_model.parameters.lr=50
smash_model.parameters.lr[9,9] = smash_model.parameters.lr[8,8] = smash_model.parameters.lr[7,7]= smash_model.parameters.lr[6,6] = 1
smash_model.parameters.cft=100
smash_model.parameters.cp=200

smash_model.run(inplace=True)
       

# ----------------------------- #
# store results
# ----------------------------- #
for id_inflow in inflows.keys():
    my_inflow = inflows[id_inflow]
    my_inflow["qsim"] = smash_model.output.qsim_domain[my_inflow["id"][0],my_inflow["id"][1],:]
    
    
plt.imshow(data["mesh"]["categorical"])
for id_inflow in inflows.keys():
    my_inflow = inflows[id_inflow]
    plt.text(x = my_inflow["id"][0],
             y = my_inflow["id"][1],
             s = id_inflow)
plt.show()
plt.savefig()
plt.close()


sum_inflows=np.zeros(shape = my_inflow["qsim"].shape)
for id_inflow in inflows.keys():
    my_inflow = inflows[id_inflow]
    sum_inflows=sum_inflows+ my_inflow["qsim"]
    if id_inflow % 2 == 0:
        plt.plot(np.arange(len(my_inflow["qsim"])),my_inflow["qsim"], label = f"Inflow: {id_inflow}")
    else:
        plt.plot(np.arange(len(my_inflow["qsim"])),
                 my_inflow["qsim"], 
                 label = f"Inflow: {id_inflow}",
                 marker = "x")
        
        
    plt.legend()
plt.show()
plt.close()

### Compare inflows to outflows
my_shape = smash_model.output.qsim_domain.shape
outflow = smash_model.output.qsim_domain[my_shape[0]-1,my_shape[0]-1] 
    
plt.plot(np.arange(len(my_inflow["qsim"])),sum_inflows, 
         label = f"sum inflows ")
plt.plot(np.arange(len(my_inflow["qsim"])),outflow, 
         label = f"outflow")
plt.legend()
plt.show()
plt.close()



inflow_cumsum_minus1 = np.zeros(shape = my_inflow["qsim"].shape)
inflow_cumsum = np.zeros(shape = my_inflow["qsim"].shape)

for id_inflow in inflows.keys():
        my_inflow = inflows[id_inflow]
        print( np.mean(my_inflow["qsim"]))
        inflow_cumsum=inflow_cumsum+ my_inflow["qsim"]
        plt.fill_between(x =np.arange(len(my_inflow["qsim"])), 
                         y1 = inflow_cumsum, y2 =inflow_cumsum_minus1)
        plt.text(x = 40, 
                 y = 0.99*(inflow_cumsum[-1] + inflow_cumsum_minus1[-1])/2, 
                 s = f"Inflow: {id_inflow}"  ) 
        plt.plot(np.arange(len(my_inflow["qsim"])), 
                 inflow_cumsum, 
                 color = "black", 
                 linewidth = 3, linestyle = "dotted")
        inflow_cumsum_minus1 = inflow_cumsum.copy()
        
plt.fill_between(x =np.arange(len(my_inflow["qsim"])), y1 = outflow, y2 = sum_inflows)
plt.text(x = 33, 
         y =0.99* (outflow[30] + sum_inflows[30])/2, 
         s = r"net rainfall on $ \Omega_{rr \backslash up}$"  ) 

plt.plot(np.arange(len(my_inflow["qsim"])), sum_inflows, color = "red", linewidth = 3, label = "sum inflows") 
plt.plot(np.arange(len(my_inflow["qsim"])), outflow, color = "blue", linewidth = 3, label = "Outflow")

plt.plot(np.arange(len(my_inflow["qsim"])), data["rain"] * 10 **-3 * 10**8 /dt_smash
         , color = "black", linewidth = 3, label = "raw rainfall", linestyle = "dashed")

plt.xlabel("time (Hour)")
plt.ylabel("Water Volume fluxes $m^3/s$")
plt.legend(loc=8,
          ncol=3, fancybox=True, shadow=True); plt.show();plt.close()



id_river_rr = np.where(data["mesh"]["categorical"]==1)
river_cells = dict()
sum_rain = np.zeros(my_shape[2])
for river_cell in range(len(id_river_rr[0])):
    res = dict()
    res["id"] = (id_river_rr[0][river_cell], id_river_rr[1][river_cell])
    res["rain"] =smash_model.output.net_prcp_domain[res["id"][0], res["id"][1]]
    res["rain_m3s"]= res["rain"] /3.6
    sum_rain = sum_rain + res["rain_m3s"]
    river_cells[river_cell]  = res
    
plt.plot(np.arange(len(sum_inflows)),sum_inflows, color= 'blue' ,
         label = f"sum inflows ")
plt.plot(np.arange(len(sum_rain)),sum_inflows + sum_rain, 
         label = f"sum inflows+rain ", color = "red", marker="1")
plt.plot(np.arange(len(my_inflow["qsim"])),outflow, 
         label = f"outflow", marker = "2", color ="purple")
plt.legend()
plt.show()
plt.close()




##########
# GEOM CASE 1: 
#    1 inflow
  
# GEOM CASE 2: 
#    2 inflow
  
# GEOM CASE 2: 
#    2 inflow
    

def write_ascii_raster_smash(smash_model, pathfile,my_rast=None,):
    
    if my_rast is None:
        rast_to_plot = smash_model.mesh.flwacc[:,:]
    else:
        rast_to_plot = my_rast
    
    xmin=smash_model.mesh.xmin
    ymax=smash_model.mesh.ymax
    nrow=smash_model.mesh.nrow
    ncol=smash_model.mesh.ncol
    cellsize=smash_model.mesh.dx
    nodata_value=-99
    
    f = open(f"{pathfile}", "w")
    f.write(f"NCOLS {ncol} \n") 
    f.write(f"NROWS {nrow} \n")
    f.write(f"XLLCORNER {xmin} \n")
    f.write(f"YLLCORNER {ymax} \n")
    f.write(f"CELLSIZE {cellsize} \n")
    f.write(f"NODATA_VALUE {nodata_value} \n")
    
    for id_row in range(nrow):
        to_write = ""
        data = rast_to_plot[:,id_row]
        for i in range(len(data)):
                    to_write = to_write +" " + str(data[i])
        to_write = to_write + ' \n'
        f.write(to_write)
    f.close()

write_ascii_raster_smash(smash_model = smash_model, pathfile = "/home/livillenave/Documents/distant/SD-FLOOD/synthetic_case_coupling/tmp.asc")




def write_df2d_hydrograph_from_smash_model(write_dir, inflows, smash_model):
    time = np.arange(start = 0, stop = smash_model.setup.dt * len(inflows[0]["qsim"]), step = smash_model.setup.dt ) 
    with open(f'{write_dir}/hydrograph.txt', 'w') as f:
        f.write("#comment\n")
        f.write("#comment\n")
        f.write("#comment\n")
        f.write(str(len(inflows))+"\n")
        for inflow_dic in inflows.values():
            print(inflow_dic)
            f.write("#comment\n")
            f.write("#comment\n")
            f.write("#comment\n")
            f.write(f"{len(time)}\n")
            for i in range(len(time)):
                f.write(f"{time[i]} {inflow_dic['qsim'][i]}\n")
    
    return("")
    

def write_df2d_hydrograph_alamano(write_dir, inflows, time):
    with open(f'{write_dir}/hydrograph.txt', 'w') as f:
        f.write("#comment\n")
        f.write("#comment\n")
        f.write("#comment\n")
        f.write(str(len(inflows))+"\n")
        for inflow_dic in inflows.values():
            print(inflow_dic)
            f.write("#comment\n")
            f.write("#comment\n")
            f.write("#comment\n")
            f.write(f"{len(time)}\n")
            for i in range(len(time)):
                f.write(f"{time[i]} {inflow_dic['qsim'][i]}\n")
    
    return("")
    
###################
# prepare inflow mono bc
##################

# mono bc
inflows_monobc = {0:dict()}
inflows_monobc[0]["qsim"] = np.zeros(shape= inflows[0]["qsim"].shape)
for my_inflow_key in inflows.keys():
    inflows_monobc[0]["qsim"] =     inflows_monobc[0]["qsim"] + inflows[my_inflow_key]["qsim"]
    
# 3 bcs
# 0 <- 0,1,2
# 6 <- 8,6,4
# 5 <- 7,5,3
    
inflows_3bc =  {0:dict(), 5:dict(),6:dict()}  
inflows_3bc[0]["qsim"] = inflows_3bc[5]["qsim"] = inflows_3bc[6]["qsim"] = np.zeros(shape= inflows[0]["qsim"].shape)
for my_inflow_key in inflows.keys():
    print(my_inflow_key)
    if my_inflow_key < 3:
        inflows_3bc[0]["qsim"] =  inflows_3bc[0]["qsim"] + inflows[my_inflow_key]["qsim"]
    if my_inflow_key == 3 or my_inflow_key == 5  or my_inflow_key == 7 :
        inflows_3bc[5]["qsim"] =  inflows_3bc[5]["qsim"]  + inflows[my_inflow_key]["qsim"]
    if my_inflow_key == 4 or my_inflow_key == 6  or my_inflow_key == 8 :
        inflows_3bc[6]["qsim"] =  inflows_3bc[6]["qsim"] + inflows[my_inflow_key]["qsim"]
    
    
###########################
# DASSFLOW EXPERIMENT 
###########################
        
# -------------------------------
# scale mesh = 1000
# -------------------------------

nb_bc = 1

base_dir = f"{case_dir}/bin_dirs/1000"
source_bin = f"{case_dir}/bin_dirs/1000/bin_default"
bin_dir = f"{base_dir}/{nb_bc}"

# copy files from source_bin to  bin_dir
# this account for: input.txt, land_uses.txt, rating_curve.txt
os.system(f"cp {source_bin}/input.txt {bin_dir}/input.txt")
os.system(f"cp {source_bin}/land_uses.txt {bin_dir}")
os.system(f"cp {source_bin}/rating_curve.txt {bin_dir}")


# write hydrograph
time = np.arange(start = 0, stop = smash_model.setup.dt * len(inflows[0]["qsim"]), 
                 step = smash_model.setup.dt )   
time = time[1:]
time_to_add= np.arange(start = 0, stop = 3600, 
                 step = 100)  
time = np.concatenate(  (time_to_add,  time) )



inflows_monobc_alamano = copy.deepcopy(inflows_monobc) 
for key,value in inflows_monobc_alamano.items():
    q0 = value["qsim"][0]
    value["qsim"]= value["qsim"][1:]
    q_to_add = np.repeat(q0, len(time_to_add))
    value["qsim"]= np.concatenate(  (q_to_add,  value["qsim"]) )

inflows_3bc_alamano = copy.deepcopy(inflows_3bc) 
for key,value in inflows_3bc_alamano.items():
    q0 = value["qsim"][0]
    value["qsim"]= value["qsim"][1:]
    q_to_add = np.repeat(q0, len(time_to_add))
    value["qsim"]= np.concatenate(  (q_to_add,  value["qsim"]) )
    
    
inflows_alamano = copy.deepcopy(inflows) 
for key,value in inflows_alamano.items():
    q0 = value["qsim"][0]
    value["qsim"]= value["qsim"][1:]
    q_to_add = np.repeat(q0, len(time_to_add))
    value["qsim"]= np.concatenate(  (q_to_add,  value["qsim"]) )
    
    
    
    
    
if nb_bc == 1:
#    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
#                                       inflows = inflows_monobc, 
#                                       smash_model = smash_model)
    write_df2d_hydrograph_alamano(write_dir=bin_dir, 
                                       inflows = inflows_monobc_alamano, 
                                       time = time)

model = df2d.dassflowmodel(bin_dir =bin_dir , hdf5_path=f"{bin_dir}/res/simu.hdf5", 
                           run_type="direct", clean = True )

my_ts = smash_model.setup.dt * len(inflows_monobc[0]["qsim"])

my_dtw = my_ts/100

my_dtp = my_dtw

config = {"ts":my_ts, "dtw":my_dtw, "dtp":my_dtp  }
model.config.set(config)

model.config.get()
model.init_all()
model.kernel.dof0.h[:]=0.1
model.run()

plt.plot(model.outputs.post.all_time, 
         model.outputs.post.sum_mass_flux_inflow,
         label = "inflows", color = "blue")
plt.plot(model.outputs.post.all_time, 
         model.outputs.post.sum_mass_flux_outflow,
         label = "outflows", color = "red")
plt.title("sum mass flux")
plt.legend()
plt.show()
plt.close()


fig, ax1 = plt.subplots()
ax1.plot(model.outputs.post.all_time, 
         model.outputs.post.sum_q_inflow,
         label = "inflows", color = "blue")
plt.legend()
plt.title("sum q")

ax2 = ax1.twinx()
ax2.plot(model.outputs.post.all_time, 
         model.outputs.post.sum_q_outflow,
         label = "outflows", color = "red")
plt.legend()
plt.show()
plt.close()

#
#for i in range(model.outputs.result.h.shape[1]):
#    model.meshing.mesh_pyvista.plot( scalars = model.outputs.result.h[:,i])
#

#for i in range(model.outputs.result.zs.shape[1]):
#    model.meshing.mesh_pyvista.plot( scalars = model.outputs.result.zs[:,i])
#    

df2d.wrapping.call_model.clean_model(model.kernel)
# -------------------------------
# scale mesh = 100
# -------------------------------

nb_bc = 1

base_dir = f"{case_dir}/bin_dirs/100"
source_bin = f"{case_dir}/bin_dirs/100/bin_default"
bin_dir = f"{base_dir}/{nb_bc}"


# copy files from source_bin to  bin_dir
# this account for: input.txt, land_uses.txt, rating_curve.txt
os.system(f"cp {source_bin}/input.txt {bin_dir}/input.txt")
os.system(f"cp {source_bin}/land_uses.txt {bin_dir}")
os.system(f"cp {source_bin}/rating_curve.txt {bin_dir}")

# write hydrograph
if nb_bc == 1:
    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
                                       inflows = inflows_monobc, 
                                       smash_model = smash_model)



model_100 = df2d.dassflowmodel(bin_dir =bin_dir , hdf5_path=f"{bin_dir}/res/simu.hdf5", run_type="direct", clean = True)
my_ts = smash_model.setup.dt * len(inflows_monobc[0]["qsim"])
my_dtw = my_ts/100
my_dtp = my_dtw
config = {"ts":my_ts, "dtw":my_dtw, "dtp":my_dtp  }
model_100.config.set(config)

model_100.config.get()
model_100.init_all()
model_100.kernel.dof0.h[:]=0.1
model_100.run()    
df2d.wrapping.call_model.clean_model(model_100.kernel)



# ==================== #
# COMPARE 1BC for different scales
# ==================== #

all_model = {"1000":model, "100":model_100}
all_color = {"1000":"blue", "100":"red"}

# compare imposed discharge

for my_scale in all_model.keys():
    print(my_scale)
    all_model[my_scale].boundary.plot()
    plt.show()
    plt.close()
    
# compare sum mass flux inflow
plt.plot(np.arange(0,48*3600, 3600), sum_inflows, label = "smash sum qin")
for my_scale in all_model.keys():    
    plt.plot(all_model[my_scale].outputs.post.all_time, 
             all_model[my_scale].outputs.post.sum_q_inflow,
             label = my_scale, color = all_color[my_scale])
plt.title("sum q inflow, for different mesh scale, for 1 unique inflow")
plt.legend(title = "mesh scale")
plt.show()
plt.close()

# compare sum mass q inflow
plt.plot(np.arange(0,48*3600, 3600), sum_inflows, label = "smash sum qin")
for my_scale in all_model.keys():    
    plt.plot(all_model[my_scale].outputs.post.all_time, 
             all_model[my_scale].outputs.post.sum_mass_flux_inflow,
             label = my_scale, color = all_color[my_scale])
plt.title("sum mass flux inflow, for different mesh scale, for 1 unique inflow")
plt.legend()
plt.show()
plt.close()



# compare sum mass flux outflow
plt.plot(np.arange(0,48*3600, 3600), outflow, label = "smash qout")
for my_scale in all_model.keys():    
    plt.plot(all_model[my_scale].outputs.post.all_time, 
             all_model[my_scale].outputs.post.sum_q_outflow,
             label = my_scale, color = all_color[my_scale])
plt.title("sum q outflow, for different mesh scale, for 1 unique inflow")
plt.legend(title = "mesh scale")
plt.show()
plt.close()

# compare sum q outflow
plt.plot(np.arange(0,48*3600, 3600), outflow, label = "smash qout")
for my_scale in all_model.keys():    
    plt.plot(all_model[my_scale].outputs.post.all_time, 
             all_model[my_scale].outputs.post.sum_mass_flux_outflow,
             label = my_scale, color = all_color[my_scale])
plt.title("sum mass flux outflow, for different mesh scale, for 1 unique inflow")
plt.legend()
plt.show()
plt.close()


# compare sum q outflow
for my_scale in all_model.keys():    
    plt.plot(all_model[my_scale].outputs.post.all_time, 
             all_model[my_scale].outputs.post.water_vol,
             label = my_scale, color = all_color[my_scale])
plt.title("sum mass flux outflow, for different mesh scale, for 1 unique inflow")
plt.legend()
plt.show()
plt.close()

#plt.plot(np.arange(0,48*3600, 3600), sum_inflows - outflow, label = "smash sum qin - outflow")
#for my_scale in all_model.keys():    
#    plt.plot(all_model[my_scale].outputs.post.all_time, 
#             all_model[my_scale].outputs.post.sum_mass_flux_inflow - 
#             all_model[my_scale].outputs.post.sum_mass_flux_outflow,
#             label = my_scale, color = all_color[my_scale])
#
#plt.title("sum mass flux inflow,- sum mass flow outflow for different mesh scale, for 1 unique inflow")
#plt.legend()
#plt.show()



fig, axs= plt.subplots(2)

axs[0].plot(np.arange(0,48*3600, 3600), sum_inflows, label = "smash ")
axs[1].plot(np.arange(0,48*3600, 3600), outflow, label = "smash ")

for my_scale in all_model.keys():    

    axs[0].plot(all_model[my_scale].outputs.post.all_time, 
             all_model[my_scale].outputs.post.sum_mass_flux_inflow,
             label = my_scale, color = all_color[my_scale])

    axs[1].plot(all_model[my_scale].outputs.post.all_time, 
             all_model[my_scale].outputs.post.sum_mass_flux_outflow,
             label = my_scale, color = all_color[my_scale])

axs[0].legend()
axs[1].legend()
axs[0].title.set_text("Inflow")
axs[1].title.set_text("Outflows")
plt.show()



fig, axs= plt.subplots(2)
axs[0].plot(np.arange(0,48*3600, 3600),  np.cumsum(sum_inflows), label = "smash ")
axs[1].plot(np.arange(0,48*3600, 3600),  np.cumsum(outflow), label = "smash ")

for my_scale in all_model.keys():    

    axs[0].plot(all_model[my_scale].outputs.post.all_time, 
             np.cumsum(all_model[my_scale].outputs.post.sum_mass_flux_inflow),
             label = my_scale, color = all_color[my_scale])

    axs[1].plot(all_model[my_scale].outputs.post.all_time, 
              np.cumsum(all_model[my_scale].outputs.post.sum_mass_flux_outflow),
             label = my_scale, color = all_color[my_scale])

axs[0].legend()
axs[1].legend()
axs[0].title.set_text("Inflow")
axs[1].title.set_text("Outflows")
plt.legend()
plt.show()




##############################
# SAME BUT WITH 3 BOUNDARIES
##############################

# -------------------------------
# scale mesh = 1000, nc_bc == n
# -------------------------------
nb_bc = 3


base_dir = f"{case_dir}/bin_dirs/1000"
source_bin = f"{case_dir}/bin_dirs/1000/bin_default"
bin_dir = f"{base_dir}/{nb_bc}"


# copy files from source_bin to  bin_dir
# this account for: input.txt, land_uses.txt, rating_curve.txt
os.system(f"cp {source_bin}/input.txt {bin_dir}/input.txt")
os.system(f"cp {source_bin}/land_uses.txt {bin_dir}")
os.system(f"cp {source_bin}/rating_curve.txt {bin_dir}")

# write hydrograph
if nb_bc == 1:
#    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
#                                       inflows = inflows_monobc, 
#                                       smash_model = smash_model)
    
    write_df2d_hydrograph_alamano(write_dir=bin_dir, 
                                       inflows = inflows_monobc_alamano, 
                                       time = time)

if nb_bc == 3:
#    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
#                                       inflows = inflows_3bc, 
#                                       smash_model = smash_model)    
    write_df2d_hydrograph_alamano(write_dir=bin_dir, 
                                       inflows = inflows_3bc_alamano, 
                                       time = time)




model_1000_3bc = df2d.dassflowmodel(bin_dir =bin_dir , hdf5_path=f"{bin_dir}/res/simu.hdf5", run_type="direct", clean = True)
my_ts = smash_model.setup.dt * len(inflows_monobc[0]["qsim"])
my_dtw = my_ts/100
my_dtp = my_ts/10000
config = {"ts":my_ts, "dtw":my_dtw, "dtp":my_dtp  }
model_1000_3bc.config.set(config)

model_1000_3bc.config.get()
model_1000_3bc.init_all()
model_1000_3bc.kernel.dof0.h[:]=0.1
model_1000_3bc.run()

plt.plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "inflows", color = "blue")
plt.plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "outflows", color = "red")
plt.title("sum mass flux")
plt.legend()
plt.show()
plt.close()


fig, ax1 = plt.subplots()
ax1.plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_q_inflow, axis = 1) ,
         label = "inflows", color = "blue")
plt.legend()
plt.title("sum q")
ax2 = ax1.twinx()
ax2.plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_q_outflow, axis = 1),
         label = "outflows", color = "red")
plt.legend()
plt.show()
plt.close()
#    
df2d.wrapping.call_model.clean_model(model_1000_3bc.kernel)




plt.plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "inflows", color = "blue")
plt.plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "outflows", color = "red")
plt.title("sum mass flux")
plt.legend()
plt.show()
plt.close()


plt.plot(model.outputs.post.all_time, 
         np.sum(model.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "inflows", color = "blue")
plt.plot(model.outputs.post.all_time, 
         np.sum(model.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "outflows", color = "red")
plt.title("sum mass flux")
plt.legend()
plt.show()
plt.close()




fig, axs= plt.subplots(2)
axs[0].plot(np.arange(0,48*3600, 3600), sum_inflows, label = "smash ")
axs[1].plot(np.arange(0,48*3600, 3600), outflow, label = "smash ")


axs[0].plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "3bc", color = "blue")

axs[0].plot(model.outputs.post.all_time, 
         np.sum(model.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "1bc", color = "red")


axs[1].plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "3bc", color = "blue")

axs[1].plot(model.outputs.post.all_time, 
         np.sum(model.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "1bc", color = "red")

axs[0].legend()
axs[1].legend()
axs[0].title.set_text("Inflow")
axs[1].title.set_text("Outflows")
plt.suptitle("mesh scale = 1000")
plt.legend()
plt.show()



df2d.wrapping.call_model.clean_model(model_1000_3bc.kernel)


# -------------------------------
# scale mesh = 1000, nc_bc == n
# -------------------------------

nb_bc = "n"
base_dir = "/home/livillenave/Documents/distant/SD-FLOOD/synthetic_case_coupling/new_0603/mesh_gen/MESHINGS/1000"
source_bin = "/home/livillenave/Documents/distant/SD-FLOOD/synthetic_case_coupling/new_0603/mesh_gen/MESHINGS/1000/bin_default"

# bin_dir already contain bc.txt and mesh.geo adapted files to nb_bc (1,3,n)
bin_dir = f"{base_dir}/{nb_bc}"


# copy files from source_bin to  bin_dir
# this account for: input.txt, land_uses.txt, rating_curve.txt
os.system(f"cp {source_bin}/input.txt {bin_dir}/input.txt")
os.system(f"cp {source_bin}/land_uses.txt {bin_dir}")
os.system(f"cp {source_bin}/rating_curve.txt {bin_dir}")

# write hydrograph
if nb_bc == 1:
#    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
#                                       inflows = inflows_monobc, 
#                                       smash_model = smash_model)
    
    write_df2d_hydrograph_alamano(write_dir=bin_dir, 
                                       inflows = inflows_monobc_alamano, 
                                       time = time)

if nb_bc == 3:
#    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
#                                       inflows = inflows_3bc, 
#                                       smash_model = smash_model)    
    write_df2d_hydrograph_alamano(write_dir=bin_dir, 
                                       inflows = inflows_3bc_alamano, 
                                       time = time)
    
if nb_bc == "n":
    write_df2d_hydrograph_alamano(write_dir=bin_dir, 
                                  inflows = inflows_alamano,
                                       time = time)


model_1000_nbc = df2d.dassflowmodel(bin_dir =bin_dir , hdf5_path=f"{bin_dir}/res/simu.hdf5", run_type="direct", clean = True)
my_ts = smash_model.setup.dt * len(inflows_monobc[0]["qsim"])
my_dtw = my_ts/100
my_dtp = my_ts/10000
config = {"ts":my_ts, "dtw":my_dtw, "dtp":my_dtp  }
model_1000_nbc.config.set(config)

model_1000_nbc.config.get()

#os.chdir(model_1000_nbc.bin_dir)
#print("call dassflow2d.wrapping.m_mesh.msh()")
#model_1000_nbc.kernel.mesh = df2d.wrapping.m_mesh.msh()
#print("call         dassflow2d.wrapping.call_model.init_solver(self.kernel)")
#df2d.wrapping.call_model.init_solver(model_1000_nbc.kernel)
#print("call   Meshing(mesh_fortran = self.kernel.mesh)")
#model_1000_nbc.meshing = Meshing(mesh_fortran = model_1000_nbc.kernel.mesh)

for i in range(model_1000_nbc.kernel.mesh.neb) :
    print(model_1000_nbc.kernel.mesh.edgeb[0].typlim)

model_1000_nbc.init_all()
model_1000_nbc.kernel.dof0.h[:]=0.1
model_1000_nbc.run()

plt.plot(model_1000_nbc.outputs.post.all_time, 
         np.sum(model_1000_nbc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "inflows", color = "blue")
plt.plot(model_1000_nbc.outputs.post.all_time, 
         np.sum(model_1000_nbc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "outflows", color = "red")
plt.title("sum mass flux")
plt.legend()
plt.show()
plt.close()


fig, ax1 = plt.subplots()
ax1.plot(model_1000_nbc.outputs.post.all_time, 
         np.sum(model_1000_nbc.outputs.post.sum_q_inflow, axis = 1) ,
         label = "inflows", color = "blue")
plt.legend()
plt.title("sum q")
ax2 = ax1.twinx()
ax2.plot(model_1000_nbc.outputs.post.all_time, 
         np.sum(model_1000_nbc.outputs.post.sum_q_outflow, axis = 1),
         label = "outflows", color = "red")
plt.legend()
plt.show()
plt.close()
#    
df2d.wrapping.call_model.clean_model(model_1000_nbc.kernel)




plt.plot(model_1000_nbc.outputs.post.all_time, 
         np.sum(model_1000_nbc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "inflows", color = "blue")
plt.plot(model_1000_nbc.outputs.post.all_time, 
         np.sum(model_1000_nbc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "outflows", color = "red")
plt.title("sum mass flux")
plt.legend()
plt.show()
plt.close()


plt.plot(model.outputs.post.all_time, 
         np.sum(model.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "inflows", color = "blue")
plt.plot(model.outputs.post.all_time, 
         np.sum(model.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "outflows", color = "red")
plt.title("sum mass flux")
plt.legend()
plt.show()
plt.close()




fig, axs= plt.subplots(2)
axs[0].plot(np.arange(0,48*3600, 3600), sum_inflows, label = "smash ")
axs[1].plot(np.arange(0,48*3600, 3600), outflow, label = "smash ")


axs[0].plot(model_1000_nbc.outputs.post.all_time, 
         np.sum(model_1000_nbc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "nbc", color = "black", marker = ".")

axs[0].plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "3bc", color = "blue")


axs[0].plot(model.outputs.post.all_time, 
         np.sum(model.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "1bc", color = "red")


axs[1].plot(model_1000_nbc.outputs.post.all_time, 
         np.sum(model_1000_nbc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "nbc", color = "black")

axs[1].plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "3bc", color = "blue")

axs[1].plot(model.outputs.post.all_time, 
         np.sum(model.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "1bc", color = "red")

axs[0].legend()
axs[1].legend()
axs[0].title.set_text("Inflow")
axs[1].title.set_text("Outflows")
plt.suptitle("mesh scale = 1000")
plt.legend()
plt.show()




df2d.wrapping.call_model.clean_model(model_1000_nbc.kernel)

####################################
# MESH SCALE == 100
####################################

# -------------------------------
# scale mesh = 100, nbc = 1
# -------------------------------
nb_bc = 1

base_dir = f"{case_dir}/bin_dirs/100"
source_bin = f"{case_dir}/bin_dirs/100/bin_default"
bin_dir = f"{base_dir}/{nb_bc}"


# copy files from source_bin to  bin_dir
# this account for: input.txt, land_uses.txt, rating_curve.txt
os.system(f"cp {source_bin}/input.txt {bin_dir}/input.txt")
os.system(f"cp {source_bin}/land_uses.txt {bin_dir}")
os.system(f"cp {source_bin}/rating_curve.txt {bin_dir}")

# write hydrograph

# write hydrograph
if nb_bc == 1:
    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
                                       inflows = inflows_monobc, 
                                       smash_model = smash_model)
if nb_bc == 3:
    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
                                       inflows = inflows_3bc, 
                                       smash_model = smash_model)
if nb_bc == "n":
    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
                                       inflows = inflows, 
                                       smash_model = smash_model)



model_100 = df2d.dassflowmodel(bin_dir =bin_dir , hdf5_path=f"{bin_dir}/res/simu.hdf5", run_type="direct", clean = True)
my_ts = smash_model.setup.dt * len(inflows_monobc[0]["qsim"])
my_dtw = my_ts/100
my_dtp = my_dtw
config = {"ts":my_ts, "dtw":my_dtw, "dtp":my_dtp  }
model_100.config.set(config)

model_100.config.get()
model_100.init_all()
model_100.kernel.dof0.h[:]=0.1
model_100.run()

df2d.wrapping.call_model.clean_model(model_100.kernel)



# -------------------------------
# scale mesh = 100, nbc = 3
# -------------------------------

nb_bc = 3

base_dir = f"{case_dir}/bin_dirs/100"
source_bin = f"{case_dir}/bin_dirs/100/bin_default"
bin_dir = f"{base_dir}/{nb_bc}"

# copy files from source_bin to  bin_dir
# this account for: input.txt, land_uses.txt, rating_curve.txt
os.system(f"cp {source_bin}/input.txt {bin_dir}/input.txt")
os.system(f"cp {source_bin}/land_uses.txt {bin_dir}")
os.system(f"cp {source_bin}/rating_curve.txt {bin_dir}")

# write hydrograph
if nb_bc == 1:
    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
                                       inflows = inflows_monobc, 
                                       smash_model = smash_model)
if nb_bc == 3:
    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
                                       inflows = inflows_3bc, 
                                       smash_model = smash_model)
if nb_bc == "n":
    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
                                       inflows = inflows, 
                                       smash_model = smash_model)




model_100_3bc = df2d.dassflowmodel(bin_dir =bin_dir , hdf5_path=f"{bin_dir}/res/simu.hdf5", run_type="direct", clean = True)
my_ts = smash_model.setup.dt * len(inflows_monobc[0]["qsim"])
my_dtw = my_ts/100
my_dtp = my_dtw
config = {"ts":my_ts, "dtw":my_dtw, "dtp":my_dtp  }
model_100_3bc.config.set(config)

model_100_3bc.config.get()
model_100_3bc.init_all()
model_100_3bc.kernel.dof0.h[:]=0.1
model_100_3bc.run()

df2d.wrapping.call_model.clean_model(model_100_3bc.kernel)


# -------------------------------
# scale mesh = 100
# -------------------------------

nb_bc = "n"

base_dir = f"{case_dir}/bin_dirs/100"
source_bin = f"{case_dir}/bin_dirs/100/bin_default"
bin_dir = f"{base_dir}/{nb_bc}"

# copy files from source_bin to  bin_dir
# this account for: input.txt, land_uses.txt, rating_curve.txt
os.system(f"cp {source_bin}/input.txt {bin_dir}/input.txt")
os.system(f"cp {source_bin}/land_uses.txt {bin_dir}")
os.system(f"cp {source_bin}/rating_curve.txt {bin_dir}")


# write hydrograph
if nb_bc == 1:
    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
                                       inflows = inflows_monobc, 
                                       smash_model = smash_model)
if nb_bc == 3:
    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
                                       inflows = inflows_3bc, 
                                       smash_model = smash_model)
if nb_bc == "n":
    write_df2d_hydrograph_from_smash_model(write_dir=bin_dir, 
                                       inflows = inflows, 
                                       smash_model = smash_model)




model_100_nbc = df2d.dassflowmodel(bin_dir =bin_dir , hdf5_path=f"{bin_dir}/res/simu.hdf5", run_type="direct", clean = True)
my_ts = smash_model.setup.dt * len(inflows_monobc[0]["qsim"])
my_dtw = my_ts/100
my_dtp = my_dtw
config = {"ts":my_ts, "dtw":my_dtw, "dtp":my_dtp  }
model_100_nbc.config.set(config)

model_100_nbc.config.get()
model_100_nbc.init_all()
model_100_nbc.kernel.dof0.h[:]=0.1
model_100_nbc.run()

df2d.wrapping.call_model.clean_model(model_100_nbc.kernel)




fig, axs= plt.subplots(2)
fig.set_size_inches(w=15,h=8)
fig.tight_layout(pad=2.0)
axs[0].plot(np.arange(0,48*3600, 3600), sum_inflows, label = "smash ")
axs[1].plot(np.arange(0,48*3600, 3600), outflow, label = "smash ")


axs[0].plot(model_100_nbc.outputs.post.all_time, 
         np.sum(model_100_nbc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "DassFlow - 9 inflows", color = "black")
axs[0].plot(model_100_3bc.outputs.post.all_time, 
         np.sum(model_100_3bc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "DassFlow - 3  inflows", color = "blue")
axs[0].plot(model_100.outputs.post.all_time, 
         np.sum(model_100.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "DassFlow - 1  inflow", color = "red")


axs[1].plot(model_100_nbc.outputs.post.all_time, 
         np.sum(model_100_nbc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "DassFlow - 9 inflows", color = "black")
axs[1].plot(model_100_3bc.outputs.post.all_time, 
         np.sum(model_100_3bc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "DassFlow - 3 inflows", color = "blue")
axs[1].plot(model_100.outputs.post.all_time, 
         np.sum(model_100.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "DassFlow - 1 inflow", color = "red")

axs[0].legend()
axs[1].legend()
axs[0].title.set_text(" Sum of inflow discharge")
axs[1].title.set_text("Outflow discharge")

axs[1].set_xlabel("time [s]")
axs[0].set_ylabel(r"$Discharge [m^3/s]$")
axs[1].set_ylabel(r"$Discharge [m^3/s]$")
#plt.suptitle("mesh scale = 100")
plt.legend()
plt.show()


fig, axs= plt.subplots(2,  figsize=(15, 15))
axs[0].plot(np.arange(0,48*3600, 3600), sum_inflows, label = "smash ")
axs[1].plot(np.arange(0,48*3600, 3600), outflow, label = "smash ")



axs[0].plot(model_100_nbc.outputs.post.all_time, 
         np.sum(model_100_nbc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "nbc", color = "black")

axs[0].plot(model_100_3bc.outputs.post.all_time, 
         np.sum(model_100_3bc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "3bc", color = "blue")


axs[0].plot(model_100.outputs.post.all_time, 
         np.sum(model_100.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "1bc", color = "red")


axs[1].plot(model_100_nbc.outputs.post.all_time, 
         np.sum(model_100_nbc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "nbc", color = "black")

axs[1].plot(model_100_3bc.outputs.post.all_time, 
         np.sum(model_100_3bc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "3bc", color = "blue")

axs[1].plot(model_100.outputs.post.all_time, 
         np.sum(model_100.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "1bc", color = "red")

axs[0].legend()
axs[1].legend()
axs[0].title.set_text("Inflow")
axs[1].title.set_text("Outflows")
#plt.suptitle("mesh scale = 100")
plt.legend()
plt.show()






plt.plot(np.arange(0,48*3600, 3600), sum_inflows, label = "sum inflows ")
plt.plot(np.arange(0,48*3600, 3600), outflow, label = "sum outflow ")
plt.plot(np.arange(0,48*3600, 3600), np.repeat(1111.1,48), label = "raw rainfall contributing to the outlet")
plt.legend()
plt.xlabel("time [s]")
plt.ylabel("Water volume fluxes $[m^3/s]$")
plt.show()
plt.close()






fig, axs= plt.subplots(2)
axs[0].plot(np.arange(0,48*3600, 3600), sum_inflows, label = "smash ")
axs[1].plot(np.arange(0,48*3600, 3600), outflow, label = "smash ")


axs[0].plot(model_100_nbc.outputs.post.all_time, 
         np.sum(model_100_nbc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "100-nbc", color = "black")

axs[0].plot(model_100_3bc.outputs.post.all_time, 
         np.sum(model_100_3bc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "100-3bc", color = "blue")


axs[0].plot(model.outputs.post.all_time, 
         np.sum(model.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "100-1bc", color = "red")

axs[0].plot(model_1000_nbc.outputs.post.all_time, 
         np.sum(model_1000_nbc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "1000-nbc", color = "black", marker = "1")
axs[0].plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "1000-3bc", color = "blue", marker = "1")
axs[0].plot(model.outputs.post.all_time, 
        np.sum(model.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "1000-1bc", color = "red", marker = "1")



axs[1].plot(model_100_nbc.outputs.post.all_time, 
         np.sum(model_100_nbc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "100-nbc", color = "black")

axs[1].plot(model_100_3bc.outputs.post.all_time, 
         np.sum(model_100_3bc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "100-3bc", color = "blue")


axs[1].plot(model.outputs.post.all_time, 
         np.sum(model.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "100-1bc", color = "red")


axs[1].plot(model_1000_nbc.outputs.post.all_time, 
         np.sum(model_1000_nbc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "1000-nbc", color = "black", marker = "1")
axs[1].plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "1000-3bc", color = "blue", marker = "1")
axs[1].plot(model.outputs.post.all_time, 
        np.sum(model.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "1000-1bc", color = "red", marker = "1")

axs[0].legend()
axs[1].legend()
axs[0].title.set_text("INFLOW")
axs[1].title.set_text("OUTFLOW")
plt.suptitle("compare mesh scales")
plt.legend()
plt.show()




fig, axs= plt.subplots(3)

axs[0].plot(np.arange(0,48*3600, 3600), outflow, label = "smash ")
axs[1].plot(np.arange(0,48*3600, 3600), outflow, label = "smash ")
axs[2].plot(np.arange(0,48*3600, 3600), outflow, label = "smash ")


axs[0].plot(model_100_nbc.outputs.post.all_time, 
         np.sum(model_100_nbc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "100-nbc", color = "blue")
axs[0].plot(model_1000_nbc.outputs.post.all_time, 
         np.sum(model_1000_nbc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "1000-nbc", color = "red", marker = "1")


axs[1].plot(model_100_3bc.outputs.post.all_time, 
         np.sum(model_100_3bc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "100-3bc", color = "blue")
axs[1].plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "1000-3bc", color = "red", marker = "1")


axs[2].plot(model.outputs.post.all_time, 
         np.sum(model.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "100-1bc", color = "blue")
axs[2].plot(model.outputs.post.all_time, 
        np.sum(model.outputs.post.sum_mass_flux_outflow, axis = 1),
         label = "1000-1bc", color = "red", marker = "1")


axs[0].legend()
axs[1].legend()
axs[0].title.set_text("n bc")
axs[1].title.set_text("3 bc")
axs[2].title.set_text("1 bc")
plt.suptitle("compare outflows")
plt.legend()
plt.show()


fig, axs= plt.subplots(3)

axs[0].plot(np.arange(0,48*3600, 3600), sum_inflows, label = "smash ", marker = ".", color = "black")
axs[1].plot(np.arange(0,48*3600, 3600), sum_inflows, label = "smash ", marker = ".",  color = "black")
axs[2].plot(np.arange(0,48*3600, 3600), sum_inflows, label = "smash ", marker = ".",  color = "black")

axs[0].plot(model_100_nbc.outputs.post.all_time, 
         np.sum(model_100_nbc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "100-nbc", color = "blue", marker = "+")
axs[0].plot(model_1000_nbc.outputs.post.all_time, 
         np.sum(model_1000_nbc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "1000-nbc", color = "red")

axs[1].plot(model_100_3bc.outputs.post.all_time, 
         np.sum(model_100_3bc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "100-3bc", color = "blue", marker = "+")

axs[1].plot(model_1000_3bc.outputs.post.all_time, 
         np.sum(model_1000_3bc.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "1000-3bc", color = "red")

axs[2].plot(model.outputs.post.all_time, 
         np.sum(model.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "100-1bc", color = "blue", marker = "+")


axs[2].plot(model.outputs.post.all_time, 
        np.sum(model.outputs.post.sum_mass_flux_inflow, axis = 1),
         label = "1000-1bc", color = "red")

axs[0].legend()
axs[1].legend()
axs[0].title.set_text("n bc")
axs[1].title.set_text("3 bc")
axs[2].title.set_text("1 bc")
plt.suptitle("compare inflows")
plt.legend()
plt.show()



##############################
# Find back smash parameters
##############################

#
#for i in range(model_100_nbc.outputs.result.h.shape[1]):
#    model_100_nbc.meshing.mesh_pyvista.plot( scalars = model_100_nbc.outputs.result.zs[:,i], cpos = "xy", show_edges = True)
##
    


smash_discharge_times = np.arange(start = 0 , 
                                  step = 3600, 
                                  stop = 3600*
                                  len(smash_model.output.qsim[0]) )

smash_observed_discharges = np.ndarray(shape = smash_model.input_data.qobs.shape)

for i in range(len(smash_discharge_times)):
    if i<len(smash_discharge_times)-1:
        time0=smash_discharge_times[i]
        time1=smash_discharge_times[i+1]
    else:
        time0=smash_discharge_times[i]
        time1=smash_discharge_times[i]*100
        
    index_dassflow = np.where(np.logical_and(model_100_nbc.outputs.post.all_time>=time0,  model_100_nbc.outputs.post.all_time<time1))[0]
    smash_observed_discharges[0,i] = np.mean(model_100_nbc.outputs.post.sum_mass_flux_outflow[index_dassflow])


smash_model_infer = smash_model.copy()
smash_model_infer.input_data.qobs = smash_observed_discharges

smash_model_infer.optimize(mapping = "uniform", inplace = True)
smash_model_infer.optimize(mapping = "distributed", inplace = True)




plt.plot(smash_discharge_times ,
smash_model_infer.output.qsim[0,:], color = "red", label = "sim")
plt.plot(smash_discharge_times , smash_model_infer.input_data.qobs[0,:], color = "blue", label = "obs" )
plt.plot(smash_discharge_times , smash_model.output.qsim[0], color = "black", label = "prior")
plt.legend()
plt.show()
plt.close()


plt.imshow(smash_model_infer.parameters.cp); plt.colorbar(); plt.show();plt.close()
plt.imshow(smash_model_infer.parameters.cft); plt.colorbar(); plt.show();plt.close()
plt.imshow(smash_model_infer.parameters.lr); plt.colorbar(); plt.show();plt.close()
plt.imshow(smash_model_infer.parameters.exc); plt.colorbar(); plt.show();plt.close()




# Dassflow has been built manually

tmp_gauge_pos = np.ndarray(shape = (len(inflows)+1, 2))
tmp_gauge_pos[0,:]=data["mesh"]["gauge_pos"]

tmp_code = np.ndarray(shape = (len(inflows)+1, 1), dtype='<U13')
tmp_code[0] = data["mesh"]["code"] 
compteur = 0
tmp_code2 = []
tmp_code2.append( str(data["mesh"]["code"][0] ))
for my_key in inflows.keys():
    compteur = compteur+1
    #tmp_code[compteur] = my_key
    tmp_code2.append(str(np.asanyarray(my_key, dtype='<U13')))
    tmp_gauge_pos[compteur,:] = np.asanyarray(inflows[my_key]["id"])
    
    
new_data = copy.deepcopy(data)
new_data["mesh"]["ng"]  = len(inflows)+1
new_data["mesh"]["gauge_pos"] = tmp_gauge_pos
new_data["mesh"]["code"] = tmp_code2

smash_model_inference = smash.Model( new_data["setup"], 
                                     new_data["mesh"])  

for i in range(smash_model_inference.input_data.qobs.shape[0]):
    if i == 0:
        smash_model_inference.input_data.qobs[i,:] = smash_observed_discharges[0]
    else:
        discharge = inflows[ int(new_data["mesh"]["code"][i])]["qsim"]
        smash_model_inference.input_data.qobs[i,:] = discharge

smash_model_inference.input_data.prcp = np.broadcast_to(data["rain"], smash_model.input_data.prcp.shape)
smash_model_inference.input_data.pet = data["pet"]
smash_model_inference.parameters.lr=50
smash_model_inference.parameters.lr[9,9] = smash_model.parameters.lr[8,8] = smash_model.parameters.lr[7,7]= smash_model.parameters.lr[6,6] = 1
smash_model_inference.parameters.cft=100
smash_model_inference.parameters.cp=200

smash_model_inference.run(inplace=True)
smash_model_inference.optimize(mapping = "uniform",gauge="all", wgauge=[0.0001,0.1,
                                                                        0.1,0.1,
                                                                        0.1,0.1,
                                                                        0.1,0.1,
                                                                        0.1, 0.1], inplace = True)#[ "0", "Practice_case"]

    
smash_model_inference.optimize(mapping = "distributed",gauge="all", wgauge=[0.0001,0.1,
                                                                            0.1,0.1,
                                                                            0.1,0.1,
                                                                            0.1,0.1,
                                                                            0.1, 0.1], inplace = True)


for i in range(smash_model_inference.output.qsim.shape[0]):
    plt.plot(smash_discharge_times ,
             smash_model_inference.output.qsim[i,:], color = "red", label = "sim")
    
    plt.plot(smash_discharge_times, 
             smash_model_inference.input_data.qobs[i,:], color = "blue", label = "obs" )
    plt.title(f'inflow {new_data["mesh"]["code"][i] } ')
    plt.legend()
    plt.show()
    plt.close()    

plt.imshow(smash_model_inference.parameters.cp); plt.colorbar(); plt.show();plt.close()
plt.imshow(smash_model_inference.parameters.cft); plt.colorbar(); plt.show();plt.close()
plt.imshow(smash_model_inference.parameters.lr); plt.colorbar(); plt.show();plt.close()
plt.imshow(smash_model_inference.parameters.exc); plt.colorbar(); plt.show();plt.close()




