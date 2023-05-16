#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 14:26:49 2023

@author: livillenave
"""

# ======> Require dassflow_model instance in the buffer

import dassflow2d as df2d
import gmsh
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import os        
import pandas as pd
from shapely import geometry
from osgeo import gdal
import smash
import pyvista as pv

# import developed librairies
os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/Tools/1_pre-treatment/2_shapefile_to_gmsh/organised_version/libs")
from GIS_libs import *
from dassflow_mesh_manipulation_libs import *
from smash_coupling_libs import *
from gmsh_libs import *
import numpy as np
import shutil

demo_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/Demo"
source_dir =  f"{demo_dir}/source_file/"
path_minor_bed = f"{source_dir}/hydraulic_shapefile/minor_bed.shp"


def get_cells_in_minor_bed(dassflow_model, path_minor_bed):
    
    x,y = get_xy_cells(dassflow_model.meshing.mesh_fortran)
    points = gpd.points_from_xy(x = x, y = y)
    test_points =gpd.geoseries.GeoSeries(points)
    raw_minor_bed = gpd.read_file( filename = f"{path_minor_bed}")  
    
    points_in_major_bed = []
    for i in range(len(test_points)):
        mytest = raw_minor_bed.contains(test_points[i])
        mytest = np.asarray(mytest, dtype = "bool")        
        if mytest:
            points_in_major_bed.append(i)
    return(points_in_major_bed)



# ---------------------- #
# Dassflow parameter
# ---------------------- #
dassflow_dtw = 10
dassflow_dtp = 10

dassflow_dtw_heating = 10
dassflow_dtp_heating = 10


# ---------------------- #
#run
# ---------------------- #
bin_dir = f"{demo_dir}/bin2"
print(bin_dir)

dassflow_model = df2d.dassflowmodel(bin_dir = bin_dir, 
                           hdf5_path=bin_dir+"/res/simu.hdf5", 
                           run_type= "direct",  
                           clean = True)

df2d.wrapping.m_common.set_mesh_name("final_mesh.geo")

df2d.wrapping.m_common.set_ts(  (len(qin_heating[0])-1) * dt_hydrology  )
df2d.wrapping.m_common.set_dtw( dassflow_dtw_heating )
df2d.wrapping.m_common.set_dtp( dassflow_dtp_heating )
dassflow_model.init_all()

# ------------------------- FORCE INITIAL DOF

points_in_major_bed = get_cells_in_minor_bed(dassflow_model = dassflow_model, 
                                             path_minor_bed = path_minor_bed)

#test_points[points_in_major_bed].plot()
nb_cell = dassflow_model.meshing.mesh_fortran.nc 
initial_h = np.zeros( nb_cell)
initial_h[points_in_major_bed] = 1

dassflow_model.kernel.dof.h[:nb_cell] = dassflow_model.kernel.dof0.h[:nb_cell] = initial_h[:]
# manually set correctly  on boundary cell
dassflow_model.kernel.dof.h[474] = dassflow_model.kernel.dof0.h[474] =  1
dassflow_model.kernel.dof0.h[758] = dassflow_model.kernel.dof.h[758] = dassflow_model.kernel.dof0.h[5842]= dassflow_model.kernel.dof.h[5842] = dassflow_model.kernel.dof.h[474]


# manually set correctly bc
#dassflow_model.kernel.dof.h[474] = dassflow_model.kernel.dof0.h[474] =   1
#dassflow_model.kernel.dof0.h[758] = dassflow_model.kernel.dof.h[758] = dassflow_model.kernel.dof0.h[5842]= dassflow_model.kernel.dof.h[5842] = dassflow_model.kernel.dof.h[474]
dassflow_model.run()
#
#for i in dassflow_model.outputs.all_res.keys():
#    
#    dassflow_model.meshing.mesh_pyvista.plot(scalars = dassflow_model.outputs.all_res[i]["h"], 
#                                    cpos = "xy", show_edges = False)
#    
    
    
##############################################################################
##############################################################################
    
    
dassflow_model.post = df2d.core.output.Post(boundary_metadata = dassflow_model.boundary.metadata, bin_dir = dassflow_model.bin_dir)




#----------------------------#
# PLOTS
#----------------------------#

# - - - - - -#
# water heigh
# - - - - - - #

best_diff = 1000 #todo
all_keys = [x for x in dassflow_model.outputs.all_res.keys()]
old_key= all_keys[10]
for i in dassflow_model.outputs.all_res.keys():
    
    previous_h = np.asanyarray(dassflow_model.outputs.all_res[old_key]["h"])
    h_array = np.asanyarray(dassflow_model.outputs.all_res[i]["h"])
    diff = previous_h[:] - h_array
    diff = np.linalg.norm(diff)
   # print(i, "diff=", diff)
    if diff<best_diff:
        print(i, "diff=", diff)
        best_diff =  diff
    dassflow_model.meshing.mesh_pyvista.plot(scalars = dassflow_model.outputs.all_res[i]["h"], 
                                    cpos = "xy", show_edges = False)
    old_key = i
    
    
# - - - - - - #
# inflow + outflow discharge
# - - - - - - #

fig = plt.figure(figsize = (12,12))
fig.suptitle("Inflow discharges (sum_q)")
ax = plt.gca()

ax.scatter( x = dassflow_model.post.sum_q_inflow["0"][:,0],
           y = dassflow_model.post.sum_q_inflow["0"][:,1], marker = "o", 
           label = "group 1", c = "red")
ax.scatter( x = dassflow_model.post.sum_q_inflow["1"][:,0],
           y = dassflow_model.post.sum_q_inflow["1"][:,1], marker = "+", 
           label = "group 2" , c = "blue")
ax.scatter( x = dassflow_model.post.sum_q_inflow["2"][:,0],
           y = dassflow_model.post.sum_q_inflow["2"][:,1], marker = "x", 
           label = "group 3", c = "green")
ax.legend()



fig = plt.figure()
fig.suptitle("outflow discharges")
ax = plt.gca()

ax.scatter( x = dassflow_model.post.sum_mass_flux_outflow["0"][:,0],
           y = dassflow_model.post.sum_mass_flux_outflow["0"][:,1], marker = "o", 
           label = "group 1", c = "red")
ax.set_yscale('log')
ax.legend()




# - - - - - - #
# Inflow
# - - - - - - #

all_id_inflow = np.asanyarray([0,1,2])
for group in all_id_inflow+1:
    #group = 3
    post_group = str(group -1)
    
    fig = plt.figure(figsize = (20,12))
    fig.suptitle(f"Inflow discharge {group}")
    ax = plt.gca()
    
    ax.plot( dassflow_model.post.sum_mass_flux_inflow[post_group][:,0],
                dassflow_model.post.sum_mass_flux_inflow[post_group][:,1], 
               label = "sum_mass_flux_inflow", c = "red", linestyle = "-")
    ax.plot( dassflow_model.post.sum_q_inflow[post_group][:,0],
                dassflow_model.post.sum_q_inflow[post_group][:,1], 
               label = "sum_q_inflow", c = "blue", linestyle = "--")
    ax.plot(np.arange(start = 0, stop = dt_hydrology * len(qin_heating[0]), step  = dt_hydrology),
                qin_heating[group-1], 
               label = "qin_hydrology", c = "black", linestyle = "dotted")
    
    id_boundary_cell = dassflow_model.boundary.corresp["discharg1"][group]["id_cell"]
    id_boundary_cell = id_boundary_cell-1 # pass to python index ?
    length_edge = dassflow_model.meshing.mesh_fortran.edge[dassflow_model.boundary.corresp["discharg1"][group]["id_edge"]].length
    
    q_inflow_sum = np.zeros(shape = (len(dassflow_model.outputs.all_res)))
    
    h = []
    u = []
    v = []
    for key, val in dassflow_model.outputs.all_res.items():
        h.append(np.asarray(val["h"][id_boundary_cell]))
        u.append(np.asarray(val["u"][id_boundary_cell]))
        v.append(np.asarray(val["v"][id_boundary_cell]))    
    h = np.asanyarray(h)
    u = np.asanyarray(u) 
    v = np.asanyarray(v)
    q = length_edge * h[:] * np.sqrt(u[:]**2 + v[:]**2)
    q_inflow_sum[:] = q_inflow_sum[:] + q[:,0] #IMPORTANT: built q_inflow_sum during the plot
    ax.plot(np.arange(start = 0, stop = dassflow_dtw_heating * len(q), step  = dassflow_dtw_heating),
                q, 
               label = "q for DOF variable (python)", c = "purple")
    
    sum1=sum( dassflow_model.post.sum_mass_flux_inflow[post_group][:,1] ) / len( dassflow_model.post.sum_mass_flux_inflow[post_group][:,1])
    sum2=sum( dassflow_model.post.sum_q_inflow[post_group][:,1]         ) / len(dassflow_model.post.sum_q_inflow[post_group][:,1]         )
    sum3=sum( new_qin[group-1] ) / len(new_qin[group-1] )
    sum4=sum( q ) / len(q)
    ax.text(0.5, 0.9 * max(q), s= f" sums/nb elements \n sum sum_mass_flux_inflow = {sum1}  \n sum sum_q_inflow  = {sum2} \n sum qin_hydrology =  {sum3} \n sum q dof =  {sum4}")
    
    ax.legend()



# - - - - - - #
# Outflow
# - - - - - - #
group = all_id_inflow[-1]+2
post_group = str(0)

fig = plt.figure(figsize = (12,12))
fig.suptitle(f"Outflow discharge {group}")
ax = plt.gca()

ax.plot( dassflow_model.post.sum_mass_flux_outflow[post_group][:,0],
            dassflow_model.post.sum_mass_flux_outflow[post_group][:,1], 
           label = "sum_mass_flux_outflow", c = "red", linestyle = "-")
ax.plot( dassflow_model.post.sum_q_outflow[post_group][:,0],
            dassflow_model.post.sum_q_outflow[post_group][:,1], 
           label = "sum_q_outflow", c = "blue", linestyle = "--")
#ax.plot(np.arange(start = 0, stop = dt_hydrology * len(new_qin[0]), step  = dt_hydrology),
#            new_qin[group-1], 
#           label = "qin_hydrology", c = "black", linestyle = "dotted")

id_boundary_cell = dassflow_model.boundary.corresp["transm"][group]["id_cell"]
id_boundary_cell = id_boundary_cell-1 # pass to python index ?
length_edge = dassflow_model.meshing.mesh_fortran.edge[dassflow_model.boundary.corresp["transm"][group]["id_edge"]].length


h = []
u = []
v = []
for key, val in dassflow_model.outputs.all_res.items():
    h.append(np.asarray(val["h"][id_boundary_cell]))
    u.append(np.asarray(val["u"][id_boundary_cell]))
    v.append(np.asarray(val["v"][id_boundary_cell]))    
h = np.asanyarray(h)
u = np.asanyarray(u)
v = np.asanyarray(v)
q = length_edge * h[:] * np.sqrt(u[:]**2 + v[:]**2)
q_outflow = q
ax.plot(np.arange(start = 0, stop = dassflow_dtw_heating * len(q_outflow), step  = dassflow_dtw_heating),
            q_outflow, 
           label = "q for DOF variable (python)", c = "purple")
ax.text(0.5, 0.9 * max(q_outflow), s= f"sum sum_mass_flux_outflow =        {sum(dassflow_model.post.sum_mass_flux_outflow[post_group][:,1])} \n sum sum_mass_flux_outflow  = {sum(dassflow_model.post.sum_q_outflow[post_group][:,1])} \n sum q DOF  =  {sum(q_outflow)}")

ax.legend()



# - - - - - - #
# Compare inflow - outflow dassflow
# - - - - - - #

qin = dassflow_model.post.sum_mass_flux_inflow["0"][:,1]+ dassflow_model.post.sum_mass_flux_inflow["1"][:,1]+ dassflow_model.post.sum_mass_flux_inflow["2"][:,1]
qout = dassflow_model.post.sum_mass_flux_outflow["0"][:,1]

fig = plt.figure(figsize = (12,12))
fig.suptitle(f"sum_mass_flux :  Inflow & outflow")
ax = plt.gca()
ax.plot(qin, label = "sum inflows")
ax.plot(qout, label = "outflow")
ax.legend()

qin = dassflow_model.post.sum_q_inflow["0"][:,1]+ dassflow_model.post.sum_q_inflow["1"][:,1]+ dassflow_model.post.sum_q_inflow["2"][:,1]
qout = dassflow_model.post.sum_q_outflow["0"][:,1]

fig = plt.figure(figsize = (12,12))
fig.suptitle(f"sum_q_flux :  Inflow & outflow")
ax = plt.gca()
ax.plot(qin, label = "sum inflows")
ax.plot(qout, label = "outflow")
ax.text(0.5, max(qin), s= f" sum qin = {sum(qin)} \n sum qout =  {sum(qout)}")
ax.legend()



fig = plt.figure(figsize = (12,12))
fig.suptitle(f"DOF Q :  Inflow & outflow")
ax = plt.gca()
ax.plot(q_inflow_sum, label = "q inflows")
ax.plot(q_outflow, label = "q outflow")
ax.text(0.5, max(qin), s= f" sum qin = {sum(q_inflow_sum)} \n sum qout =  {sum(qout)}")
ax.legend()



#--------------------     #
# water volume evolution
#------------------------ #
time_post = dassflow_model.post.water_vol[:,0]
water_volume = dassflow_model.post.water_vol[:,1]

time_res = np.arange(start = 0, stop = dassflow_dtw_heating * len(q), step  = dassflow_dtw_heating)

fig = plt.figure(figsize = (12,12))
fig.suptitle(f"water volume from post file 'water_vol.dat' ")
ax = plt.gca()
ax.plot(time_post,water_volume, label = "post.water_vol")
ax.text(0.5, 0.9* max(water_volume), s= "sum water volume ")
ax.legend()


fig = plt.figure(figsize = (12,12))
fig.suptitle(f"water volume diff (red curve) compared to DOF inflow & outflow")
ax1 = plt.subplot()
l1, = ax1.plot(time_post[1:],water_volume[1:]-water_volume[:-1], label = "post.water_vol diff", color='red')
ax2 = ax1.twinx()
l2, = ax2.plot(time_res,q_inflow_sum, label = "q inflows", color='blue')
l3, = ax2.plot(time_res,q_outflow, label = "q outflow", color='green')
ax1.legend(loc=0);ax2.legend(loc=2)
ax1.set_ylabel("Water volume (m3?)");ax2.set_ylabel("dicharge [m3/s]")



qdiff = q_inflow_sum[:]-q_outflow[:].squeeze()[:]

all_dof_volume = []

for my_key in   dassflow_model.outputs.all_res.keys() :
    dof_volume=[]
    for id_cell in range(dassflow_model.meshing.mesh_fortran.nc):
        h = dassflow_model.outputs.all_res[my_key]["h"][id_cell]
        surface =  dassflow_model.meshing.mesh_fortran.cell[id_cell].surf
        volume = h * surface
        dof_volume.append(volume)
    all_dof_volume.append(dof_volume)
    
sum_dof_volume = [np.sum(x) for x in all_dof_volume]
sum_dof_volume= np.asanyarray(sum_dof_volume)

fig = plt.figure(figsize = (12,12))
fig.suptitle(f"Delta water volume at each timestep \n SOURCE [ blue:DOF-fluxes, black:DOF-volumes, red:post/water_vol ]")
ax1 = plt.subplot()
l1, = ax1.plot(time_post[1:],water_volume[1:]-water_volume[:-1], label = "post.water_vol diff", color='red')
ax2 = ax1.twinx()
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("axes", 1.1)) 
l2, = ax2.plot(time_res,qdiff, label = "q inflows-outflow", color='blue')
l3, = ax3.plot(time_res[1:],sum_dof_volume[1:]-sum_dof_volume[:-1],
               label = "diff dof volume", color='black', marker='o', linestyle='None')


qdiff = q_inflow_sum[:]-q_outflow[:].squeeze()[:]

all_dof_volume = []

for my_key in   dassflow_model.outputs.all_res.keys() :
    dof_volume=[]
    for id_cell in range(dassflow_model.meshing.mesh_fortran.nc):
        h = dassflow_model.outputs.all_res[my_key]["h"][id_cell]
        surface =  dassflow_model.meshing.mesh_fortran.cell[id_cell].surf
        volume = h * surface
        dof_volume.append(volume)
    all_dof_volume.append(dof_volume)
    
sum_dof_volume = [np.sum(x) for x in all_dof_volume]
sum_dof_volume= np.asanyarray(sum_dof_volume)

fig = plt.figure(figsize = (12,12))
fig.suptitle(f"Delta water volume at each timestep \n SOURCE [ blue:DOF-fluxes, black:DOF-volumes, red:post/water_vol ]")
ax1 = plt.subplot()
l1, = ax1.plot(time_post[1:],water_volume[1:]-water_volume[:-1], label = "post.water_vol diff", color='red')
ax2 = ax1.twinx()
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("axes", 1.1)) 
l2, = ax2.plot(time_res,qdiff, label = "q inflows-outflow", color='blue')
l3, = ax3.plot(time_res[1:],sum_dof_volume[1:]-sum_dof_volume[:-1], 
               label = "diff dof volume", color='black', marker='o', linestyle='None')
ax1.legend(loc=0);ax2.legend(loc=2);ax3.legend(loc=3)
ax1.set_ylabel("Water volume (m3?)");ax2.set_ylabel("dicharge [m3/s]");
ax3.set_ylabel("DOF water volume")

ax1.legend(loc=0);ax2.legend(loc=2);ax3.legend(loc=3)
ax1.set_ylabel("Water volume (m3?)");ax2.set_ylabel("dicharge [m3/s]");
ax3.set_ylabel("DOF water volume")


#plt.plot(dassflow_model.post.water_vol_num_add[:,0],dassflow_model.post.water_vol_num_add[:,1], label = "water vol num add")



    