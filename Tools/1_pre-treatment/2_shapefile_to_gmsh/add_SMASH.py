#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 09:17:50 2022

@author: livillenave
"""

import smash
import matplotlib.pyplot as plt
import multiprocessing
import os
from matplotlib.colors import LogNorm, SymLogNorm
import numpy as np


import pyvista as pv
import geopandas as gpd


def get_hydrological_coupling(mesh, 
                              treshold_drained_area):
    # id_cell, tupple of (x,y) if it was a matrix and not a numpy ndarray
    
    def get_inflowing_id(id_cell, flow_dir):
                    
            direction = flow_dir[id_cell[0], id_cell[1]].copy()
            direction = int(direction)
            
            if direction ==1:
                res =  (id_cell[0]-1,id_cell[1])
            elif direction ==2:
                res =  (id_cell[0]-1,id_cell[1]+1)
            elif direction ==3:
                res =  (id_cell[0],id_cell[1]+1)
            elif direction ==4:
                res =  (id_cell[0]+1,id_cell[1]+1)
            elif direction ==5:
                res =  (id_cell[0]+1,id_cell[1])
            elif direction ==6:
                res =  (id_cell[0]+1,id_cell[1]-1)
            elif direction ==7:
                res =  (id_cell[0],id_cell[1]-1)
            elif direction ==8:
                res =  (id_cell[0]-1,id_cell[1]-1)
            else:
                print("direction =", direction)
                print("flow directions must be between 1 and 8 (SMASH convention)")
                return
            return(res)
            
    drained_area =  mesh["drained_area"]
    flow_dir = mesh["flwdir"]
    
    # >>> define treshold 
    
    # tmp = drained_area.ravel().copy()
    # tmp.sort()
    # treshold_drained_area = tmp[int(len(tmp) * 0.93)]
    
    
    # >>> identify index of "hydraulic cells"
    
    id_hydraulic = np.where(drained_area > treshold_drained_area)
    
    new = []
    for i in range(len(id_hydraulic[0])) :
                        # row            "col
        new.append( (id_hydraulic[0][i],id_hydraulic[1][i]) )    
    id_hydraulic = new
    
    # -------------------- PLOT
    tmp = mesh["drained_area"].copy()
    tmp[:,:] = 0
    for my_id in id_hydraulic:
        tmp[my_id[0],my_id[1]] = 1
    plt.imshow(tmp)
    # -------------------- OUT PLOT

    
    all_inflows=dict()
    all_adj=[]
    nb_inflow=0
    # >>> identify   inflowing cells in hydraulic domain
        # --- all first
    for i in range(len(id_hydraulic)) :
        idrow_max =  flow_dir.shape[0]-1
        idcol_max =   flow_dir.shape[1]-1
        
        idrow = id_hydraulic[i][0]      
        idcol = id_hydraulic[i][1]  
                
        adjacent = []
        # get id adjacent cells
        for row_modif in [-1,0,1]:  # i for row index modif
             for col_modif in [-1,0,1]: # j for column index modif
                if not (row_modif == 0 and col_modif == 0): # if not the cell itself
                 print("idrow + row_modif,  idcol  +col_modif = ", idrow + row_modif,  idcol  +col_modif)
                 print("not in id_hydraulic ??? ", ((idrow + row_modif,  idcol  +col_modif) not in id_hydraulic) )
                 if (idrow + row_modif,  idcol  +col_modif) not in id_hydraulic:  # if the adjacent cell is not of hydraulic domain
                     check1 = ( idrow+row_modif<= idrow_max)
                     check2=(idrow+row_modif>=0)
                     check3 = (idcol+col_modif<= idcol_max)
                     check4 =( idcol+col_modif>=0)
                     check =   np.all( [check1,check2,check3,check4])
                     print("is in domain ??? ", ((idrow + row_modif,  idcol  +col_modif) not in id_hydraulic) )                 
                     if check :  # if the adjacent cell is in domain
                         my_id_cell = (idrow + row_modif,  idcol  +col_modif)
                         if not isinstance(flow_dir[my_id_cell[0], my_id_cell[1]], np.ma.core.MaskedConstant) :
                             if get_inflowing_id(id_cell = my_id_cell , flow_dir = flow_dir) ==  (idrow,idcol):     # if adjacent cell inflow the target cell
                                 adjacent.append( (idrow + row_modif,  idcol  +col_modif)  )
                    
        for k in range(len(adjacent)):
            all_inflows[nb_inflow] = dict()
            all_inflows[nb_inflow]["id"]= adjacent[k]
            all_inflows[nb_inflow]["inflowed"]= (idrow,idcol)
            all_inflows[nb_inflow]["self_drained_area"]= drained_area[adjacent[k][0], adjacent[k][1]]
            all_inflows[nb_inflow]["inflowed_drained_area"]= drained_area[idrow, idcol]
          
            nb_inflow=nb_inflow+1
            
        all_adj.append(adjacent)
    categorical = drained_area.copy()    
    categorical[:,:]    = 0
                    
    for i in range(len(all_inflows)):
        id_inflow = all_inflows[i]["id"]                
        categorical[id_inflow[0],id_inflow[1]] =2
                
    for i in range(len(id_hydraulic)):
        id_inflow = id_hydraulic[i]               
        categorical[id_inflow[0],id_inflow[1]] =1
                
#    plt.imshow(categorical)
#    plt.colorbar()
#    plt.show()    
        
    mesh["categorical"] = categorical
    
    return(mesh, all_inflows)
   




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 14:43:47 2022

@author: livillenave
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 14:40:36 2022

@author: livillenave
"""



# ===================================================================== #
# SCRIPT
# ===================================================================== #

#---------- Y1422030  
# x_gauge = 654_521 
# y_gauge = 6_234_639
# area =    3_035_579_999
#code_gauge = "Y1422030"

#---------- Y1422020  
x_gauge = 662_594 
y_gauge = 6_233_537
area =    3_192_579_999
code_gauge = "Y1422030"

# >>> Build MESHING + find Outlet find  outlet for SMASH MODEL
meshing = smash.mesh.meshing.generate_mesh(
        path =f'/home/livillenave/Documents/data/DONNEES/LEBLOIS_DATA/10/flow_dir.asc',#f'/home/livillenave/Documents/data/DONNEES/LEBLOIS_DATA/10/flow_dir.asc',#f'{coupled_dir}/{my_scale}/SMASH_files/FLOW.asc', 
        epsg=2154,
        x=x_gauge,
        y=y_gauge,
        area= area, #3_192_100_000, 
        #bbox = (601_435, 663_000 ,6_159_000,6_258_977 ),
        code =code_gauge)
       

new_mesh, all_inflow = get_hydrological_coupling(mesh=meshing , treshold_drained_area=1500)



# =========================================== #
# VISUALISATION 1
# =========================================== #

x = new_mesh["active_cell"].copy()
y = new_mesh["active_cell"].copy()
z = new_mesh["active_cell"].copy()

xmin = new_mesh["xmin"] 
ymax = new_mesh["ymax"] 
xres = new_mesh["dx"]

compteur = 0
for i in range(x.shape[0]): #rows
    for j in range(x.shape[1]):       # cols 
        compteur = compteur+1
        y[i,j]= ymax -  i * xres     # rows
        x[i,j]= xmin +  j * xres     # cols
        z[i,j]= new_mesh["categorical"][i,j]  
#
test = pv.StructuredGrid(x, y, z)
#test.plot(cpos = "xy", show_edges=True, notebook = True, show_bounds=True)
points = test.points

raw_major_bed = gpd.read_file( filename = "/home/livillenave/Documents/distant/SD-FLOOD/real_case-AUDE/DATA/DASSFLOW/V2/hydraulic/tmp/contour1_v2.shp" )  
raw_minor_bed = gpd.read_file( filename = "/home/livillenave/Documents/distant/SD-FLOOD/real_case-AUDE/DATA/DASSFLOW/V2/hydraulic/tmp/contour2.shp" )   
cropped_riverline = gpd.read_file( filename = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/mesh_generation/river_network.shp" )    
raw_riverline =  gpd.read_file("/home/livillenave/Documents/data/DONNEES/spatial/COURS_D_EAU_NATUREL_mono_oriente_1a5.shp", 
                            bbox = cropped_riverline)#(np.min(x[:,:]), np.max(x[:,:]), np.min(y[:,:]), np.max(y[:,:])))


fig,ax = plt.subplots(2,2, figsize=(15,8))

fig.suptitle("Source GIS FILES", fontsize=16)
raw_major_bed.plot(ax = ax[0,0])
ax[0,0].set_title("raw major bed")
raw_minor_bed.plot(ax = ax[1,0])
ax[1,0].set_title("raw Minor bed")
raw_riverline.plot(ax = ax[0,1])
ax[0,1].set_title("raw river network")
cropped_riverline.plot(ax = ax[1,1])
ax[1,1].set_title("Cropped river network")
plt.show()
plt.close()



zoom_coord = raw_major_bed.bounds
raw_major_bed.plot(figsize=(20, 20))
plt.scatter(points[:, 0],points[:, 1], c=points[:, 2], marker="o",  s=100)
plt.scatter(x_gauge,y_gauge, marker = "*", c= "r", s = 8)
plt.axis("image")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.xlim(np.asarray(zoom_coord["minx"])[0], np.asarray(zoom_coord["maxx"])[0])
plt.ylim(np.asarray(zoom_coord["miny"])[0], np.asarray(zoom_coord["maxy"])[0])
plt.savefig("/home/livillenave/Images/o.png", dpi = 500)
plt.show()



fig,ax = plt.subplots(1,2, figsize=(15,8))
ax[0].imshow(new_mesh["drained_area"], extent = [min(points[:,0]), max(points[:,0]), min(points[:,1]), max(points[:,1])])
ax[1].imshow(new_mesh["categorical"], extent = [min(points[:,0]), max(points[:,0]), min(points[:,1]), max(points[:,1])])
plt.show()
plt.close()


# ========================================================== #
# GENERATE MAPPING BETWEEN dassflow_mesh and SMASH_mesh      #
# ========================================================== #


smash_mesh = meshing
import dassflow2d as df2d
# >>> Add xy coords of smash inflows
def get_xy_smash_grid(smash_mesh, id_row, id_col):
    # geometrical properties
    dcol = smash_mesh["dx"]
    drow = -smash_mesh["dx"]
    upper_left_x =  smash_mesh["xmin"]
    upper_left_y =  smash_mesh["ymax"]
    
    x = upper_left_x + dcol * id_col
    y = upper_left_y + drow * id_row
    
    return((x,y))
    
for i in range(len(all_inflow)):
    all_inflow[i]["self_xy"]= get_xy_smash_grid(smash_mesh=smash_mesh, id_row = all_inflow[i]["id"][0], id_col = all_inflow[i]["id"][1])
    all_inflow[i]["inflowed_xy"]= get_xy_smash_grid(smash_mesh=smash_mesh, id_row = all_inflow[i]["inflowed"][0], id_col = all_inflow[i]["inflowed"][1])




# >>> Source xy coords of 
df2d_mesh =  read_mesh_from_textfile(path = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/final_mesh.geo",
                               read_boundary = True)

model = df2d.dassflowmodel(bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/", 
                           hdf5_path="/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/res/simu.hdf5", 
                           run_type= "direct",  
                           clean = True)

df2d.wrapping.m_common.set_mesh_name("final_mesh.geo")
model.init_all()

#  get xy cells of dassflow
df2d_x, df2d_y = get_xy_cells(model.meshing.mesh_fortran)


id_fortran_cell_inlets  = []
xy_fortran_cell_inlets = []
# get xy of inlet cells dassflow
for key, value in df2d_mesh["boundaries"]["inlet"].items():
    index_fortran_cell = value["index_fortran_cell"]
    index_python_cell  =  index_fortran_cell-1
    id_fortran_cell_inlets.append(index_fortran_cell)
    xy_fortran_cell_inlets.append( (df2d_x[index_python_cell], df2d_y[index_fortran_cell]) )
    
    
    
# >>> Identify closest bc to main_inflow

tmp = []
for key, value in all_inflow.items():
    xy_inflow_smash = value["self_xy"]    
    tmp.append(xy_inflow_smash)
    id_closest_inlet = -1
    best_dist = 100000000000000000000
     
    #print("loop")
    for id_inlet in range(len(xy_fortran_cell_inlets)):  
   #     print("id_inlet")
        xy_inlet_dassflow =  xy_fortran_cell_inlets[id_inlet]
    #    print("xy_inlet_dassflow", xy_inlet_dassflow)
        dist = np.sqrt( (xy_inflow_smash[0] - xy_inlet_dassflow[0] )**2 + ( xy_inflow_smash[1] - xy_inlet_dassflow[1] )**2  )
        
        if dist<best_dist:
            best_dist = dist
            id_closest_inlet = id_inlet
    
    all_inflow[key]["bc_group_dassflow"] = id_closest_inlet
    
    print(best_dist)    
    print(id_inlet)
    
#
x1 = [x[0] for x in xy_fortran_cell_inlets]
y1 = [x[1] for x in xy_fortran_cell_inlets]
#
x2 = [x[0] for x in tmp]
y2 = [x[1] for x in tmp]
colorvalue = [ all_inflow[x]["bc_group_dassflow"] for x in all_inflow.keys()]

plt.figure(figsize = (20,20))
plt.scatter(x1, y1, marker = '>', c = np.arange(len(x1)));
plt.scatter(x2, y2, c = colorvalue)
for i in range(len(xy_fortran_cell_inlets)):
    plt.annotate(f"{i+1}_dassflow", xy_fortran_cell_inlets[i])
for i in range(len(tmp)):
    plt.annotate(i+1, tmp[i])
#plt.savefig("/home/livillenave/Images/o2.png", dpi = 500)

plt.show()
plt.close()

    


all_inflow
plt.imshow(meshing["categorical"])



