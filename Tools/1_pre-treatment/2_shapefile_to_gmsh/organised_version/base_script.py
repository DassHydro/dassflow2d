#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 09:54:00 2022

@author: livillenave
"""
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
from gmsh_libs import *
from smash_coupling_libs import *
from dassflow_mesh_manipulation_libs import *
import numpy as np
import shutil


#####################
# parameters
#####################

# Dassflow mesh manipulation
demo_dir = "/home/livillenave/Documents/distant/Demo"
    
# default bin dir, containing necessary file for dassflow to compile
default_bin_dir  = f"{demo_dir}/default_bin" 
# were to invoke dassflow2d    
bin_dir = f"{demo_dir}/bin7"
# source file location
source_dir =  f"{demo_dir}/source_file/"
# were most of intermediate files (mostly GIS are written)
mesh_generation_dir = f"{bin_dir}/mesh_generation"
# DEL dirctory
mnt_dir =  f"{source_dir}/DEM"
# raw dem to source
source_mnt_name = "aude_mean.tif"
# dem cropped  on extent  major bed defined below 
new_mnt_name    = "aude_cropped.tif"

path_raw_river_network = f"{source_dir}/river_network/COURS_D_EAU_NATUREL_mono_oriente_1a5.shp"# path to shapefiles
path_major_bed = f"{source_dir}/hydraulic_shapefile/major_bed_aude.shp"
path_minor_bed = f"{source_dir}/hydraulic_shapefile/minor_bed.shp"
#path_major_bed = f"{source_dir}/hydraulic_shapefile/lit_majeur_opti.shp"
#path_minor_bed = "/home/livillenave/Documents/distant/Demo/source_file/tosort/minor_bed_2.shp"#f"{source_dir}/hydraulic_shapefile/lit_mineur_opti.shp"
# path flow diretion
path_leblois_data = f"{source_dir}/DEM/LEBLOIS_DATA/10/flow_dir.asc"

# name mesh
mesh_name = "aude"

# mesh size in meters
mesh_size_major = 100 # major bed
mesh_size_minor = 10  # minor bed

# set False  to disable minor bed
include_minor_bed =  True


# ================================================================ #
#   Script 
# 
# 0 - generate empty directories
# 1 - Build Gmsh mesh
# 2- traduce it to dassflow_mesh.geo without boundary conditions
# 3- Add boundary condition to mesh
# 4- Add bathymetry to mesh
# 5- generate SMASH meshing () 
# 6- Associate  (map) dassflow  and smash mesh
# ================================================================ #

#----------------------------------------------------#
# 0- prepare directories and default files
#----------------------------------------------------# 

# >>> Directories

if not os.path.exists(f'{bin_dir}'):
    os.mkdir(f'{bin_dir}')
else:
    shutil.rmtree(f'{bin_dir}')
    os.mkdir(f'{bin_dir}')
if not os.path.exists(f'{mesh_generation_dir}'):
    os.mkdir(f'{mesh_generation_dir}')
os.chdir(f'{mesh_generation_dir}')

# >>> Files
for my_filename in os.listdir(default_bin_dir):
    shutil.copyfile(src=f"{default_bin_dir}/{my_filename}", dst=f"{bin_dir}/{my_filename}")

#----------------------------------------------------#
# 1- Build Gmsh mesh
#----------------------------------------------------# 
nodes, cells = build_gmsh_mesh(path_major_bed = path_major_bed, 
                                mesh_size_major = mesh_size_major, 
                                include_minor_bed = include_minor_bed, 
                                path_minor_bed = path_minor_bed,  
                                mesh_size_minor = mesh_size_minor,
                                write_dir = mesh_generation_dir,
                                mesh_name = mesh_name)

#----------------------------------------------------# 
# 2- Traduce gmsh to  dassflow mesh
#----------------------------------------------------#

gen_empty_mesh_dassflow(nodes = nodes, cells = cells, path = f"{mesh_generation_dir}/raw_{mesh_name}.geo")
shutil.copyfile(src=f"{mesh_generation_dir}/raw_{mesh_name}.geo", dst=f"{bin_dir}/dassflow_mesh_default.geo")
#gen_empty_mesh_dassflow(nodes = nodes, cells = cells, path = f"{bin_dir}/dassflow_mesh_default.geo") # dassflow_mesh_default is basic sourced name by dassflow2d


#----------------------------------------------------#
#  3-  Add boundary condition to mesh
# note that we will only use dassflow meshing 
#----------------------------------------------------#

# initialise fortran instance, and python corrponding data
model = df2d.dassflowmodel(bin_dir = bin_dir, 
                           hdf5_path=bin_dir+"/res/simu.hdf5", 
                           run_type= "direct",  
                           clean = True)
df2d.wrapping.m_common.set_mesh_name("dassflow_mesh_default.geo")
model.init_all()
    
#  get + write contour hy shapefile
multi_line = get_countour_hy(model = model, save_shapefile= True, write_dir = mesh_generation_dir, file_name = "contour_hy.shp")


# get_river_network()
# source river network
# bbox allow to crop source shapefile
river_network = gpd.read_file(filename = f"{path_raw_river_network}",
                              bbox = multi_line,
                              crs = 2154)
river_network.to_file(f'{mesh_generation_dir}/river_network.shp')

##############

# get values
(x_intersect, y_intersect) = get_xy_intersect_shapely_lines(river_network=river_network, contour_hy=multi_line)

# save
df = pd.DataFrame({'x': x_intersect, 'y': y_intersect})
gdf = gpd.GeoDataFrame(
    df, geometry=gpd.points_from_xy(df['x'], df['y']))
gdf.to_file(f'{mesh_generation_dir}/intersect_point.shp')


index_wall_edges =  get_index_wall_edges(model = model)
# get all id of edges intersecting hydrographic network
all_edge_intersect = []
for i in range(len(x_intersect)): # loop on each intersection point
    # extract xy for this intersection point:
        x = x_intersect[i]
        y = y_intersect[i]
        best_edge = 1
        best_dist = 10**10
        
        best_edge = get_index_closest_edge(mesh = model.meshing.mesh_fortran, x_ref = x, y_ref = y, 
                               edge_index_to_test = index_wall_edges)
        all_edge_intersect.append(best_edge)


boundaries_to_add = dict()
for i in range(len(all_edge_intersect)):
    my_edge = all_edge_intersect[i]
    boundaries_to_add[my_edge] = get_boundary_metadata_for_textfile_mesh(mesh = model.meshing.mesh_fortran,
                                        index_edge_intersect = my_edge)
    
# outlet gauge coordinate
xy_gauge = (652874.12,6235496.94)

closest_id = -1
best_dist = 10000000000000000
for key, value in boundaries_to_add.items():
    x,y = get_xy_center_edge(mesh = model.meshing.mesh_fortran, index_wall_edges=key)
    dist = np.sqrt( (x-xy_gauge[0])**2 + (y-xy_gauge[1])**2 )
    
    if dist<=best_dist:
        best_dist = dist
        closest_id = key

for key, value in boundaries_to_add.items():
    if int(key) == int(closest_id):
        boundaries_to_add[key]["bc_type"] =  "transm"
    else:
        boundaries_to_add[key]["bc_type"] =  "discharg1"
        
        
mesh = read_mesh_from_textfile(f'{mesh_generation_dir}/raw_{mesh_name}.geo')
# add boundaries to mesh file
mesh["boundaries"] = dict()
mesh["boundaries"]["header"] = dict()
mesh["boundaries"]["inlet"] = dict()
mesh["boundaries"]["outlet"]= dict()

compteur_inflow = 0
compteur_outflow=0

for key, value in boundaries_to_add.items():
    if value["bc_type"] == "discharg1":
        compteur_inflow = compteur_inflow+1
        mesh["boundaries"]["inlet"][key]= {"index_fortran_cell":value["index_cell"],                     # id cell
         "index_relative_edge":value["edge_index_in_cell_structure"]+1,  # relative edge
         "ghost_bathy":value["edge_index_in_cell_structure"],  # ghost cell bathy
         "group":compteur_inflow }
        
    if value["bc_type"] == "transm":
        compteur_outflow=compteur_outflow+1
        mesh["boundaries"]["outlet"][key]= {"index_fortran_cell":value["index_cell"],                     # id cell
         "index_relative_edge":value["edge_index_in_cell_structure"]+1,  # relative edge
         "ghost_bathy":value["edge_index_in_cell_structure"],  # ghost cell bathy
         "group":compteur_outflow }


mesh["boundaries"]["header"]["nb_inflow"] = len(mesh["boundaries"]["inlet"]) 
mesh["boundaries"]["header"]["nb_cell_inflow"] = len(mesh["boundaries"]["inlet"])   # to update  for multi-cell bc
mesh["boundaries"]["header"]["nb_outflow"] = len(mesh["boundaries"]["outlet"]) 
mesh["boundaries"]["header"]["nb_cell_outflow"] = len(mesh["boundaries"]["outlet"])    # to update for multi-cell bc


gen_mesh_dassflow(mesh_dict = mesh, path = f"{mesh_generation_dir}/bc_{mesh_name}.geo")
df2d.wrapping.call_model.clean_model(model.kernel)
del model


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

model = df2d.dassflowmodel(bin_dir = bin_dir, 
                           hdf5_path=bin_dir+"/res/simu.hdf5", 
                           run_type= "direct",  
                           clean = True)
df2d.wrapping.m_common.set_mesh_name("dassflow_mesh_default.geo")
model.init_all()
    
    


# -----------------------#
# PREPARE MNT 
# ----------------------- # 
# >>> crop MNT on area
xcell, ycell = get_xy_cells(model.meshing.mesh_fortran)
ulx=min(xcell)-1000
uly=max(ycell)+1000
lrx=max(xcell)+1000
lry=min(ycell)-1000
os.chdir(mnt_dir)
os.system(f'gdal_translate -projwin {ulx} {uly} {lrx} {lry}  -of GTiff -co "TILED=YES" {source_mnt_name} {new_mnt_name}')
os.chdir(bin_dir)


# -----------------------#
# SOURCE MNT + store xy cell coordinates
# ----------------------- # 
# >>> LOAD CROPPED MNT
ds = gdal.Open(f"{mnt_dir}/{new_mnt_name}")
ulx, xres, _, uly, _, yres  = ds.GetGeoTransform()

upper_left_x = ulx
upper_left_y = uly
x_resolution = xres
y_resolution = yres


# >>> Source value
data_array = ds.GetRasterBand(1).ReadAsArray()

x_ndarray = data_array.copy()
y_ndarray = data_array.copy()

# >>> calc zcell
for row_id in range(x_ndarray.shape[0]):
    for col_id in range(x_ndarray.shape[1]):
        # x change à chaque colomne
        # y change à chaque ligne
        x_ndarray[row_id, col_id]= upper_left_x + x_resolution * col_id
        y_ndarray[row_id, col_id]= upper_left_y + y_resolution * row_id


plt.imshow(x_ndarray); plt.colorbar();plt.show();plt.close()
plt.imshow(y_ndarray); plt.colorbar();plt.show();plt.close()



# -----------------------#
# ATTRIBUTE Z value to each cell
# ----------------------- # 

rowcol_corresp = []
for i in range(len(xcell)): #
    x = xcell[i]
    y = ycell[i]
    
    best_row = -1
    best_col = -1
    best_dist = 1000000000000
    
    raw_best_col = np.where(abs(x_ndarray[0,:] - x) == min(abs(x_ndarray[0,:] - x)))[0][0]
    raw_best_row = np.where(abs(y_ndarray[:,0] - y) == min(abs(y_ndarray[:,0] - y)))[0][0]
    
    candidates_rows = np.arange(raw_best_row-5,raw_best_row+5, dtype = "int")
    candidates_cols = np.arange(raw_best_col-5,raw_best_col+5, dtype = "int")
    
    # clean and remove index outside raster (too big because of the +5)
    to_keep_rows = np.where(candidates_rows < x_ndarray.shape[0])[0]
    to_keep_cols = np.where(candidates_cols < x_ndarray.shape[1])[0]
    
    candidates_rows = candidates_rows[to_keep_rows]
    candidates_cols = candidates_cols[to_keep_cols]
    
    
    for row_id in candidates_rows:
            for col_id in candidates_cols:
            
                dist = np.sqrt((x- x_ndarray[row_id, col_id])**2 + (y- y_ndarray[row_id, col_id])**2)
                
                if dist < best_dist:
                    best_dist=dist
                    best_row = row_id
                    best_col = col_id
    rowcol_corresp.append( (best_row, best_col) )




zcell = []
for rowcol in rowcol_corresp:
    zcell.append(data_array[rowcol[0],rowcol[1]])
    
zcell = np.asanyarray(zcell)


# >>> UPDATE MESH

my_mesh = read_mesh_from_textfile(path = f"{mesh_generation_dir}/bc_{mesh_name}.geo",
                               read_boundary = True)
my_mesh["cell"]["bathy"] = zcell
gen_mesh_dassflow(mesh_dict = my_mesh, path = f"{mesh_generation_dir}/final_{mesh_name}.geo",)

gen_mesh_dassflow(mesh_dict = my_mesh, path = f"{bin_dir}/final_mesh.geo",)
plt.scatter(xcell, ycell, c = zcell, cmap = "magma")
plt.colorbar()

df2d.wrapping.call_model.clean_model(model.kernel)
del model

model = df2d.dassflowmodel(bin_dir = bin_dir, 
                           hdf5_path=bin_dir+"/res/simu.hdf5", 
                           run_type= "direct",  
                           clean = True)

df2d.wrapping.m_common.set_mesh_name("final_mesh.geo")

model.init_all()

#model.boundary.plot(what = "meshing")




#---------- Y1422020  
x_gauge = 662_594 
y_gauge = 6_233_537
area =    3_192_579_999
code_gauge = "Y1422030"
treshold_coupling= 1500

# >>> Build MESHING + find Outlet find  outlet for SMASH MODEL
meshing = smash.mesh.meshing.generate_mesh(
        path =f'{path_leblois_data}',#f'/home/livillenave/Documents/data/DONNEES/LEBLOIS_DATA/10/flow_dir.asc',#f'{coupled_dir}/{my_scale}/SMASH_files/FLOW.asc', 
        epsg=2154,
        x=x_gauge,
        y=y_gauge,
        area= area, #3_192_100_000, 
        #bbox = (601_435, 663_000 ,6_159_000,6_258_977 ),
        code =code_gauge)
       

new_mesh, all_inflow = get_hydrological_coupling(mesh=meshing , treshold_drained_area=treshold_coupling)



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

raw_major_bed = gpd.read_file( filename = f"{path_major_bed}" )  
raw_minor_bed = gpd.read_file( filename = f"{path_minor_bed}" )   
cropped_riverline = gpd.read_file( filename = f"{mesh_generation_dir}/river_network.shp" )    
raw_riverline =  gpd.read_file(f"{source_dir}/river_network/COURS_D_EAU_NATUREL_mono_oriente_1a5.shp" , 
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
ax[0].imshow(new_mesh["flwacc"], extent = [min(points[:,0]), max(points[:,0]), min(points[:,1]), max(points[:,1])])
ax[1].imshow(new_mesh["categorical"], extent = [min(points[:,0]), max(points[:,0]), min(points[:,1]), max(points[:,1])])
plt.show()
plt.close()

df2d.wrapping.call_model.clean_model(model.kernel)
del model
# ========================================================== #
# GENERATE MAPPING BETWEEN dassflow_mesh and SMASH_mesh      #
# ========================================================== #

smash_mesh = meshing
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
df2d_mesh =  read_mesh_from_textfile(path = f"{bin_dir}/final_mesh.geo",
                               read_boundary = True)

dassflow_model = df2d.dassflowmodel(bin_dir =  bin_dir, 
                           hdf5_path= f"{bin_dir}/res/simu.hdf5", 
                           run_type= "direct",  
                           clean = True)

df2d.wrapping.m_common.set_mesh_name("final_mesh.geo")
dassflow_model.init_all()

#  get xy cells of dassflow
df2d_x, df2d_y = get_xy_cells(dassflow_model.meshing.mesh_fortran)


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
    
#    print(best_dist)    
#    print(id_inlet)
    
# PIcKLE SAVE
import pickle

a = all_inflow
with open(f'{bin_dir}/filename.pickle', 'wb') as handle:
    pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open(f'{bin_dir}/filename.pickle', 'rb') as handle:
    b = pickle.load(handle)

print(a == b)

#-------------
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

plt.imshow(meshing["categorical"])


import matplotlib.patches as mpatches
categorical_as_an_array = np.unique(meshing["categorical"].data)

## first you need to define your color map and value name as a dic
t = 1 ## alpha value
cmap = {0:[1.0,0.5,0.1,t],1:[0.1,0.1,1.0,t],2:[1.0,0.1,0.1,t]}
labels = {0:'None', 1:'river network',2:'main inflows'}
arrayShow = np.array([[cmap[i] for i in j] for j in meshing["categorical"]])    
## create patches as legend
patches =[mpatches.Patch(color=cmap[i],label=labels[i]) for i in cmap.keys()]

plt.figure(figsize=(20,20))
plt.title("Determination of SMASH outlets inflowing Dassflow (red cells)  ")
plt.imshow(arrayShow)
plt.legend(handles=patches, loc=4, borderaxespad=0.)
plt.show()



df2d.wrapping.call_model.clean_model(dassflow_model.kernel)
del dassflow_model

#################################
#################################
# perform smash simulation and extract inflows
################################"

x_gauge = 662_594 
y_gauge = 6_233_537
area =    3_192_579_999
code_gauge = "Y1422030"

calib_options={
    'structure':'gr-a',
    'dt':3600,
		'start_time':        '2018-10-14 12:00',
		'end_time':        '2018-10-20 20:00',
    
    'read_qobs':True,
    'qobs_directory':'/home/livillenave/Documents/data/FORCING/DEBIT/60',
    
    'read_prcp':True,
    'prcp_directory':'/home/livillenave/Documents/data/FORCING/PLUIE/J+1/1H',
    'prcp_conversion_factor':0.1,
    
    
    'read_pet':True,
    'pet_conversion_factor':1,
    'pet_directory':'/home/livillenave/Documents/data/FORCING/ETP-SFR-FRA-INTERA_L93',
    
    'daily_interannual_pet':True,
    'save_qsim_domain':True,
    "save_production_level":'TRUE',
    "save_others_variables":'TRUE',
    "save_net_rainfall":'TRUE' }


meshing = smash.mesh.meshing.generate_mesh(
        path =f'{path_leblois_data}',#f'/home/livillenave/Documents/data/DONNEES/LEBLOIS_DATA/10/flow_dir.asc',#f'{coupled_dir}/{my_scale}/SMASH_files/FLOW.asc', 
        epsg=2154,
        x=x_gauge,
        y=y_gauge,
        area= area, #3_192_100_000, 
        #bbox = (601_435, 663_000 ,6_159_000,6_258_977 ),
        code =code_gauge)
       
                    
model = smash.Model(calib_options, mesh = meshing)


model.optimize(mapping="uniform",
              #algorithm='sbs', 
                                   gauge=code_gauge,
                                   control_vector=["cp", "cft", "exc", "lr"], 
                                   jobs_fun='nse', 
                                   inplace = True, 
                                   options={'maxiter': 4})
# launch calibration
model.optimize(mapping="distributed",
                                   gauge=code_gauge,
                                   control_vector=["cp", "cft", "exc", "lr"], 
                                   jobs_fun='nse', 
                                   inplace = True, 
                                   options={'maxiter': 4})
                                   #ost)
                            
	
plt.plot(model.output.qsim[0,:], c= "red", label ="inf");plt.plot(model.input_data.qobs[0,:], c="blue", label ="obs");plt.legend();plt.show();plt.close()


all_q = np.ndarray(shape = (len(all_inflow),model.output.qsim_domain.shape[2]))
for key, value in all_inflow.items():
    my_id  = value["id"]
    all_q[key,:] = model.output.qsim_domain[my_id[0],my_id[1]]

sum_q = np.zeros(model.output.qsim_domain.shape[2])



for i in range(all_q.shape[0]):
    sum_q[:] =  sum_q[:] + all_q[i,:]
    

plt.plot(model.output.qsim[0,:], c= "red", label ="sim")
plt.plot(sum_q, c="blue", label ="sum inflow")
plt.legend()
plt.show()
plt.close()


# subdivide in 3 inflows according to mapping performed earlier in all_inflow (param bc_group_dassflow)

# --> get necessary bc parameter
# -->  

all_id_inflow = np.arange(df2d_mesh["boundaries"]["header"]["nb_inflow"])
qin = dict()
for mykey in all_id_inflow :
    qin[mykey] = np.zeros((model.output.qsim_domain.shape[2]))


for key, value in all_inflow.items():  
    my_id  = value["id"]
    if value["bc_group_dassflow"] in all_id_inflow:
         qin[value["bc_group_dassflow"]][:] =  qin[value["bc_group_dassflow"]][:] + model.output.qsim_domain[my_id[0],my_id[1]][:]


for key, val in qin.items():
    plt.plot(val, label =f"discharge {key}")

plt.plot(qin[0]+qin[1]+qin[2], label =f"sum_inflow")
plt.plot(model.output.qsim[0,:], c= "black", label ="inf");
#plt.plot(model.input_data.qobs[0,:], c="black", label ="obs");
plt.legend();
plt.show()
plt.close()

##### write bc.txt  and hydrograph.txt files

# qin_dict, dictionary, the key is the id of the hydrogram
# write_dir: (absolute) path where to write the file 
# return; write the hydrograph.txt file where precised 
def write_hydrograph(qin_dict, write_dir, dt ):
      # dt = calib_options["dt"]
    time = np.arange(start = 0, stop = dt * len(qin_dict[0]), step = dt) 
    with open(f'{write_dir}/hydrograph.txt', 'w') as f:
        f.write("#comment\n")
        f.write("#comment\n")
        f.write("#comment\n")
        f.write(str(len(qin_dict))+"\n")
        for hyd in qin_dict.values():
            f.write("#comment\n")
            f.write("#comment\n")
            f.write("#comment\n")
            f.write(f"{len(hyd)}\n")
            for i in range(len(time)):
                f.write(f"{time[i]} {hyd[i]}\n")
                
            print(hyd)
        f.write("#comment\n")
                
                
# dirty lilian
#qin[0][:] = 1000
qin[1][:] = 1000
qin[2][:] = 1000
write_hydrograph(qin_dict = qin, write_dir = bin_dir , dt = calib_options["dt"] )

            


####################################################
## perform dassflow simulation
####################################################
#print(bin_dir)
#model = df2d.dassflowmodel(bin_dir = bin_dir, 
#                           hdf5_path=bin_dir+"/res/simu.hdf5", 
#                           run_type= "direct",  
#                           clean = True)
#
#df2d.wrapping.m_common.set_mesh_name("final_mesh.geo")
#
#df2d.wrapping.m_common.set_ts( calib_options["dt"] * len(qin[0]))
#df2d.wrapping.m_common.set_dtw( 120)
#df2d.wrapping.m_common.set_dtp( 60)
#model.init_all()
#
#model.kernel.dof.h[:] = model.kernel.dof0.h[:] = 0.2
#model.run()
#
#model.save_all()
#
#model.post = df2d.core.output.Post(boundary_metadata = model.boundary.metadata, bin_dir = model.bin_dir)
#
#model.outputs.all_res[0]
#
#
#best_diff = 1000 #todo
#all_keys = [x for x in model.outputs.all_res.keys()]
#old_key= all_keys[10]
#for i in model.outputs.all_res.keys():
#    
#    previous_h = np.asanyarray(model.outputs.all_res[old_key]["h"])
#    h_array = np.asanyarray(model.outputs.all_res[i]["h"])
#    diff = previous_h[:] - h_array
#    diff = np.linalg.norm(diff)
#   # print(i, "diff=", diff)
#    if diff<best_diff:
#        print(i, "diff=", diff)
#        best_diff =  diff
#    model.meshing.mesh_pyvista.plot(scalars = model.outputs.all_res[i]["h"], 
#                                    cpos = "xy", show_edges = False)
#    old_key = i
#    
#    #for i in model.outputs.all_res.keys():
##    model.meshing.mesh_pyvista.plot(scalars = model.outputs.all_res[i]["v"], 
##                                    cpos = "xy", show_edges = False)
##    
##for i in model.outputs.all_res.keys():
##    model.meshing.mesh_pyvista.plot(scalars = model.outputs.all_res[i]["u"], 
##                                    cpos = "xy", show_edges = False)    
#    
#    # velocity
#for i in model.outputs.all_res.keys():
#    model.meshing.mesh_pyvista.plot(scalars = np.sqrt(model.outputs.all_res[i]["u"]**2 +model.outputs.all_res[i]["v"]**2), 
#                                    cpos = "xy", show_edges = False)    
#    
#    # Froude
#for i in model.outputs.all_res.keys():
#    model.meshing.mesh_pyvista.plot(scalars = np.sqrt(model.outputs.all_res[i]["u"]**2 +model.outputs.all_res[i]["v"]**2) / np.sqrt(10 * model.outputs.all_res[i]["h"]**2), 
#                                    cpos = "xy", show_edges = False)
#    
#model.meshing.mesh_pyvista.plot(scalars = model.outputs.all_res[0]["bathy"], 
#                                    cpos = "xy", show_edges = False)
#
#model.meshing.mesh_fortran.nc
##
#best_diff = 1000 #todo
#all_keys = [x for x in model.outputs.all_res.keys()]
#old_key= all_keys[10]
#for i in model.outputs.all_res.keys():
#    
#    previous_h = np.asanyarray(model.outputs.all_res[old_key]["h"])
#    h_array = np.asanyarray(model.outputs.all_res[i]["h"])
#    diff = previous_h[:] - h_array
#    diff = np.linalg.norm(diff)
#   # print(i, "diff=", diff)
#    if diff<best_diff:
#        print(i, "diff=", diff)
#        best_diff =  diff
#    model.meshing.mesh_pyvista.plot(scalars = model.outputs.all_res[i]["h"], 
#                                    cpos = "xy", show_edges = False)
#    old_key = i
##df2d.wrapping.call_model.clean_model(model.kernel)
##del model
#
#



