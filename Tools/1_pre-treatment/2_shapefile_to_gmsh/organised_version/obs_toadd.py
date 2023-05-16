#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 10:05:40 2023

@author: livillenave
"""

#--------------------------- #

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


#

from osgeo import gdal, ogr
import os

import geopandas as gpd
import pandas as pd
import numpy as np
import pickle
import shapely
from shapely.geometry import Point
import dassflow2d as df2d

import pyvista as pv
pv.global_theme.color = 'white'
pv.global_theme.show_edges = True

#--------------------------------#
# BIN directory containing actual mesh file
#--------------------------------#
demo_dir = "/home/livillenave/Documents/distant/Demo"
source_bin_dir = f"{demo_dir}/built_bin3" # built_bin
mesh_name = "final_mesh.geo"

#--------------------------------#
# Laisse crue
#--------------------------------#


# Build polygon of laisse de crue

cropper_path = "/home/livillenave/Documents/distant/Demo/source_file/hydraulic_shapefile/major_bed_aude.shp"
zone_crue1 = "/home/livillenave/Documents/data/select_ftp_pics/laisse_crue/aude_amont.shp"
zone_crue2 = "/home/livillenave/Documents/data/select_ftp_pics/laisse_crue/aude_aval.shp"
tempdir = "/home/livillenave/Documents/trash_perso/"
res_dir = "/home/livillenave/Documents/data/select_ftp_pics/prepared/"

# clip zone crue 1;zone crue 2
os.system(f"qgis_process run native:clip --distance_units=meters --area_units=m2 --ellipsoid=EPSG:2154 --INPUT={zone_crue1} --OVERLAY={cropper_path} --OUTPUT={tempdir}/tempres1.shp")
os.system(f"qgis_process run native:clip --distance_units=meters --area_units=m2 --ellipsoid=EPSG:2154 --INPUT={zone_crue2} --OVERLAY={cropper_path} --OUTPUT={tempdir}/tempres2.shp")
# merge zones1 AND 2
os.system(f"qgis_process run native:union --distance_units=meters --area_units=m2 --ellipsoid=EPSG:2154 --INPUT={tempdir}/tempres1.shp --OVERLAY={tempdir}/tempres2.shp --OVERLAY_FIELDS_PREFIX= --OUTPUT={res_dir}/merged_laisse_crue.shp")



dassflow_model = df2d.dassflowmodel(bin_dir = source_bin_dir, 
                           hdf5_path=source_bin_dir+"/res/simu.hdf5", 
                           run_type= "direct",  
                           clean = True)

df2d.wrapping.m_common.set_mesh_name(mesh_name)

dassflow_model.init_all()

laisse_crue = gpd.read_file(f"{res_dir}/merged_laisse_crue.shp")
x=y=point=[]
is_flooded=[]
for idpython_cell in range(dassflow_model.meshing.mesh_fortran.nc):
    dassflow_model.meshing.mesh_fortran.cell[idpython_cell]
    local_x = dassflow_model.meshing.mesh_fortran.cell[idpython_cell].grav.x
    local_y = dassflow_model.meshing.mesh_fortran.cell[idpython_cell].grav.y  
    x.append(local_x)
    y.append(local_y)
    local_point= Point(local_x, local_y)
    point.append(local_point)
    is_flooded.append(np.any(laisse_crue.contains(local_point))*1)
    
is_flooded = np.asanyarray(is_flooded, dtype = "int")
dassflow_model.meshing.mesh_pyvista.plot(cpos="xy", show_edges = True, scalars=is_flooded)

xcell_dassflow = np.asanyarray(x)
ycell_dassflow = np.asanyarray(y)
#--------------------------------#
# PHE
#--------------------------------#

phe =  gpd.read_file(f"/home/livillenave/Documents/data/select_ftp_pics/PHE_Hocini_et_al/phe2018_ID.shp")
laisse_crue.contains(phe["geometry"])


phe_indomain=[]
for id_phe in range(np.shape(phe)[0]):
    phe_indomain.append(np.any(laisse_crue.contains(phe["geometry"][id_phe])))

id_phe_indomain= np.where( np.asanyarray(phe_indomain, dtype = "int") ==1 )

phe_cropped = phe.iloc[id_phe_indomain]


phe_cropped["geometry"].plot(figsize = (12,20))
for i in range(len(laisse_crue["geometry"])):
    if isinstance(laisse_crue["geometry"][i], shapely.geometry.polygon.Polygon):
        x,y = laisse_crue["geometry"][i].exterior.xy
        plt.plot(x, y, c="red")
    elif isinstance(laisse_crue["geometry"][i], shapely.geometry.multipolygon.MultiPolygon):
        for j in range(len(laisse_crue["geometry"][i])):
            x,y = laisse_crue["geometry"][i][j].exterior.xy
            plt.plot(x, y, c="red")
    else:
        print("bug")
    
plt.show()
    


phe_z = np.asanyarray(phe_cropped.Z_ref)
phe_x = np.asanyarray(phe_cropped.geometry.x)
phe_y = np.asanyarray(phe_cropped.geometry.y)

phe_points = dict()
for id_phe in range(len(phe_x)):
    phe_points[id_phe] = Point(phe_x[id_phe],phe_y[id_phe])
    phe_points[id_phe]


phe_cropped["geometry"].keys()

def calc_dist(x1, y1, x2, y2):
    return(np.sqrt((x1-x2)**2+(y1-y2)**2))
    
    
    
def get_id_dassflow_cell_corresponding(dict_points, dassflow_mesh):
    id_dassflow_cell = []
    xcell,ycell = get_xy_cells(dassflow_mesh)
    for my_key in dict_points.keys():
        x =  dict_points[my_key].x
        y =  dict_points[my_key].y
        dist = calc_dist(xcell,ycell, x, y)
        id_closest_cell = np.argmin(dist)
        id_dassflow_cell.append(id_closest_cell)
        
    return(id_dassflow_cell)   
    
    
id_phe_dassflow_cell = get_id_dassflow_cell_corresponding(dict_points=phe_cropped["geometry"], dassflow_mesh = dassflow_model.meshing.mesh_fortran)
id_phe_dassflow_cell = np.asanyarray(id_phe_dassflow_cell)

phe_ed_cell = np.zeros(dassflow_model.meshing.mesh_fortran.nc) # goal is to creaty 
phe_ed_cell[id_phe_dassflow_cell] = 1


dassflow_model.meshing.mesh_pyvista.plot(cpos="xy", scalars=phe_ed_cell, background="white")
#--------------------------------#
# add MNT INFO TO PHE
#--------------------------------#
# source file location
source_dir =  f"{demo_dir}/source_file/"
mnt_dir =  f"{source_dir}/DEM"
new_mnt_name    = "aude_cropped.tif"


# >>> LOAD CROPPED MNT
ds = gdal.Open(f"{mnt_dir}/{new_mnt_name}")
ulx, xres, _, uly, _, yres  = ds.GetGeoTransform()
# gives more explicit name
upper_left_x = ulx
upper_left_y = uly
x_resolution = xres
y_resolution = yres

# >>> Source value
data_array = ds.GetRasterBand(1).ReadAsArray()
x_ndarray = data_array.copy()
y_ndarray = data_array.copy()
# >>> get xy of each raster cell
for row_id in range(x_ndarray.shape[0]):
    for col_id in range(x_ndarray.shape[1]):
        # x change à chaque colomne
        # y change à chaque ligne
        x_ndarray[row_id, col_id]= upper_left_x + x_resolution * col_id
        y_ndarray[row_id, col_id]= upper_left_y + y_resolution * row_id


plt.imshow(x_ndarray); plt.colorbar();plt.show();plt.close()
plt.imshow(y_ndarray); plt.colorbar();plt.show();plt.close()



rowcol_corresp = []
for i in range(len(phe_x)): #
    x = phe_x[i]
    y = phe_y[i]
    
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



phe_z_mnt = []
for rowcol in rowcol_corresp:
    phe_z_mnt.append(data_array[rowcol[0],rowcol[1]])
    
phe_z_mnt = np.asanyarray(phe_z_mnt)
phe_zref = phe_z


# create figure and axis objects with subplots()
fig,ax = plt.subplots()

plt.title("closest DEM value vs Hocini PHE")
ax.plot(phe_z_mnt, label = "mnt  source (closest raster cell value)")
ax.plot(phe_zref, label = "data  source (raw file hocini)")
ax.legend()
# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
ax2.plot(phe_zref -  phe_z_mnt, label = "hocini - mnt", color = "red", marker = "o")
ax2.legend()
plt.show()
plt.close()




#=================================#
# PYVISTA JOB
#==================================#
#-----------------
# plot river network
#-----------------
river_network_path = "/home/livillenave/Documents/distant/Demo/built_bin/mesh_generation/river_network.shp"
river_network = gpd.read_file(river_network_path)
print(type(river_network["geometry"][0]))

river_lines = []
# convert each edge into a line
for _, row in river_network.iterrows():
    type(row['geometry'])
    x_pts = row['geometry'].xy[0]
    y_pts = row['geometry'].xy[1]
    z_pts = np.zeros(len(x_pts))
    pts = np.column_stack((x_pts, y_pts, z_pts))
    line = pv.lines_from_points(pts)
    river_lines.append(line)

combined_river_lines = river_lines[0].merge(river_lines[1:])
combined_river_lines.plot(line_width=3, cpos='xy')

pv.global_theme.color = 'white'
plotter = pv.Plotter(off_screen = False, notebook = False)
actor1 = plotter.add_mesh(mesh = combined_river_lines, color = "blue")
actor2 = plotter.add_mesh(mesh = dassflow_model.meshing.mesh_pyvista, show_edges = True)
plotter.show(cpos = "xy")

#-----------------
# plot lit mineur
#-----------------
minor_bed_path = "/home/livillenave/Documents/distant/Demo/source_file/hydraulic_shapefile/minor_bed.shp"
minor_bed = gpd.read_file(minor_bed_path)
print(type(minor_bed["geometry"][0]))
minor_bed_lines = []
# convert each edge into a line
for _, row in minor_bed.iterrows():
    type(row['geometry'])
    x_pts, y_pts = row['geometry'].exterior.xy
    z_pts = np.zeros(len(x_pts))
    pts = np.column_stack((x_pts, y_pts, z_pts))
    line = pv.lines_from_points(pts)
    minor_bed_lines.append(line)
combined_minor_bed_lines = minor_bed_lines[0].merge(minor_bed_lines[1:])

pv.global_theme.color = 'white'
plotter = pv.Plotter(off_screen = False, notebook = False)
actor1 = plotter.add_mesh(mesh = combined_minor_bed_lines, color = "blue")
actor2 = plotter.add_mesh(mesh = dassflow_model.meshing.mesh_pyvista, show_edges = True)
plotter.show(cpos = "xy")



#-----------------
# plot laisse crue
#-----------------
print(type(laisse_crue["geometry"][0]))
minor_bed_lines = []
# convert each edge into a line
for _, row in laisse_crue.iterrows():
    type(row['geometry'])
    if(isinstance(row['geometry'], shapely.geometry.polygon.Polygon)):
        x_pts, y_pts = row['geometry'].exterior.xy
        z_pts = np.zeros(len(x_pts))
        pts = np.column_stack((x_pts, y_pts, z_pts))
        line = pv.lines_from_points(pts)
        minor_bed_lines.append(line)
    else:
        for j in range(len(row['geometry'])):
            x_pts, y_pts = row['geometry'][j].exterior.xy
            z_pts = np.zeros(len(x_pts))
            pts = np.column_stack((x_pts, y_pts, z_pts))
            line = pv.lines_from_points(pts)
            minor_bed_lines.append(line)
combined_minor_bed_lines = minor_bed_lines[0].merge(minor_bed_lines[1:])

plotter = pv.Plotter(off_screen = False, notebook = False)
plotter.set_background( 'white', top=None)
actor1 = plotter.add_mesh(mesh = combined_minor_bed_lines, color = "red", line_width=5)
actor2 = plotter.add_mesh(mesh = dassflow_model.meshing.mesh_pyvista, scalars=is_flooded, show_edges = True, color =pv.Color("white", opacity=0))
plotter.show(cpos = "xy")



#-----------------
# plot PHE localisation
#-----------------

points = phe_x, phe_y
import numpy as np
import pyvista

rng = np.random.default_rng()

points = rng.random((len(phe_x), 3))
points[:,0]=phe_x  #x
points[:,1]=phe_y  #y
points[:,2]= 0     #z
pset = pyvista.PointSet(points)

plotter = pv.Plotter(off_screen = False, notebook = False)
plotter.set_background( 'white', top=None)
actor1 = plotter.add_mesh(mesh = pset, color = "red", line_width=5)
actor2 = plotter.add_mesh(mesh = dassflow_model.meshing.mesh_pyvista, scalars=is_flooded, show_edges = True, color =pv.Color("white", opacity=0))
plotter.show(cpos = "xy")

del plotter
# ---------------------
# Plot DEM
# -------------
x_ndarray
y_ndarray


xpoint=[]
ypoint=[]
zpoint=[]
for row in range(x_ndarray.shape[0]):
    for col in range(x_ndarray.shape[1]):
        xpoint.append(x_ndarray[row,col])
        ypoint.append(y_ndarray[row,col])
        zpoint.append(data_array[row,col])
finalx = np.asanyarray(xpoint)
finaly = np.asanyarray(ypoint)
finalz = np.asanyarray(zpoint)
points = rng.random((len(finalx), 3))
points[:,0]=finalx  #x
points[:,1]=finaly #y
points[:,2]= 0     #z
pset = pyvista.PointSet(points)



plotter = pv.Plotter(off_screen = False, notebook = False)
plotter.set_background( 'white', top=None)
actor3 = plotter.add_mesh(mesh =  pset, point_size = 5, opacity = 0.01)
actor1 = plotter.add_mesh(mesh = pset, color = "red", line_width=5)
actor2 = plotter.add_mesh(mesh = dassflow_model.meshing.mesh_pyvista, scalars=is_flooded, show_edges = True, color =pv.Color("white", opacity=0))
plotter.show(cpos = "xy")


test = pset.cast_to_unstructured_grid()