#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 10:19:56 2022

@author: livillenave
"""
def get_xy_cells(mesh_fortran):

    xcell = []
    ycell = []
    
    
    for i in range(mesh_fortran.nc ):
        xcell.append(mesh_fortran.cell[i].grav.x)
        ycell.append(mesh_fortran.cell[i].grav.y)
        
    xcell = np.asanyarray(xcell)
    ycell = np.asanyarray(ycell)
    return( xcell, ycell)
    


# ======================================================= #
#  SCRIPT                                                 #
# ======================================================= #
from osgeo import gdal
from osgeo import ogr
import matplotlib.pyplot as plt
import os
import dassflow2d as df2d

bin_dir =  "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A"
mnt_dir =  "/home/livillenave/Documents/data/DONNEES/prepared_data/MNT/"
source_mnt_name = "aude_mean.tif"
new_mnt_name    = "aude_cropped.tif"
output_dir =  f"{bin_dir}/mesh_generation"



model = df2d.dassflowmodel(bin_dir = bin_dir, 
                           hdf5_path=bin_dir+"/res/simu.hdf5", 
                           run_type= "direct",  
                           clean = True)
df2d.wrapping.m_common.set_mesh_name("dassflow_mesh_default.geo")
model.init_all()
    
    

# >>> crop MNT on area
xcell, ycell = get_xy_cells(model.meshing.mesh_fortran)
ulx=min(xcell)-1000
uly=max(ycell)+1000
lrx=max(xcell)+1000
lry=min(ycell)-1000
os.chdir(mnt_dir)
os.system(f'gdal_translate -projwin {ulx} {uly} {lrx} {lry}  -of GTiff -co "TILED=YES" {source_mnt_name} {new_mnt_name}')
os.chdir(bin_dir)


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

import numpy as np
rowcol_corresp = []
for i in range(len(xcell)): #
    print(i)
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


for i in range(len(rowcol_corresp)):
    print(rowcol_corresp[i])


zcell = []
for rowcol in rowcol_corresp:
    zcell.append(data_array[rowcol[0],rowcol[1]])
    
zcell = np.asanyarray(zcell)


# >>> UPDATE MESH

my_mesh = read_mesh_from_textfile(path = f"{output_dir}/dassflow_mesh2.geo",
                               read_boundary = True)
my_mesh["cell"]["bathy"] = zcell
gen_mesh_dassflow(mesh_dict = my_mesh, path = f"{output_dir}/dassflow_mesh3.geo",)


gen_mesh_dassflow(mesh_dict = my_mesh, path = f"{bin_dir}/final_mesh.geo",)


plt.scatter(xcell, ycell, c = zcell, cmap = "magma")
plt.colorbar()


df2d.wrapping.call_model.clean_model(model.kernel)
del model