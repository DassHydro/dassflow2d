#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 20:04:07 2023

@author: livillenave
"""
import osr
import gdal 

def array2raster(newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array, path):
   """ This function works fine, until the raster mask. It does the work whatsoever. 
        Not to be used unless necessary. """
   reversed_arr = array[::-1]    
   cols = reversed_arr.shape[1]
   rows = reversed_arr.shape[0]
   originX = rasterOrigin[0]
   originY = rasterOrigin[1]

   driver = gdal.GetDriverByName('GTiff')
   outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float64)
   outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
   outband = outRaster.GetRasterBand(1)                                                                                                                              

   outband.WriteArray(reversed_arr)
   outband.SetNoDataValue(-999)
   outRasterSRS = osr.SpatialReference()
   outRasterSRS.ImportFromEPSG(2154) # http://spatialreference.org/ref/epsg/32614/
   outRaster.SetProjection(outRasterSRS.ExportToWkt())
   outband.FlushCache()
   
   
cell_size = 1000
array2raster('jan1950' + '.tif', [smash_model.mesh.xmin,smash_model.mesh.ymax], cell_size, cell_size, meshing["drained_area"], "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/")

array= meshing["active_cell"]
array[np.where(meshing["active_cell"]==0)] = -999
array2raster('jan1950' + '.tif', [smash_model.mesh.xmin,smash_model.mesh.ymax], cell_size, cell_size,array, "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/")


import rasterio
from rasterio.features import shapes
import geopandas as gp

mask = None
with rasterio.Env():
    with rasterio.open('/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/'+'jan1950' + '.tif') as src:
        image = src.read(1) # first band
        image = np.asanyarray(image, dtype = "float32")
        results = (
        {'properties': {'raster_val': v}, 'geometry': s}
        for i, (s, v) 
        in enumerate(
            shapes(image, mask=mask, transform=src.transform)))
geoms = list(results)

gpd_polygonized_raster  = gp.GeoDataFrame.from_features(geoms)

gpd_polygonized_raster.to_file('/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/toshape')


mypoly=[]
for vec in rasterio.features.shapes(myarray):
    mypoly.append(shape(vec))