#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 09:58:03 2022

@author: livillenave
"""

import geopandas as gpd
import numpy as np



# get xy numpy ndarray from path using goepandas package
# from shapefile.shp
def get_xy_shapefile(path):
    
    shapefile = gpd.read_file( path  )    
    x_coord = []
    y_coord = []    
    for index, row in shapefile.iterrows():
         for pt in list(row['geometry'].exterior.coords): 
            x_coord.append(pt[0])
            y_coord.append(pt[1])
            
    res = np.ndarray(shape=(len(x_coord),2) )
    res[:,0] = x_coord
    res[:,1] = y_coord
    return(res)
    
    
# mainly return xy numpy ndarray from get_xy_shapedfile(path) 
# clean points if there are repetitions
def source_shapefile(path):
        
    # load shapefile
    xy_points =  get_xy_shapefile(path)  
    # clean if last point and first are the same (remove the last one)
    if    (xy_points[-1,0], xy_points[-1,1]) == (xy_points[0,0], xy_points[0,1])  :
        xy_points = xy_points[:-1,:]        
    return(xy_points)
    

#############
# a la mano pour intersection river_network / hydraulic contour
# INPUT
# river_network =  geopandas.geodataframe.GeoDataFrame
# contour_hy =  shapely.geometry.multilinestring.MultiLineString  (with 1 linestring for each edge)
# OUTPUT
# x_intersect, y_intersect: lists of arrays of coordinate
def get_xy_intersect_shapely_lines(river_network, contour_hy):

    points_intersect =  river_network.intersection(contour_hy)
    
    all_intersect_point = []
    compteur = 0
    for elem in points_intersect:
        
        if isinstance(elem, geometry.point.Point) :
            compteur = compteur + 1
            all_intersect_point.append([elem.x, elem.y])
            
        elif isinstance(elem, geometry.MultiPoint) :
            for i in range(len(elem)):
                compteur = compteur + 1
                all_intersect_point.append([elem[i].x, elem[i].y])               
    
    x_intersect = [point[0] for point in all_intersect_point]
    y_intersect = [point[1] for point in all_intersect_point]
    return(x_intersect, y_intersect)
  