#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 10:24:36 2022

@author: livillenave
"""
import numpy as np
import gmsh 

import geopandas as gpd
import matplotlib.pyplot as plt
# USES
# >>>  get_xy_shapefile, source_shapefile from GIS_libs.py
    

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
#----------------------------------------------------#
# Build_gmesh mesh
# uses:
#  - path_major_bed
#  - path_minor_bed
#  - mesh_size_major
#  - mesh_size_minor
#  - include_minor_bed
# write .gmsh file  at f"{mesh_generation_dir}/{mesh_name}.msh"
# th gmsh model is in the instance via gmsh API and directly accessible
#----------------------------------------------------#
def build_gmsh_mesh(path_major_bed, path_minor_bed, mesh_size_major, mesh_size_minor, include_minor_bed, write_dir, mesh_name):
    
    # intialise information linked to the major bed
    point_majeur = source_shapefile(path_major_bed)    
    all_points_majeur= []
    all_lines_majeur= []
    lc_majeur = mesh_size_major
    
    
    # intialise information linked to the minor bed bed
    if include_minor_bed:
        point_mineur = source_shapefile(path_minor_bed)
        all_points_mineur= []
        all_lines_mineur= []
        lc_mineur = mesh_size_minor
    
    
    
    # INITIALIZE
    gmsh.initialize()
    # name model
    model = gmsh.model.add(f"{mesh_name}")
    
    # add nodes to mesh
    print(lc_majeur)
    for i in range(len(point_majeur)):
        a = gmsh.model.geo.addPoint(x = point_majeur[i,0], y = point_majeur[i,1], z = 0,  meshSize = lc_majeur, tag = i+1)
        all_points_majeur.append(a)
        
    # define lines from nodes
    for i  in range(len(point_majeur)-1) :
        a = gmsh.model.geo.addLine(startTag =  all_points_majeur[i],  endTag = all_points_majeur[i+1], tag = i+1)
        all_lines_majeur.append(i+1) 
        
    a = gmsh.model.geo.addLine(startTag = all_points_majeur[-1],  endTag =  all_points_majeur[0], tag = i+2)
    all_lines_majeur.append(a)
    
    # plot points  
    plt.scatter(point_majeur[:,0],point_majeur[:,1])
    
    if include_minor_bed:
        plt.scatter(point_mineur[:,0],point_mineur[:,1])
    
    if include_minor_bed:
        # id total 
        id_start = len(all_points_majeur)+1
        # add points
        print(lc_mineur)
        for i in range(len(point_mineur)):
            b = gmsh.model.geo.addPoint(x = point_mineur[i,0], y = point_mineur[i,1], z = 0, meshSize = lc_mineur)
            all_points_mineur.append(b)
            
        # define lines from nodes
        for i  in range(len(point_mineur)-1) :
            a = gmsh.model.geo.addLine(startTag =  all_points_mineur[i],  endTag = all_points_mineur[i+1], tag = id_start + i)
            all_lines_mineur.append(a) 
            
        a = gmsh.model.geo.addLine(startTag = all_points_mineur[-1],  endTag =  all_points_mineur[0], tag = id_start + i +1)
        all_lines_mineur.append(a)
        
    
    
    gmsh.model.geo.addCurveLoop(all_lines_majeur , 1)
    
    if include_minor_bed:
        gmsh.model.geo.addCurveLoop(all_lines_mineur , 2)
        
        
    # add plane surface
    
    if include_minor_bed:
        gmsh.model.geo.addPlaneSurface(wireTags =[1,2], tag=1)
        gmsh.model.geo.synchronize()
        gmsh.model.geo.addPlaneSurface(wireTags =[2], tag=2)
        gmsh.model.geo.synchronize()
    else : 
        gmsh.model.geo.addPlaneSurface(wireTags =[1], tag=1)
        
    
    gmsh.model.mesh.Coherence = True
    gmsh.model.mesh.generate(2)
    gmsh.write(f"{write_dir}/{mesh_name}.msh")
    
    
    
    #---------------------------
    # Extract nodes and cells for dassflow
    #---------------------------
    
    nodes = gmsh.model.mesh.getNodes()
    
    raw_cells =   gmsh.model.mesh.get_elements(dim = 2, tag = 1)
    id_elem =  raw_cells[1][0]
    
    
    if include_minor_bed:
        raw_cells2 =   gmsh.model.mesh.get_elements(dim = 2, tag = 2)
        id_elem2 =  raw_cells2[1][0]
        id_elem = np.concatenate((id_elem, id_elem2))
    
    cells = dict()
    compteur=0
    full_nodes = []
    for i in id_elem:
        compteur = compteur +1 
        my_cell = gmsh.model.mesh.get_element(i)
        cells[compteur] = dict()
        cells[compteur]["type"] = my_cell[0]
        cells[compteur]["nodeTags"] = my_cell[1]
        cells[compteur]["dim"] = my_cell[2]
        cells[compteur]["tag"] = my_cell[3]
        # for check below
        full_nodes.append( my_cell[1])
    
    # check all nodes are considered
    tmp = np.concatenate(full_nodes)
    tmp = np.unique(tmp)
    if np.all(tmp == nodes[0]) :
        print("all right")
    else:
        print("WARNING ALL NODES ARE NOT IN CELL MAYBE CHECK --- this will generate issue while dassflow meshing")
        
            
    # end gmsh instance
    gmsh.finalize()
    return(nodes, cells)