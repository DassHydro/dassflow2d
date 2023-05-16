#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 15:51:49 2023

@author: livillenave
"""

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
bin_dir = f"{demo_dir}/bin6"
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
#path_major_bed = f"{source_dir}/hydraulic_shapefile/lit_majeur_opti.shp"
#path_minor_bed = "/home/livillenave/Documents/distant/Demo/source_file/tosort/minor_bed_2.shp"#f"{source_dir}/hydraulic_shapefile/lit_mineur_opti.shp"
# path flow diretion
path_leblois_data = f"{source_dir}/DEM/LEBLOIS_DATA/10/flow_dir.asc"

# name mesh
mesh_name = "aude"

# mesh size in meters
mesh_size_major = 100 # major bed
mesh_size_minor = 50  # minor bed




path_major_bed = f"/home/livillenave/Documents/distant/SD-FLOOD/synthetic_case_coupling/new_0603/major_bed.shp"
path_minor_bed = f"/home/livillenave/Documents/distant/SD-FLOOD/synthetic_case_coupling/new_0603/minor_bed.shp"
mesh_generation_dir= "/home/livillenave/Documents/distant/SD-FLOOD/synthetic_case_coupling/new_0603/mesh_gen"

# set False  to disable minor bed
include_minor_bed =  False


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

if not os.path.exists(f'{mesh_generation_dir}'):
    os.mkdir(f'{mesh_generation_dir}')
os.chdir(f'{mesh_generation_dir}')


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
        gmsh.model.geo.synchronize()
        
    
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
    
    
#----------------------------------------------------#
# 1- Build Gmsh mesh
#----------------------------------------------------# 
    
all_size  = [1000, 500, 100, 50, 5, 1]
for my_size in all_size:
    my_mesh_name = f"mesh_{my_size}"
    
    nodes, cells = build_gmsh_mesh(path_major_bed = path_major_bed, 
                                    mesh_size_major = my_size, 
                                    include_minor_bed = include_minor_bed, 
                                    path_minor_bed = path_minor_bed,  
                                    mesh_size_minor = mesh_size_minor,
                                    write_dir = mesh_generation_dir,
                                    mesh_name = my_mesh_name)
    
    
    gen_empty_mesh_dassflow(nodes = nodes, 
                            cells = cells, 
                            path = f"{mesh_generation_dir}/raw_{my_mesh_name}.geo")



