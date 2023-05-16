#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 15:30:15 2022

@author: livillenave
"""
# -------------------------
# Import librairies
# -------------------------
import gmsh
import shapely
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import os


# -------------------------
# define librairies
# -------------------------
# get xy numpy ndarray from path using goepandas package
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
    

# 
def gen_empty_mesh_dassflow(nodes, cells, path = 'dassflow_mesh.geo' ):
    
    all_id_node = nodes[0]
    
    nb_node = len(all_id_node)
    nb_cell = len(cells)
    
    with open(path, 'w') as f:
# Write HEADER
        f.write('# mesh generated from gmsh using gen_empty_mesh_dassflow() \n')  # python will convert \n to os.linesep
        f.write(f'{nb_node} {nb_cell} 0. \n') 
# Write Nodes block
        f.write('#Nodes||| id node, x coord, y coord, bathymetry\n ') 
        for id_node in all_id_node:
              python_id_node = int(id_node -1)
              x = nodes[1][python_id_node * 3 + 0] 
              y = nodes[1][python_id_node * 3 + 1] 
              z = nodes[1][python_id_node * 3 + 2] 
              f.write(f'{id_node} {x} {y} {z} \n') 
# Write cell block
        f.write('#cells||| id cell, id_node1, id_node2, id_node3, id_node4, patch_manning, bathymetry \n ') 
        for id_cell, my_cell in cells.items():
              nb_node = len(my_cell['nodeTags'])
              
              id_node1=my_cell['nodeTags'][0]
              id_node2=my_cell['nodeTags'][1]
              id_node3=my_cell['nodeTags'][2]
              if nb_node == 3:
                  id_node4=0#my_cell['nodeTags'][0]
              elif nb_node == 4:
                  id_node4=my_cell['nodeTags'][3]
              else:
                  print("WARNING: nb node not 3 or 4")
                  
              f.write(f'{id_cell} {id_node1} {id_node2} {id_node3} {id_node4} 1 0.\n')
# write empty boundary block
        f.write("# boundaries \n")
        f.write("INLET 0 0 \n")
        f.write("OUTLET 0 0 \n") 

# -------------------------
# parameters
# -------------------------
bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/"
write_dir = f"{bin_dir}/mesh_generation"
shapefile_dir = "/home/livillenave/Documents/distant/SD-FLOOD/real_case-AUDE/DATA/DASSFLOW/V2/hydraulic/tmp/"
os.chdir(f'{write_dir}')

# name mesh
mesh_name = "aude"

# path to shapefiles
path_major_bed = f"{shapefile_dir}/contour1_v2.shp"
path_minor_bed = f"{shapefile_dir}/contour2.shp"


# mesh size in meters
mesh_size_major = 200 # major bed
mesh_size_minor = 100 # minor bed

# set False  to disable minor bed
include_minor_bed =  True

# -------------------------
# Script
# -------------------------

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




gmsh.initialize()
model = gmsh.model.add(f"{mesh_name}")

# add nodes to mesh
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
plt.scatter(point_mineur[:,0],point_mineur[:,1])

if include_minor_bed:
    id_start = len(all_points_majeur)+1
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
    
gmsh.model.geo.addPlaneSurface(wireTags =[1,2], tag=1)
gmsh.model.geo.synchronize()
gmsh.model.geo.addPlaneSurface(wireTags =[2], tag=2)
gmsh.model.geo.synchronize()
gmsh.model.mesh.Coherence = True
gmsh.model.mesh.generate(2)
gmsh.write(f"{mesh_name}.msh")





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
    print("WARNING ALL NODES ARE NOT IN CELL MAYBE CHECK")


gen_empty_mesh_dassflow(nodes = nodes, cells = cells)
gen_empty_mesh_dassflow(nodes = nodes, cells = cells, path = f"{bin_dir}/dassflow_mesh_default.geo")

gmsh.finalize()




