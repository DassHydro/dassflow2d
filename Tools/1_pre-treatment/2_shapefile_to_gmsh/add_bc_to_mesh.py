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

# get index of wall edges in the mesh_fortran of model object
# input:
# mesh_fortran: df2d.core.meshing.Meshing.mesh_fortran
# output: 
# index_wall_edges: numpy array of PYTHON (from 0) index of "wall" edges  in FORTRAN.
def get_index_wall_edges(model):
    
    index_wall_edges = []
    for ie in range(model.meshing.mesh_fortran.ne):
        if model.meshing.mesh_fortran.edge[ie].boundary == 1:
            index_boundary_edge = model.meshing.mesh_fortran.edge[ie].lim -1
            if model.meshing.mesh_fortran.edgeb[index_boundary_edge].typlim.decode("utf-8") == "wall"  or model.meshing.mesh_fortran.edgeb[index_boundary_edge].typlim.decode("utf-8") == "discharg1":
                index_wall_edges.append(ie)
            
    index_wall_edges = np.asanyarray( index_wall_edges)
    return(index_wall_edges)
    
    
# get nodes coordinates for each edge_index in mesh_fortran
# INPUT: 
# mesh_fortran :   df2d.core.meshing.Meshing.mesh_fortran
# index_edge: index from wich extract points
# OUTPUT
# xy_by_edge: dictionary (keys are the index_edge)
# contains xy_by_edge[key_edge] = dict()  containg {"x1", "y1", "x2", "y2"}
def get_xy_nodes_from_edge(model, index_wall_edges):
    all_nodes = []
    nodes_ordered = dict()
    xy_by_edge = dict()
    
    for ie in index_wall_edges:
        nodes_ordered[ie] = dict()
        nodes_ordered[ie][1]= model.meshing.mesh_fortran.edge[ie].node[0]
        nodes_ordered[ie][2]= model.meshing.mesh_fortran.edge[ie].node[1]
        
    for key, value in nodes_ordered.items():
        xy_by_edge[key] = dict()
        index_point1= value[1]
        index_point2= value[2]
        index_python_point1= index_point1 - 1 
        index_python_point2= index_point2 - 1
        
        x1= model.meshing.mesh_fortran.node[index_python_point1].coord.x
        y1= model.meshing.mesh_fortran.node[index_python_point1].coord.y
        
        x2= model.meshing.mesh_fortran.node[index_python_point2].coord.x
        y2= model.meshing.mesh_fortran.node[index_python_point2].coord.y
        
        xy_by_edge[key]["x1"] = x1
        xy_by_edge[key]["y1"] = y1
        xy_by_edge[key]["x2"] = x2
        xy_by_edge[key]["y2"] = y2

    return(xy_by_edge)
        
# Build multi line strings from all nodes gotten by  get_nodes_from_edge
def build_multi_line_string(all_nodes):
    
    all_lines = []
    for key, value in all_nodes.items():
        point1 = geometry.Point(value["x1"], value["y1"])
        point2 = geometry.Point(value["x2"], value["y2"])
        line1 = geometry.LineString(point1.coords[:] + point2.coords[:])    
        all_lines.append(line1)
        del line1 
    multi_line = geometry.MultiLineString(all_lines)
    return(multi_line)
    
    
  

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
  
    
    
  
# 
def get_xy_center_edge(mesh, index_wall_edges):
    
    all_x = []
    all_y = []
    
    if isinstance(index_wall_edges, np.ndarray ):
        for index in index_wall_edges:
            x = mesh.edge[index].center.x
            y = mesh.edge[index].center.y
            all_x.append(x)
            all_y.append(y)
    elif isinstance(index_wall_edges, int ) or isinstance(index_wall_edges, np.int64 ):
        all_x.append(mesh.edge[index_wall_edges].center.x)
        all_y.append(mesh.edge[index_wall_edges].center.y)
        
    res_x = np.asanyarray(all_x)
    res_y = np.asanyarray(all_y) 
    
    return(res_x, res_y)  
    
# mesh = mesh_fortran
# x_ref: float : x coordinate reference
# y_ref: float : y coordinate reference
# edge_index_to_test : float : 
def get_index_closest_edge(mesh, x_ref, y_ref, edge_index_to_test):
    
    all_x,all_y = get_xy_center_edge(mesh, index_wall_edges=index_wall_edges)
    best_edge = 1
    best_dist = 10**10
    for j in range(len(all_x)):
        x_edge = all_x[j]
        y_edge = all_y[j]            
        dist = np.sqrt((x_ref-x_edge)**2 + (y_ref-y_edge)**2 )            
        if dist < best_dist:
            best_dist = dist
            best_edge = index_wall_edges[j]
    return(best_edge)


        
def get_boundary_metadata_for_textfile_mesh(mesh, index_edge_intersect):
    
    index_python_edge_intersect =  index_edge_intersect 
    
    index_cell = mesh.edge[index_python_edge_intersect].cell[0]
    index_python_cell = index_cell -1
    edge_index_in_cell_structure = int(np.where(mesh.cell[index_python_cell].edge -1 ==  index_edge_intersect)[0])
    
    res = dict()
    res["index_cell"]= index_cell
    res["index_python_cell"]= index_cell -1
    res["edge_index_in_cell_structure"]= edge_index_in_cell_structure
    return(res)


    
    

def read_mesh_from_textfile(path, read_boundary = False ):

# >>> Source File
    with open(path, 'r') as f:
        a = f.readlines()    
    
# >>> extract "metadata" + build list index
    nb_node, nb_cell, scale = np.asanyarray(a[1].split(), dtype = "float")

    read_node_indexes = np.arange(nb_node, dtype = "int") +3
    read_cell_indexes = np.arange(nb_cell, dtype = "int") + read_node_indexes[-1] + 2

# same for boundary
    if read_boundary:
      # inlet
        inlet_header_index = read_cell_indexes[-1]+2
        nb_inlet = int(a[inlet_header_index].split()[1])
        nb_cells_inlet = int(a[inlet_header_index].split()[2]    )
        read_inlet_indexes = np.arange(nb_cells_inlet, dtype = "int") +inlet_header_index +1
      # outlet
        outlet_header_index = read_inlet_indexes[-1] +1
        nb_outlet = int(a[outlet_header_index].split()[1])
        nb_cells_outlet = int(a[outlet_header_index].split()[2])
        read_outlet_indexes = np.arange(nb_cells_outlet, dtype = "int") +outlet_header_index +1
        
    
# >>> Define variables to fill
    nodes_x = []
    nodes_y = []
    nodes_z = []
    nodes_index = []
    
    
    all_cell_index = []
    node_1 = []
    node_2 = []
    node_3 = []
    node_4 = []
    all_land_use =[]
    all_bathy = []
    
    if read_boundary:        
        #  
        INLET_index_fortran_cell =  []
        INLET_index_relative_edge = []
        INLET_ghost_bathy = []
        INLET_group = []        
        # 
        OUTLET_index_fortran_cell =  []
        OUTLET_index_relative_edge = []
        OUTLET_ghost_bathy = []
        OUTLET_group = []
        
# Nodes section ---------------------------------------
    
    for node_index in read_node_indexes:
        index_node, x_coord_node, y_coord_node, z_coord_node =  np.asanyarray(a[node_index].split(), dtype = "float") 
        nodes_x.append(x_coord_node)
        nodes_y.append(y_coord_node)
        nodes_z.append(z_coord_node)        
        nodes_index.append(int(index_node))
        
# Cells section ---------------------------------------
        
    for cell_index in read_cell_indexes: 
        index_cell, node1, node2, node3, node4, land_use, bathy =  np.asanyarray(a[cell_index].split(), dtype = "float")      
        all_cell_index.append(int(index_cell))
        node_1.append(int(node1))
        node_2.append(int(node2))
        node_3.append(int(node3))
        node_4.append(int(node4))
        all_land_use.append( int(land_use ))
        all_bathy.append(bathy)
        
        
# Boundary section ---------------------------------------
    if read_boundary:
        # read inlet
        for inlet_index in read_inlet_indexes :
            index_fortran_cell, index_relative_edge , useless, ghost_bathy, group =  np.asanyarray(a[inlet_index].split(), dtype = "int")   
            
            INLET_index_fortran_cell.append(index_fortran_cell)
            INLET_index_relative_edge.append(index_relative_edge)
            INLET_ghost_bathy.append(ghost_bathy)
            INLET_group.append(group)   
            
        for outlet_index in read_outlet_indexes :
            index_fortran_cell, index_relative_edge , useless, ghost_bathy, group =  np.asanyarray(a[outlet_index].split(), dtype = "int")   
            
            OUTLET_index_fortran_cell.append(index_fortran_cell)
            OUTLET_index_relative_edge.append(index_relative_edge)
            OUTLET_ghost_bathy.append(ghost_bathy)
            OUTLET_group.append(group)   
        
        
        
        
    res = dict()
    
    res["header"] = dict()
    res["header"]["nb_node"]=nb_node
    res["header"]["nb_cell"]=nb_cell
    
    res["node"] = dict()
    res["node"]["x"]=nodes_x
    res["node"]["y"]=nodes_y
    res["node"]["z"]=nodes_z    
    res["node"]["index"]=nodes_index
    
    res["cell"] = dict()
    res["cell"]["index"]    = all_cell_index
    res["cell"]["node_1"]   = node_1
    res["cell"]["node_2"]   = node_2
    res["cell"]["node_3"]   = node_3
    res["cell"]["node_4"]   = node_4
    res["cell"]["land_use"] = all_land_use
    res["cell"]["bathy"]    = all_bathy

    if read_boundary:
        
        res["boundaries"] = dict()
        res["boundaries"]["header"] = dict()
        res["boundaries"]["header"]["nb_inflow"] = nb_inlet
        res["boundaries"]["header"]["nb_cell_inflow"] = nb_cells_inlet
        res["boundaries"]["header"]["nb_outflow"] = nb_outlet
        res["boundaries"]["header"]["nb_cell_outflow"] = nb_cells_outlet
        
        res["boundaries"]["inlet"] = dict()
        for key in range(nb_inlet):
             res["boundaries"]["inlet"][key]  = dict()
             
             res["boundaries"]["inlet"][key]["index_fortran_cell"] = INLET_index_fortran_cell[key]
             res["boundaries"]["inlet"][key]["index_relative_edge"] = INLET_index_relative_edge[key]
             res["boundaries"]["inlet"][key]["ghost_bathy"] = INLET_ghost_bathy[key]
             res["boundaries"]["inlet"][key]["group"] = INLET_group[key]
        
        res["boundaries"]["outlet"] = dict()
        for key in range(nb_outlet):
            res["boundaries"]["outlet"][key]  = dict()
            res["boundaries"]["outlet"][key]["index_fortran_cell"]  = OUTLET_index_fortran_cell[key]
            res["boundaries"]["outlet"][key]["index_relative_edge"] = OUTLET_index_relative_edge[key]
            res["boundaries"]["outlet"][key]["ghost_bathy"]         = OUTLET_ghost_bathy[key]
            res["boundaries"]["outlet"][key]["group"]               = OUTLET_group[key]

        
    return(res)
    
    
    


def gen_mesh_dassflow(mesh_dict, path = 'dassflow_mesh.geo'):
       
    nb_node = int(mesh_dict["header"]["nb_node"])
    nb_cell = int(mesh_dict["header"]["nb_cell"])
    
    with open(path, 'w') as f:
# Write HEADER
        f.write('# mesh generated from gmsh using gen_empty_mesh_dassflow() \n')  # python will convert \n to os.linesep
        f.write(f'{nb_node} {nb_cell} 0. \n') 
# Write Nodes block
        f.write('#Nodes||| id node, x coord, y coord, bathymetry\n ') 
        for id_node in range(nb_node):
              towrite_id_node = mesh_dict["node"]["index"][id_node]
              x = mesh_dict["node"]["x"][id_node]
              y = mesh_dict["node"]["y"][id_node]
              z = mesh_dict["node"]["z"][id_node]
              f.write(f'{towrite_id_node} {x} {y} {z} \n') 
# Write cell block
        f.write('#cells||| id cell, id_node1, id_node2, id_node3, id_node4, patch_manning, bathymetry \n ') 
                
        for id_cell in range(nb_cell):              
              cell_index = mesh_dict["cell"]["index"][id_cell]
              id_node1 =  mesh_dict["cell"]["node_1"][id_cell]
              id_node2 =  mesh_dict["cell"]["node_2"][id_cell]
              id_node3 = mesh_dict["cell"]["node_3"][id_cell]
              id_node4 =  mesh_dict["cell"]["node_4"][id_cell]
              land_use =   mesh_dict["cell"]["land_use"][id_cell]
              bathy =   mesh_dict["cell"]["bathy"][id_cell]                  
              f.write(f'{cell_index} {id_node1} {id_node2} {id_node3} {id_node4} {land_use} {bathy}\n') 


        if "boundaries" in mesh_dict.keys():
              f.write("# boundaries \n")
                      
              f.write(f'INLET {mesh_dict["boundaries"]["header"]["nb_inflow"]} {mesh_dict["boundaries"]["header"]["nb_cell_inflow"]}\n')
              for key, values in mesh_dict["boundaries"]["inlet"].items():
                  f.write(f'{values["index_fortran_cell"]} {values["index_relative_edge"]} {values["index_relative_edge"]} {values["ghost_bathy"]}  {values["group"]}\n')
                  
              f.write(f'OUTLET {mesh_dict["boundaries"]["header"]["nb_outflow"]} {mesh_dict["boundaries"]["header"]["nb_cell_outflow"]}\n')
              for key, values in mesh_dict["boundaries"]["outlet"].items():
                  f.write(f'{values["index_fortran_cell"]} {values["index_relative_edge"]} {values["index_relative_edge"]} {values["ghost_bathy"]}  {values["group"]}\n')
        f.close()
        return()


# ======================================================= #
#  SCRIPT                                                #
# ======================================================= #
        
import dassflow2d as df2d

import numpy as np  
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt      
import dassflow2d as df2d
from shapely import geometry


bin_dir =  "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A"
output_dir =  f"{bin_dir}/mesh_generation"




# initialise fortran instance, and python corrponding data
model = df2d.dassflowmodel(bin_dir = bin_dir, 
                           hdf5_path=bin_dir+"/res/simu.hdf5", 
                           run_type= "direct",  
                           clean = True)

df2d.wrapping.m_common.set_mesh_name("dassflow_mesh_default.geo")

model.init_all()
    
    


index_wall_edges = get_index_wall_edges(model = model)
all_nodes = get_xy_nodes_from_edge(model = model, index_wall_edges=index_wall_edges)
multi_line = build_multi_line_string(all_nodes = all_nodes)


###########
# write contour hy shapefile

contour_hy = gpd.GeoDataFrame({"geometry": [multi_line]}, crs=2154)
contour_hy.to_file(f'{output_dir}/contour_hy.shp')



##########
# source river network

river_network = gpd.read_file(filename = "/home/livillenave/Documents/data/DONNEES/spatial/COURS_D_EAU_NATUREL_mono_oriente_1a5.shp",
                              bbox = multi_line,
                              crs = 2154)

river_network.to_file(f'{output_dir}/river_network.shp')

##############

# get values
(x_intersect, y_intersect) = get_xy_intersect_shapely_lines(river_network=river_network, contour_hy=multi_line)

# save
df = pd.DataFrame({'x': x_intersect, 'y': y_intersect})
gdf = gpd.GeoDataFrame(
    df, geometry=gpd.points_from_xy(df['x'], df['y']))
gdf.to_file(f'{output_dir}/intersect_point.shp')



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
                                        index_edge_intersect = my_edge
                                        )
    
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
        
        
mesh = read_mesh_from_textfile(f'{output_dir}/dassflow_mesh.geo')
# add boundaries to mesh file
mesh["boundaries"] = dict()
mesh["boundaries"]["header"] = dict()
mesh["boundaries"]["inlet"] = dict()
mesh["boundaries"]["outlet"]= dict()

compteur_inflow = 0

for key, value in boundaries_to_add.items():
    if value["bc_type"] == "discharg1":
        compteur_inflow = compteur_inflow+1
        mesh["boundaries"]["inlet"][key]= {"index_fortran_cell":value["index_cell"],                     # id cell
         "index_relative_edge":value["edge_index_in_cell_structure"]+1,  # relative edge
         "ghost_bathy":value["edge_index_in_cell_structure"],  # ghost cell bathy
         "group":compteur_inflow }
        
    if value["bc_type"] == "transm":
        mesh["boundaries"]["outlet"][key]= {"index_fortran_cell":value["index_cell"],                     # id cell
         "index_relative_edge":value["edge_index_in_cell_structure"]+1,  # relative edge
         "ghost_bathy":value["edge_index_in_cell_structure"],  # ghost cell bathy
         "group":compteur_inflow }

mesh["boundaries"]["header"]["nb_inflow"] = len(mesh["boundaries"]["inlet"]) 
mesh["boundaries"]["header"]["nb_cell_inflow"] = len(mesh["boundaries"]["inlet"])   # to update  for multi-cell bc
mesh["boundaries"]["header"]["nb_outflow"] = len(mesh["boundaries"]["outlet"]) 
mesh["boundaries"]["header"]["nb_cell_outflow"] = len(mesh["boundaries"]["outlet"])    # to update for multi-cell bc



gen_mesh_dassflow(mesh_dict = mesh, path = f"{output_dir}/dassflow_mesh2.geo")


df2d.wrapping.call_model.clean_model(model.kernel)
del model