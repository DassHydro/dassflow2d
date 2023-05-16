#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 17:03:43 2022

@author: livillenave
"""


# -------------------------
# Import librairies
# -------------------------
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

# -------------------------
# define librairies
# -------------------------
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
    

# Write 'dassflow_mesh.geo'  without bc (inlet 0 0 and outlet 0 0)
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
  
    
    
# mesh : dassflow2d.core.Meshing.mesh.mesh_fortran
# index_wall_edge:  np.ndarray , containing index of edge from which you want to extract center 
# return 2 arrays x and y which correspond to the center of edge
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


# as the title says
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


    
    

# as the title says
# path : path to textfile
# if readboundary is true, read and get boundary data, else do not try to read this part of the file
# return dictionary containing all mesh information
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
    
    
    

# write mesh from  mesh_dict, wich has the same structure as the one gotten with read_mesh
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
                  values["ghost_bathy"]= 300 # DIRTY lilian
                  f.write(f'{values["index_fortran_cell"]} {values["index_relative_edge"]} {values["index_relative_edge"]} {values["ghost_bathy"]}  {values["group"]}\n')
                  
              f.write(f'OUTLET {mesh_dict["boundaries"]["header"]["nb_outflow"]} {mesh_dict["boundaries"]["header"]["nb_cell_outflow"]}\n')
              for key, values in mesh_dict["boundaries"]["outlet"].items():
                  f.write(f'{values["index_fortran_cell"]} {values["index_relative_edge"]} {values["index_relative_edge"]} {values["ghost_bathy"]}  {values["group"]}\n')
        f.close()
        return()



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
                if not (row_modif == 0 and col_modif == 0): # if not the cell itselfprint("not in id_hydraulic ??? ", ((idrow + row_modif,  idcol  +col_modif) not in id_hydraulic) )
                 if (idrow + row_modif,  idcol  +col_modif) not in id_hydraulic:  # if the adjacent cell is not of hydraulic domain
                     check1 = ( idrow+row_modif<= idrow_max)
                     check2=(idrow+row_modif>=0)
                     check3 = (idcol+col_modif<= idcol_max)
                     check4 =( idcol+col_modif>=0)
                     check =   np.all( [check1,check2,check3,check4])               
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
   

# -------------------------
# parameters
# -------------------------
    
demo_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/Demo"
    
    
    
    
    
# were to invoke dassflow2d    
bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/Demo/bin"
# were most of intermediate files (mostly GIS are written)
mesh_generation_dir = f"{bin_dir}/mesh_generation"

# source file location
source_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/Demo/source_file/"


mnt_dir =  f"{source_dir}/DEM"

source_mnt_name = "aude_mean.tif"
new_mnt_name    = "aude_cropped.tif"
path_raw_river_network = f"{source_dir}/river_network/COURS_D_EAU_NATUREL_mono_oriente_1a5.shp"# path to shapefiles
path_major_bed = f"{source_dir}/hydraulic_shapefile/major_bed_aude.shp"
path_minor_bed = f"{source_dir}/hydraulic_shapefile/minor_bed.shp"
path_leblois_data = f"{source_dir}/DEM/LEBLOIS_DATA/10/flow_dir.asc"

if not os.path.exists(f'{mesh_generation_dir}'):
    os.mkdir(f'{mesh_generation_dir}')
os.chdir(f'{mesh_generation_dir}')

# name mesh
mesh_name = "aude"


# mesh size in meters
mesh_size_major = 300 # major bed
mesh_size_minor = 50 # minor bed

# set False  to disable minor bed
include_minor_bed =  True

# ================================================================ #
# Script
# 
# 1 - Build Gmsh mesh
# 2- traduce it to dassflow_mesh.geo without boundary conditions
# 3- Add boundary condition to mesh
# 4- Add bathymetry to mesh
# 5- generate SMASH meshing () 
# 6- Associate  (map) dassflow  and smash mesh
# ================================================================  #


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
gmsh.write(f"{mesh_generation_dir}/{mesh_name}.msh")





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



#----------------------------------------------------#
# 2- Traduce gmsh to  dassflow mesh
#----------------------------------------------------#
gen_empty_mesh_dassflow(nodes = nodes, cells = cells, path = f"{mesh_generation_dir}/raw_{mesh_name}.geo")
gen_empty_mesh_dassflow(nodes = nodes, cells = cells, path = f"{bin_dir}/dassflow_mesh_default.geo")

# end gmsh instance
gmsh.finalize()



#----------------------------------------------------#
#  3-  Add boundary condition to mesh
#----------------------------------------------------#



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
contour_hy.to_file(f'{mesh_generation_dir}/contour_hy.shp')



##########
# source river network

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
    'dt':900,
		'start_time':        '2018-10-14 12:00',
		'end_time':        '2018-10-15 20:00',
    
    'read_qobs':True,
    'qobs_directory':'/home/livillenave/Documents/data/FORCING/DEBIT/15',
    
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
                                   options={'maxiter': 10})
# launch calibration
model.optimize(mapping="distributed",
                                   gauge=code_gauge,
                                   control_vector=["cp", "cft", "exc", "lr"], 
                                   jobs_fun='nse', 
                                   inplace = True, 
                                   options={'maxiter': 10})
                                   #ost)
                            

plt.figure(figsize = (20,20))
plt.plot(model.output.qsim[0,:], c= "red", label ="inf")
plt.plot(model.input_data.qobs[0,:], c="blue", label ="obs")
plt.scatter(x = np.arange(0,150,4), c="blue", label ="obs")
plt.legend()
plt.show()
plt.close()


model.run(inplace = True)

plt.figure(figsize = (20,20))
plt.plot(model.output.qsim[0,:], c= "red", label ="inf")
plt.plot(model.input_data.qobs[0,:], c="blue", label ="obs")
plt.scatter(x = np.arange(0,150,4), y = np.zeros(int(150/4)), c="blue", label ="obs")
plt.legend()
plt.show()
plt.close()

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
#qin[1][:] = 1000
#qin[2][:] = 1000
write_hydrograph(qin_dict = qin, write_dir = bin_dir , dt = calib_options["dt"] )
                
    
###################################################
# perform dassflow simulation
###################################################
print(bin_dir)
model = df2d.dassflowmodel(bin_dir = bin_dir, 
                           hdf5_path=bin_dir+"/res/simu.hdf5", 
                           run_type= "direct",  
                           clean = True)

df2d.wrapping.m_common.set_mesh_name("final_mesh.geo")

df2d.wrapping.m_common.set_ts( calib_options["dt"] * len(qin[0]))
df2d.wrapping.m_common.set_dtw( 60)
df2d.wrapping.m_common.set_dtp( 10)
model.init_all()

model.kernel.dof.h[:] = model.kernel.dof0.h[:] = 0.2
model.run()

model.save_all()


model.outputs.all_res[0]


best_diff = 1000 #todo
all_keys = [x for x in model.outputs.all_res.keys()]
old_key= all_keys[10]
for i in model.outputs.all_res.keys():
    
    previous_h = np.asanyarray(model.outputs.all_res[old_key]["h"])
    h_array = np.asanyarray(model.outputs.all_res[i]["h"])
    diff = previous_h[:] - h_array
    diff = np.linalg.norm(diff)
   # print(i, "diff=", diff)
    if diff<best_diff:
        print(i, "diff=", diff)
        best_diff =  diff
    model.meshing.mesh_pyvista.plot(scalars = model.outputs.all_res[i]["h"], 
                                    cpos = "xy", show_edges = False)
    old_key = i
#for i in model.outputs.all_res.keys():
#    model.meshing.mesh_pyvista.plot(scalars = model.outputs.all_res[i]["v"], 
#                                    cpos = "xy", show_edges = False)
#    
#for i in model.outputs.all_res.keys():
#    model.meshing.mesh_pyvista.plot(scalars = model.outputs.all_res[i]["u"], 
#                                    cpos = "xy", show_edges = False)    
    
    # velocity
for i in model.outputs.all_res.keys():
    model.meshing.mesh_pyvista.plot(scalars = np.sqrt(model.outputs.all_res[i]["u"]**2 +model.outputs.all_res[i]["v"]**2), 
                                    cpos = "xy", show_edges = False)    
    
    # Froude
for i in model.outputs.all_res.keys():
    model.meshing.mesh_pyvista.plot(scalars = np.sqrt(model.outputs.all_res[i]["u"]**2 +model.outputs.all_res[i]["v"]**2) / np.sqrt(10 * model.outputs.all_res[i]["h"]**2), 
                                    cpos = "xy", show_edges = False)
    
model.meshing.mesh_pyvista.plot(scalars = model.outputs.all_res[0]["bathy"], 
                                    cpos = "xy", show_edges = False)

model.meshing.mesh_fortran.nc
#
best_diff = 1000 #todo
all_keys = [x for x in model.outputs.all_res.keys()]
old_key= all_keys[10]
for i in model.outputs.all_res.keys():
    
    previous_h = np.asanyarray(model.outputs.all_res[old_key]["h"])
    h_array = np.asanyarray(model.outputs.all_res[i]["h"])
    diff = previous_h[:] - h_array
    diff = np.linalg.norm(diff)
   # print(i, "diff=", diff)
    if diff<best_diff:
        print(i, "diff=", diff)
        best_diff =  diff
    model.meshing.mesh_pyvista.plot(scalars = model.outputs.all_res[i]["h"], 
                                    cpos = "xy", show_edges = False)
    old_key = i
#df2d.wrapping.call_model.clean_model(model.kernel)
#del model





