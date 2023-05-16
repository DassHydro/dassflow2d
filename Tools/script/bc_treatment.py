#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 10:27:54 2022

@author: livillenave
"""

import matplotlib.pyplot as plt    
from matplotlib.colors import ListedColormap
import pyvista as pv
import numpy as np


def plot_hyd(hydrograph, title_plot= "Inflow discharge", subtitle = "group 0"):
    
    plt.plot(hydrograph.t, hydrograph.q)
    plt.xlabel("time [s]")
    plt.ylabel("Discharge [m3/s]")
    plt.suptitle(title_plot)
    plt.title(subtitle)
    plt.show()
    plt.close()
    
def plot_ratcurve(rat, title_plot= "Rating curve"):
    
    plt.plot(rat.q, rat.h)
    plt.xlabel("Discharge [m3/s]")
    plt.ylabel("Surface heigth [m]")
    plt.suptitle(title_plot)
    plt.title(f"group {rat.group}")
    plt.show()
    plt.close()
    
    
    
def plot_zspresc(zspresc, title_plot= "Rating curve"):
    
    plt.plot(zspresc.t, zspresc.z)
    plt.xlabel("time [s]")
    plt.ylabel("Imposed surface height [m]")
    plt.suptitle(title_plot)
    plt.title(f"group {zspresc.group}")
    plt.show()
    plt.close()


def plot_hpresc(hpresc, title_plot= "Rating curve"):
    
    plt.plot(hpresc.t, hpresc.h)
    plt.xlabel("time [s]")
    plt.ylabel("Water elevation [m]")
    plt.suptitle(title_plot)
    plt.title(f"group {hpresc.group}")
    plt.show()
    plt.close()
    
bc = df2d.wrapping.m_model.get_bc()


id_hyd = list(bc.hyd.indices)
id_zspresc = list(bc.zspresc.indices)
id_hpresc = list(bc.hpresc.indices)
id_rat = list(bc.rat.indices)

if len(id_hyd)>0:    
   for i in range(len(id_hyd)):
       plot_hyd(bc.hyd[i], title_plot= "Inflow discharge ", subtitle = f'group {i+1}')
       

if len(id_rat)>0:
   for i in range(len(id_rat)):
       plot_ratcurve(bc.rat[i])
       

if len(id_zspresc)>0:    
   for i in range(len(id_zspresc)):
       plot_zspresc(bc.zspresc[i])


if len(id_hpresc)>0:    
   for i in range(len(id_hpresc)):
       plot_hpresc(bc.hpresc[i])
                     
       
       
       
# example pour subplot
#id = 411
#plt.subplot(id)
#plot_hyd(hyd)  
#id = 412
#plt.subplot(id)
#plot_hyd(hyd)
#id = 413
#plt.subplot(id)
#plot_hyd(hyd)
#id = 414
#plt.subplot(id)
#plot_hyd(hyd)
    
    


grid = python_interface.grid


#a = python_interface.plot_meshing("edge")

typlim  = [] # typlim of the edge
group = [] # gorup of the edge
index = [] # global edge index
for i in list(python_interface.kernel.mesh.edgeb.indices):
    typlim.append(python_interface.kernel.mesh.edgeb[i].typlim.decode("utf8"))
    group.append( python_interface.kernel.mesh.edgeb[i].group)
    index.append( python_interface.kernel.mesh.edgeb[i].ind)
    
    
def find_indices(lst, condition):
       return np.array([i for i, elem in enumerate(lst) if elem == condition], dtype = "int")
   
    
def find_indices_double(lst, condition1, condition2):
       return np.array([i for i, elem in enumerate(lst) if elem == condition1 and elem == condition2], dtype = "int")



all_bc = dict()
<<<<<<< HEAD
conditions = [        "wall" , "discharg1", "discharg2", "neumann", "rat"    ,"hpresc" ,"zspresc"]
colors =     ["white","black", "red"      , "orange"   , "yellow" , "#2A495C"   ,"#779231", "#329C46" ]
=======
conditions = [        "wall" , "discharg1", "discharg2", "neumann", "rat"    ,"hpresc" ,"zpresc"]
colors =     ["black","white", "red"      , "orange"   , "yellow" , "#2A495C"   ,"#779231", "#329C46" ]
>>>>>>> dev_lilian
             # no bc,  wall, discharg1, etc...


res = dict()

for my_condition in conditions:
    
    nb = typlim.count(my_condition)   
    if my_condition == "wall":
        tofind = find_indices(typlim, my_condition )
        ind_final = []
        for i in tofind :
            ind_final.append(index[i]-1)
        res[my_condition] = {0:ind_final}
    else:
        if nb>0:
            my_ind = find_indices(typlim, my_condition )
            all_group =[]
            
            for tmp_ind in my_ind:                
                all_group.append(group[tmp_ind])  
            
            if len(np.unique(all_group)) == len(all_group):
                if my_condition not in res.keys():
                            ind_final = []
                            for i in my_ind :
                                ind_final.append(index[i]-1)            
                            res[my_condition] = {all_group[0]:ind_final}
            else:
                for i in all_group:        
                    ind_group = find_indices(group, i)
                    # remind that my_ind correspond to index with correct bc type
                    
                    # intersection of group and type give unical bc
                    tmp = np.intersect1d(my_ind, ind_group)
                    ind_final = []
                    for j in tmp :
                                ind_final.append(index[j]-1)   
                    
                    if my_condition not in res.keys():
                        res[my_condition] = {f"{group[ind_group[0]]}":ind_final}
                    else:
                        res[my_condition][f"{group[ind_group[0]]}"] = ind_final
              


def plot_bc(grid, bc_corresp):
    
		    
    def find_indices(lst, condition):
		       return np.array([i for i, elem in enumerate(lst) if elem == condition], dtype = "int")
		   
		    
    def find_indices_double(lst, condition1, condition2):
		       return np.array([i for i, elem in enumerate(lst) if elem == condition1 and elem == condition2], dtype = "int")
			
    conditions = [        "wall" , "discharg1", "discharg2", "neumann", "rat"    ,"hpresc" ,"zpresc"]
    colors =     ["black","white", "red"      , "orange"   , "yellow" , "#2A495C"   ,"#779231", "#329C46" ]
                  
    newcol=["black"]
    new_cond=["none"]
    for k in bc_corresp.keys():
        id = find_indices(conditions, k)[0]+1
        newcol.append(colors[id])
        new_cond.append(k)
    
        
    legend_entries = []
    # not a bc condition given first
    legend_entries.append(['none', 'black'])
    
    edges = grid.extract_all_edges()            
    edge_color= np.zeros(edges.n_faces) 
    
    
    for i, k in enumerate(bc_corresp.keys()):
        color_index= i+1# color index correspond also to the color calue given to edge_color
        for g in bc_corresp[k].keys():
                edge_color[bc_corresp[k][g]] = color_index
                print(edge_color[bc_corresp[k][g]])
        # add legend entry
        col = newcol[color_index]
        print(i,k,  col)
        legend_entries.append([f'{k}', f'{col}'])
    
    # create scalar typlim (within pyvista polydata object )
    edges["typlim"] = edge_color
    
    tmp = pv.Plotter(notebook = False)
    
    tmp.add_mesh(edges,scalars = "typlim", cmap=ListedColormap(newcol), show_scalar_bar=True)
    tmp.add_legend(legend_entries, face = "line")
    tmp.show_bounds(
        grid='front',
        location='outer',
        all_edges=False)
    #tmp.show_bounds(padding=0.05, location = 'outer')
    tmp.show(cpos = "xy")
    
    
    
    
    
    
    
    
    
    #############################################
    # SIMILAR PLOT WITH MAIN INFLOW POSITION
    #############################################
    
    
    
    
    #ind_str = [str(i+1) for i, el in enumerate(points[:,0])]
    #ind_str = [str(i+1) for i in range(points.shape[0])]
    #dict_ind_str = {str(i+1): f"main_inflow={i+1} \n goto" for i in range(points.shape[0])}
    #ind_str = [str(i+1) for i in range(points.shape[0])]
    
   # lab =  [f'main_inflow {i+1} \n goto hydraulic_bc {corresp2[i]}' for i in range(points.shape[0])]
    
    
    legend_entries=[
            ["main inflow hydrology", "blue"],
            ["inflow hydraulic", "red"]
                   ]
    tmp = pv.Plotter(notebook = False)
    
    tmp.add_mesh(edges,scalars = "typlim", cmap=ListedColormap(newcol), show_scalar_bar=True)
    tmp.add_legend(legend_entries, face = "rectangle")
    tmp.show_bounds(
        grid='front',
        location='outer',
        all_edges=False)
    
    tmp.add_point_labels(points, lab, point_size = 10, text_color = "blue", point_color = "blue")
    tmp.show(cpos = "xy")



#############################################
# SIMILAR PLOT only MAIN INFLOW POSITION
#############################################

def dict_to_array(dictionary):
    result = dictionary.items()
    data = list(result)
    numpyArray = np.array(data)
    return(numpyArray)
  
def SPECIALdict_to_array(dictionary):
    result = dictionary.items()
    data = list(result)
    numpyArray = np.array(data)
    
    res = np.zeros( (len(numpyArray),3) )
    for i in range(len(numpyArray)):
        res[i,0]=numpyArray[i,1][0]
        res[i,1]=numpyArray[i,1][1]
    return(res)
    
    
point_bcin = dict_to_array(center_inflow)
point_bcinv2 = SPECIALdict_to_array(center_inflow)
#ind_str = [str(i+1) for i, el in enumerate(points[:,0])]

lab_point =  [f'main_inflow {i+1} \n goto hydraulic_bc {corresp2[i]}' for i in range(points.shape[0])]
lab_hyd  =  [f'bc {point_bcin[i][0]} \n from m_inflow {[x+1 for x in corresp1[point_bcin[i][0]]]}' for i in range(len(point_bcin))]
lab_hyd  =  [f'bc {point_bcin[i][0]}' for i in range(len(point_bcin))]



legend_entries=[
        ["main inflow hydrology", "blue"],
        ["inflow hydraulic", "red"]
               ]



tmp = pv.Plotter(notebook = False)

tmp.add_mesh(edges,scalars = "typlim", cmap=ListedColormap(newcol), show_scalar_bar=True)
tmp.add_legend(legend_entries, face = "rectangle")
#tmp.show_bounds(
#    grid='front',
#    location='outer',
#    all_edges=False)
tmp.add_point_labels(point_bcinv2, lab_hyd, point_size = 10, text_color = "red", point_color = "red")
tmp.show(cpos = "xy")





#############################################
# SIMILAR PLOT WITH MAIN INFLOW POSITION
#############################################

def dict_to_array(dictionary):
    result = dictionary.items()
    data = list(result)
    numpyArray = np.array(data)
    return(numpyArray)
  
def SPECIALdict_to_array(dictionary):
    result = dictionary.items()
    data = list(result)
    numpyArray = np.array(data)
    
    res = np.zeros( (len(numpyArray),3) )
    for i in range(len(numpyArray)):
        res[i,0]=numpyArray[i,1][0]
        res[i,1]=numpyArray[i,1][1]
    return(res)
    
    
point_bcin = dict_to_array(center_inflow)
point_bcinv2 = SPECIALdict_to_array(center_inflow)
#ind_str = [str(i+1) for i, el in enumerate(points[:,0])]

lab_point =  [f'main_inflow {i+1} \n goto hydraulic_bc {corresp2[i]}' for i in range(points.shape[0])]
lab_hyd  =  [f'bc {point_bcin[i][0]} \n from m_inflow {[x+1 for x in corresp1[point_bcin[i][0]]]}' for i in range(len(point_bcin))]
#lab_hyd  =  [f'bc {point_bcin[i][0]}' for i in range(len(point_bcin))]



legend_entries=[
        ["main inflow hydrology", "blue"],
        ["inflow hydraulic", "red"]
               ]



tmp = pv.Plotter(notebook = False)

tmp.add_mesh(edges,scalars = "typlim", cmap=ListedColormap(newcol), show_scalar_bar=True)
tmp.add_legend(legend_entries, face = "rectangle")
tmp.show_bounds(
    grid='front',
    location='outer',
    all_edges=False)

tmp.add_point_labels(points, lab_point, point_size = 10, text_color = "blue", point_color = "blue")
tmp.add_point_labels(point_bcinv2, lab_hyd, point_size = 10, text_color = "red", point_color = "red")
tmp.show(cpos = "xy")






  
# print the numpy array



# ILUMINE CELL

a = np.zeros(grid.n_cells)
a[:] = 0
a[1207-1] = 1
a[1133-1] = 1
a[1067-1] = 1

grid["a"] = a
tmp = pv.Plotter(notebook = False)
tmp.add_mesh(grid,scalars = "a")
tmp.show(cpos = "xy")
