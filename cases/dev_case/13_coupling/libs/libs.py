#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 13:05:32 2022

@author: livillenave
"""
# add_bc
#----
# 1 - add to bc.txt  the add boundary
# 2 - if do not exist, add create text file associated to the boundary
# 3 - if .txt file, update it to invlue the given table value
#
# warnings
#----
# works only if inflow boundary
#
# inputs
# bin_dir : bin directory where text files are append
#----
import re
#re.sub("[^0-9]", "", "sdkjh987978asd098as0980a98sd")


# generate time value table for dassflow from smash output qsimgrid
# param
# ----
# tuple_index, tuple of x,y index,  it correspond to x,y coordin
import numpy as np
# qin, nb_timestep, dt
def gen_table_fromqsim(qin, nb_timestep, dt):
    time_serie = np.zeros(nb_timestep)
    for i in range(nb_timestep):
        time_serie[i] = dt*i    
    return(np.row_stack( (time_serie, qin) ))
    

def add_bc(bin_dir, bc_type, table_values):    
    
#============================#
# Param                       #
#============================#
    path_bctxt = f"{bin_dir}/bc.txt"
    # path tabletxt defined in if condition bellow
    
    # define if it is a in, out, internal, --wall--
    if bc_type in ["discharg1","discharg2"]:
        bc_class = "in" 
        path_tabletxt = f"{bin_dir}/hydrograph.txt"
    else:
        warning("only dicharg implemented")
        return()
        
#============================#
# bc.txt treatment
#============================#
    f = open(path_bctxt)
    data = f.readlines()
    for i, line in enumerate(data):
        print("i=", i)
       # print("line=",line)
        if i ==3 :
            # get actual number of boundary
            nb_bc = data[i].replace('\n', '')
            nb_bc = np.asarray(nb_bc, dtype='int64')
            
            # add 1 boundary
            nb_bc = nb_bc+1
            #replace line of text file
            line = f"{nb_bc}\n"
            data[i] = line
    last_line = data[-1]
    x = last_line.split()
    out_type = x[1]
    
    # update  last as added inflow bc
    data[-1] = f"{x[0]} {bc_type} file\n" 
    
    # add  new line as outflow bc (updated index)
    data.append(f"{nb_bc} {out_type} file\n")
     
        # write file
    with open(path_bctxt, 'w', encoding='utf-8') as file:
        file.writelines(data)
        
        
        
#============================#
# hydrograph.txt treatment
#============================#
    f = open(path_tabletxt)
    data = f.readlines()
    for i, line in enumerate(data):
        print("i=", i)
       # print("line=",line)
        if i ==3 :
            # get actual number of hydrograph
            nb_hydrograph = data[i].replace('\n', '')
            nb_hydrograph = np.asarray(nb_hydrograph, dtype='int64')            
            # add 1 hydrograph in count
            nb_hydrograph = nb_hydrograph+1
            
            #replace line of text file
            line = f"{nb_hydrograph}\n"
            # update data structure
            data[i] = line
            
            
    # add necessary 3 comment rows hydrograph
    data.append("!===========================!\n")
    data.append("!===========================!\n")
    data.append("!===========================!\n")
    # add table values
    data.append(f"{np.shape(table_values)[1]}\n")    
    for i in range(np.shape(table_values)[1]):
        data.append(f"{table_values[0,i]} {table_values[1,i].round(6)}\n")   
        
    # wrinte file
    with open(path_tabletxt, 'w', encoding='utf-8') as file:
        file.writelines(data)
    
    
            
            
        
        # rotation matrix to align dassflow case with hydrological case
def rot_mat(deg):
    theta = deg * np.pi / 180 
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c,-s], [s,c]])


# id_cell, tupple of (x,y) if it was a matrix and not a numpy ndarray
def add_xy_coord(mesh):
        shape = mesh["drained_area"].shape
        x_coord = np.zeros(shape)
        y_coord = np.zeros(shape)

        for i in range(shape[0]):
             for j in reversed(range(shape[1])):
                 if i == 0:
                     x_coord[:,i] =mesh["dx"]/2
                 else :
                     x_coord[:,i] = x_coord[:,i-1] + mesh["dx"]                     
                 if j == shape[1]-1:
                     y_coord[j,:] =mesh["dx"]/2
                 else:
                     y_coord[j,:] = y_coord[j+1,:] + mesh["dx"]
        mesh["x_coord"] = x_coord 
        mesh["y_coord"] = y_coord
                 
                 
                 
        return(mesh)