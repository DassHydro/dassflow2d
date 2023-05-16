#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 10:17:13 2022

@author: livillenave
"""
import numpy as np
import matplotlib.pyplot as plt


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
            
    drained_area =  mesh["flwacc"]
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
    tmp = mesh["flwacc"].copy()
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
   