#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 10:12:56 2022

@author: livillenave
"""


import geopandas as gpd
import numpy as np


#  main inflows hydrology

shapefile = gpd.read_file("/home/livillenave/Téléchargements/aude_new/hydrology/coordinates_main_inflow.shp")

points = np.zeros((shapefile.geometry.shape[0], 3))

for i in range(shapefile.geometry.shape[0]):
    points[i,0] = shapefile.geometry.x[i]
    points[i,1] = shapefile.geometry.y[i]
    points[i,2] = 0
    
    
#  inflows hydraulic
         # res produced in bc_treatment

center_inflow = dict()

for group in  ["discharg1"].keys():
    if len(res["discharg1"][group]) ==1 :
        id_edge = res["discharg1"][group][0]
        to_fill = np.zeros((3))
        to_fill[0] =  direct_model.model.mesh.edge[id_edge].center.x
        to_fill[1] =  direct_model.model.mesh.edge[id_edge].center.y
        to_fill[2] =  0
        coords_group = to_fill
    elif  len(res["discharg1"][group]) ==0 :
        print("ERROR, group has no edge linked")
    else :
        id_edge = res["discharg1"][group]
        
        to_fill = np.zeros((len(id_edge), 3))
        for ie in range(len(id_edge)):
            to_fill[ie, 0] =  direct_model.model.mesh.edge[id_edge[ie]].center.x
            to_fill[ie, 1] =  direct_model.model.mesh.edge[id_edge[ie]].center.y
            to_fill[ie, 2] =  0
            coords_group = np.mean(to_fill, axis = 0)
    center_inflow[group] = coords_group
            
        



#####################
# Correspondance V1
# find closest main inflow to existing hydraulic bc    
#####################"
#corresp1 = dict()
#
#for group in  center_inflow.keys():
#   reference = center_inflow[group]
#   distances = np.zeros(len(points))
#   for i in range(len(points)):
#       distances[i] = np.sqrt( (points[i,0]-reference[0])**2 +  (points[i,1]-reference[1])**2  )
#       id_main_inflow = distances.argmin()
#   corresp[group] = id_main_inflow    
#   
#numpy_corresp = np.zeros((len(corresp1),2))
#
#for i in range(len(corresp1)):
#    id_hydraulic = list(corresp1.keys())[i]
#    id_hydrology = corresp[id_hydraulic]
#    numpy_corresp[i,0] = int(id_hydraulic)
#    numpy_corresp[i,1] = id_hydrology
#   
#
#test = numpy_corresp[:,1]
#id_ordered = np.argsort(test)
#numpy_corresp[:,:] = numpy_corresp[id_ordered,:]









###############################################################
# Correspondance V2
# find closest hydraulic bc    to any  main inflow
###############################################################
corresp2 = dict()

# i = id main inflow
# reference = coordinate of associated main inflow :[x,y, 0] 
for i,reference in  enumerate(points):
   print(i,reference)
   distances = np.zeros(len(center_inflow))
   # k 
   for id_distance, k in enumerate(center_inflow.keys()):
       distances[id_distance] = np.sqrt( (center_inflow[k][0]-reference[0])**2 +  (center_inflow[k][1]-reference[1])**2  )
   id_key = distances.argmin()
   id_group = list(center_inflow.keys())[id_key]
   corresp2[i] = id_group    
   
   


###############################################################
# --> return dictionary for each bc of which closest 
###############################################################
    
def getKeysByValue(dictOfElements, valueToFind):
    listOfKeys = list()
    listOfItems = dictOfElements.items()
    for item  in listOfItems:
        if item[1] == valueToFind:
            listOfKeys.append(item[0])
    return  listOfKeys

corresp1 = dict()

for group in center_inflow.keys():
    corresp1[group] =getKeysByValue(corresp2, group)
    
    