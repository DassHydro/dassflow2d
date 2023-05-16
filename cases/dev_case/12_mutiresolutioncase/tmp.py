#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 11:26:42 2022

@author: livillenave
"""


model.meshing.mesh_fortran.nc

xcoord=[]
ycoord=[]

for i in range(model.meshing.mesh_fortran.nc):
    xcoord.append(model.meshing.mesh_fortran.cell[i].grav.x)
    ycoord.append(model.meshing.mesh_fortran.cell[i].grav.y)
    
    
res = np.ndarray( (model.meshing.mesh_fortran.nc, 4) )

res[:,0]  = xcoord
res[:,1]  = ycoord
res[:,2]  = 60
res[:,3]  = 1

f = open("new_obs.txt", "w")
for i in range(np.shape(res)[0]):
    f.write(str(res[i,0])+ " "+ str(res[i,1])+ " "+ str(res[i,2])+ " "+ str(res[i,3])+ " \n")
    
f.close()

