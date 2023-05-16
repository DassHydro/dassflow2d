#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 09:00:39 2022

@author: livillenave
"""



import dassflow2d as df2d
import numpy as np
import pandas as pd
import sys
import os
import random
from matplotlib.colors import ListedColormap
import h5py



all_manning = dict()

for key in ["uniform", "1", "2", "3", "4" ]:
    print(key)
    tmp = h5py.File(f"/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/hdf5_file/inf{key}.hdf5")
    all_manning[key] = dict()
    all_manning[key]["mesh_correspondance"] = tmp["input/param/manning/mesh_correspondance"][...]
    all_manning[key]["patch_value"] = tmp["input/param/manning/patch_value"][:]
    tmp.close()

os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A")
run_type = "direct"
model = df2d.dassflowmodel(bin_dir = os.getcwd(), 
                           hdf5_path=os.getcwd()+"/res/simu.hdf5", 
                           run_type = run_type,
                           clean=True)
model.init_all()


for key in ["uniform", "1", "2", "3", "4"]:
    # get platch value
    values = np.array(
            [all_manning[key]["patch_value"][int(x-1)]     for x in     all_manning[key]["mesh_correspondance"]])
    model.meshing.plot(my_scalar =values, title_plot=f"key = {key}")


for key in ["uniform", "1", "2", "3", "4"]:
    # get platch value
    values = np.array(
            [all_manning[key]["patch_value"][int(x-1)]     for x in     all_manning[key]["mesh_correspondance"]])
    model.meshing.plot(my_scalar =values, title_plot=f"key = {key}")


key = "true"
tmp = h5py.File(f"/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/hdf5_file/true_simu.hdf5")
all_manning["true"] = dict()
all_manning["true"]["mesh_correspondance"] = tmp["input/param/manning/mesh_correspondance"][...]
all_manning["true"]["patch_value"] = tmp["input/param/manning/patch_value"][:]
ref_values = np.array(
            [all_manning["true"]["patch_value"][int(x-1)]     for x in     all_manning[key]["mesh_correspondance"]])

model.meshing.plot(my_scalar =ref_values, title_plot=f"key = {key}")


all_values =[]
for key in ["uniform", "1", "2", "3", "4"]:
    # get platch value
    values = np.array(
            [all_manning[key]["patch_value"][int(x-1)]     for x in     all_manning[key]["mesh_correspondance"]])
    all_values.append(values)
    model.meshing.plot(my_scalar =values, title_plot=f"key = {key}")


import matplotlib.pyplot as plt
plt.plot(all_values[0]-ref_values, c= "red")
plt.plot(all_values[1]-ref_values, c= "yellow")
plt.plot(all_values[2]-ref_values, c = "purple")
plt.plot(all_values[3]-ref_values, c = "green")
plt.plot(all_values[4]-ref_values, c = "orange")



import matplotlib.pyplot as plt
plt.figure(figsize=(15, 10), dpi=80)
plt.plot(ref_values, c= "blue")
plt.plot(all_values[0], c= "red")
plt.plot(all_values[1], c= "yellow")
plt.plot(all_values[2], c = "purple")
plt.plot(all_values[3], c = "green")
plt.plot(all_values[4], c = "orange")


init_val = ref_values.copy()
init_val[:]  =  0.001
for i in range(len(all_values)-1):
    j=i+1
    if i==0 :        
        plt.figure(figsize=(15, 10), dpi=80)
        plt.plot(ref_values, c = "blue")
        plt.plot(init_val, c = "green")
        plt.plot(all_values[j], c = "red")
        plt.show()
    else:
        plt.figure(figsize=(15, 10), dpi=80)
        plt.plot(ref_values, c = "blue")
        plt.plot(all_values[i], c = "green")
        plt.plot(all_values[j], c = "red")
        plt.show()



init_val = ref_values.copy()
init_val[:]  =  0.001
for i in range(len(all_values)-1):
    j=i+1
    if i==0 :        
        plt.figure(figsize=(8, 5), dpi=80)
        plt.plot(ref_values[0:100], c = "blue")
        plt.plot(init_val[0:100], c = "green")
        plt.plot(all_values[i][0:100], c = "red")
        plt.show()
    else:
        plt.figure(figsize=(8, 5), dpi=80)
        plt.plot(ref_values[0:100], c = "blue")
        plt.plot(all_values[i-1][0:100], c = "green")
        plt.plot(all_values[i][0:100], c = "red")
        plt.show()
        
        
init_val = ref_values.copy()
init_val[:]  =  0.001
for i in range(len(all_values)-1):
    j=i+1
    if i==0 :        
        plt.figure(figsize=(8, 5), dpi=80)
        plt.plot(ref_values[400:500], c = "blue", label = "true")
        plt.plot(init_val[400:500], c = "green", label = "prior")
        plt.plot(all_values[i][400:500], c = "red", label = "inf")
        plt.legend()
        plt.show()
    else:
        plt.figure(figsize=(8, 5), dpi=80)
        plt.plot(ref_values[400:500], c = "blue", label = "true")
        plt.plot(all_values[i-1][400:500], c = "green", label = "prior")
        plt.plot(all_values[i][400:500], c = "red", label = "inf")
        plt.legend()
        plt.show()