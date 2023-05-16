#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:31:25 2023

@author: livillenave
"""

import dassflow2d as df2d
import numpy as np
import pandas as pd
import sys
import os
import random
from matplotlib.colors import ListedColormap

os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/cases/dev_case/12_mutiresolutioncase/organised_script/")
import libs


bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A"
case_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/cases/dev_case/12_mutiresolutioncase/bin_A"
os.chdir(bin_dir)
os.system(f"cp -a {case_dir}/* {bin_dir}/")



#-------------------------------------------
# THIRD : PREPARE PATCHES FOR INFERENCE
#-------------------------------------------
def half_cut(my_square):
    
    xmin, xmax = my_square["xmin"],my_square["xmax"]
    ymin, ymax = my_square["ymin"],my_square["ymax"]
    
    xmin1=xmin
    xmax1=(xmax+xmin)/2    
    xmin2=xmax1 
    xmax2=xmax
    
    ymin1=ymin
    ymax1=(ymax+ymin)/2    
    ymin2=ymax1 
    ymax2=ymax
    
    all_squares=dict()
    
    all_squares[1] =dict()
    all_squares[1]["xmin"],all_squares[1]["xmax"]=xmin1, xmax1
    all_squares[1]["ymin"],all_squares[1]["ymax"]=ymin1, ymax1
    
    
    all_squares[2] =dict()
    all_squares[2]["xmin"],all_squares[2]["xmax"]=xmin2, xmax2
    all_squares[2]["ymin"],all_squares[2]["ymax"]=ymin2, ymax2
    
    all_squares[3] =dict()
    all_squares[3]["xmin"],all_squares[3]["xmax"]=xmin1, xmax1
    all_squares[3]["ymin"],all_squares[3]["ymax"]=ymin2, ymax2
    
    all_squares[4] =dict()
    all_squares[4]["xmin"],all_squares[4]["xmax"]=xmin2, xmax2
    all_squares[4]["ymin"],all_squares[4]["ymax"]=ymin1, ymax1
    
    return( all_squares[1],  all_squares[2], all_squares[3], all_squares[4])


def cell_in_patch(x_coord, y_coord, id_cell, patch):
    x_coord=np.array(x_coord)
    y_coord=np.array(y_coord)
    
    x=x_coord[id_cell]
    y=y_coord[id_cell]
    
    for key, my_square in patch.items():
        if x< my_square["xmax"] and x >=my_square["xmin"] and y < my_square["ymax"] and y >= my_square["ymin"]:
            return(key)


###########################################
# subroutines for simpler plots
def random_color():
    levels = range(32,256,32)
    my_tuple = tuple(random.choice(levels)/256 for _ in range(3))
    my_tuple +=(1,)
    return my_tuple

def get_colormap(id_patch):
    newcolors = np.empty((len(id_patch), 4))
    for i in np.unique(id_patch):
        id_modif=np.where(id_patch == i)
        newcolors[id_modif] = np.array(random_color())
    
    my_colormap = ListedColormap(newcolors)
    return(my_colormap)
###########################################


def build_patches(x_coord, y_coord, xmin,xmax,ymin,ymax): 
        
    my_square = dict()
    my_square["xmin"],my_square["xmax"]=xmin, xmax
    my_square["ymin"],my_square["ymax"]=ymin, ymax
    
    patch_uniform = dict()
    patch_uniform[0] = my_square
    
    patch_lvl2    = dict()
    patch_lvl2[0] = {'xmin': 0,  'xmax': 100, 'ymin': 0, 'ymax': 100}
    patch_lvl2[1] = {'xmin': 100, 'xmax': 200, 'ymin': 0, 'ymax': 100}
    patch_lvl2[2] = {'xmin': 200, 'xmax':300, 'ymin': 0, 'ymax': 100}
    patch_lvl2[3] = {'xmin': 300, 'xmax': 400, 'ymin': 0, 'ymax': 100}
    patch_lvl2[4] = {'xmin': 400, 'xmax': 500, 'ymin': 0, 'ymax': 100}
    patch_lvl2[5] = {'xmin': 500, 'xmax': 600, 'ymin': 0, 'ymax': 100}
    patch_lvl2[6] = {'xmin': 600, 'xmax': 700, 'ymin': 0, 'ymax': 100}
    patch_lvl2[7] = {'xmin': 700, 'xmax': 800, 'ymin': 0, 'ymax': 100}
    patch_lvl2[8] = {'xmin': 800, 'xmax': 900, 'ymin': 0, 'ymax': 100}
    patch_lvl2[9] = {'xmin': 900, 'xmax': 1000, 'ymin': 0, 'ymax': 100}
    
    patch_lvl3 = dict()
    count=0
    for i in range(len(patch_lvl2)):
        sq1,sq2,sq3,sq4=half_cut(patch_lvl2[i])
        patch_lvl3[count]=sq1
        patch_lvl3[count+1]=sq2
        patch_lvl3[count+2]=sq3
        patch_lvl3[count+3]=sq4
        count = count+4
    
    
    patch_lvl4 = dict()
    count=0
    for i in range(len(patch_lvl3)):
        sq1,sq2,sq3,sq4=half_cut(patch_lvl3[i])
        patch_lvl4[count]=sq1
        patch_lvl4[count+1]=sq2
        patch_lvl4[count+2]=sq3
        patch_lvl4[count+3]=sq4
        count = count+4
    
    # store ids and plot
    all_patches = dict()
    
    # ----------
    # check uniform
    # ------------
    
    id_patch = []
    for id_cell in range(len(x_coord)):            
        id_patch.append(cell_in_patch(x_coord=x_coord, y_coord=y_coord, id_cell = id_cell, 
                                      patch = patch_uniform))
    
    id_patch = np.array(id_patch)
    model.meshing.plot(my_scalar = id_patch)
    
    model.meshing.mesh_pyvista.plot(scalars=id_patch+1, cpos = "xy", 
                                    show_edges=True, 
                                    notebook = True, show_bounds=True)
    
    all_patches[0] = id_patch
    # ----------
    # check lvl2
    # ------------
    
    id_patch = []
    for id_cell in range(len(x_coord)):            
        id_patch.append(cell_in_patch(x_coord=x_coord, y_coord=y_coord, id_cell = id_cell, 
                                      patch = patch_lvl2))
    
    id_patch = np.array(id_patch)
    #pv.plotter()
    model.meshing.mesh_pyvista.plot(scalars = id_patch, cpos = "xy", 
                                    show_edges=True, 
                                    notebook = True)
    
    all_patches[1] = id_patch
    
    # ------------
    # check lvl3
    # ------------
    id_patch = []
    for id_cell in range(len(x_coord)):            
        tmp_id_patch = cell_in_patch(x_coord=x_coord, y_coord=y_coord, id_cell = id_cell, 
                                      patch = patch_lvl3)
        print(tmp_id_patch)
        id_patch.append(tmp_id_patch)
    
    id_patch = np.array(id_patch)
    
    
    my_colormap = get_colormap(id_patch)
    
    model.meshing.mesh_pyvista.plot(scalars=id_patch,  
                                    cpos="xy", 
                                    cmap = my_colormap,
                                    show_edges=True, 
                                    notebook=True)
    
    all_patches[2] = id_patch
    # ----------
    # check lvl4
    # ------------
    
    id_patch = []
    for id_cell in range(len(x_coord)):            
        id_patch.append(cell_in_patch(x_coord=x_coord, y_coord=y_coord, id_cell = id_cell, 
                                      patch = patch_lvl4))
    
    id_patch = np.array(id_patch)
    
    
    my_colormap = get_colormap(id_patch)
    
    
    model.meshing.mesh_pyvista.plot(scalars=id_patch, 
                                    show_edges=True, 
                                    cpos="xy", 
                                    cmap = my_colormap, 
                                    notebook=True)
    
    
    
    all_patches[3] = id_patch
    
    
    # ----------
    # cell level
    # ----------
    all_patches[4] = np.arange(len(x_coord), step =1)
    my_colormap = get_colormap(all_patches[4])
    #model.meshing.mesh_pyvista.plot(scalars=all_patches[4], 
    #                                show_edges=True, 
    #                                cpos="xy", 
    #                                cmap = my_colormap, 
    #                                notebook=False)

    # set 1 as minimum value
    for i,patch in all_patches.items():
        patch = patch+1
        all_patches[i] = patch
    return(all_patches)


def all_patches_dic_2_ndarray(all_patches):
    
    tmp = np.ndarray(shape = (len(all_patches[0]), 
                              len(all_patches)+1), 
                              dtype = 'int')
    
    tmp[:,0] = np.arange(len(all_patches[0]))+1
    
    for i in range(len(all_patches)):
            tmp[:,i+1] = all_patches[i][:]
    
    reorder=  np.argsort(tmp[:,4])
    tmp = tmp[reorder,:]
    return(tmp)
    
    

#-------------------------------------------
# FIRST : CHECK AVAILABLE CASE
# Load python data
#-------------------------------------------
run_type = "direct"
# initialise fortran instance, and python corrponding data
model = df2d.dassflowmodel(bin_dir = bin_dir, 
                           hdf5_path= bin_dir +"/res/simu.hdf5", 
                           run_type = run_type,
                           clean=True)

# SET UP MANNING VALUE
new_config = model.config.get()
new_config["use_obs"] = 0
new_config["w_obs"] = 1
model.config.set(new_config)
model.init_all()
model.run()
model.save_all()


df2d.wrapping.call_model.clean_model(self = model.kernel)
reference_manning = model.param.friction.manning.copy()

df2d.wrapping.call_model.clean_model(self = model.kernel)



#-------------------------------------------
# SECOND : GENERATE  TRUE RUN
#-------------------------------------------

x_coord = model.outputs.result.x.copy()
y_coord = model.outputs.result.y.copy()

# INPUTS
run_type = "direct"
true_mesh_corresp = np.arange(1, len(reference_manning["mesh_correspondance"])+1,
                              1, dtype= "int32" )
true_values = 0.01 + np.array(x_coord)/20_000 #+ abs(50 - np.array(y_coord)) / 10_000
# random.random() generate a number between 0 and 1


def manning_de_x(x_value):
    return(0.01 + x_value/20_000)

# initialise fortran instance, and python corrponding data
model = df2d.dassflowmodel(bin_dir = os.getcwd(), 
                           hdf5_path=os.getcwd()+"/res/simu.hdf5", 
                           run_type = run_type,
                           clean=True)

# SET UP MANNING VALUE
model.config["use_obs"] = 0
model.config["w_obs"] = 1
model.config.set()

model.init_all()

model.kernel.dof0.h[...]= model.kernel.dof.h[...] = 1
model.kernel.dof0.u[...]= model.kernel.dof.u[...] = 1
model.kernel.dof0.v[...]= model.kernel.dof.v[...] = 1

model.param.friction.manning["mesh_correspondance"] = true_mesh_corresp
model.param.friction.manning["patch_correspondance"] = true_mesh_corresp
model.param.friction.manning["patch_value"] = true_values
model.param.set()
model.run()

model.save_all()

df2d.wrapping.call_model.clean_model(self = model.kernel)

#
os.chdir(model.bin_dir)
os.system(f"mkdir obs")
os.system(f"cp ./res/obs/* ./obs")

os.system(f"mkdir hdf5_file")
os.system(f"mv ./res/simu.hdf5 ./hdf5_file/true_simu.hdf5")






x_coord = model.outputs.result.x.copy()
y_coord = model.outputs.result.y.copy()

xmin, xmax=0,1000
ymin, ymax=0,100


all_patches = build_patches(x_coord,y_coord,xmin,xmax,ymin,ymax)
all_patch_as_ndarray = all_patches_dic_2_ndarray(all_patches)


reorder=  np.argsort(all_patch_as_ndarray[:,0])
all_patch_as_ndarray[reorder,:]
matrixes = dict()

for i in range(len(all_patches)):
    if( i == 0):
        # i = 0 is column of id cell --> do nothing
        print("")
    elif( i == 1):
        # i = 1 is fist patch (uniforme)
        matrixes[i-1] = libs.build_p_matrix(
                tmp= all_patch_as_ndarray, 
                previous_patch = i, 
                new_patch = i+1)
    else:
        matrixes[i-1] = libs.build_p_matrix(
                tmp= all_patch_as_ndarray, 
                previous_patch = i, 
                new_patch = i+1)
        
#with np.printoptions(threshold=np.inf):
#    print(a)

import matplotlib.pyplot as plt
for i in range(1, len(matrixes)):
    plt.figure( figsize=(15,15))
    plt.imshow(matrixes[i]);plt.show(); plt.close()
    
def flatten(l):
    return [item for sublist in l for item in sublist]

id_reorder = []
for i in range(matrixes[3].shape[1]):
    array_to_match = np.zeros(shape =matrixes[3].shape[1])
    array_to_match[i]=1
    
    id_to_add=[]
    for j in range(matrixes[3].shape[0]):
        if sum(matrixes[3][j,:]) !=  1:
            print("problem")
        if np.all(matrixes[3][j,:] ==array_to_match):
            id_to_add.append(j)
    id_reorder.append(id_to_add)

test = matrixes[3][np.array(flatten(id_reorder)),:]


cmap1=ListedColormap(["gold"])
cmap = ListedColormap(["darkblue","gold"])
fig, axs = plt.subplots(1, 4,figsize=(8*2, 6*2))
axs[0].imshow(matrixes[0], cmap =cmap1, aspect = "auto")
axs[1].imshow(matrixes[1], cmap =cmap, aspect = "auto")
axs[2].imshow(matrixes[2], cmap =cmap, aspect = "auto")
axs[3].imshow(test, cmap =cmap, aspect = "auto")
axs[3].legend()

axs[0].set_title('uniform to 10 patches ')
axs[1].set_title('10 to 40 patches ')
axs[2].set_title('40 to 160 patches ')
axs[3].set_title('160 to 891 patches')

#plt.yticks(np.array(flatten(id_reorder))[::6], rotation = 45, labelsize=3)


def flatten(l):
    return [item for sublist in l for item in sublist]

# ======================================================================================= #
#  ID PATCH == 1
# ======================================================================================= #

bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A"


all_min = dict()
all_best_patch= dict()

for my_id_patch in range(len(all_patches)):
    
    if my_id_patch==0:        
        #-------------------------------------------
        # FORTH : Perform inference uniform
        #-------------------------------------------
        run_type = "min"
        
        patch_value =  np.array(0.001,ndmin = 1)
        
        #
        # >>> Generate corresponding case
        #
        mesh_dict = libs.read_mesh_from_textfile(path = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/channel.geo", read_boundary = True )

        
        libs.rewrite_friction_dassflow(bin_dir = bin_dir, 
                                  mesh_name =  "channel.geo", 
                                  land_use_name = "land_uses.txt", 
                                  mesh_dict = mesh_dict, 
                                  mesh_corresp = all_patches[0], 
                                  land_use_value = patch_value
                                  )
        
        
        
        
        # initialise fortran instance, and python corrponding data
        model = df2d.dassflowmodel(bin_dir =bin_dir, 
                                   hdf5_path=bin_dir+"/res/simu.hdf5", 
                                   run_type = run_type,
                                   clean=True)
        
        # SET UP MANNING VALUE
        model.config["use_obs"]   = 1
        model.config["w_obs"]     = 1
        model.config["c_manning"] = 1
        model.config["c_hydrograph"] = 0
        model.config["eps_min"] = 0.001
        model.config.set()
        
        model.init_all()
        
        
        model.kernel.dof0.h[...] =  1
        model.kernel.dof0.u[...] =  1
        model.kernel.dof0.v[...] =  1
        model.kernel.dof.h[...]  =  1
        model.kernel.dof.u[...]  =  1
        model.kernel.dof.v[...]  =  1
        
        model.run()
        
        (optimized_patch_value, optimized_patch_corresp) = libs.source_min_res(min_dir = model.bin_dir + "/min")   
        
        
        all_min[my_id_patch] = model.min
        all_best_patch[my_id_patch] = optimized_patch_value
        
        df2d.wrapping.call_model.clean_model(self = model.kernel)
        
        
        # ----------------------------------------------
        # perform direct run with optimized values
        # -----------------------------------------
        run_type = "direct"
        
        # initialise fortran instance, and python corrponding data
        model = df2d.dassflowmodel(bin_dir = os.getcwd(), 
                                   hdf5_path=os.getcwd()+"/res/simu.hdf5", 
                                   run_type = run_type,
                                   clean=True)
        
        # SET UP MANNING VALUE
        new_config = model.config.get()
        new_config["use_obs"]   = 0
        new_config["w_obs"]     = 1
        # SET UP MANNING VALUE
        model.config.set(new_config)        
        
        libs.rewrite_friction_dassflow(bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A", 
                                  mesh_name =  "channel.geo", 
                                  land_use_name = "land_uses.txt", 
                                  mesh_dict = mesh_dict, 
                                  mesh_corresp = all_patches[0], 
                                  land_use_value =  optimized_patch_value)
        model.init_all()
        
        
        model.kernel.dof0.h[...] =  1
        model.kernel.dof0.u[...] =  1
        model.kernel.dof0.v[...] =  1
        model.kernel.dof.h[...]  =  1
        model.kernel.dof.u[...]  =  1
        model.kernel.dof.v[...]  =  1
        
        model.run()        
        model.save_all()        
        ref_manning =  model.param.friction.manning.copy()        
        df2d.wrapping.call_model.clean_model(self = model.kernel)
        
        os.chdir(model.bin_dir)
        os.system(f"mv ./res/simu.hdf5 ./hdf5_file/inf_uniform.hdf5")
    else:        
        
        del model
        print("---------------------------------------------------------------------")       
        print( my_id_patch )
        print("---------------------------------------------------------------------")            
        
        run_type = "min"        
        optimized_patch_values = optimized_patch_value       
        
        libs.rewrite_friction_dassflow(bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A", 
                                  mesh_name =  "channel.geo", 
                                  land_use_name = "land_uses.txt", 
                                  mesh_dict = mesh_dict, 
                                  mesh_corresp = all_patches[my_id_patch], 
                                  land_use_value = matrixes[my_id_patch-1] @ optimized_patch_value                               
                                  )       
        
        # initialise fortran instance, and python corrponding data
        model = df2d.dassflowmodel(bin_dir = os.getcwd(), 
                                   hdf5_path=os.getcwd()+"/res/simu.hdf5", 
                                   run_type = run_type,
                                   clean=True)        
        
        # SET UP MANNING VALUE
        new_config = model.config.get()
        new_config["use_obs"]   = 1
        new_config["w_obs"]     = 1
        new_config["c_manning"] = 1
        new_config["eps_min"] = 0.0001
        model.config.set(new_config)
         
        model.init_all()      
          
        model.kernel.dof0.h[...] =  1
        model.kernel.dof0.u[...] =  1
        model.kernel.dof0.v[...] =  1
        model.kernel.dof.h[...]  =  1
        model.kernel.dof.u[...]  =  1
        model.kernel.dof.v[...]  =  1        
        
        model.run()        
        (optimized_patch_value, optimized_patch_corresp) = libs.source_min_res(min_dir = model.bin_dir + "/min")        
        
        all_min[my_id_patch] = model.min
        all_best_patch[my_id_patch] = optimized_patch_value
        
        # ----------------------------------------------
        # perform direct run with optimized values
        # -----------------------------------------
        
        run_type = "direct"
        
        libs.rewrite_friction_dassflow(bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A", 
                                  mesh_name =  "channel.geo", 
                                  land_use_name = "land_uses.txt", 
                                  mesh_dict = mesh_dict, 
                                  mesh_corresp = all_patches[my_id_patch], 
                                  land_use_value = optimized_patch_value                               
                                  )
        
        # initialise fortran instance, and python corrponding data
        model = df2d.dassflowmodel(bin_dir = os.getcwd(), 
                                   hdf5_path=os.getcwd()+"/res/simu.hdf5", 
                                   run_type = run_type,
                                   clean=True)
        
        # SET UP MANNING VALUE
        new_config = model.config.get()
        new_config["use_obs"]   = 0
        new_config["w_obs"]     = 1
        new_config["c_manning"] = 1
        new_config["eps_min"] = 0.001
        model.config.set(new_config)
        
        model.init_all()
        
        
        model.kernel.dof0.h[...] =  1
        model.kernel.dof0.u[...] =  1
        model.kernel.dof0.v[...] =  1
        model.kernel.dof.h[...]  =  1
        model.kernel.dof.u[...]  =  1
        model.kernel.dof.v[...]  =  1
        
        
        model.run()
        model.save_all()
        
        df2d.wrapping.call_model.clean_model(self = model.kernel)
        #
        os.chdir(model.bin_dir)
        os.system(f"mv ./res/simu.hdf5 ./hdf5_file/inf{my_id_patch}.hdf5")

      
fig, axs = plt.subplots(1, 5, figsize = (15,5))

axs[0].plot(all_min[0].all_ite, all_min[0].j     ,'b',linewidth=1,label=r'$J_{hy}$')
my_colors = ["r-", "y-", "o-", "g-"]
for i in range(all_min[0].gradj.shape[1]) :
    print(i)
    tex = r"$|| \nabla_{" +fr"{all_min[0].gradj_name[i]}" + "} J_{hy}||$"
    tex = fr"$|| \nabla {{{all_min[0].gradj_name[i]}}} ||$"
    # name update necesary soon
    axs[0].plot(all_min[0].all_ite,all_min[0].gradj[:,i],my_colors[i],linewidth=0.5, label= rf"{tex}")

#plt.plot(ite,grad_cost2,'r',linewidth=1,label=r'$||\nabla_{Q_2} J_{hy}||$')

#plt.xscale('log')
axs[0].set_yscale('log')
#           


axs[1].plot(all_min[1].all_ite, all_min[1].j     ,'b',linewidth=1,label=r'$J_{hy}$')
my_colors = ["r-", "y-", "o-", "g-"]
for i in range(all_min[1].gradj.shape[1]) :
    print(i)
    tex = r"$|| \nabla_{" +fr"{all_min[1].gradj_name[i]}" + "} J_{hy}||$"
    tex = fr"$|| \nabla {{{all_min[1].gradj_name[i]}}} ||$"
    # name update necesary soon
    axs[1].plot(all_min[1].all_ite,all_min[1].gradj[:,i],my_colors[i],linewidth=0.5, label= rf"{tex}")
#axs[1].legend()
#plt.plot(ite,grad_cost2,'r',linewidth=1,label=r'$||\nabla_{Q_2} J_{hy}||$')

#plt.xscale('log')
axs[1].set_yscale('log')
#           


axs[2].plot(all_min[2].all_ite, all_min[2].j     ,'b',linewidth=1,label=r'$J_{hy}$')
my_colors = ["r-", "y-", "o-", "g-"]
for i in range(all_min[2].gradj.shape[1]) :
    print(i)
    tex = r"$|| \nabla_{" +fr"{all_min[2].gradj_name[i]}" + "} J_{hy}||$"
    tex = fr"$|| \nabla {{{all_min[2].gradj_name[i]}}} ||$"
    # name update necesary soon
    axs[2].plot(all_min[2].all_ite,all_min[2].gradj[:,i],my_colors[i],linewidth=0.5, label= rf"{tex}")
#axs[2].legend()
#plt.plot(ite,grad_cost2,'r',linewidth=1,label=r'$||\nabla_{Q_2} J_{hy}||$')
axs[2].set_xlabel('Iterations')
#plt.xscale('log')
axs[2].set_yscale('log')
#           



axs[3].plot(all_min[3].all_ite, all_min[3].j     ,'b',linewidth=1,label=r'$J_{hy}$')
my_colors = ["r-", "y-", "o-", "g-"]
for i in range(all_min[3].gradj.shape[1]) :
    print(i)
    tex = r"$|| \nabla_{" +fr"{all_min[3].gradj_name[i]}" + "} J_{hy}||$"
    tex = fr"$|| \nabla {{{all_min[3].gradj_name[i]}}} ||$"
    # name update necesary soon
    axs[3].plot(all_min[3].all_ite,all_min[3].gradj[:,i],my_colors[i],linewidth=0.5, label= rf"{tex}")
#axs[3].legend()
#plt.plot(ite,grad_cost2,'r',linewidth=1,label=r'$||\nabla_{Q_2} J_{hy}||$')
axs[3].set_xlabel('Iterations')
#plt.xscale('log')
axs[3].set_yscale('log')
#         

axs[3].legend(loc='upper center', 
             bbox_to_anchor=(0.0, -0.4),fancybox=False, shadow=False, ncol=3)  




axs[4].plot(all_min[4].all_ite, all_min[4].j     ,'b',linewidth=1,label=r'$J_{hy}$')
my_colors = ["r-", "y-", "o-", "g-"]
for i in range(all_min[4].gradj.shape[1]) :
    print(i)
    tex = r"$|| \nabla_{" +fr"{all_min[4].gradj_name[i]}" + "} J_{hy}||$"
    tex = fr"$|| \nabla {{{all_min[4].gradj_name[i]}}} ||$"
    # name update necesary soon
    axs[4].plot(all_min[4].all_ite,all_min[4].gradj[:,i],my_colors[i],linewidth=0.5, label= rf"{tex}")
#axs[4].legend()
#plt.plot(ite,grad_cost2,'r',linewidth=1,label=r'$||\nabla_{Q_2} J_{hy}||$')
axs[4].set_xlabel('Iterations')
#plt.xscale('log')
axs[4].set_yscale('log')
#         

axs[4].legend(loc='upper center', 
             bbox_to_anchor=(0.0, -0.4),fancybox=False, shadow=False, ncol=3)  


axs[0].set_title('Patch uniform')
axs[1].set_title('Patch 10 values')
axs[2].set_title('Patch 40 values')
axs[3].set_title('Patch 160 values')
axs[4].set_title('Patch 891 values')
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.8,
                    hspace=0.4)


plt.show()
            




        # initialise fortran instance, and python corrponding data
model = df2d.dassflowmodel(bin_dir = os.getcwd(), 
                           hdf5_path=os.getcwd()+"/res/simu.hdf5", 
                           run_type = run_type,
                           clean=False)        
# SET UP MANNING VALUE
new_config = model.config.get()
new_config["use_obs"]   = 0
new_config["w_obs"]     = 1
new_config["c_manning"] = 1
new_config["eps_min"] = 0.001
model.config.set(new_config)        
model.init_all()



from matplotlib import cm
n_colors = 160
colours = cm.rainbow(np.linspace(0, 1, n_colors))
rng = np.random.default_rng()
colours = rng.permuted(colours, axis=0)
colours = rng.permuted(colours, axis=0)
mapping = np.linspace(all_patches[3].min(), all_patches[3].max(), 256)
newcolors = np.empty((256, 4))
for i in range(len(np.unique(all_patches[3]))):
    newcolors[ (mapping > i-1) ] = colours[i-1]
# Make the colormap from the listed colors
my_colormap = ListedColormap(newcolors)
model.meshing.mesh_pyvista.plot(scalars=all_patches[3], cpos = "xy", 
                                    show_edges=True, cmap =my_colormap,
                                    notebook = True, show_bounds=True)
import matplotlib.pyplot as plt
for i in range(1, len(matrixes)):
    plt.figure( figsize=(15,15))
    plt.imshow(matrixes[i]);plt.show(); plt.close()   
    

xcell = np.zeros(shape = (model.meshing.mesh_fortran.nc))
ycell = np.zeros(shape = (model.meshing.mesh_fortran.nc))
for i in range(model.meshing.mesh_fortran.nc):
    xcell[i]=model.meshing.mesh_fortran.cell[i].grav.x
    ycell[i]=model.meshing.mesh_fortran.cell[i].grav.y

def get_center_patch(mesh_corresp, xcell):
    id_patch = np.unique(mesh_corresp)
    x_center_patch = np.zeros(shape = (len(id_patch)))
    for i in id_patch:
        id_cell_in_patch = np.where(mesh_corresp == i)
        x_center_patch[i-1] = np.mean(xcell[id_cell_in_patch])
    return(x_center_patch)
    
 
def get_categorical_center_patch(mesh_corresp, ycell):
    id_patch = np.unique(mesh_corresp)
    y_center_patch = np.zeros(shape = (len(id_patch)))
    for i in id_patch:
        id_cell_in_patch = np.where(mesh_corresp == i)
        print(i, id_cell_in_patch)
        y_center_patch[i-1] =  np.mean(ycell[id_cell_in_patch])
        
    patch_line_value = np.unique(y_center_patch.round(decimals=4))
    line_id = dict()
    
    np.where(y_center_patch.round(decimals=4) - patch_line_value[1] < 0.0001)
    return(y_center_patch)   
    
    
all_x_center = dict()
all_y_center = dict()
for my_key in all_patches.keys():
    all_x_center[my_key] = get_center_patch(all_patches[my_key], xcell)
    all_y_center[my_key] =  get_center_patch(all_patches[my_key], ycell)
    
    
    
fig, axs = plt.subplots(5, 1,figsize=(8*2, 6*2))

my_key=0
xtrue=all_x_center[4]
true=true_values
inf = 0
prior =  np.array(0.001,ndmin = 1)
axs[my_key].plot(xtrue, true_values, color = "blue")
axs[my_key].plot( [0,1000],np.repeat(prior, 2),  color = "green", linestyle='dashed' )
axs[my_key].plot([0,1000],  np.repeat(all_best_patch[my_key],2), color = "red")

my_key=1
prior = matrixes[my_key-1] @ all_best_patch[my_key-1]
plt.figure(figsize=(4,12))
axs[my_key].plot(xtrue, true_values, color = "blue")
axs[my_key].plot(all_x_center[my_key],prior,  color = "green", linestyle='dashed' )
axs[my_key].plot(all_x_center[my_key],  all_best_patch[my_key],color = "red")


my_key=2
to_keep = np.where(  abs(all_y_center[my_key] - np.unique(all_y_center[my_key])[0] ) < 0.01)[0]
prior = matrixes[my_key-1] @ all_best_patch[my_key-1]


new_x =all_x_center[my_key][to_keep]
new_prior =  prior[to_keep]
new_patch = all_best_patch[my_key][to_keep]

id_reorder = np.argsort(new_x)
new_x = new_x[id_reorder]
new_prior =  new_prior[id_reorder]
new_patch = new_patch[id_reorder]

axs[my_key].plot(xtrue, true_values, color = "blue")
axs[my_key].plot(new_x, 
                 new_prior, 
                 color = "green", linestyle='dashed' )
axs[my_key].plot( new_x, 
                 new_patch, 
                 color = "red", label = "first guess")

my_key=3
to_keep = np.where( abs(all_y_center[my_key] - 33.333) < 0.001)[0]
prior = matrixes[my_key-1] @ all_best_patch[my_key-1]

new_x =all_x_center[my_key][to_keep]
new_prior =  prior[to_keep]
new_patch = all_best_patch[my_key][to_keep]

id_reorder = np.argsort(new_x)
new_x = new_x[id_reorder]
new_prior =  new_prior[id_reorder]
new_patch = new_patch[id_reorder]

axs[my_key].plot(xtrue, true_values, color = "blue")
axs[my_key].plot(new_x, 
                 new_prior, 
                 color = "green", linestyle='dashed' )
axs[my_key].plot( new_x, 
                 new_patch, 
                 color = "red", label = "first guess")

my_key=4
to_keep = np.where( abs(all_y_center[my_key] - 50.) < 0.01)[0]
prior = matrixes[my_key-1] @ all_best_patch[my_key-1]


new_x =all_x_center[my_key][to_keep]
new_prior =  prior[to_keep]
new_patch = all_best_patch[my_key][to_keep]

id_reorder = np.argsort(new_x)
new_x = new_x[id_reorder]
new_prior =  new_prior[id_reorder]
new_patch = new_patch[id_reorder]

axs[my_key].plot(xtrue, true_values, color = "blue")
axs[my_key].plot(new_x, 
                 new_prior, 
                 color = "green", linestyle='dashed' )
axs[my_key].plot( new_x, 
                 new_patch, 
                 color = "red", label = "first guess")

axs[0].set_title('Patch uniform')
axs[1].set_title('Patch 10 values')
axs[2].set_title('Patch 40 values')
axs[3].set_title('Patch 160 values')
axs[4].set_title('Patch 891 values')
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

axs[2].legend(loc='upper center', 
             bbox_to_anchor=(0.1, -0.1),fancybox=False, shadow=False, ncol=3)  



########




fig, axs = plt.subplots(5, 1,figsize=(8*2, 6*2))

my_key=0
xtrue=all_x_center[4]
true=true_values
inf = 0
new_true = np.zeros(shape = (2))
compteur = -1
for my_x in [0,1000]:
    compteur = compteur +1
    new_true[compteur] = manning_de_x(my_x)
prior =  np.array(0.001,ndmin = 1)
axs[my_key].plot(xtrue, true_values-true_values, color = "blue")
axs[my_key].scatter( [0,1000],np.repeat(prior, 2)-new_true,  color = "green", linestyle='dashed' )
axs[my_key].plot([0,1000],  np.repeat(all_best_patch[my_key],2)-new_true, color = "red")

my_key=1
new_x = all_x_center[my_key]
prior = matrixes[my_key-1] @ all_best_patch[my_key-1]


new_true = np.zeros(shape = new_x.shape)
compteur = -1
for my_x in new_x:
    compteur = compteur +1
    new_true[compteur] = manning_de_x(my_x)

axs[my_key].plot(xtrue, true_values - true_values, color = "blue")
axs[my_key].scatter(all_x_center[my_key],prior - new_true,  color = "green", linestyle='dashed' )
axs[my_key].plot(all_x_center[my_key],  all_best_patch[my_key]- new_true,
   color = "red")


my_key=2
to_keep = np.where(  abs(all_y_center[my_key] - np.unique(all_y_center[my_key])[0] ) < 0.01)[0]
prior = matrixes[my_key-1] @ all_best_patch[my_key-1]


new_x =all_x_center[my_key][to_keep]
new_prior =  prior[to_keep]
new_patch = all_best_patch[my_key][to_keep]

id_reorder = np.argsort(new_x)
new_x = new_x[id_reorder]
new_prior =  new_prior[id_reorder]
new_patch = new_patch[id_reorder]


new_true = np.zeros(shape = new_x.shape)
compteur = -1
for my_x in new_x:
    compteur = compteur +1
    new_true[compteur] = manning_de_x(my_x)


axs[my_key].plot(xtrue, true_values-true_values, color = "blue")
axs[my_key].scatter(new_x, 
                 new_prior - new_true, 
                 color = "green", linestyle='dashed' )
axs[my_key].plot( new_x, 
                 new_patch - new_true, 
                 color = "red", label = "first guess")

my_key=3
to_keep = np.where( abs(all_y_center[my_key] - 33.333) < 0.001)[0]
prior = matrixes[my_key-1] @ all_best_patch[my_key-1]

new_x =all_x_center[my_key][to_keep]
new_prior =  prior[to_keep]
new_patch = all_best_patch[my_key][to_keep]

id_reorder = np.argsort(new_x)
new_x = new_x[id_reorder]
new_prior =  new_prior[id_reorder]
new_patch = new_patch[id_reorder]


new_true = np.zeros(shape = new_x.shape)
compteur = -1
for my_x in new_x:
    compteur = compteur +1
    new_true[compteur] = manning_de_x(my_x)


axs[my_key].plot(xtrue, true_values-true_values, color = "blue")
axs[my_key].scatter(new_x, 
                 new_prior - new_true, 
                 color = "green", linestyle='dashed' )
axs[my_key].plot( new_x, 
                 new_patch - new_true, 
                 color = "red", label = "first guess")

my_key=4
to_keep = np.where( abs(all_y_center[my_key] - 50.) < 0.01)[0]
prior = matrixes[my_key-1] @ all_best_patch[my_key-1]


new_x =all_x_center[my_key][to_keep]
new_prior =  prior[to_keep]
new_patch = all_best_patch[my_key][to_keep]

id_reorder = np.argsort(new_x)
new_x = new_x[id_reorder]
new_prior =  new_prior[id_reorder]
new_patch = new_patch[id_reorder]


new_true = np.zeros(shape = new_x.shape)
compteur = -1
for my_x in new_x:
    compteur = compteur +1
    new_true[compteur] = manning_de_x(my_x)


axs[my_key].plot(xtrue, true_values-true_values, color = "blue")
axs[my_key].scatter(new_x, 
                 new_prior - new_true, 
                 color = "green", linestyle='dashed' )
axs[my_key].plot( new_x, 
                 new_patch - new_true, 
                 color = "red", label = "first guess")

axs[0].set_title('Patch uniform')
axs[1].set_title('Patch 10 values')
axs[2].set_title('Patch 40 values')
axs[3].set_title('Patch 160 values')
axs[4].set_title('Patch 891 values')
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

axs[2].legend(loc='upper center', 
             bbox_to_anchor=(0.1, -0.1),fancybox=False, shadow=False, ncol=3)  






fig, axs = plt.subplots(5, 1,figsize=(8*2, 6*2))

my_key=0
xtrue=all_x_center[4]
true=true_values
inf = 0
new_true = np.zeros(shape = (2))

prior =  np.array(0.001,ndmin = 1)
axs[my_key].plot(xtrue, true_values, color = "blue")
axs[my_key].scatter( [0,1000],np.repeat(prior, 2),  color = "green", linestyle='dashed' )
axs[my_key].plot([0,1000],  np.repeat(all_best_patch[my_key],2), color = "red")

my_key=1
new_x = all_x_center[my_key]
prior = matrixes[my_key-1] @ all_best_patch[my_key-1]

axs[my_key].plot(xtrue, true_values , color = "blue")
axs[my_key].scatter(all_x_center[my_key],prior ,  color = "green", linestyle='dashed' )
axs[my_key].plot(all_x_center[my_key],  all_best_patch[my_key],
   color = "red")


my_key=2
to_keep = np.where(  abs(all_y_center[my_key] - np.unique(all_y_center[my_key])[0] ) < 0.01)[0]
prior = matrixes[my_key-1] @ all_best_patch[my_key-1]

new_x =all_x_center[my_key][to_keep]
new_prior =  prior[to_keep]
new_patch = all_best_patch[my_key][to_keep]

id_reorder = np.argsort(new_x)
new_x = new_x[id_reorder]
new_prior =  new_prior[id_reorder]
new_patch = new_patch[id_reorder]

axs[my_key].plot(xtrue, true_values, color = "blue")
axs[my_key].scatter(new_x, 
                 new_prior, 
                 color = "green", linestyle='dashed' )
axs[my_key].plot( new_x, 
                 new_patch , 
                 color = "red", label = "first guess")

my_key=3
to_keep = np.where( abs(all_y_center[my_key] - 33.333) < 0.001)[0]
prior = matrixes[my_key-1] @ all_best_patch[my_key-1]

new_x =all_x_center[my_key][to_keep]
new_prior =  prior[to_keep]
new_patch = all_best_patch[my_key][to_keep]

id_reorder = np.argsort(new_x)
new_x = new_x[id_reorder]
new_prior =  new_prior[id_reorder]
new_patch = new_patch[id_reorder]

axs[my_key].plot(xtrue, true_values, color = "blue")
axs[my_key].scatter(new_x, 
                 new_prior , 
                 color = "green", linestyle='dashed' )
axs[my_key].plot( new_x, 
                 new_patch, 
                 color = "red", label = "first guess")

my_key=4
to_keep = np.where( abs(all_y_center[my_key] - 50.) < 0.01)[0]
prior = matrixes[my_key-1] @ all_best_patch[my_key-1]


new_x =all_x_center[my_key][to_keep]
new_prior =  prior[to_keep]
new_patch = all_best_patch[my_key][to_keep]

id_reorder = np.argsort(new_x)
new_x = new_x[id_reorder]
new_prior =  new_prior[id_reorder]
new_patch = new_patch[id_reorder]


axs[my_key].plot(xtrue, true_values, color = "blue")
axs[my_key].scatter(new_x, 
                 new_prior, 
                 color = "green", linestyle='dashed' )
axs[my_key].plot( new_x, 
                 new_patch, 
                 color = "red", label = "first guess")

axs[0].set_title('Patch uniform')
axs[1].set_title('Patch 10 values')
axs[2].set_title('Patch 40 values')
axs[3].set_title('Patch 160 values')
axs[4].set_title('Patch 891 values')
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

axs[2].legend(loc='upper center', 
             bbox_to_anchor=(0.1, -0.1),fancybox=False, shadow=False, ncol=3)  

########


#######


model.meshing.mesh_pyvista.plot(scalars = true_values, cpos = "xy", 
                                show_edges=True, 
                                notebook = True)


    
black = np.array([11 / 256, 11 / 256, 11 / 256, 1.0])
red = np.array([1.0, 0.0, 0.0, 1.0])
orange = np.array([255/256, 153/256, 51/256, 1.0])
yellow = np.array([255 / 256, 247 / 256, 0 / 256, 1.0])
grey = np.array([189 / 256, 189 / 256, 189 / 256, 1.0])
green = np.array([2 / 256, 255 / 256, 0 / 256, 1.0])
blue = np.array([0, 0, 204/256, 1.0])

mapping = np.array( [0,1,5,10,15,20,30], dtype = "float64")

newcolors = np.empty((7, 4))
newcolors[np.where(mapping >= 30)] = black
newcolors[np.where(mapping < 30)] = red
newcolors[mapping < 20] = orange
newcolors[mapping < 15] = yellow
newcolors[mapping < 10] = grey
newcolors[mapping < 5] = blue
newcolors[mapping < 1] = green

# Make the colormap from the listed colors
my_colormap = ListedColormap(newcolors)
    
    
for i in [0,1,2,3,4]:
    cell_value = all_patches[i].copy()
    cell_value = np.array(cell_value, dtype = "float64")
    
    for j in range(len(np.unique( all_patches[i]))):
        cell_value[np.where(all_patches[i]==j+1)[0]] = all_best_patch[i][j]
    
    model.meshing.mesh_pyvista.plot(scalars = cell_value-true_values, cpos = "xy", 
                                    show_edges=True, 
                                    notebook = True)



for i in [0,1,2,3,4]:
    cell_value = all_patches[i].copy()
    cell_value = np.array(cell_value, dtype = "float64")
    
    for j in range(len(np.unique( all_patches[i]))):
        cell_value[np.where(all_patches[i]==j+1)[0]] = all_best_patch[i][j]
    
    categorical = 100 * (cell_value-true_values)/true_values
    abs_cate = abs(categorical)
    new_scalar = abs_cate.copy()
    new_scalar[np.where(abs_cate >= 30)] = 30
    new_scalar[np.where(abs_cate < 30)] = 20
    new_scalar[np.where(abs_cate < 20)] = 15
    new_scalar[np.where(abs_cate < 15)] = 10
    new_scalar[np.where(abs_cate < 10)] = 5
    new_scalar[np.where(abs_cate < 5)] = 1
    new_scalar[np.where(abs_cate < 1)] = 0
    
    cmap_colors = []
    category=[0.,1.,5.,10.,15.,20.,30.]
    cols =  ["green", "blue", "grey","yellow", "orange", "red", "black"]
    compteur = -1
    for my_test in category:
        compteur = compteur+1
        if my_test in np.unique(new_scalar):
            cmap_colors.append(cols[compteur])
        

    
    
    model.meshing.mesh_pyvista.plot(scalars = new_scalar, 
                                    cpos = "xy", 
                                    cmap = cmap_colors,
                                    show_edges=True, 
                                    notebook = True)






for i in [0,1,2,3,4]:
    cell_value = all_patches[i].copy()
    cell_value = np.array(cell_value, dtype = "float64")
    
    for j in range(len(np.unique( all_patches[i]))):
        cell_value[np.where(all_patches[i]==j+1)[0]] = all_best_patch[i][j]
    
    categorical = 100 * (cell_value-true_values)/true_values
    abs_cate = categorical
    new_scalar = abs_cate.copy()
    new_scalar[np.where(abs_cate >= 15)] = 15
    new_scalar[np.where(abs_cate < 15)] = 10
    new_scalar[np.where(abs_cate < 10)] = 5
    new_scalar[np.where(abs_cate < 5)] = 1
    new_scalar[np.where(abs_cate < 1)] = 0
    new_scalar[np.where(abs_cate < -1)] = -1
    new_scalar[np.where(abs_cate < -5)] = -5
    new_scalar[np.where(abs_cate < -10)] = -10
    new_scalar[np.where(abs_cate < -15)] = -15
    
    cmap_colors = []
    category=[-15   ,-10    ,-5     ,-1      ,0      ,1        ,5        ,10       ,15]
    cols =  ["green", "blue", "grey", "brown", "pink", "purple", "yellow", "orange", "red", "black"]
    compteur = -1
    for my_test in category:
        compteur = compteur+1
        if my_test in np.unique(new_scalar):
            cmap_colors.append(cols[compteur])
        

    
    
    model.meshing.mesh_pyvista.plot(scalars = new_scalar, 
                                    cpos = "xy", 
                                    cmap = cmap_colors,
                                    show_edges=True, 
                                    notebook = True)
