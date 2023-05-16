#==============================================================#
#==============================================================#
# set initial conditions
#==============================================================#
#==============================================================#

import numpy as np
# 
#
# @param : hdf5_path path to hdf5 file
def ic_to_hdf5(hdf5_path, my_model):
        
    simu = h5py.File(hdf5_path, "a")
    ic = simu.create_group("initial_condition")
    
    ic.create_dataset("h", data = my_model.model.dof0.h[:]  )
    ic.create_dataset("u", data = my_model.model.dof0.u[:]  )
    ic.create_dataset("v", data = my_model.model.dof0.v[:]  )
            
    simu.close()
    
    
# @param : reference_hdf5_path path to hdf5 file
#            --> This file contains reference initial conditions to load
# my_model: Dassflowmodel instance, on which we apply the values
def hdf5_to_initial_condition(reference_hdf5_path, my_model):

    simu = h5py.File(reference_hdf5_path, "r")
    ic = simu["initial_condition"]
    
    my_model.model.dof0.h[:]  = ic["h"][:]
    my_model.model.dof0.u[:]  = ic["u"][:]
    my_model.model.dof0.v[:]  = ic["v"][:]
            
    simu.close()
    dof0 = my_model.model.dof0
    return(dof0)
    
    
def get_xy(my_model, what = "cell"):        
    if what !="cell" and what != "ghostcell":
        print('argument what wrong')
    elif what == "cell":
        cell = my_model.plot_meshing("cell")
        table = np.ndarray(shape = (2,len(cell)))
        for id in range(len(cell)):
            table[0,id] =  cell[id]["grav_x"]
            table[1,id] =  cell[id]["grav_y"]
        return(table)
    elif what == "ghostcell":
        table = np.ndarray(shape = (2,my_model.model.mesh.ncb))
        for i in range( my_model.model.mesh.ncb):
            table[0,i] =  my_model.model.mesh.cellb[i].grav.x
            table[1,i] = my_model.model.mesh.cellb[i].grav.y
        return(table)
        
        
def set_xyarray_to_cell(func, my_model, what = "h"):
    
    xy = get_xy(my_model, "cell")
    xyb =  get_xy(my_model, "ghostcell")

    if what =="h" :
        var = my_model.model.dof0.h
    if what =="u" :
        var = my_model.model.dof0.u
    if what =="v" :
        var = my_model.model.dof0.v
    
    for k in range(my_model.model.mesh.nc) :
        var[k] = func(x = xy[0,k], y = xy[1,k])
        
    for k in range(my_model.model.mesh.ncb) :
        var[my_model.model.mesh.nc + k] = func(x = xyb[0,k], y = xyb[1,k])   
        
    return(var)








#=======================================================#
#  EXAMPLE Direct run plus save resut (initial conditions)
#=======================================================#
    
#import dassflow2d as df2d
#import h5py
#import os
#
#bin_dir =  "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A"
#        
#os.chdir(bin_dir)
#os.chdir("../.")
#os.system("make cleanres cleanmin")
#os.system(f"rm {bin_dir}/restart.bin")
#os.chdir(bin_dir)    
## initialise fortran instance, and python corrponding data
#my_model = df2d.DassFlowModel(bin_dir = bin_dir, 
#                              run_type = "direct") # run_type can be min or direct (grad ?)
#my_model.build_grid()
#my_model.update_fortran()
#
#
#
#
#def h_init(x,y):
#    
#    if(y>50):
#        sym_y = y -50
#    if(y<=50):
#        sym_y = abs(y - 50)
#    if(y==50):
#        sym_y = 0
#        
#    array  =  2  - (x /1000) + (sym_y/100)
#    return(array)
#    
#    
#    
#my_model.model.dof0.h[:] = set_xyarray_to_cell(func = h_init, my_model = my_model, what ="h")
#my_model.model.dof0.u[:] = 1
#my_model.model.dof0.v[:] = 0.0000000
#
#my_model.run()
#
#my_model.save_res()
#os.chdir("./res/")
#ic_to_hdf5("simu.hdf5", my_model)
#os.chdir("../")
#df2d.call_model.clean_model(my_model.model)
#
#import shutil
#shutil.copyfile("./res/simu.hdf5", "./ic.hdf5")
#
##==============================================================#
## use ic.hdf5
##==============================================================#
#
#os.system("make cleanres cleanmin")
#os.system(f"rm {bin_dir}/restart.bin")
#        
#os.chdir(bin_dir)
## initialise fortran instance, and python corrponding data
#my_model = df2d.DassFlowModel(bin_dir = bin_dir, 
#                              run_type = "direct") # run_type can be min or direct (grad ?)
#my_model.build_grid()
#my_model.update_fortran()
#
#tmp = hdf5_to_initial_condition("ic.hdf5", my_model)
#    
#my_model.run()
#my_model.save_res()
#df2d.call_model.clean_model(my_model.model)