#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 10:54:29 2022

@author: livillenave
"""

#=======================================================#
# Source librairies
#=======================================================#
import dassflow2d as df2d
import numpy as np
import os
import matplotlib.pyplot as plt

from subprocess import call # terminal comand 
import h5py


### channel macdonnald
g= 10
lx= 1000
q=2
n=0.033

def h_exact(x, g, lx):
    res = (4/g)**(1/3) * ( 1 + 0.5 * np.exp(-16 *(x/lx - 0.5 )**2 ) )
    return(res)
    


#=======================================================#
# copy existing case files
#=======================================================#

dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"

os.chdir( f"{dassflow_dir}/code/")
os.system("make cleanres cleanmin")

#print("Would you like to paste case file")
#args = input("Press Y or N to continue.") # ajouter exit si pas Y ou N
#
#if args == "Y":
#    # delete all files in your simulation directory before starting
#    os.system(f"rm -r {dassflow_dir}/code/bin_A/*")
#    # Copy recursively the files provided in DassFlow case repository into your own simulation directory **code/bin_A/**.
#    os.system(f"cp -r {dassflow_dir}/cases/tuto_case/10_demo_wrap/* {dassflow_dir}/code/bin_A")
#    os.chdir( f"{dassflow_dir}/code/")
#    os.system("make cleanres cleanmin")



#=======================================================#
#=======================================================#
#              DIRECT RUN --- SAVE OBS   
#    
#=======================================================#
#=======================================================#
    
    
#=======================================================#
# initialise 
#=======================================================#

# initialise fortran instance, and python corrponding data (see DassFlowModel in sphinx documentation)
my_model = df2d.DassFlowModel(bin_dir =  f"{dassflow_dir}/code/bin_A" , arg = "direct")
my_model.build_grid()

# ----------------- PLOT MESH ----------------- #
print("would you like to plot mesh data")
args = input("Press Y or N to continue.") # ajouter exit si pas Y ou N
if args == "Y":
    my_model.plot_meshing(what = "cell")
    my_model.plot_meshing(what = "node")
    my_model.plot_meshing(what = "edge")    
        
    for i in range( my_model.model.mesh.ne):
        print("------------------------------------")
        print("i =", i)
        print("------------------------------------")
        tmp = my_model.model.mesh.edge[i].lim  
        if tmp>0:
            res = my_model.model.mesh.edgeb[tmp-1].ind
            print("i=", i+1)
            print("res=",res)
        
# ----------------- PRINT MESH INFO ----------------- #    
    for i in range( my_model.model.mesh.neb):
        
        if my_model.model.mesh.edgeb[i].typlim.decode("utf-8")  != 'wall':
             print(my_model.model.mesh.edgeb[i].ind  )
             print(my_model.model.mesh.edgeb[i].typlim.decode("utf-8")  )
             print( my_model.model.mesh.edgeb[i].group)
             print(my_model.model.mesh.edgeb[i].perio)
    
    
    
    for i in range( my_model.model.mesh.ncb):
             print(my_model.model.mesh.cellb[i].ind, 
                   my_model.model.mesh.cellb[i].typlim.decode("utf-8"),  
                   my_model.model.mesh.cellb[i].group,
                   my_model.model.mesh.cellb[i].cell,
                   my_model.model.mesh.cellb[i].grav.x,
                   my_model.model.mesh.cellb[i].grav.y)
             
    for i in range( my_model.model.mesh.nnb):
             print(my_model.model.mesh.nodeb[i].ind, 
                   my_model.model.mesh.nodeb[i].typlim.decode("utf-8"),  
                   my_model.model.mesh.nodeb[i].group)





# ----------------- Define initial conditions + get x coordinates ----------------- #    
x = np.zeros(my_model.model.mesh.nc)

for i in range(my_model.model.mesh.nc):
    x[i] = my_model.model.mesh.cell[i].grav.x
    my_model.model.dof0.h[i]=  h_exact(x = my_model.model.mesh.cell[i].grav.x, g = g, lx=lx)
    my_model.model.dof0.u[i]= q/ my_model.model.dof0.h[i]
    my_model.model.dof0.v[i]=0

dx=x[1] - x[0]
lx = x[len(x)-1]  +  dx/2

# ----------------- Define initial conditions for GHOST CELLS ----------------- #    
    # specific to my test case, the code to generalise has to be produced
for i in range(my_model.model.mesh.neb):
    ind_edge_global = my_model.model.mesh.edgeb[i].ind -1
    ind_ghost_cell = my_model.model.mesh.edge[ind_edge_global].cell[1]-1
    ind_true_cell = my_model.model.mesh.edge[ind_edge_global].cell[0]-1
    
    if my_model.model.mesh.edgeb[i].typlim.decode("utf-8") == "wall" :
#        my_model.model.my_param_model.bathy_cell[ind_ghost_cell] = my_model.model.my_param_model.bathy_cell[ind_true_cell]
        my_model.model.dof0.h[ind_ghost_cell] = my_model.model.dof0.h[ind_true_cell]
        my_model.model.dof0.u[ind_ghost_cell] = - my_model.model.dof0.u[ind_true_cell]
        my_model.model.dof0.v[ind_ghost_cell] = my_model.model.dof0.v[ind_true_cell]
        
    if my_model.model.mesh.edgeb[i].typlim.decode("utf-8") == "discharg1" :
        print("----------- discharg1-------")
        print("dof0 :", my_model.model.dof0.h[ind_ghost_cell], my_model.model.dof0.u[ind_ghost_cell])
        print("dof final :",my_model.model.dof.h[ind_ghost_cell], my_model.model.dof.u[ind_ghost_cell])
        #print("bathy 0 :", my_model.model.my_param_model.bathy_cell[ind_ghost_cell])        
        #print("bathy final :",df2d.m_model.get_array_bathy_cell()[ind_ghost_cell])
#        my_model.model.my_param_model.bathy_cell[ind_ghost_cell] = my_model.model.my_param_model.bathy_cell[ind_true_cell]
        my_model.model.dof0.h[ind_ghost_cell] = my_model.model.dof0.h[ind_true_cell]
        my_model.model.dof0.u[ind_ghost_cell] = my_model.model.dof0.u[ind_true_cell]
        my_model.model.dof0.v[ind_ghost_cell] = my_model.model.dof0.v[ind_true_cell]
        
    if my_model.model.mesh.edgeb[i].typlim.decode("utf-8") == "hpresc" :
        print("----------- hpresc---------------")
        print("dof0 :", my_model.model.dof0.h[ind_ghost_cell], my_model.model.dof0.u[ind_ghost_cell])
        print("dof final :",my_model.model.dof.h[ind_ghost_cell], my_model.model.dof.u[ind_ghost_cell])
        #print("bathy 0 :", my_model.model.my_param_model.bathy_cell[ind_ghost_cell])        
        #print("bathy final :",df2d.m_model.get_array_bathy_cell()[ind_ghost_cell])
#        my_model.model.my_param_model.bathy_cell[ind_ghost_cell] = my_model.model.my_param_model.bathy_cell[ind_true_cell]
        my_model.model.dof0.h[ind_ghost_cell] = h_exact(x = lx+dx/2, g = g, lx=lx)
        my_model.model.dof0.u[ind_ghost_cell] = q / h_exact(x = lx+dx/2, g = g, lx=lx)
        my_model.model.dof0.v[ind_ghost_cell] = my_model.model.dof0.v[ind_true_cell]       
        
my_model.model.dof = my_model.model.dof0
 

# ----------------- Parameter definition----------------- #   
        
# >>> Manning parameter
my_model.model.my_friction.manning[:] = n

# ----------------- FORTRAN CALL ----------------- #  
#   update python variables in fortran (bathy_cell, dof, manning)
my_model.update_fortran()


# ----------------- Boundary condtion definition ----------------- #    
# must be done after my_model.update_fortran() call, that initialise (allocate 
# and set values) the arrays of boundary conditions, based on furnished txt files

my_bc = df2d.wrapping.m_model.get_bc()
#my_bc.hyd[0].q[:] = 20
#my_bc.hpresc[0].h[:] =   h_exact(x = lx+dx/2, g = g, lx=lx)
#df2d.wrapping.m_model.set_bc(my_bc)




# ----------------- FORTRAN CALL ----------------- #  
#  run and save result in hdf5 format
# run model
my_model.run()
# save simulation results in hdf5 files
my_model.save_res()

################# check bc
for i in range(my_model.model.mesh.neb):
    ind_edge_global = my_model.model.mesh.edgeb[i].ind -1
    ind_ghost_cell = my_model.model.mesh.edge[ind_edge_global].cell[1]-1
    ind_true_cell = my_model.model.mesh.edge[ind_edge_global].cell[0]-1

        
    if my_model.model.mesh.edgeb[i].typlim.decode("utf-8") == "discharg1" :
        print("----------- discharg1-------")
        print("dof0 :", my_model.model.dof0.h[ind_ghost_cell], my_model.model.dof0.u[ind_ghost_cell])
        print("dof final :",my_model.model.dof.h[ind_ghost_cell], my_model.model.dof.u[ind_ghost_cell])
        print("bathy 0 :", my_model.model.my_param_model.bathy_cell[ind_ghost_cell])        
        print("bathy final :",df2d.m_model.get_array_bathy_cell()[ind_ghost_cell])
        
    if my_model.model.mesh.edgeb[i].typlim.decode("utf-8") == "hpresc" :
        print("----------- hpresc---------------")
        print("dof0 :", my_model.model.dof0.h[ind_ghost_cell], my_model.model.dof0.u[ind_ghost_cell])
        print("dof final :",my_model.model.dof.h[ind_ghost_cell], my_model.model.dof.u[ind_ghost_cell])
        print("bathy 0 :", my_model.model.my_param_model.bathy_cell[ind_ghost_cell])        
        print("bathy final :",df2d.m_model.get_array_bathy_cell()[ind_ghost_cell])

# necessary for plots # builds callable objects
my_model.build_grid() 





print("Would you like to plot initial conditions")
args = input("Press Y or N to continue.") # ajouter exit si pas Y ou N
if args == "Y":    
    my_model.plot_grid(my_model.model.my_param_model.bathy_cell[0:my_model.model.mesh.nc], title_plot = "bathymetry")
    my_model.plot_grid(my_model.model.dof0.h[0:my_model.model.mesh.nc], title_plot = "initial h")
    my_model.plot_grid(my_model.model.my_param_model.bathy_cell[0:my_model.model.mesh.nc]+
                       my_model.model.dof0.h[0:my_model.model.mesh.nc], title_plot = "initial zs")
    my_model.plot_grid(my_model.model.dof0.u[0:my_model.model.mesh.nc], title_plot = "initial u")
    my_model.plot_grid(my_model.model.dof0.v[0:my_model.model.mesh.nc], title_plot = "initial v")
        
    
    
print("would you like to plot hydrograph and rating curve")
args = input("Press Y or N to continue.") # ajouter exit si pas Y ou N
if args == "Y":

    i= 0
    plt.plot( my_bc.hyd[i].t, my_bc.hyd[i].q )
    plt.xlabel("time [s]")
    plt.ylabel("Q [m3/s]")
    plt.title(f"hydrograph {i} file")
    plt.show()
    
    plt.plot(my_bc.rat[i].q, my_bc.rat[i].h )
    plt.xlabel("Discharge, Q [m3/s]")
    plt.ylabel("Water Heigth, H [M]")
    plt.title(f"Rating curve {i} file")
    plt.show()
    
    plt.plot(my_bc.hpresc[i].t, my_bc.hpresc[i].h )
    plt.xlabel("time [s]")
    plt.ylabel("Water Heigth, H [M]")
    plt.title(f"hpresc {i} file")
    plt.show()



print("Would you like to plot some model outputs - from hdf5-file")
args = input("Press Y or N to continue.") # ajouter exit si pas Y ou N

if args == "Y":
    my_model.plot_var(what = "bathy", when = "initial", title_plot = "bahtymetry", save_plot=True, filename = "./res/bathy")
    my_model.plot_var(what = "h", when = "initial", title_plot = "INITIAL h", save_plot=True, filename = "./res/h_0")
    my_model.plot_var(what = "zs", when = "initial", title_plot = "INITIAL zs", save_plot=True, filename = "./res/zs_0")
    my_model.plot_var(what = "u", when = "initial", title_plot = "INITIAL u", save_plot=True, filename = "./res/u_0")
    my_model.plot_var(what = "v", when = "initial", title_plot = "INITIAL v", save_plot=True, filename = "./res/v_0")
    
    # etc ...
    my_model.plot_var(what = "h", when = 0, title_plot = "INITIAL h")
    my_model.plot_var(what = "h", when = 1, title_plot = "h at second time step", save_plot=True,filename = "./res/h_fin")
    
    # result at the end of the simulation:
    my_model.plot_var(what = "vel", when = "final", title_plot = "norm(u,v) at final time step", save_plot=True,filename = "./res/velocity_fin")
    my_model.plot_var(what = "zs", when = "final", title_plot = "Zs(m) at final time step", save_plot=True,filename = "./res/Zs_fin")
    

print("Would you like to plot temporal evolution of the water height")
args = input("Press Y or N to continue.") # ajouter exit si pas Y ou N

if args == "Y":
    
    for i in range(11):
        my_model.plot_var(what = "h", when = i, title_plot = f" h  at write timestep= {i}")
        


print("Would you like to plot temporal evolution of the water heigth 1D")
args = input("Press Y or N to continue.") # ajouter exit si pas Y ou N

if args == "Y":
    a = h5py.File(name = f"{dassflow_dir}/code/bin_A/res/simu.hdf5")

    for i in range(a["spatial_var/bathy"].shape[1]):
        b=a["spatial_var/bathy"][:,i] 
        zs = a["spatial_var/h"][:,i] +  a["spatial_var/bathy"][:,0] 
        h = a["spatial_var/h"][:,i]
        fig, ax1 = plt.subplots()
        mytime = a["time"][i].round(2)
        ax1.set_title(f"time = {mytime} s")
        ax2 = ax1.twinx()
        ax1.scatter(x = np.arange(len(zs )),
                            y = zs, 
                            linewidths = 1,
                            label = "zs")
        ax1.scatter(x = np.arange(len(b )),
                            y = b, 
                            c="red", 
                            linewidths = 1,
                            label = "bathymetry" )
        ax2.scatter(x = np.arange(len(h)),
                            y = h, 
                            c="purple", 
                            linewidths = 1,
                            label = "h" ) 
        ax2.scatter(x = np.arange(len(h)),
                            y = h_exact(x=x, g=g, lx=lx), 
                            c="green", 
                            linewidths = 1,
                            label = "h_TRUE" )
        
        ax1.set_xlabel('X [m]')
        ax1.set_ylabel('elevation [m]', color='r')
        ax2.set_ylabel('Water height [m]', color='g')
#        ax1.legend()
#        ax2.legend()
        fig.legend()

for i in range(a["spatial_var/bathy"].shape[1]):
    plt.scatter(x = np.arange(len(b )),
                        y = a["spatial_var/h"][:,i] * a["spatial_var/u"][:,i], 
                        c="red", 
                        linewidths = 1,
                        label = "h*u" )
    
    plt.scatter(x = np.arange(len(b )),
                        y = a["spatial_var/h"][:,i], 
                        c="green", 
                        linewidths = 1,
                        label = "h" )
    
    plt.scatter(x = np.arange(len(b )),
                        y = a["spatial_var/u"][:,i], 
                        c="blue", 
                        linewidths = 1,
                        label = "u" )
    
    
    plt.legend()
    plt.show()
    

    plt.scatter(x = np.arange(len(a["spatial_var/h"][99,:])),
                        y = a["spatial_var/h"][99,:] * a["spatial_var/u"][99,:], 
                        c="red", 
                        linewidths = 1,
                        label = "h*u" )
    
    plt.scatter(x = np.arange(len(a["spatial_var/h"][99,:])),
                        y = a["spatial_var/h"][99,:], 
                        c="green", 
                        linewidths = 1,
                        label = "h" )
    
    plt.scatter(x = np.arange(len(a["spatial_var/h"][99,:])),
                        y = a["spatial_var/u"][9.9,:], 
                        c="blue", 
                        linewidths = 1,
                        label = "u" )
    
    
    plt.legend()
    plt.show()

    plt.scatter(x = np.arange(len(a["spatial_var/h"][0,:])),
                        y = a["spatial_var/h"][0,:] * a["spatial_var/u"][0,:], 
                        c="red", 
                        linewidths = 1,
                        label = "h*u" )
    
    plt.scatter(x = np.arange(len(a["spatial_var/h"][0,:])),
                        y = a["spatial_var/h"][0,:], 
                        c="green", 
                        linewidths = 1,
                        label = "h" )
    
    plt.scatter(x = np.arange(len(a["spatial_var/h"][0,:])),
                        y = a["spatial_var/u"][0,:], 
                        c="blue", 
                        linewidths = 1,
                        label = "u" )
    
    
    plt.legend()
    plt.show()

bathy = a["spatial_var/bathy"][:,i]
dzdx = np.zeros(len(bathy))
dzdx[1:len(bathy)]= bathy[1:len(bathy)] -  bathy[0:len(bathy)-1]

h = a["spatial_var/h"][:,i]
dhdx = np.zeros(len(bathy))
dhdx[1:len(bathy)]= h[1:len(h)] -  h[0:len(h)-1]

dx= 10
res = dzdx -( dhdx * (q**2 /(g*h**3)) - dx * (n**2 * q**2) / (h ** (10/3)))
plt.scatter(x = np.arange(len(res )),
                    y = res, 
                    c="blue", 
                    linewidths = 1,
                    label = "u" )
            
def cp_dir(source, target):
  call(['cp', '-a', source, target]) # Linux


def cp_hdf5_file(source, target):
    source = f'res/{source}'
    target = f'hdf5/{target}'
    call(['mkdir', "hdf5"])
    call(['rm', target]) 
    call(['cp', source, target]) # Linux

if df2d.wrapping.m_common.get_w_obs() == 1 : 
    cp_dir('./res/obs', '.')
    cp_hdf5_file(source = "simu.hdf5", target = "true.hdf5")
    # cp_dir copy files from source directory to target directory.
# note that the directory is copied (not only the files within)

                

df2d.wrapping.call_model.clean_model(my_model.model)         # deallocate correctly (necessary action)







#=======================================================#
#=======================================================#
#              INFERENCE RUN --- SAVE OBS   
#    
#=======================================================#
#=======================================================#

# initialise fortran instance, and python corrponding data
inf_model = df2d.DassFlowModel(bin_dir =  f"{dassflow_dir}/code/bin_A" , run_type = "min")
df2d.wrapping.m_common.set_use_obs(1)
inf_model.build_grid()



my_model.model.my_friction.manning[:] = 0.04
#my_model.model.my_param_model.bathy_cell[:] = 1
inf_model.model.dof0.h[:]=1
inf_model.model.dof0.u[:]=0
inf_model.model.dof0.v[:]=0


print("CHOOSE INFERENCE TYPE (1 hydrograph, 2 land_use, 3 = bathy, 4 = initial conditions)")
inference_type = input("Enter 1,2 or 3 \n") # ajouter exit si pas Y ou N



if inference_type == "1":
    df2d.m_model.set_c_hydrograph(1)
    df2d.m_model.set_c_manning(0)
    df2d.m_model.set_c_bathy(0)
    inf_model.update_fortran()
    # -- infer hydrograph --
    my_bc = df2d.wrapping.m_model.get_bc()
    my_bc.hyd[0].q[:] = 20
    df2d.wrapping.m_model.set_bc(my_bc)
    
elif inference_type == "2":
    df2d.m_model.set_c_manning(1)
    # -- same for manning -- #
    inf_model.model.my_friction.manning[:] = 1
    inf_model.update_fortran()
    
elif inference_type == "3":
    df2d.m_model.set_c_bathy(1)
    # -- infer hydrograph --#   *
    inf_model.model.my_param_model.bathy_cell[:] = 1    
    inf_model.update_fortran()

elif inference_type == "4":
    inf_model.model.dof0.h[:]=2
    inf_model.model.dof0.u[:]=-1
    inf_model.model.dof0.v[:]=0
    inf_model.update_fortran()

my_bc = df2d.wrapping.m_model.get_bc()
my_bc.hpresc[0].h[:] = 1
if inference_type != "1":
    my_bc.hyd[0].q[:] = 10
    df2d.wrapping.m_model.set_bc(my_bc)    

inf_model.run()
inf_model.save_res()               

df2d.wrapping.call_model.clean_model(inf_model.model)         # deallocate correctly (necessary action)
