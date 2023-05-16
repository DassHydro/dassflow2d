#==============================================================#
#==============================================================#
# INTERFACING configuration parameters (python / fortran / hdf5 file)
#==============================================================#
#==============================================================#

import dassflow2d as df2d
import h5py
import copy




# Load input param for fortran kernel
# @param: input_param :: dictionarry with keys corresponging to fields to be filled
# @return: dictionary input_param filled with fortran kernel values
# df2d istance must be on
def get_input_param(input_param):
    
    res = copy.deepcopy(input_param)
#    all_keys = input_param.keys()
    for k in res:        
       print(k)
        # in module m_common
       if k == 'mesh_name':      
             print(res[k])
             res[k] = df2d.m_common.get_mesh_name().decode('utf-8')
             print(res[k])
       elif k == 'ts':
             print(res[k])
             res[k] = df2d.m_common.get_ts()
             print(res[k])
               # DTA IN ASSIMILATION FILE (OBSOLETE HERE ?)
#       elif k == 'dta':
#               res[k] = df2d.m_common.get_dta()
       elif k == 'dtw':
               print("dtw")
               res[k] = df2d.m_common.get_dtw()
       elif k == 'dtp':
               res[k] = df2d.m_common.get_dtp()
       elif k == 'adapt_dt':
               res[k] = df2d.m_common.get_adapt_dt()
       elif k == 'dt':
               res[k] = df2d.m_common.get_dt()
       elif k == 'cfl':
               res[k] = df2d.m_common.get_cfl()
       elif k == 'temp_scheme':
               res[k] = df2d.m_common.get_temp_scheme().decode('utf-8')
       elif k == 'spatial_scheme':
               res[k] = df2d.m_common.get_spatial_scheme().decode('utf-8')               
       elif k == 'w_tecplot':
               res[k] = df2d.m_common.get_w_tecplot()
       elif k == 'w_gnuplot':
               res[k] = df2d.m_common.get_w_gnuplot()
       elif k == 'w_vtk':
               res[k] = df2d.m_common.get_w_vtk()
       elif k == 'w_obs':
               res[k] = df2d.m_common.get_w_obs()
       elif k == 'use_obs':
               res[k] = df2d.m_common.get_use_obs()
       elif k == 'max_nt_for_adjoint':
               res[k] = df2d.m_common.get_max_nt_for_adjoint()
       elif k == 'max_nt_for_direct':
               res[k] = df2d.m_common.get_max_nt_for_direct()
       elif k == 'restart_min':
               res[k] = df2d.m_common.get_restart_min()
       elif k == 'get_eps_min':
               res[k] = df2d.m_common.get_eps_min()

           #in module m_model
       elif k == 'feedback_inflow':
               res[k] = df2d.m_model.get_feedback_inflow()
       elif k == 'coef_feedback':              
               res[k] = df2d.m_model.get_coef_feedback()
       elif k == 'friction':
               res[k] = df2d.m_model.get_friction()
       elif k == 'g':           
               res[k] = df2d.m_model.get_g()
       elif k == 'c_manning':           
              res[k] =  df2d.m_model.get_c_manning()
       elif k == 'c_manning_beta':           
               res[k] = df2d.m_model.get_c_manning_beta()
       elif k == 'c_bathy':           
               df2d.m_model.get_c_bathy()
       elif k == 'c_hydrograph':           
               res[k] = df2d.m_model.get_c_hydrograph()               
       elif k == 'c_rain':           
               res[k] = df2d.m_model.get_c_rain()             
       elif k == 'c_ic':           
               res[k] = df2d.m_model.get_c_ic()    
       
        
    return(res)   
    

# update values in fortran kernel based on provided input
# @param: input_param :: dictionary with keys corresponging to fields to be filled in fortran + values set
# @return: nothing : set the values in fortran kernel
# df2d istance must be on
def set_input_param(input_param):
#    all_keys = input_param.keys()
    for k in input_param:        
       print(k)
        # in module m_common
       if k == 'mesh_name':
                df2d.m_common.set_mesh_name(input_param[k])
       elif k == 'ts':
                df2d.m_common.set_ts(input_param[k])
               # DTA IN ASSIMILATION FILE (OBSOLETE HERE ?)
#       if k == 'dta':
#               input_param[k] = df2d.m_common.set_dta()
       elif k == 'dtw':
               df2d.m_common.set_dtw(input_param[k])
       elif k == 'dtp':
               df2d.m_common.set_dtp(input_param[k])
       elif k == 'adapt_dt':
               df2d.m_common.set_adapt_dt(input_param[k])
       elif k == 'dt':
               df2d.m_common.set_dt(input_param[k])
       elif k == 'cfl':
               df2d.m_common.set_cfl(input_param[k])
       elif k == 'temp_scheme':
               df2d.m_common.set_temp_scheme(input_param[k])
       elif k == 'spatial_scheme':
               df2d.m_common.set_spatial_scheme(input_param[k])           
       elif k == 'w_tecplot':
               df2d.m_common.set_w_tecplot(input_param[k])
       elif k == 'w_gnuplot':
               df2d.m_common.set_w_gnuplot(input_param[k])
       elif k == 'w_vtk':
               df2d.m_common.set_w_vtk(input_param[k])
       elif k == 'w_obs':
               df2d.m_common.set_w_obs(input_param[k])
       elif k == 'use_obs':
               df2d.m_common.set_use_obs(input_param[k])
       elif k == 'max_nt_for_adjoint':
               df2d.m_common.set_max_nt_for_adjoint(input_param[k])
       elif k == 'max_nt_for_direct':
               df2d.m_common.set_max_nt_for_direct(input_param[k])
       elif k == 'restart_min':
               df2d.m_common.set_restart_min(input_param[k])
       elif k == 'set_eps_min':
               df2d.m_common.set_eps_min(input_param[k])

           #in module m_model
       elif k == 'feedback_inflow':
               df2d.m_model.set_feedback_inflow(input_param[k]) 
       elif k == 'coef_feedback':              
               df2d.m_model.set_coef_feedback(input_param[k]) 
       elif k == 'friction':
               df2d.m_model.set_friction(input_param[k]) 
       elif k == 'g':           
               df2d.m_model.set_g(input_param[k]) 
       elif k == 'c_manning':           
               df2d.m_model.set_c_manning(input_param[k]) 
       elif k == 'c_manning_beta':           
               df2d.m_model.set_c_manning_beta(input_param[k]) 
       elif k == 'c_bathy':           
               df2d.m_model.set_c_bathy(input_param[k]) 
       elif k == 'c_hydrograph':           
               df2d.m_model.set_c_hydrograph(input_param[k])           
    #   elif k == 'c_rain':           
    #           df2d.m_model.set_c_rain()              
       #elif k == 'c_ic':           
               #df2d.m_model.set_c_ic()    
    print("values set")
    return()   
    
    
    
def print_input_param(input_param):
    for key, value in input_param.items():
        print(key, ' : ', value)
           


# @param hdf5_file_path: path to hdf5 file
# @¶eturn set from fortran kernel values to hdf5_file_path file
def input_param_to_hdf5(hdf5_file_path):
    
    input_param = get_input_param( { "mesh_name":'channel.geo',
                                            "ts":0, 
                                            "dta":0, 
                                            "dtw":0, 
                                            "dtp":0,
                                            "dt":0, 
                                            "temp_scheme":'euler', 
                                            "spatial_scheme":'first_b1', 
                                            "adapt_dt":0, 
                                            "cfl":0., 
                                            "feedback_inflow":0,                                    
                                            "coef_feedback":0.,
                                            "heps":0, 
                                            "friction":0, 
                                            "g":10, 
                                            "w_tecplot":0,
                                            "w_vtk":0, 
                                            "w_gnuplot":0, 
                                            "w_obs":0, 
                                            "use_obs":0, 
                                            "max_nt_for_adjoint":0,
                                            "c_manning":0,
                                            "c_manning_beta":0,
                                            "c_bathy":0,
                                            "c_hydrograph":0,
                                            "c_ratcurve":0,
                                            "c_rain":0,
                                            #"c_infil":0,
                                            "c_ic":0,
                                            "restart_min":0,
                                            "eps_min":0} )
        
    
    h = h5py.File(hdf5_file_path, "a")
    input_data = h.create_group("input_data")
    
    for k, v in input_param.items():
            input_data.attrs[k] = v
            
    h.close()


# hdf5 to dictionary input_param
# @return : dictionary of input_param
def hdf5_to_input_param(hdf5_file_path):
    

    new_input_param = {}
    file = h5py.File(hdf5_file_path, 'r')
    input_data = file['input_data']
    
    for k in input_data.attrs.keys():
        new_input_param[k]= input_data.attrs[k]    
 
    
    return(new_input_param)



#=======================================================#
#  EXAMPLE
#=======================================================#
#input_param= { "ts":2000 } 
#bin_dir =  "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A"
#        
#import os
#os.chdir(bin_dir)
## initialise fortran instance, and python corrponding data
#my_model = df2d.DassFlowModel(bin_dir = bin_dir, 
#                              run_type = "direct") # run_type can be min or direct (grad ?)
#set_input_param(input_param)
#
#my_model.update_fortran()
#my_model.run()
#
#my_model.save_res()
#input_param = get_input_param(input_param)
#
#os.chdir("./res/")
#input_param_to_hdf5("simu.hdf5")
#os.chdir("../")
#
#df2d.call_model.clean_model(my_model.model)
