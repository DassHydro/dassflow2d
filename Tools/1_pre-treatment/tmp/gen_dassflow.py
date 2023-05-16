from __future__ import annotations

import errno
import os
import numpy as np
import pandas as pd
import gen_channel_case

__all__ = ["_gen_channel_case"]

lx = 1000 # length of channel in x direction
nx = 101 # number of nodes in x direction

ly = 10 # length of channel in y direction
ny = 11  # number of nodes in y direction

# ---------- Inflow boudary condition
in_typ="discharg1"            # discharge imposed  (approximate backwater curve)
# table of imposed dischage is created (stationary values so 2 rows are enough)
qin = np.ndarray( (2,2) )
qin[:,0] = np.arange(start=0, stop = 10000, step =10000)
qin[:,1] = 20

# ---------- Outflow boudary condition
out_typ="hpresc" 		# water eigth is imposed following analytical solution
# table of imposed water eigth is created (stationary values so 2 rows are enough)
hpresc= np.ndarray( (2,2) )
hpresc[:,0] = np.arange(start=0, stop = 10000, step =10000)
hpresc[:,1] = gen_channel_case.h_true_macdo(lx + (lx / (nx - 1))/2, lx = lx, g = 10)

n = 0.033 # manning-strickler coefficient value

def gen_dassflow_files(mesh_param 	 = {"nx": nx, "ny": ny, "lx" : lx , "ly" :ly},
			friction_param  = {"manning_alpha":0.005, "manning_beta":0},
			bc_param        = {"in_typ":"discharg1", "table_in":qin, 
										"out_typ":"hpresc", "table_out":hpresc},
			obs_param       = {"nx_obs":2, "ny_obs":1, "xmax_obs":600, "ymax_obs":5, "xmin_obs":400, "ymin_obs":5, "dt_obs":60.},
			
			input_param     = { "mesh_name":'channel.geo',
                                        "ts":100, 
                                        "dta":100, 
                                        "dtw":10, 
                                        "dtp":10,
                                        "dt":1, 
                                        "temp_scheme":'euler', 
                                        "spatial_scheme":'first_b1', 
                                        "adapt_dt":1, 
                                        "cfl":0.8, 
                                        "feedback_inflow":1,                                    
                                        "coef_feedback":0.8,
                                        "heps":0, 
                                        "friction":1, 
                                        "g":10, 
                                        "w_tecplot":0,
                                        "w_vtk":0, 
                                        "w_gnuplot":1, 
                                        "w_obs":1, 
                                        "use_obs":0, 
                                        "max_nt_for_adjoint":2500,
                                        "c_manning":0,
                                        "c_manning_beta":0,
                                        "c_bathy":0,
                                        "c_hydrograph":0,
                                        "c_ratcurve":0,
                                        "c_rain":0,
                                        "c_infil":0,
                                        "c_ic":0,
                                        "restart_min":0,
                                        "eps_min":0.0001} ):
	
	# TO DO ADD CHECK
	# obs bien placés, bc cohérentes, nx,ny>2
	# input param vs le reste
	gen_input(input_param = input_param)
	
	gen_channel_case.gen_basic_channel(					nx = mesh_param["nx"],   #number of nodes in x direction
										ny = mesh_param["ny"],   #number of nodes in y direction
										lx = mesh_param["lx"],   #length of channel in x direction
										ly = mesh_param["ly"])   #length of channel in y direction
	
	gen_channel_case.gen_land_use(manning_alpha = friction_param["manning_alpha"], #uniform manning value
					manning_beta =  friction_param["manning_beta"]) # coefficient for other friction law
	
	gen_channel_case.gen_bc(			 in_type	 = bc_param["in_typ"],  # discharg1 /discharg2
							 out_type 	 = bc_param["out_typ"])  # ratcurve /hpresc/zpresc
	
	gen_channel_case.gen_bc_data(	bc_typ = bc_param["in_typ"], 
									nrow=bc_param["table_in"].shape[0], 
									var1 =bc_param["table_in"][:,0], 
									var2= bc_param["table_in"][:,1]  )	
	
	gen_channel_case.gen_bc_data(	bc_typ = bc_param["out_typ"], 
									nrow=bc_param["table_out"].shape[0], 
									var1 =bc_param["table_out"][:,0], 
									var2= bc_param["table_out"][:,1])
									
	gen_channel_case.gen_obs(	nx_obs 	= obs_param["nx_obs"], 
								ny_obs	=obs_param["ny_obs"], 
								xmax_obs=obs_param["xmax_obs"], 
								ymax_obs=obs_param["ymax_obs"], 
								xmin_obs=obs_param["xmin_obs"], 
								ymin_obs=obs_param["ymin_obs"], 
								dt_obs=obs_param["dt_obs"])
								
	print("files generated in", os.getcwd(), " /files")
	
	
	
	

def gen_input(input_param= { "mesh_name":'channel.geo',
                                        "ts":100, 
                                        "dta":100, 
                                        "dtw":10, 
                                        "dtp":10,
                                        "dt":1, 
                                        "temp_scheme":'euler', 
                                        "spatial_scheme":'first_b1', 
                                        "adapt_dt":1, 
                                        "cfl":0.8, 
                                        "feedback_inflow":1,                                    
                                        "coef_feedback":0.8,
                                        "heps":0, 
                                        "friction":1, 
                                        "g":10, 
                                        "w_tecplot":0,
                                        "w_vtk":0, 
                                        "w_gnuplot":1, 
                                        "w_obs":0, 
                                        "use_obs":0, 
                                        "max_nt_for_adjoint":2500,
                                        "c_manning":0,
                                        "c_manning_beta":0,
                                        "c_bathy":0,
                                        "c_hydrograph":0,
                                        "c_ratcurve":0,
                                        "c_rain":0,
                                        #"c_infil":0,
                                        "c_ic":0,
                                        "restart_min":0,
                                        "eps_min":0.0001} ):
        
    with open('files/input.txt', 'w') as f:
        f.write('!======================================================================================================================!')
        f.write('\n!Input File for Shallow-Water Model')
        f.write('\n!======================================================================================================================!')
        f.write('\n\n    &list_input')
        
        f.write('\n\n!======================================================================================================================!')
        f.write('\n!Mesh Type')
        f.write('\n!======================================================================================================================!')
        f.write("\n\n	mesh_type	=	'dassflow',")
        tmp_name = f"""\n    mesh_name    =   '{input_param["mesh_name"]}', """
        f.write(tmp_name)
        
        f.write('\n\n!======================================================================================================================!')
        f.write('\n!Simulation parameters')
        f.write('\n!======================================================================================================================!')
        tmp_name = f"""\n\n  ts                 =   {input_param["ts"]},      ! Simulation Time """
        f.write(tmp_name)
        tmp_name = f"""\n    dtw                =   {input_param["dtw"]},       ! data assimilation"""
        f.write(tmp_name)
        tmp_name = f"""\n    dtp                =   {input_param["dtp"]}, """
        f.write(tmp_name)    
        tmp_name = f"""\n    dta                =   {input_param["dta"]}, """
        f.write(tmp_name)    
        tmp_name = f"""\n    temp_scheme        =   '{input_param["temp_scheme"]}', """
        f.write(tmp_name)    
        tmp_name = f"""\n    spatial_scheme     =   '{input_param["spatial_scheme"]}', """
        f.write(tmp_name)    
        tmp_name = f"""\n    friction           =   {input_param["friction"]}, """
        f.write(tmp_name)    
        tmp_name = f"""\n    adapt_dt           =   {input_param["adapt_dt"]}, """
        f.write(tmp_name)    
        tmp_name = f"""\n    dt                 =   {input_param["dt"]}, """
        f.write(tmp_name)      
        tmp_name = f"""\n    cfl                =   {input_param["cfl"]}, """
        f.write(tmp_name)     
        tmp_name = f"""\n    feedback_inflow    =   {input_param["feedback_inflow"]}, """
        f.write(tmp_name)     
        tmp_name = f"""\n    coef_feedback      =   {input_param["coef_feedback"]}, """
        f.write(tmp_name)     
        f.write('\n\n!======================================================================================================================!')
        f.write('\n!PHYSICAL PARAMETERS')
        f.write('\n!======================================================================================================================!')
        tmp_name = f"""\n\n    g                  =   {input_param["g"]}, """
        f.write(tmp_name) 
        f.write('\n\n!======================================================================================================================!')
        f.write('\n!OUTPUT RESULTS')
        f.write('\n!======================================================================================================================!')     
        tmp_name = f"""\n\n    w_tecplot          =   {input_param["w_tecplot"]}, """
        f.write(tmp_name)   
        tmp_name = f"""\n      w_gnuplot          =   {input_param["w_gnuplot"]}, """
        f.write(tmp_name)   
        tmp_name = f"""\n      w_vtk              =   {input_param["w_vtk"]}, """
        f.write(tmp_name)   
        f.write('\n\n!======================================================================================================================!')
        f.write('\n!ASSIMILATION PARAMETER')
        f.write('\n!======================================================================================================================!') 
        tmp_name = f"""\n\n    w_obs                 =   {input_param["w_obs"]}, """
        f.write(tmp_name)   
        tmp_name = f"""\n    use_obs               =   {input_param["use_obs"]}, """
        f.write(tmp_name)   
        tmp_name = f"""\n    max_nt_for_adjoint    =   {input_param["max_nt_for_adjoint"]}, """
        f.write(tmp_name) 
       
        tmp_name = f"""\n\n    c_hydrograph          =   {input_param["c_hydrograph"]}, """
        f.write(tmp_name)    
        tmp_name = f"""\n    c_ratcurve            =   {input_param["c_ratcurve"]}, """
        f.write(tmp_name) 
        tmp_name = f"""\n    c_manning             =   {input_param["c_manning"]}, """
        f.write(tmp_name) 
        tmp_name = f"""\n    c_manning_beta        =   {input_param["c_manning_beta"]}, """
        f.write(tmp_name)
        tmp_name = f"""\n    c_bathy               =   {input_param["c_bathy"]}, """
        f.write(tmp_name)
        tmp_name = f"""\n    c_rain                =   {input_param["c_rain"]}, """
        f.write(tmp_name)
       # tmp_name = f"""\n    c_infil               =   {input_param["c_infil"]}, """
       # f.write(tmp_name)
        tmp_name = f"""\n    c_ic                  =   {input_param["c_ic"]}, """
        f.write(tmp_name)
        
        tmp_name = f"""\n\n  restart_min         =   {input_param["restart_min"]}, """
        f.write(tmp_name)                
        tmp_name = f"""\n    eps_min               =   {input_param["eps_min"]}, """
        f.write(tmp_name)    
        f.write("\n //  ")         
        

