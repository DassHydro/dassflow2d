#==================================#
# define librairies                #
#==================================#

import os
os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/Tools/2_gen_basic_channel")


import gen_channel_case
import gen_dassflow

import errno
import numpy as np
import pandas as pd


#==================================#
# Param        #
#==================================#

lx = 1000 # length of channel in x direction
nx = 1000 # number of nodes in x direction

ly = 10 # length of channel in y direction
ny = 2  # number of nodes in y direction


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

#==================================#
# in rough params
#==================================#
ts = 1000 # time of simulation
dtw = 100 # timestep to write result
cfl = 0.8 # cfl criteria 

# ... for more see default parameters set in input_param dicitonary below

#==================================#
# gen files
#==================================#

gen_dassflow.gen_dassflow_files(
			mesh_param 	 = {"nx": nx, "ny": ny, "lx" : lx , "ly" :ly}, 
			friction_param  = {"manning_alpha":n, "manning_beta":0},
			bc_param        = {"in_typ":in_typ, "table_in":qin, 
					    "out_typ":out_typ, "table_out":hpresc},
			obs_param       = {"nx_obs":2, "ny_obs":1, "xmax_obs":600, "ymax_obs":50, "xmin_obs":400, "ymin_obs":50, "dt_obs":60.},
			
			input_param     = { "mesh_name":'channel.geo',
                                        "ts":ts, 
                                        "dta":1000, 
                                        "dtw":dtw, 
                                        "dtp":10,
                                        "dt":1, 
                                        "temp_scheme":'euler', 
                                        "spatial_scheme":'first_b1', 
                                        "adapt_dt":1, 
                                        "cfl":cfl,
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
                                        "c_infil":0,
                                        "c_ic":0,
                                        "restart_min":0,
                                        "eps_min":0.0001})



