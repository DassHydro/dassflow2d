#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 11:20:57 2022

@author: livillenave
"""


os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/Tools/2_gen_basic_channel")


import gen_channel_case
import gen_dassflow


input_params = { "mesh_name":'channel.geo',
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
										"w_gnuplot":0, 
										"w_obs":0, 
										"use_obs":0, 
                                        "max_nt_for_adjoint":2500,
										"c_manning":0,
										"c_manning_beta":0,
										"c_bathy":0,
										"c_hydrograph":0,
										"c_ratcurve":0,
										"c_rain":0,
#										"c_infil":0,
										"c_ic":0,
                                        "restart_min":0,
                                        "eps_min":0.0001}

gen_dassflow.gen_input(input_param = input_params)