#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 11:51:10 2022

@author: livillenave
"""


import dassflow2d as df2d
import numpy as np
import sys
import os



# initialise fortran instance, and python corrponding data
my_model = df2d.DassFlowModel(bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/", run_type = "direct") # run_type can be min or direct (grad ?)

#=========================
# ACESS input values
# ========================
df2d.m_common.get_mesh_name().decode('utf-8')
df2d.m_common.get_ts()
df2d.m_common.get_dtw()
df2d.m_common.get_dtp()
df2d.m_common.get_temp_scheme().decode('utf-8')
df2d.m_common.get_spatial_scheme().decode('utf-8')
df2d.m_common.get_adapt_dt()
df2d.m_common.get_dt()
df2d.m_common.get_cfl()

df2d.m_model.get_feedback_inflow()
df2d.m_model.get_coef_feedback()
df2d.m_model.get_heps()
df2d.m_model.get_friction()
df2d.m_model.get_g()

df2d.m_common.get_w_tecplot()
df2d.m_common.get_w_vtk()
df2d.m_common.get_w_gnuplot()

df2d.m_common.get_w_obs()
df2d.m_common.get_use_obs()
# df2d.xxxxx.get_max_nt_for_adjoint() ---> didn't found how to access it



df2d.m_model.get_c_manning()
df2d.m_model.get_c_manning_beta()
df2d.m_model.get_c_bathy()
df2d.m_model.get_c_hydrograph()
df2d.m_model.get_c_ratcurve()
df2d.m_model.get_c_rain()
df2d.m_model.get_c_ic()


df2d.m_common.get_restart_min()
df2d.m_common.get_eps_min()


  #=========================
  # Set input values (replace get by set)
  # ========================
# set 

df2d.m_common.set_ts(1000)

  #=========================
  # run model 
  # ========================
my_model.update_fortran()
my_model.run()
my_model.save_res()









