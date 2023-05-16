#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 11:38:18 2022

@author: livillenave
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 11:51:10 2022

@author: livillenave
"""


import dassflow2d as df2d



# initialise fortran instance, and python corrponding data
my_model = df2d.DassFlowModel(bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/", run_type = "direct") # run_type can be min or direct (grad ?)


my_model.update_fortran()

# =========================
# get bc values
# ========================
my_bc = df2d.wrapping.m_model.get_bc()
print(my_bc.hyd[0].q)
print(my_bc.hyd[0].t)

print(my_bc.hpresc[0].h)
print(my_bc.hpresc[0].t)


# =========================
# set bc values
# ========================
my_bc.hyd[0].q[:] = 500
my_bc.hpresc[0].h[:] = 10



# =========================
# feed fortran kernel with new information
# =========================

df2d.wrapping.m_model.set_bc(my_bc)


# =========================
# run model 
# ========================
my_model.run()
my_model.save_res()









