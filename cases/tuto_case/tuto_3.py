#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 11:19:20 2022

@author: livillenave
"""

import dassflow2d as df2d



# main directory
dassflow_dir="/home/livillenave/Documents/distant/dassflow2d-wrap/"
# code directory, where compilation happens
code_dir =  f"{dassflow_dir}/code/"
# bin directory, where simulation happens
bin_dir = f"{code_dir}/bin_A/"


# open interface with fortran kernel
model = df2d.dassflowmodel(dassflow_dir=dassflow_dir, bin_dir=bin_dir )
help(df2d.dassflowmodel)

# main path are stored and accesible
print(model.bin_dir)
print(model.dassflow_dir)
# type of run is stored, and applied ar run() function
print(model.run_type)


# we gathered the configuration data, using previsously seen method:
print(model.config)

# to access fortran kernel value, use the get method:
model.config.get()
# to modify, update the dictionary config, and then use the set method:
model.config["ts"] = 1000
model.config.set()



# then  classical code can be performed
model.kernel.mesh = df2d.wrapping.m_mesh.msh()
model.kernel.dof  = df2d.wrapping.m_model.unk(model.kernel.mesh)
model.kernel.dof0 = model.kernel.dof 
df2d.wrapping.call_model.init_solver(model.kernel)
df2d.wrapping.call_model.init_fortran(model.kernel)
df2d.wrapping.call_model.run(model.kernel, arg = "direct")
df2d.wrapping.call_model.clean_model(model.kernel)
