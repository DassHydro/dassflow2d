#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:12:52 2023

@author: livillenave
"""

res = dassflow_model.outputs.all_res
keys = dassflow_model.outputs.all_times


best_len = 0
for my_key in all_keys:
    new_len = np.count_nonzero(res[my_key]["h"])
    if new_len > best_len:
        best_len = new_len
        best_key = my_key
        
        
best_h = np.asanyarray(res[best_key]["h"])


pv.global_theme.color = 'white'
plotter = pv.Plotter(off_screen = False, notebook = False)
actor2 = plotter.add_mesh(mesh = dassflow_model.meshing.mesh_pyvista, show_edges = False, 
                          scalars = best_h)
plotter.show(cpos = "xy")


#
#pv.global_theme.color = 'white'
#plotter = pv.Plotter(off_screen = False, notebook = False)
#actor2 = plotter.add_mesh(mesh = dassflow_model.meshing.mesh_pyvista, show_edges = True, 
#                          scalars = res[best_key]["bathy"])
#plotter.show(cpos = "xy")


# 0 pas d'eau
# 1 que de l'eau simulée
# 2 que de l'eau observée
# 3 les 2 sont d'accord

id_simu_flooded = np.where(best_h >0)[0]
id_obs_flooded = np.where(is_flooded >0)[0]

common_flooded = best_h.copy()
common_list = [c for c in id_simu_flooded if c in id_obs_flooded]
res = np.asanyarray(common_list)


final = np.zeros(shape = len(best_h))
final[id_simu_flooded] = 1
final[id_obs_flooded] = 2
final[res] = 3

pv.global_theme.color = 'white'
plotter = pv.Plotter(off_screen = False, notebook = False)
actor2 = plotter.add_mesh(mesh = dassflow_model.meshing.mesh_pyvista, show_edges = False, 
                          scalars = final)
plotter.show(cpos = "xy")

