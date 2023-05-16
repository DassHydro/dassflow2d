#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 13:14:14 2023

@author: livillenave
"""

import numpy as np
import smash
import dassflow2d as df2d
import pickle
import datetime


import shutil
import os

def yyyymmdd_hhmm_TO_datetime(yyyymmdd_hhmm) :
    
    res_date = datetime.date(year = int(yyyymmdd_hhmm[0:4]), month = int(yyyymmdd_hhmm[5:7]), day =int(yyyymmdd_hhmm[8:10] ) )
    res_time = datetime.time(hour = int(yyyymmdd_hhmm[11:13]), minute = int(yyyymmdd_hhmm[14:16]), second =00 )
    
    res = datetime.datetime.combine(res_date, res_time)
    return(res)
    
    
    
# ============================ #
# PARAMETERS
# ============================ #
    
# ------------------- #
# Paths  
# ------------------- #

demo_dir = "/home/livillenave/Documents/distant/Demo"
source_bin_dir = f"{demo_dir}/bin6"
bin_dir = f"{demo_dir}/built_bin3"
forcing_dir = "/home/livillenave/Documents/data/FORCING"


#  -----------------------------------------   #
#  Import source bin into  calculation bin_dir #
#  -----------------------------------------   #

if  os.path.exists(f'{bin_dir}'):
    shutil.rmtree(bin_dir)
else:
    os.mkdir(bin_dir)
    
shutil.copytree(src = source_bin_dir, dst = bin_dir)


# ----------------------------------------- #
# time parameter
# ----------------------------------------- #

# >>> Hydrology ----
start_hydrology =   '2018-10-10 12:00'
stop_hydrology =    '2018-10-20 20:00'
dt_hydrology   =    3600

# >>> Hydraulic ----
start_hydraulic =  '2018-10-14 22:30'
stop_hydraulic =  '2018-10-16 22:00'

start_heating_hydraulic= '2018-10-14 20:00'
stop_heating_hydraulic =  '2018-10-14 22:30'
dt_hydraulic = "adapt_dt" # NOT USED YET

# to datetime
start_hydrology_datetime = yyyymmdd_hhmm_TO_datetime(yyyymmdd_hhmm=start_hydrology)
stop_hydrology_datetime = yyyymmdd_hhmm_TO_datetime(yyyymmdd_hhmm=stop_hydrology)
start_hydraulic_datetime = yyyymmdd_hhmm_TO_datetime(yyyymmdd_hhmm=start_hydraulic)
stop_hydraulic_datetime = yyyymmdd_hhmm_TO_datetime(yyyymmdd_hhmm=stop_hydraulic)

start_heating_hydraulic_datetime = yyyymmdd_hhmm_TO_datetime(yyyymmdd_hhmm = start_heating_hydraulic )
stop_heating_hydraulic_datetime  = yyyymmdd_hhmm_TO_datetime(yyyymmdd_hhmm = stop_heating_hydraulic )

# store all in 1 big variable
all_date = dict()
all_date["rr"] = dict()
all_date["rr"]["start"] = start_hydrology_datetime
all_date["rr"]["stop"] = stop_hydrology_datetime
all_date["rr"]["dt"] = dt_hydrology
all_date["hy"] = dict()
all_date["hy"]["start"] = start_hydraulic_datetime
all_date["hy"]["stop"] = stop_hydraulic_datetime
all_date["hy"]["dt"] = dt_hydraulic
all_date["hy"]["start_heating"] = start_heating_hydraulic_datetime
all_date["hy"]["stop_heating"] = stop_heating_hydraulic_datetime


# ----------------------------------------- #
# Spatial parameter
# ----------------------------------------- #
# smash_scale     ---> defined due to provided all_inflows that give which smash cells are coupled to dafflow
# dassflow_scale  ---> defined by provided meshing

path_flow_direction =  f"{demo_dir}/source_file/DEM/LEBLOIS_DATA/10/flow_dir.asc"



# ----------------------------------------- #
# SMASH ADDITIONAL PARAMETER
# ----------------------------------------- #
# Structure (among "gr-a", ...)
structure = 'gr-a'
# observed discharge
read_qobs = True

if dt_hydrology == 900:
    qobs_directory = f"{forcing_dir}/DEBIT/15"
elif dt_hydrology == 3600:
    qobs_directory = f"{forcing_dir}/DEBIT/60"
else:
    print("dt_hydrology is not 900 or 3600, i cant work")
# precipitation
read_prcp = True
prcp_directory = f"{forcing_dir}/PLUIE/J+1/1H"
prcp_conversion_factor = 0.1
# Evapotranspiration
read_pet=True
pet_directory=f"{forcing_dir}/ETP-SFR-FRA-INTERA_L93"
pet_conversion_factor=1
daily_interannual_pet=1
# various variable to save
save_qsim_domain = True
save_production_level  = 'TRUE'
save_others_variables  = 'TRUE'
save_net_rainfall  = 'TRUE'

# ----------------------------------------- #
# DASSFLOW ADDITIONAL PARAMETER
# ----------------------------------------- #
dassflow_dtw = 360
dassflow_dtp = 360

dassflow_dtw_heating = 3600
dassflow_dtp_heating = 3600
# ts automaticaly calculated from start_hydrolic to stop_hydraulic

#+++++++++++++++++#
# Gauge parameter
##+++++++++++++++++#
#x_gauge = [662_594 , 647_966, 653_777, 639_767,]# 636_866, 625_654, 649_088]
## 639_767, 636_866, 625_654, 649_088
#y_gauge = [6_233_537, 6_235_059,6_237_208, 6_210_365 ]# , 6_221_767 , 6_242_738,6_237_846]
## 6_210_365  , 6_221_767 , 6_242_738,6_237_846
#area =    [3_192_579_999, 1_828_800_000,    242_400_000 , 15_300_000 ]#, 197_000_000 , 211_200_000, 935_500_000 ]
## , 15_300_000 , 197_000_000 , 211_200_000, 935_500_000
#code_gauge = ["Y1422030", "Y1232010",   "Y1415020","Y1141150"]#,"Y1205010","Y1314010", "Y1364010",]
#                                       # "Y1422020", aude a marseillete
#                                        #"Y1422030", aude a trebes
#                                       # "Y1524020",
#                                      #  "Y1605050", la cesse a mirepassset
#                                      #  "Y1612020" l'aude a moussan
                                      #]

# "Y1141150","Y1205010","Y1314010", "Y1364010",


x_gauge = [662_594 , 647_966, 653_777, 649_088,]# 636_866, 625_654, ]
# 639_767, 636_866, 625_654, 649_088
y_gauge = [6_233_537, 6_235_059,6_237_208, 6_237_846 ]# , 6_221_767 , 6_242_738,]
# 6_210_365  , 6_221_767 , 6_242_738,6_237_846
area =    [3_192_579_999, 1_828_800_000,    242_400_000 , 935_500_000 ]#, 197_000_000 , 211_200_000,  ]
# , 15_300_000 , 197_000_000 , 211_200_000, 935_500_000
code_gauge = ["Y1422030", "Y1232010",   "Y1415020","Y1364010"]#,"Y1205010","Y1314010", "",]
                                       # "Y1422020", aude a marseillete
                                        #"Y1422030", aude a trebes
                                       # "Y1524020",
                                      #  "Y1605050", la cesse a mirepassset
                                      #  "Y1612020" l'aude a moussan
# ============================ #
# PARAMETERS
# ============================ #



#--------------------------------
# Perform smash simulation
#--------------------------------
smash_options={
        'structure' : structure,
        'dt'        : dt_hydrology,
		'start_time': start_hydrology,
		'end_time'  : stop_hydrology,
        'read_qobs':read_qobs,
        'qobs_directory':qobs_directory,
        'read_prcp':read_prcp,
        'prcp_directory':prcp_directory,
        'prcp_conversion_factor':prcp_conversion_factor,
        'read_pet':read_pet,
        'pet_conversion_factor':pet_conversion_factor,
        'pet_directory':pet_directory,    
        'daily_interannual_pet':daily_interannual_pet,
        'save_qsim_domain':save_qsim_domain,
        "save_production_level":save_production_level,
        "save_others_variables":save_others_variables,
        "save_net_rainfall":save_net_rainfall }


meshing = smash.mesh.meshing.generate_mesh(
        path =f'{path_flow_direction}',#f'/home/livillenave/Documents/data/DONNEES/LEBLOIS_DATA/10/flow_dir.asc',#f'{coupled_dir}/{my_scale}/SMASH_files/FLOW.asc', 
        epsg=2154,
        x=x_gauge,
        y=y_gauge,
        area= area, #3_192_100_000, 
        #bbox = (601_435, 663_000 ,6_159_000,6_258_977 ),
        code =code_gauge)
       
                    
smash_model = smash.Model(smash_options, mesh = meshing)

smash_model.optimize(mapping="uniform",
              #algorithm='sbs', 
                                   gauge=code_gauge,
                                   control_vector=["cp", "cft", "exc", "lr"], 
                                   jobs_fun='nse', 
                                   inplace = True, 
                                   options={'maxiter': 2})
# launch calibration
smash_model.optimize(mapping="distributed",
                                   gauge=code_gauge,
                                   control_vector=["cp", "cft", "exc", "lr"], 
                                   jobs_fun='nse', 
                                   inplace = True, 
                                   options={'maxiter': 2})
smash_model.run()

# ----------------------------
# Transfert hydrograph (via textfile)
# ----------------------
                                   
# qin_dict, dictionary, the key is the id of the hydrogram
# write_dir: (absolute) path where to write the file 
# return; write the hydrograph.txt file where precised 
def write_hydrograph(qin_dict, write_dir, dt ):
    time = np.arange(start = 0, stop = dt * len(qin_dict[0]), step = dt) 
    with open(f'{write_dir}/hydrograph.txt', 'w') as f:
        f.write("#comment\n")
        f.write("#comment\n")
        f.write("#comment\n")
        f.write(str(len(qin_dict))+"\n")
        for hyd in qin_dict.values():
            f.write("#comment\n")
            f.write("#comment\n")
            f.write("#comment\n")
            f.write(f"{len(hyd)}\n")
            for i in range(len(time)):
                f.write(f"{time[i]} {hyd[i]}\n")
                
        f.write("#comment\n")
                
 


#+++++++++#
# Dassflow 2d mesh (for having bc information)
#+++++++++#
with open(f'{source_bin_dir}/filename.pickle', 'rb') as handle:
    all_inflow= pickle.load(handle)
df2d_mesh =  read_mesh_from_textfile(path = f"{source_bin_dir}/final_mesh.geo",
                               read_boundary = True)




#+++++++++#
# get coupled data 
#+++++++++#

all_id_inflow = np.arange(df2d_mesh["boundaries"]["header"]["nb_inflow"])
all_id_outflow = np.arange(df2d_mesh["boundaries"]["header"]["nb_outflow"])
qin = dict()
for mykey in all_id_inflow :
    qin[mykey] = np.zeros((smash_model.output.qsim_domain.shape[2]))


for key, value in all_inflow.items():  
    my_id  = value["id"]
    if value["bc_group_dassflow"] in all_id_inflow:
         qin[value["bc_group_dassflow"]][:] =  qin[value["bc_group_dassflow"]][:] + smash_model.output.qsim_domain[my_id[0],my_id[1]][:]
         
         
#-------------#
# plot
#-------------#

# >>> plot main_flow decomposition
all_fig = []
all_axes= []

for i in range(len(all_id_inflow)):
    all_fig.append(plt.figure(figsize = (12,12)))
    all_fig[i].suptitle(f"Inflow discharges (all smash pixels), group: {i}")
    all_axes.append( plt.gca())

for key, value in all_inflow.items():  
    my_id  = value["id"]
    if value["bc_group_dassflow"] in all_id_inflow:         
         for i in range(len(all_id_inflow)):    
             group = all_id_inflow[i]
             if value["bc_group_dassflow"] ==group:
                 all_axes[i].plot(  smash_model.output.qsim_domain[my_id[0],my_id[1]][:], label = f"ID: ({my_id[0]},{my_id[1]})")
#ax1.plot(qin[0][:], label = f"Sum of all", c = "black" )
[all_axes[i].plot(qin[i][:], label = f"Sum of all", c = "black" ) for i in range(len(all_id_inflow))]        
[x.legend() for x in all_axes]
all_fig[0]
all_fig[1]
all_fig[2]
plt.close()


# >>> compare sum main inflow to output
for key, val in qin.items():
    plt.plot(val, label =f"discharge {key}")
plt.plot(qin[0]+qin[1]+qin[2], label =f"sum inflows")
plt.plot(smash_model.output.qsim[0,:], c= "black", label ="inf");
plt.legend()
plt.title("Repartition of main inflows")
plt.show()
plt.close()             
 
#-------------#
# >>> END plot
#-------------#
# Transfer hydrograph to dassflow
#--------------#

range_date_rr = []
for x in range (0, np.shape(qin[0])[0]):
    range_date_rr.append(all_date["rr"]["start"] + datetime.timedelta(seconds = (x+1) * dt_hydrology))
print('check dateList[-1] == all_date["rr"]["stop"] :', 
      range_date_rr[-1] == all_date["rr"]["stop"])
def  subset_qin_dict(qin_dict, range_date_rr, 
                   start_subset, stop_subset):
        to_keep_date = [(x>=start_subset and x<=stop_subset) for x in range_date_rr]
        subset_dates = [x  for x in range_date_rr if (x>=start_subset and x<=stop_subset)] 
        new_qin = qin_dict.copy()
        for key, val in new_qin.items():
            new_qin[key] = val[to_keep_date]          
        return(new_qin, subset_dates)

(new_qin, subset_dates) = subset_qin_dict(qin_dict = qin, range_date_rr = range_date_rr, 
                   start_subset = all_date["hy"]["start"], stop_subset = all_date["hy"]["stop"])

#write_hydrograph(qin_dict = new_qin, write_dir = bin_dir , dt = dt_hydrology )
                


(qin_heating, subset_dates_heating) = subset_qin_dict(qin_dict = qin, range_date_rr = range_date_rr, 
                   start_subset = all_date["hy"]["start_heating"], stop_subset = all_date["hy"]["stop_heating"])


# -----------------------------
# Perform dassflow simulation
# -----------------------------

# disrt generate appropriate hydrograph for prepared variable
write_hydrograph(qin_dict = qin_heating, write_dir = bin_dir , dt = dt_hydrology )
# run on calm period so that we can get rid off initialization issues
print(bin_dir)
dassflow_model = df2d.dassflowmodel(bin_dir = bin_dir, 
                           hdf5_path=bin_dir+"/res/simu.hdf5", 
                           run_type= "direct",  
                           clean = True)

df2d.wrapping.m_common.set_mesh_name("final_mesh.geo")

df2d.wrapping.m_common.set_ts( dt_hydrology * ( len(qin_heating[0]) -1) )
df2d.wrapping.m_common.set_dtw( dassflow_dtw_heating )
df2d.wrapping.m_common.set_dtp( dassflow_dtp_heating )
dassflow_model.init_all()

dassflow_model.kernel.dof.h[:] = 0
# manually set correctly bc
dassflow_model.kernel.dof.h[3777] =   1
dassflow_model.kernel.dof.h[11185] =  dassflow_model.kernel.dof.h[10609] = dassflow_model.kernel.dof.h[3777]

dassflow_model.kernel.dof0 = dassflow_model.kernel.dof
dassflow_model.run()

#dassflow_model.save_all()
#dassflow_model.post = df2d.core.output.Post(boundary_metadata = dassflow_model.boundary.metadata, bin_dir = dassflow_model.bin_dir)




# Manualy clean directories that need to be removed (because clean = False in next call of dassflowmodel to keep the restart.bin file)
#

res_dir = bin_dir+"/res/"
shutil.rmtree(res_dir)
os.mkdir(res_dir)

#--------------------------------
# Perform dassflow simulation
#--------------------------------


write_hydrograph(qin_dict = new_qin, write_dir = bin_dir , dt = dt_hydrology )

print(bin_dir)
dassflow_model = df2d.dassflowmodel(bin_dir = bin_dir, 
                           hdf5_path=bin_dir+"/res/simu.hdf5", 
                           run_type= "direct",  
                           clean = True)  # WARNING: clean = true

df2d.wrapping.m_common.set_mesh_name("final_mesh.geo")

df2d.wrapping.m_common.set_ts(dt_hydrology * len(new_qin[0]) )
df2d.wrapping.m_common.set_dtw( dassflow_dtw )
df2d.wrapping.m_common.set_dtp( dassflow_dtp )
dassflow_model.init_all()
dassflow_model.kernel.dof.h[:] =  0

id_inflow_cell= [3778, 11186, 10610 ]
for my_id in id_inflow_cell:
    dassflow_model.kernel.dof.h[my_id-1] = 1
    
# manually set correctly bc
dassflow_model.kernel.dof.h[100] = 1
dassflow_model.kernel.dof.h[13752]  = 1
dassflow_model.kernel.dof.h[3332] = 1

dassflow_model.kernel.dof0.h[:]  = dassflow_model.kernel.dof.h[:] 
#dassflow_model.kernel.dof.h[:] = dassflow_model.kernel.dof0.h[:] = 1
dassflow_model.run()
dassflow_model.save_all()

dassflow_model.post = df2d.core.output.Post(boundary_metadata = dassflow_model.boundary.metadata, bin_dir = dassflow_model.bin_dir)




#----------------------------#
# PLOTS
#----------------------------#

# - - - - - -#
# water heigh
# - - - - - - #

best_diff = 1000 #todo
all_keys = [x for x in dassflow_model.outputs.all_res.keys()]
old_key= all_keys[-1]
for i in dassflow_model.outputs.all_res.keys():
    
    previous_h = np.asanyarray(dassflow_model.outputs.all_res[old_key]["h"])
    h_array = np.asanyarray(dassflow_model.outputs.all_res[i]["h"])
    diff = previous_h[:] - h_array
    diff = np.linalg.norm(diff)
   # print(i, "diff=", diff)
    if diff<best_diff:
        print(i, "diff=", diff)
        best_diff =  diff
    dassflow_model.meshing.mesh_pyvista.plot(scalars = dassflow_model.outputs.all_res[i]["h"], 
                                    cpos = "xy", show_edges = False)
    old_key = i
    
    
# - - - - - - #
# inflow + outflow discharge
# - - - - - - #

fig = plt.figure(figsize = (12,12))
fig.suptitle("Inflow discharges (sum_q)")
ax = plt.gca()

ax.scatter( x = dassflow_model.post.sum_q_inflow["0"][:,0],
           y = dassflow_model.post.sum_q_inflow["0"][:,1], marker = "o", 
           label = "group 1", c = "red")
ax.scatter( x = dassflow_model.post.sum_q_inflow["1"][:,0],
           y = dassflow_model.post.sum_q_inflow["1"][:,1], marker = "+", 
           label = "group 2" , c = "blue")
ax.scatter( x = dassflow_model.post.sum_q_inflow["2"][:,0],
           y = dassflow_model.post.sum_q_inflow["2"][:,1], marker = "x", 
           label = "group 3", c = "green")
ax.legend()



fig = plt.figure()
fig.suptitle("outflow discharges")
ax = plt.gca()

ax.scatter( x = dassflow_model.post.sum_mass_flux_outflow["0"][:,0],
           y = dassflow_model.post.sum_mass_flux_outflow["0"][:,1], marker = "o", 
           label = "group 1", c = "red")
ax.set_yscale('log')
ax.legend()




# - - - - - - #
# Inflow
# - - - - - - #

for group in all_id_inflow+1:
    #group = 3
    post_group = str(group -1)
    
    fig = plt.figure(figsize = (20,12))
    fig.suptitle(f"Inflow discharge {group}")
    ax = plt.gca()
    
    ax.plot( dassflow_model.post.sum_mass_flux_inflow[post_group][:,0],
                dassflow_model.post.sum_mass_flux_inflow[post_group][:,1], 
               label = "sum_mass_flux_inflow", c = "red", linestyle = "-")
    ax.plot( dassflow_model.post.sum_q_inflow[post_group][:,0],
                dassflow_model.post.sum_q_inflow[post_group][:,1], 
               label = "sum_q_inflow", c = "blue", linestyle = "--")
    ax.plot(np.arange(start = 0, stop = dt_hydrology * len(new_qin[0]), step  = dt_hydrology),
                new_qin[group-1], 
               label = "qin_hydrology", c = "black", linestyle = "dotted")
    
    id_boundary_cell = dassflow_model.boundary.corresp["discharg1"][group]["id_cell"]
    id_boundary_cell = id_boundary_cell-1 # pass to python index ?
    length_edge = dassflow_model.meshing.mesh_fortran.edge[dassflow_model.boundary.corresp["discharg1"][group]["id_edge"]].length
    
    q_inflow_sum = np.zeros(shape = (len(dassflow_model.outputs.all_res)))
    
    h = []
    u = []
    v = []
    for key, val in dassflow_model.outputs.all_res.items():
        h.append(np.asarray(val["h"][id_boundary_cell]))
        u.append(np.asarray(val["u"][id_boundary_cell]))
        v.append(np.asarray(val["v"][id_boundary_cell]))    
    h = np.asanyarray(h)
    u = np.asanyarray(u) 
    v = np.asanyarray(v)
    q = length_edge * h[:] * np.sqrt(u[:]**2 + v[:]**2)
    q_inflow_sum[:] = q_inflow_sum[:] + q[:,0] #IMPORTANT: built q_inflow_sum during the plot
    ax.plot(np.arange(start = 0, stop = dassflow_dtw * len(q), step  = dassflow_dtw),
                q, 
               label = "q for DOF variable (python)", c = "purple")
    
    sum1=sum( dassflow_model.post.sum_mass_flux_inflow[post_group][:,1] ) / len( dassflow_model.post.sum_mass_flux_inflow[post_group][:,1])
    sum2=sum( dassflow_model.post.sum_q_inflow[post_group][:,1]         ) / len(dassflow_model.post.sum_q_inflow[post_group][:,1]         )
    sum3=sum( new_qin[group-1] ) / len(new_qin[group-1] )
    sum4=sum( q ) / len(q)
    ax.text(0.5, 0.9 * max(q), s= f" sums/nb elements \n sum sum_mass_flux_inflow = {sum1}  \n sum sum_q_inflow  = {sum2} \n sum qin_hydrology =  {sum3} \n sum q dof =  {sum4}")
    
    ax.legend()



# - - - - - - #
# Outflow
# - - - - - - #
group = all_id_inflow[-1]+2
post_group = str(0)

fig = plt.figure(figsize = (12,12))
fig.suptitle(f"Outflow discharge {group}")
ax = plt.gca()

ax.plot( dassflow_model.post.sum_mass_flux_outflow[post_group][:,0],
            dassflow_model.post.sum_mass_flux_outflow[post_group][:,1], 
           label = "sum_mass_flux_outflow", c = "red", linestyle = "-")
ax.plot( dassflow_model.post.sum_q_outflow[post_group][:,0],
            dassflow_model.post.sum_q_outflow[post_group][:,1], 
           label = "sum_q_outflow", c = "blue", linestyle = "--")
#ax.plot(np.arange(start = 0, stop = dt_hydrology * len(new_qin[0]), step  = dt_hydrology),
#            new_qin[group-1], 
#           label = "qin_hydrology", c = "black", linestyle = "dotted")

id_boundary_cell = dassflow_model.boundary.corresp["transm"][group]["id_cell"]
id_boundary_cell = id_boundary_cell-1 # pass to python index ?
length_edge = dassflow_model.meshing.mesh_fortran.edge[dassflow_model.boundary.corresp["transm"][group]["id_edge"]].length


h = []
u = []
v = []
for key, val in dassflow_model.outputs.all_res.items():
    h.append(np.asarray(val["h"][id_boundary_cell]))
    u.append(np.asarray(val["u"][id_boundary_cell]))
    v.append(np.asarray(val["v"][id_boundary_cell]))    
h = np.asanyarray(h)
u = np.asanyarray(u)
v = np.asanyarray(v)
q = length_edge * h[:] * np.sqrt(u[:]**2 + v[:]**2)
q_outflow = q
ax.plot(np.arange(start = 0, stop = max(dassflow_model.post.sum_q_outflow[post_group][:,0]), step  =  max(dassflow_model.post.sum_q_outflow[post_group][:,0])/ len(q_outflow)  )[0:-1],  # WARNING DIRTY DEBUX (le x[-1])
            q_outflow, 
           label = "q for DOF variable (python)", c = "purple")
ax.text(0.5, 0.9 * max(q_outflow), 
        s= f"sum sum_mass_flux_outflow = {sum(dassflow_model.post.sum_mass_flux_outflow[post_group][:,1])} \n sum sum_mass_flux_outflow  = {sum(dassflow_model.post.sum_q_outflow[post_group][:,1])} \n sum q DOF  =  {sum(q_outflow)}")

ax.legend()



# - - - - - - #
# Compare inflow - outflow dassflow
# - - - - - - #

qin = dassflow_model.post.sum_mass_flux_inflow["0"][:,1]+ dassflow_model.post.sum_mass_flux_inflow["1"][:,1]+ dassflow_model.post.sum_mass_flux_inflow["2"][:,1]
qout = dassflow_model.post.sum_mass_flux_outflow["0"][:,1]

fig = plt.figure(figsize = (12,12))
fig.suptitle(f"sum_mass_flux :  Inflow & outflow")
ax = plt.gca()
ax.plot(qin, label = "sum inflows")
ax.plot(qout, label = "outflow")
ax.legend()

qin = dassflow_model.post.sum_q_inflow["0"][:,1]+ dassflow_model.post.sum_q_inflow["1"][:,1]+ dassflow_model.post.sum_q_inflow["2"][:,1]
qout = dassflow_model.post.sum_q_outflow["0"][:,1]

fig = plt.figure(figsize = (12,12))
fig.suptitle(f"sum_q_flux :  Inflow & outflow")
ax = plt.gca()
ax.plot(qin, label = "sum inflows")
ax.plot(qout, label = "outflow")
ax.text(0.5, max(qin), s= f" sum qin = {sum(qin)} \n sum qout =  {sum(qout)}")
ax.legend()



fig = plt.figure(figsize = (12,12))
fig.suptitle(f"DOF Q :  Inflow & outflow")
ax = plt.gca()
ax.plot(q_inflow_sum, label = "q inflows")
ax.plot(q_outflow, label = "q outflow")
ax.text(0.5, max(qin), s= f" sum qin = {sum(q_inflow_sum)} \n sum qout =  {sum(qout)}")
ax.legend()



#--------------------     #
# water volume evolution
#------------------------ #
time_post = dassflow_model.post.water_vol[:,0]
water_volume = dassflow_model.post.water_vol[:,1]

time_res = np.arange(start = 0, stop = dassflow_dtw * len(q), step  = dassflow_dtw)

fig = plt.figure(figsize = (12,12))
fig.suptitle(f"water volume from post file 'water_vol.dat' ")
ax = plt.gca()
ax.plot(time_post,water_volume, label = "post.water_vol")
ax.text(0.5, 0.9* max(water_volume), s= "sum water volume ")
ax.legend()


fig = plt.figure(figsize = (12,12))
fig.suptitle(f"water volume diff (red curve) compared to DOF inflow & outflow")
ax1 = plt.subplot()
l1, = ax1.plot(time_post[1:],water_volume[1:]-water_volume[:-1], label = "post.water_vol diff", color='red')
ax2 = ax1.twinx()
l2, = ax2.plot(time_res,q_inflow_sum, label = "q inflows", color='blue')
l3, = ax2.plot(time_res,q_outflow, label = "q outflow", color='green')
ax1.legend(loc=0);ax2.legend(loc=2)
ax1.set_ylabel("Water volume (m3?)");ax2.set_ylabel("dicharge [m3/s]")



qdiff = q_inflow_sum[:]-q_outflow[:].squeeze()[:]

fig = plt.figure(figsize = (12,12))
fig.suptitle(f"water volume calculated from dof inflow, outflow diff")
ax1 = plt.subplot()
l1, = ax1.plot(time_post[1:],water_volume[1:]-water_volume[:-1], label = "post.water_vol diff", color='red')
ax2 = ax1.twinx()
l2, = ax2.plot(time_res,qdiff, label = "q inflows-outflow", color='blue')
ax1.legend(loc=0);ax2.legend(loc=2)
ax1.set_ylabel("Water volume (m3?)");ax2.set_ylabel("dicharge [m3/s]")


all_dof_volume = []

for my_key in   dassflow_model.outputs.all_res.keys() :
    dof_volume=[]
    for id_cell in range(dassflow_model.meshing.mesh_fortran.nc):
        h = dassflow_model.outputs.all_res[my_key]["h"][id_cell]
        surface =  dassflow_model.meshing.mesh_fortran.cell[id_cell].surf
        volume = h * surface
        dof_volume.append(volume)
    all_dof_volume.append(dof_volume)
    
sum_dof_volume = [np.sum(x) for x in all_dof_volume]
sum_dof_volume= np.asanyarray(sum_dof_volume)

fig = plt.figure(figsize = (12,12))
fig.suptitle(f"Delta water volume at each timestep")
ax1 = plt.subplot()
l1, = ax1.plot(time_post[1:],water_volume[1:]-water_volume[:-1], label = "post.water_vol diff", color='red')
ax2 = ax1.twinx()
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("axes", 1.1)) 
l2, = ax2.plot(time_res,qdiff, label = "q inflows-outflow", color='blue')
l3, = ax3.plot(time_res[1:],sum_dof_volume[1:]-sum_dof_volume[:-1], 
               label = "diff dof volume", color='black', marker='o', linestyle='None')
ax1.legend(loc=0);ax2.legend(loc=2);ax3.legend(loc=3)
ax1.set_ylabel("Water volume (m3?)");ax2.set_ylabel("dicharge [m3/s]");
ax3.set_ylabel("DOF water volume")


#plt.plot(dassflow_model.post.water_vol_num_add[:,0],dassflow_model.post.water_vol_num_add[:,1], label = "water vol num add")


