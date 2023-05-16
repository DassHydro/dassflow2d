#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 15:06:49 2022

@author: livillenave
"""


#-------------------------------------------#
# import required librairies
#-------------------------------------------#
import smash
import numpy as np
import matplotlib.pyplot as plt
import numpy as np

def build_case(nb_timestep=72,                          # number of timestep
               dt = 3600,                               # timestep value
               start_time="2020-01-01 00:00", 
               end_time="2020-01-02 00:00",
               max_rain=12):     # rain amplitude maximum (triangle rain is generated)
    
    print("always check that nb_timstep = [start_time-end_time]/dt")
    
    # define setup
    setup = {
        "dt": dt,
        "start_time": start_time,
        "end_time": end_time,
    }
    
    
    # mesh definition
    dx = 1_000
    (nrow, ncol) = (10, 10)
    
    mesh = {
    "dx": dx,
    "nrow": nrow,
    "ncol": ncol,
    "ng": 1,
    "nac": nrow * ncol,
    "area": nrow * ncol * (dx ** 2),
    "gauge_pos": np.array([9, 9], dtype=np.int32),
    "code": np.array(["Practice_case"])}
    
        
    mesh["flwdir"] = np.array(
    
        [    
        [4, 5, 5, 5, 5, 5, 5, 5, 5, 5],    
        [3, 4, 5, 5, 5, 5, 5, 5, 5, 5],    
        [3, 3, 4, 5, 5, 5, 5, 5, 5, 5],    
        [3, 3, 3, 4, 5, 5, 5, 5, 5, 5],    
        [3, 3, 3, 3, 4, 5, 5, 5, 5, 5],    
        [3, 3, 3, 3, 3, 4, 5, 5, 5, 5],
        [3, 3, 3, 3, 3, 3, 4, 5, 5, 5],    
        [3, 3, 3, 3, 3, 3, 3, 4, 5, 5],    
        [3, 3, 3, 3, 3, 3, 3, 3, 4, 5],    
        [3, 3, 3, 3, 3, 3, 3, 3, 3, 4],    
        ],
    
        dtype=np.int32,
    
    )
    
    
    
    mesh["drained_area"] = np.array(
    
        [    
               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],    
               [1, 4, 2, 2, 2, 2, 2, 2, 2, 2],    
               [1, 2, 9, 3, 3, 3, 3, 3, 3, 3],    
               [1, 2, 3, 16, 4, 4, 4, 4, 4, 4],    
               [1, 2, 3, 4, 25, 5, 5, 5, 5, 5],    
               [1, 2, 3, 4, 5, 36, 6, 6, 6, 6],    
               [1, 2, 3, 4, 5, 6, 49, 7, 7, 7],    
               [1, 2, 3, 4, 5, 6, 7, 64, 8, 8],    
               [1, 2, 3, 4, 5, 6, 7, 8, 81, 9],    
               [1, 2, 3, 4, 5, 6, 7, 8, 9, 100],    
            ],
            dtype=np.int32,    
        )

    
    
    
    ind_path = np.unravel_index(np.argsort(mesh["drained_area"], axis=None),    
         mesh["drained_area"].shape)
    
    
    
    mesh["path"] = np.zeros(shape=(2, mesh["drained_area"].size),    
        dtype=np.int32)
    
    
    
    mesh["path"][0, :] = ind_path[0]    
    mesh["path"][1, :] = ind_path[1]
    
    # nb_timestep = model.input_data.prcp.shape[2]
        # rainfall generation
    prcp = np.zeros(shape=nb_timestep, dtype=np.float32)
    tri = np.linspace(0, 24, 10)
    prcp[0:10] = tri
    prcp[9:19] = np.flip(tri)
    
        # potential evapotranspiration generation
    pet=0
    
    data_model = {"setup":setup, "mesh":mesh, "rain":prcp, "pet":pet }
    return(data_model)

    
def run_true_smash(data):
    model = smash.Model(data["setup"], data["mesh"])
    model.input_data.prcp = np.broadcast_to(data["rain"], model.input_data.prcp.shape)
    model.input_data.pet = data["pet"]
    model.run(inplace=True)
    return(model)
    
    
def plot_direct_case(model):
    # plot drained area
    plt.imshow(model.mesh.drained_area, cmap="Spectral");
    plt.colorbar(label="Number of cells");
    plt.title("Practice case - Drained Area");
    plt.show()
    
    
    # ---------------- parameters ---------------
    
    
    # ---------------- forcings ---------------
    # plot rain
    plt.plot(model.input_data.prcp[0,0,:]);
    plt.grid(alpha=.7, ls="--");
    plt.xlabel("Time step");
    plt.ylabel("Precipitation $(mm/h)$");
    plt.show()
    
    # plot potential evapotranspiration
    plt.plot(model.input_data.pet[0,0,:]);
    plt.grid(alpha=.7, ls="--");
    plt.xlabel("Time step");
    plt.ylabel("Evapotranspiration $(mm/h)$");
    plt.show()
    
    
    # ---------------- forcings ---------------
    
    plt.plot(model.output.qsim[0,:]);
    plt.grid(alpha=.7, ls="--");
    plt.xlabel("Time step");
    plt.ylabel("Simulated discharge $(m^3/s)$")
    plt.show()





def optimize_q_dassflow(model1, model2, 
                        dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap",
                        ts = 40000):    
    
    # infer streamflows using dassflow
    os.chdir( f"{dassflow_dir}/code/")
    os.system("make cleanres cleanmin")
    
    os.system(f"rm -r {dassflow_dir}/code/bin_A/msh/*")
    inf_model = df2d.dassflowmodel(bin_dir =  f"{dassflow_dir}/code/bin_A", hdf5_path = f"{dassflow_dir}/code/bin_A/res/simu.hdf5" , run_type = "min") # initialise fortran/python instance
    
    df2d.wrapping.m_common.set_ts(ts)
    df2d.wrapping.m_common.set_use_obs(1)
    df2d.wrapping.m_model.set_c_hydrograph(1)
    #df2d.wrapping.m_common.set_max_nt_for_adjoint(2)
    
    inf_model.init_all()
    
    
    # force bc values from smash simulation
    my_bc = df2d.wrapping.m_model.get_bc()
    my_bc.hyd[0].q = model1.output.qsim.copy()
    my_bc.hyd[1].q = model2.output.qsim.copy()
    df2d.wrapping.m_model.set_bc(my_bc)
    df2d.wrapping.m_common.set_restart_min(100)
    
    inf_model.run()
    inf_model.save_res()
    print("inference done")
    
    
def get_inference(dassflow_dir):
    bin_dir = f"{dassflow_dir}/code/bin_A"
    min_dir = f"{bin_dir}/min"
    
    files = os.listdir(min_dir)
    q1_file = [x for x in files if  "hydrograph_001" in x]
    q2_file = [x for x in files if  "hydrograph_002" in x]
    
    id_q1 = [int(x[15:18]) for x in q1_file]
    id_q2 = [int(x[15:18]) for x in q2_file]
    
    
    res = dict()
    
    res["q1"] =  np.loadtxt(f"{min_dir}/{q1_file[-1]}")
    res["q2"] =  np.loadtxt(f"{min_dir}/{q2_file[-1]}")
    
    return(res)


