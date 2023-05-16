  #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 11:24:12 2022

@author: livillenave
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 18:59:18 2022

@author: livillenave
"""
import numpy as np
import dassflow2d
import os 
import h5py
import matplotlib.pyplot as plt
import h5py
import copy # to deal with 
            
class Obs(object):
        """
        Obs Class
        
        --------------------
        contains
        --------------------
        self.bin_dir: str:  bin directory: base directory of simulation
        self.obs_dir_true: str: directory where true observation are stored
        self.obs_dir_sim: str: directory where model generated observation are stored        
        self.obs_stationtxt_path: str: path to "bin_dir/obs.txt" file, which contains obs metadata
        self.station_true: str: dictionary of Station classes, of true data, where keys are integer corresponding to the id of the station and the value is a Station object
        self.station_sim str: dictionary of Station classes, of model generated data, where keys are integer corresponding to the id of the station and the value is a Station object
        
        Note:
            self.station_true exist only if use_obs == 1 in configuration
            self.station_sim exist only if w_obs == 1 in configuration
            
        --------------------
        Methods
        --------------------
        
        __init__ : initialize obs object and source values automaticaly from bin directory 
        source: source files from bin_dir and fills values of self.station_true and self.station_sim, depending on use_obs, w_obs values
        save: save in hdf5 file the value sourced in self
        sim_to_true: copy model generated observation as true observation, can be used in case of twin experiment
        """

            
        def __source_metadata(self):
            """
            __source_metadata: extract station information from {bin_dir}/obs.txt file
        
            --------------------
            CONTAINS
            --------------------
            self: Obs object
            
            --------------------
            Return
            --------------------
            dictionary of : {id_station1:Station object, id_station1:Station object, etc...}
            fills: station.x, station.y, station.dt variables for each station object
            """
        
            all_station = dict()
            
            with open(f"{  self.obs_stationtxt_path}") as f:
                for i in range(10):
                    print(i)
                    lines = f.readline()
                    if "stations" in lines:
                        break
                    
                print(lines)
                nb_station = int(lines.split()[-1])
                
                compteur = 0
                to_store = dict()
                for i in range(nb_station*10):
                    my_line = f.readline()
                    if my_line =="\n":
                        print(f"")
                    else:
                        to_store[compteur] = my_line
                        compteur=compteur+1
                    if compteur == nb_station:
                        break
                        
                for i in range(len(to_store)):
                    my_split = to_store[0].split()
                    my_station = Station(dt=np.asarray(my_split[2], dtype = np.float64), ts = self.ts) # np.asarray(my_split[2]) is dt
                    my_station.x =  np.asarray(my_split[0], dtype = np.float64).copy()
                    my_station.y =  np.asarray(my_split[1], dtype = np.float64).copy()
                    my_station.dt =  np.asarray(my_split[2], dtype = np.float64).copy()
                    all_station[i] = my_station             
                    del my_station
                return(all_station)



                    
        def __source_values(self, station_dictionary_to_fill, base_directory):
            """
            __source_values: extract time series of station from {bin_dir}/obs/obs_station_xxxx.dat file
        
            fills self.station[i] by its  time, h,u,v  parameters
            
            --------------------
            parameter
            --------------------
            base_directory, can be:
                - self.simuobs_dir for observation generated with simulation
                - self.obs_dir for observation generated with simulation
            
            """                        
            for i in station_dictionary_to_fill.keys():   
                j=int(i)+1 # compteur pour le chemin de fichier f"{base_directory}/obs_station_{j:04d}.dat"
                station_dictionary_to_fill[i].source_station_file(pathfile =  f"{base_directory}/obs_station_{j:04d}.dat")
                
            return(station_dictionary_to_fill)
            
            
        def source(self):
            """
            Source all information about observations
            
            Source true observation data if use_obs ==1 and Source model generated observation during the simulation if w_obs ==1
            
            source datas from:
                    - source metadata from "bin_dir/obs.txt" file
                    - true observation data from "bin_dir/obs/*"
                    - model generated observation data from "bin_dir/res/obs/*".
                    
            
            --------------------
            return
            --------------------
            Obs oject is filled from textfiles sourced
            
            """
                
            stations_template = copy.deepcopy(self.__source_metadata()     )        
            
            if (dassflow2d.wrapping.m_common.get_w_obs()  == 1):
                self.station_sim  = self.__source_values(station_dictionary_to_fill = copy.deepcopy(stations_template)  , base_directory =  self.obs_dir_sim)
            if (dassflow2d.wrapping.m_common.get_use_obs()  == 1):
                self.station_true  = self.__source_values(station_dictionary_to_fill =   copy.deepcopy(stations_template)  , base_directory =  self.obs_dir_true)

        def __init__(self, bin_dir,ts=3600):
            
            # define basic param
            self.bin_dir = bin_dir
            self.obs_dir_true = f"{bin_dir}/obs/"                 # true obs
            self.obs_dir_sim = f"{bin_dir}/res/obs/"         # generated obs  with simulation            
            self.obs_stationtxt_path = f"{bin_dir}/obs.txt"       # path to file  obs.txt containing obs metadata
            self.ts = ts


            # check configuration coherence (use_obs, w_obs, and which file should be sourced/used)
            if (dassflow2d.wrapping.m_common.get_w_obs()  == 0 and  dassflow2d.wrapping.m_common.get_use_obs()== 0) :
                print("w_obs == 0 and use_obs = 0 in configuration, no observation considered for this case")
                return()
            else:                               
                if (dassflow2d.wrapping.m_common.get_w_obs()  == 1 and  dassflow2d.wrapping.m_common.get_use_obs() == 0) :
                    print("w_obs == 1 and use_obs = 0 in configuration --> observation written from results are sourced")                                   
                elif (dassflow2d.wrapping.m_common.get_w_obs()  == 1 and  dassflow2d.wrapping.m_common.get_use_obs() == 1) :
                    print("w_obs == 1 and use_obs = 1 in configuration --> observation written from results are sourced")                                  
                elif (dassflow2d.wrapping.m_common.get_w_obs()  == 0 and  dassflow2d.wrapping.m_common.get_use_obs() == 1) :
                    print("w_obs == 0 and use_obs = 1 in configuration --> NOT CHECKED IF FUNCTIONAL")
                    #return()
            self.source()
            
               
        def save(self, hdf5_path):
            """
            Save  Obs object to pre-existing hdf5 file.

            -----------
            parameter
            -----------
            self:    Obs class
            hdf5_path : path to hdf5 file to save

            -----------
            return
            -----------
            hdf5 file filled with data stored in Obs class, this account for: my_hdf5["obs"]["station_true"] // my_hdf5["obs"]["station_sim"]
             my_hdf5["obs"]["station_xxx"][id_station] is a dataset of size (nb timestep, 4), it stores (in order) the variables time,h,u,v.
             the dataset is filled of attributes dt, x,y, they are stored as strings
            """
            
            # create dataset containing table in hdf5 file
            
            f = h5py.File(hdf5_path, "a")
            obs =f.create_group("obs")
            
            if 'station_sim' in dir(self):
                obs_sim = obs.create_group("station_sim")
                for key, value in self.station_sim.items():
                    my_strkey= str(key)
                    station = obs_sim.create_dataset(my_strkey, (len(value.time),4) ) 
                    station[:,0] = value.time
                    station[:,1] = value.h
                    station[:,2] = value.u
                    station[:,3] = value.v
                    station.attrs["dt"] = str(value.dt)
                    station.attrs["x"] = str(value.x)
                    station.attrs["y"] = str(value.y)
                
            if 'station_true' in dir(self):
                obs_true = obs.create_group("station_true")
                for key, value in self.station_true.items():
                    my_strkey= str(key)
                    station = obs_true.create_dataset(my_strkey, (len(value.time),4) ) 
                    station[:,0] = value.time
                    station[:,1] = value.h
                    station[:,2] = value.u
                    station[:,3] = value.v
                    station.attrs["dt"] = str(value.dt)
                    station.attrs["x"] = str(value.x)
                    station.attrs["y"] = str(value.y)
                
            f.close()
       
        def sim_to_true(self, update_self=False):
            """
            copy model generated observation as true observation, can be used in case of twin experiment
            
            ----------
            parameter
            ----------
            update_self: if True, update self.station_true from self.station_sim. Warning, self.station_true must exist and use_obs ==1    
            
            ----------
            return
            ----------
            if update_self is True, update self object,
            always copy files from bin_dir/res/obs/ to bin_dir/obs/
            """
             # replace in self
            if update_self:
                 self.station_true = self.station_sim
             
             # replace in direcories
            os.system(f"cp -r {self.obs_dir_sim}/* {self.obs_dir_true}/")


class Station(object):        
        """
        Station Class, contains all data associated to one measure station.
        
        --------------------
        contains
        --------------------        
        self.time: np.array( nb_writing_timestep ): writting time (s) corresponding to each timestep
        self.h:    np.array( nb_writing_timestep ):  water eight (m) observed at each timestep
        self.u:    np.array( nb_writing_timestep ):  u-velocity (m/s) at each timestep
        self.v:    np.array( nb_writing_timestep ):  v-velocity (m/s) at each timestep
        self.x:    int 1  :  x_coordinates at each timestep (m)
        self.y:    int 1  :  y_coordinates at each timestep (m)
        self.dt:   int 1:  writing timestep (s). Sourced from obs.txt file by default.
        
        --------------------
        methods
        --------------------
        __init__: initialise station class
        source_station_file: source  fortran generated files of observation at a measure station
        """
    
        def __init__(self, ts,  dt=60.): # model.config.ts
            """
            
            __init__: initialise station class
            
            -----------
            parameter
            -----------
            self: Station object
            ts: total time of simulation, can be sourced from dassflowmodel.config["ts"]
            dt: timestep of observation data, can be sourced from obs.txt file or deduced from result files
            
            -----------
            return
            -----------
            initialized station object: the object are created with default sizes (depending on dt and ts values)
            
            """
            self.ts = ts
            self.dt = dt
            self.x = int()
            self.y = int()
            
            #evaluate size
            self.nb_dt = int(self.ts/self.dt)+1            
            self.time = np.zeros(self.nb_dt)
            self.h = np.zeros(self.nb_dt)
            self.u = np.zeros(self.nb_dt)
            self.v = np.zeros(self.nb_dt)
        
        def source_station_file(self, pathfile):
            """
            source station file, it must be file of format: "obs_station_xxx.dat", they are fortran generated files
            in bin_dir/obs or bin_dir/res/obs
            
            -----------
            parameter
            -----------
            self: Station object
            pathfile: path to the file to source
            
            -----------
            return
            -----------
            Fills self of values for time,h,u,v variables
            """
            res = np.loadtxt(pathfile, skiprows = 1, usecols= (0,1,2,3) ) # ignore w column
            self.time  = res[:,0]
            self.h     = res[:,1]
            self.u     = res[:,2]
            self.v     = res[:,3]