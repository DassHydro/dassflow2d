import dassflow2d
import os
import re
import pandas as pd
import numpy as np
import h5py
from importlib.metadata import version  # get package version














class Output(object):
    """
    define output from fortran kernel simulation results


    ===============
    CONTAINS
    ===============
    self.results: df2d.core.output.Result object: Result from simulation at each writting timestep
    self.post:    df2d.core.output.Post   object: Post result (flux at boundary edges) at each post writting timestep (to check dassflow_model.config["dtp"] influence)

    ==============
    methods
    ==============

    __init__ : fills values from text files, path are automaticaly generated from bin_dir root and boundary_metadata objects
    save :    save in hdf5 file the data stored in Output class
    source: source textfiles from bin_dir and fills values of self.result and self.post

    """
    def __init__(self, bin_dir,ts, boundary_metadata,  verbose=2) :
        """
        Initialize Output class:
            - initialize Result class
            - initialize Post class
        """
        self.bin_dir = bin_dir
        self.res_dir  = self.bin_dir + "/res/"

        self.result = Result(bin_dir = self.bin_dir, ts =  ts)
        self.post = Post(bin_dir = self.bin_dir, boundary_metadata =  boundary_metadata)

        if verbose >= 2:
        	print(">>> Outputs initialized ")
        	print(f" - Simulation Result files were loaded from : {self.res_dir} ")



    def save(self, hdf5_path):
        """
        Save in hdf5 format in hdf5_file["output"] group, both result and post create_dataset
        """

        if hdf5_path is not None:
            print("")
        else:
            hdf5_path = self.hdf5_path

        self.result.save(hdf5_path = hdf5_path)
        self.post.save(hdf5_path = hdf5_path)

    def source(self, boundary_metadata):
        """
        Source from text file both Result and Post objects
        """
        self.result.source()
        self.post.source(boundary_metadata = boundary_metadata)




class Result(object):
    """
    Result class

    contains result data produced as result_xxx.dat files in "bin_dir/res/" directory

    ===============
    CONTAINS
    ===============
    self.bin_dir: str : path to bin directory
    self.res_dir: str : path to bin directory

    self.all_time:  list of num : time values
    self.all_files: list of str : absolute path to each result file sourced (ordered)

    self.h: np.ndarray(nb_cell, nb_writting_timestep): h value on each cell at each timestep [m]
    self.u: np.ndarray(nb_cell, nb_writting_timestep): u-velocity value on each cell at each timestep [m/s]
    self.v: np.ndarray(nb_cell, nb_writting_timestep): v-velocity value on each cell at each timestep [m/s]

    self.id_cell_fortran: np.array(nb_cell): FORTRAN index of correponding cell (same ordering for all objects)
    self.x: np.array(nb_cell): x coordinate of cell
    self.y: np.array(nb_cell): y coordinate of cell

    self.bathy: np.array(nb_cell): bathymemtry value
    self.zs: np.array(nb_cell): water surface height
    self.manning_alpha: np.array(nb_cell): manning coefficient parameter



    ===============
    Methods
    ===============
    __init__ : initialize object and
    list_dat_files(): that list .dat result files in result directory (bin_a/res)
    source() :  that source the necessary files (listed by list dat file)
    save():     that save to hdf5 file the data
    """


    def __list_dat_files(self):
        """
        List files availaible in self.res_dir,

        res_dir: dir where to scrap from

        RETURN: list of files (absolute path) of result to source (that are the files that respect the ".dat" result file format).
        Files are ordered in ascending time
        """
        def _my_grep(test_list, pattern):
            """
            classicall grep function
            """
            return [i for i in test_list if re.search(pattern, i)]

        onlyfiles = [f for f in os.listdir(self.res_dir) if os.path.isfile(os.path.join(self.res_dir, f))]
        all_files = _my_grep(test_list= onlyfiles, pattern = ".dat")
        all_times = np.empty(shape = len(all_files))

        # replace initial and final times by appropriate numeric values
        for id_file in range(len(all_files)):
            tmp = all_files[id_file].split("_")[1]
            tmp = re.sub(".dat", "", tmp)
            if tmp == "initial":
                tmp = 0
            if tmp == "final":
                tmp = self.ts
            all_times[id_file] = tmp

        # reorder in increasing time
        id_reorder = np.argsort(all_times)
        all_times = all_times[id_reorder]
        all_files = [self.res_dir + "/" + all_files[i] for i in id_reorder]

        return(all_times, all_files)

    def __source_for_initialisation(self, all_files):
        """
        source all result file, and return np.ndarray(nb_cell, nb_writted_variables)
        variable are ordered like: id_cell_fortran x y bathy h zs Manning u v

        self: Result object
        all_files: list produced by list_result_files methods

        """

        for id_timestep in range(len(all_files)):
            res = np.loadtxt((all_files[id_timestep] ))
            # order of columns in textfile:
            # 0:i 1:x 2:y 3:bathy 4:h 5:zs 6:Manning 7:u 8:v
            self.h[:,id_timestep] = res[:,4]
            self.zs[:,id_timestep] = res[:,5]
            self.u[:,id_timestep] = res[:,7]
            self.v[:,id_timestep] = res[:,8]
        # use last file to fill variables that do not vary in time
        self.id_cell_fortran = res[:,0]
        self.x = res[:,1]
        self.y = res[:,2]
        self.bathy = res[:,3]
        self.manning_alpha = res[:,6]


    def source(self):
        # list result files (.dat) and corresponding times
        self.all_time, self.all_files  = self.__list_dat_files()
        # fill them
        self.__source_for_initialisation(all_files =  self.all_files)

    def __init__(self, bin_dir, ts =  9999999999999999):
        """
        Initialize Result structure, fills all arrays self.(h,u,v, id_cell,x,y,bathy,zs,manning_alpha),
        call source() method
        """

        self.bin_dir = bin_dir
        self.res_dir  = self.bin_dir + "/res/"
        self.ts = ts

        if self.ts == 9999999999999999: # some check
            print("warning, default ts is set for output results, this means ts = 9999999999999999")


        # list result files (.dat) and corresponding times
        self.all_time, self.all_files  = self.__list_dat_files()
        # temporary variables to determine shapes
        nb_cell =  np.loadtxt(self.all_files[0]).shape[0]
        nb_writting_timestep = len(self.all_time)

        # Generate empty variables
        # dof
        self.h = np.ndarray( shape = (nb_cell, nb_writting_timestep))
        self.u = np.ndarray( shape = (nb_cell, nb_writting_timestep))
        self.v = np.ndarray( shape = (nb_cell, nb_writting_timestep))
        # geometry
        self.id_cell_fortran = np.ndarray( shape = (nb_cell))
        self.x = np.ndarray( shape = (nb_cell))
        self.y = np.ndarray( shape = (nb_cell))
        # parameters
        self.bathy =  np.zeros(nb_cell)
        self.manning_alpha =  np.zeros(nb_cell)
        # by-products:
        self.zs = np.ndarray( shape = (nb_cell, nb_writting_timestep))

        # fill them

        self.__source_for_initialisation(all_files =  self.all_files)

    def save(self, hdf5_path ):
        """
        Save  Result class to pre-existing hdf5 file.


        -----------
        parameter
        -----------

        self:    Result class
        hdf5_path : path to hdf5 file to save

        -----------
        return
        -----------
        hdf5 file filled with data stored in Result class, this account for: my_hdf5["output"]["result"]["i"]
        with i =h,u,v,zs,x,y, id_cell_fortran,manning_aplha, bathy, all_times
        """
        # Open hdf5 file to append
        f = h5py.File(hdf5_path, "a")

        # Identify dimentions
        nb_cell     =  self.h.shape[0]
        nb_timestep =  self.h.shape[1]


        # create groups to be filled
        output = f.create_group("output")
        result = output.create_group("result")

        # Fill HDF5 file with self.xxxx values

        # dof
        h_new = result.create_dataset("h", (nb_cell,nb_timestep), 'f')
        result["h"][:,:]= self.h
        u_new = result.create_dataset("u", (nb_cell,nb_timestep), 'f')
        result["u"][:,:]= self.u
        v_new = result.create_dataset("v", (nb_cell,nb_timestep), 'f')
        result["v"][:,:] = self.v

        # geometry
        id_cell_fortran_new = result.create_dataset("id_cell_fortran", (nb_cell), 'f')
        result["id_cell_fortran"][:] = self.id_cell_fortran
        x_new = result.create_dataset("x", (nb_cell), 'f')
        result["x"][:] = self.x
        y_new = result.create_dataset("y", (nb_cell), 'f')
        result["y"][:] = self.y


        # parameters
        bathy_new = result.create_dataset("bathy", (nb_cell), 'f')
        result["bathy"][:] = self.bathy
        manning_alpha_new = result.create_dataset("manning_alpha", (nb_cell), 'f')
        result["manning_alpha"][:] = self.manning_alpha

        #by-products --> ignored (surface height = bathy + h)

        # time/file values
        time_new = result.create_dataset("all_time",  (nb_timestep), 'f')
        result["all_time"][:] = self.all_time

        f.close()

        print(f"Results.save(): File {hdf5_path} updated with Result class data")

class Post(object):
        """
        Post Class:

        contains result data produced  in "bin_dir/res/post" directory
        they mainly contains  data about flows at the edges of interest


        inflow.metadata      : tofiil
        inflow.sum_mass_flux :
        inflow.sum_q         :

        outflow.metadata      :
        outflow.sum_mass_flux :
        outflow.sum_q         :

        TODO: internals

        """

        def __count_nb_bc(self, boundary_metadata):
            """
            Use boundary_metadata dictionary to identify the number of inflow, the number of outflow and the number total of boundaries


            return: tuple (nb_q_inflow, nb_q_outflow, nb_total_bc)
            nb_q_inflow: number of inflow
            nb_q_inflow: number of outlfow
            nb_total_bc: total number of bc, accounting for internals also (not walls)
            """

            nb_q_inflow = 0
            nb_q_outflow = 0
            nb_internal_2d = 0

            # count nb_inflow
            if "discharg1" in boundary_metadata.keys():
                nb_q_inflow = nb_q_inflow + boundary_metadata["discharg1"]
            if "discharg2" in boundary_metadata.keys():
                nb_q_inflow = nb_q_inflow + boundary_metadata["discharg2"]

            # count nb_outflow
            if "transm" in boundary_metadata.keys():
                nb_q_outflow = boundary_metadata["transm"]
            elif "ratcurve" in boundary_metadata.keys():
                nb_q_outflow = boundary_metadata["ratcurve"]

           # counter total number of bc (consider internals), this is used to access the correct outflow id file for the file sum_mass_flux_outflow_ID_MAX_BC
            nb_total_bc = 0
            for my_key in boundary_metadata.keys():
                if my_key in ["discharg1", "discharg2", "transm", "ratcurve", "internal_1D", "internal_2D" ]:
                    nb_total_bc = nb_total_bc +boundary_metadata[my_key]

            return(nb_q_inflow, nb_q_outflow, nb_total_bc)

        def __source_inflow_values(self):

            template_file =  np.loadtxt(self.filenames["sum_mass_flux_inflow"][0], comments="#")
            # inflow
            self.sum_mass_flux_inflow = np.ndarray(shape = (template_file.shape[0], self.nb_q_inflow))
            self.sum_q_inflow =  np.ndarray(shape = (template_file.shape[0], self.nb_q_inflow))
            for my_id in range(self.nb_q_inflow):
                self.sum_mass_flux_inflow[:, my_id] 	= np.loadtxt(self.filenames["sum_mass_flux_inflow"][my_id], comments="#")[:,1]
                self.sum_q_inflow[:, my_id]		    = np.loadtxt(self.filenames["sum_q_inflow"][my_id], comments="#")[:,1]


        def __source_outflow_values(self):

            template_file =  np.loadtxt(self.filenames["sum_mass_flux_outflow"][0], comments="#")
            # outflow
            self.sum_mass_flux_outflow =  np.ndarray(shape = (template_file.shape[0], self.nb_q_outflow))
            self.sum_q_outflow  =  np.ndarray(shape = (template_file.shape[0], self.nb_q_outflow))
            for my_id in range(self.nb_q_outflow):
                self.sum_mass_flux_outflow[:, my_id] 	= np.loadtxt(self.filenames["sum_mass_flux_outflow"][my_id], comments="#")[:,1]
                self.sum_q_outflow[:, my_id]   		= np.loadtxt(self.filenames["sum_q_outflow"][my_id], comments="#")[:,1]


        def source(self, boundary_metadata):

            self.nb_q_inflow, self.nb_q_outflow, self.nb_total_bc = self.__count_nb_bc(boundary_metadata = boundary_metadata)

            # list all file path
            self.filenames =  dict()
            self.filenames["sum_mass_flux_inflow"] =  [f"{self.post_dir}/sum_mass_flux_inflow_{'{:03}'.format(x+1)}.dat" for x in range(self.nb_q_inflow)]
            self.filenames["sum_q_inflow"] = [f"{self.post_dir}/sum_q_inflow_{'{:03}'.format(x+1)}.dat" for x in range(self.nb_q_inflow)]
            self.filenames["sum_mass_flux_outflow"] = [f"{self.post_dir}/sum_mass_flux_outflow_{'{:03}'.format(x+1)}.dat" for x in range(self.nb_total_bc-self.nb_q_outflow, self.nb_total_bc)]
            self.filenames["sum_q_outflow"] = [f"{self.post_dir}/sum_q_outflow_{'{:03}'.format(x+1)}.dat" for x in range(self.nb_q_outflow) ] # WARNING : TO DEBUG
            self.filenames["time_step"] = f"{self.post_dir}/time_step.dat"
            self.filenames["water_vol"] = f"{self.post_dir}/water_vol.dat"
            self.filenames["water_vol_num_add"] = f"{self.post_dir}/water_vol_num_add.dat"
            
            
            
            if np.any([x=="transm" for x in boundary_metadata.keys()]):
            	self.filenames["sum_q_outflow"] = [f"{self.post_dir}/sum_q_outflow_{'{:03}'.format(x)}.dat" for x in range(self.nb_q_outflow) ]
            	self.filenames["sum_mass_flux_outflow"] = [f"{self.post_dir}/sum_mass_flux_outflow_{'{:03}'.format(self.nb_total_bc)}.dat"]
            

            # source values from files
            self.__source_inflow_values()
            self.__source_outflow_values()
            self.timestep 					= np.loadtxt(self.filenames["time_step"] , comments="#")[:,1]
            self.all_time 					=  np.loadtxt(self.filenames["time_step"], comments="#")[:,0]
            self.water_vol 					= np.loadtxt(    self.filenames["water_vol"] , comments="#")[:,1]
            self.water_vol_num_add 				=  np.loadtxt(self.filenames["water_vol_num_add"], comments="#")[:,1]


        def __init__(self, bin_dir, boundary_metadata  ):

            self.post_dir = f"{bin_dir}/res/post/"
            self.source(boundary_metadata = boundary_metadata)
            

        def save(self, hdf5_path):
            """
            Save  Result class to pre-existing hdf5 file.


            -----------
            parameter
            -----------

            self:    Result class
            hdf5_path : path to hdf5 file to save

            -----------
            return
            -----------
            hdf5 file filled with data stored in Posr class, this account for: my_hdf5["output"]["post"]["i"]
            with i = sum_q_inflow, sum_mass_flux_inflow, ..._outflow, time_step, water_vol, water_vol_num_add
            """
            # Open hdf5 file to append
            f = h5py.File(hdf5_path, "a")

            output = f["output"]
            # create groups to be filled
            post = output.create_group("post")

            # Fill HDF5 file with self.xxxx values

            # inflow
            sum_mass_flux_inflow = post.create_dataset("sum_mass_flux_inflow", (self.sum_mass_flux_inflow.shape[0],self.sum_mass_flux_inflow.shape[1]), 'f')
            post["sum_mass_flux_inflow"][:,:]= self.sum_mass_flux_inflow
            sum_q_inflow = post.create_dataset("sum_q_inflow", (self.sum_q_inflow.shape[0],self.sum_q_inflow.shape[1]), 'f')
            post["sum_q_inflow"][:,:]= self.sum_mass_flux_inflow

            # inflow
            sum_mass_flux_outflow = post.create_dataset("sum_mass_flux_outflow", (self.sum_mass_flux_outflow.shape[0],self.sum_mass_flux_outflow.shape[1]), 'f')
            post["sum_mass_flux_outflow"][:,:]= self.sum_mass_flux_outflow
            sum_q_outflow = post.create_dataset("sum_q_outflow", (self.sum_q_outflow.shape[0],self.sum_q_outflow.shape[1]), 'f')
            post["sum_q_outflow"][:,:]= self.sum_mass_flux_outflow

            # time
            time_step = post.create_dataset("time_step",  (len(self.timestep)), 'f')
            post["time_step"][:] = self.timestep
            all_time = post.create_dataset("all_time",  (len(self.all_time)), 'f')
            post["all_time"][:] = self.all_time

            # water volumes
            water_vol = post.create_dataset("water_vol",  (len( self.water_vol)), 'f')
            post["water_vol"][:] = self.water_vol
            water_vol_num_add = post.create_dataset("water_vol_num_add",  (len(self.water_vol_num_add)), 'f')
            post["water_vol_num_add"][:] = self.water_vol_num_add

            f.close()

            print(f"Post.save(): File {hdf5_path} updated with Result class data")
