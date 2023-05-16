from __future__ import annotations

import dassflow2d         # package
import numpy as np         # numpy
import matplotlib.pyplot as plt 
import os                # for system discussion (get path, etc...)
import re                 # for grep
import h5py                # for save file
import pandas as pd        # for read tables
from importlib.metadata import version  # get package version
import pyvista as pv    # plot grids



from dassflow2d.core.config import Config
from dassflow2d.core.meshing import Meshing
from dassflow2d.core.boundary import Boundary
from dassflow2d.core.output import Output
from dassflow2d.core.param import Param
from dassflow2d.core.min import Min
from dassflow2d.core.obs import Obs




#####################
# Main class
#######################

__all__ = ["dassflowmodel"]


class dassflowmodel(object):
    
    """
    Main object of ``DassFlow-2d``.

    
    -------------

    Parameters
    ----------
    dassflow_dir: str,   default None
                    Path to the root of your dassflow2d-wrap directory 
                    
    bin_dir: str, default None
            Path to your bin_directory, the configuration file input.txt must be at the root of bin_dir.

    run_type: str, default "direct"
             The type of simulation you want to run, can be either
            - 'direct' : if you want to run the forward code
            - 'min'    : if you want to perform an inference
            (Variational Data Assimilation using the adjoint code)
            - 'grad'   : if you want to perform a gradient based optimization
    clean: if True,  keeps only input files in bin directoy (apply all command "make cleanxxx"  existing)

    **kwarg: not used, bonus argument : pass any arguments to a unique function of your choice within this function

    """


#=================================================================#
# INITIALISATION                                                  #
#=================================================================#
#dassflow_dir
#       dassflow_dir: str, 
#             absolute path to the root of your dassflow2d-wrap directory   
    def __init__(self, bin_dir,hdf5_path, run_type="direct",clean = False, **kwargs): 

        """
        This initialize **fortran** instance and some **python** objects


       **Parameters**
                  
       bin_dir: str, default None
             absolute path to your bin_directory, the configuration file input.txt must be at the root of bin_dir             
       hdf5_path: str, default None
             absolute path to your hdf5 file (where you will save results)       
       run_type: str, default "direct"
             The type of simulation you want to run, can be either
             - 'direct' : if you want to run the forward code
             - 'min'    : if you want to perform an inference
             (Variational Data Assimilation using the adjoint code)
             - 'grad'   : if you want to perform a gradient based optimization
             
             
       clean: boolean,
             if True,  keeps only input files in bin directoy (apply all command "make cleanxxx"  existing)
             
             
       **Return**
       
       
       python variable:
          
             - self.bin_dir
             - self.dassflow_dir
             - self.run_type
             - self.config: dict configuration data stored in a dictionary
             
       fortran:
                
             - "self.kernel.mpi": mpi setup information
             - "config": global variables in m_common and m_model
             - "kernel": "model" in smash --> this correspond to the elements that are directly and continuously interfaced to fortran kernel

       
        """
        
        # ------------ FIRST mpi treatment------------------
        # this is necessary for both mpi and non mpi cases in order to initialise some variables
        dassflow2d.wrapping.m_mpi.init_mpi()
        rank = dassflow2d.wrapping.m_mpi.get_proc() # get the rank and number of processors
        nproc = dassflow2d.wrapping.m_mpi.get_np()
        self.mpi = [rank, nproc]
        
        
        
        
        if run_type not in ["direct", "min", "grad"]:
            raise ValueError("argument run_type incorect, must be either 'min', 'direct' or 'grad' ")

        # ------------ store inputs as data within DassFlowModel------------------
        # store main path
        self.bin_dir        = bin_dir
        self.hdf5_path      = hdf5_path
        # run type
        self.run_type     = run_type
        
        self.clean = clean
        # kwar
        self.kwargs        = kwargs

        # ------------ Define Default values------------------
        self.config = Config()

        # ------------ Clean files ------------------
        print("clean=", clean)
        if clean:
                os.chdir(f"{self.bin_dir}")
                if os.path.exists(f"{self.bin_dir}/res/"):
                    os.system(f'rm -r  {self.bin_dir}/res/*')               # removes all in bin_dir/res directory
                if os.path.exists(f"{self.bin_dir}/min/"):
                    os.system(f'rm -r  {self.bin_dir}/min/*')               # removes all in bin_dir/res directory
                if os.path.exists(f"{self.bin_dir}/msh/"):
                    os.system(f'rm -r  {self.bin_dir}/msh/')               # removes all in bin_dir/msh directory
                if os.path.exists(f"{self.bin_dir}/restart.bin"):
                    os.system(f"rm {self.bin_dir}/restart.bin")   # removes all in bin_dir/msh directory


        # ------------ input file treatment------------------
        os.chdir(bin_dir)                     # open bin directory
        dassflow2d.wrapping.read_input(f"{bin_dir}/input.txt") # read input    and update in fortrnan kernel
        self.config.get()                     # update configuration in python from fortran kerne value


        # ------------ Model = mdl in fortran, kernel in python ------------------
        self.kernel = dassflow2d.wrapping.call_model.Model()
        

    def init_mesh(self):
        """
        Initialize mesh instance, coth in fortran kernel and python.
        Also initialise fortran kernel solver
        
        parameters:
            self: dassflow2d object (with initialised self.bin_dir, dans self.config.mesh_name)
            
        return:
            self.meshing: dassflow2d.core.Meshing object initialised with fortran kernel values
        """
        os.chdir(self.bin_dir)
        print("call dassflow2d.wrapping.m_mesh.msh()")
        self.kernel.mesh = dassflow2d.wrapping.m_mesh.msh()
        print("call         dassflow2d.wrapping.call_model.init_solver(self.kernel)")
        dassflow2d.wrapping.call_model.init_solver(self.kernel)
        print("call   Meshing(mesh_fortran = self.kernel.mesh)")
        self.meshing = Meshing(mesh_fortran = self.kernel.mesh)
        
    def init_dof(self, **kwargs):
        """
        Initialize dof (h,u,v) fortran kernel instance.
        
        parameters:
            self: dassflow2d object (with initialised self.bin_dir, dans self.config.mesh_name)
            require **self.bin_dir**
            
        return:
            self.meshing: dassflow2d.core.Meshing object initialised with fortran kernel values
        """
        os.chdir(self.bin_dir)
        self.kernel.dof = dassflow2d.wrapping.m_model.unk(self.meshing.mesh_fortran)
        self.kernel.dof0 = dassflow2d.wrapping.m_model.unk(self.meshing.mesh_fortran)        
        
    def init_fortran(self, **kwargs):       
        """
        Perform all last fortran kernel initialisatins
        
        it accounts for:
            - initial conditions (from ic.bin file / restart.bin / dof.init)
            - Boundary condition table of values
            - some rain treatment
            - observation station treatment if self.config.use_obs = 1
        
        parameters:
            self: dassflow2d object (with initialised self.bin_dir, dans self.config.mesh_name)
            requires dassflowmodel.__init__() dassflowmodel.init_mesh(), dassflowmodel.init_dof() already performed
            
        return:
            self.boundary: dassflow2d.core.Boundary object initialised with fortran kernel values
            and initialiased fortran kernel   valeus
        """
        dassflow2d.wrapping.call_model.init_fortran(self.kernel)
        self.boundary  =  Boundary(mesh_pyvista = self.meshing.mesh_pyvista, mesh_fortran=self.meshing.mesh_fortran)
        self.boundary = self.boundary.get_mesh_corresp()
        
    def init_fortran_mpi(self):
        """
        light init fortran so that mpi is ok
        (no python interaction)
        """
        dassflow2d.wrapping.call_model.init_fortran(self.kernel)
    

    def init_param(self, **kwargs):       
        """
        Get parameters values from fortran kernel and put them in dassflowmodel class
        
        it accounts for:
            - friction.manning
          
        return:
            self.param: dassflow2d.core.param object
        """
        self.param = Param()
        self.param.get()
        
    def init_all(self):
        """
        Perform in correct order all fortran initialisation
        build ready to run model
        """
        self.init_mesh()
        # initialise dof structure
        self.init_dof() 
        # initialise remaining structures
        self.init_fortran()
        # source param kernel data to python
        self.init_param()
#=================================================================#
# RUN                                                              #
#=================================================================#
    def run(self, **kwargs):
        """
        Run Dassflow model
        
        

        when all values are correctly set up perform wished simumation
        

        **Parameters**
             
        self: dassflowmodel
        
        

        **Requirement**
        
        kernel: dassflow.kernel initialised object
        
        run_type: str, default "direct"
             - 'direct' : if you want to run the forward code
             - 'min'    : if you want to perform an inference
             (Variational Data Assimilation using the adjoint code)
             - 'grad'   : if you want to perform a gradient based optimization
        """
        
        dof_init = dassflow2d.wrapping.call_model.dof_copy(self.kernel.dof)
        dassflow2d.wrapping.call_model.run(self = self.kernel, arg = self.run_type)
        
        if self.run_type =="min":
            # source minimization result
            self.min = Min(bin_dir = self.bin_dir)
            self.min.source_all()
            
            # if hydrograph infered, replace the hydrograph.txt file, so that we can perform direct simulation with infered parameter
            if dassflow2d.wrapping.m_model.get_c_hydrograph() ==1:
                os.chdir(self.bin_dir)
                os.system('cp ./hydrograph.txt ./hydrograph_prior.txt')
                key_last_ite = len(self.min.param["hydrograph"])
                with open('./hydrograph.txt', 'w') as f:
                    f.write("#comment\n")
                    f.write("#comment\n")
                    f.write("#comment\n")
                    f.write(f"{len(self.min.param['hydrograph'][key_last_ite])}\n")
                    for key_id_hydrograph, hyd in  self.min.param["hydrograph"][key_last_ite].items():
                        nb_timstep = np.shape(hyd)[0]
                        f.write("#comment\n")
                        f.write("#comment\n")
                        f.write("#comment\n")
                        f.write(f"{str(nb_timstep)}\n")
                        for id_timestep in range(nb_timstep):
                            f.write(f"{hyd[id_timestep, 0]} {hyd[id_timestep, 1]}\n")
                        
                # hydrograph.txt file has been rewriten, lets clean the kernel, reinitilize the model, and perform direct simulation with appropriate parameter                
#                dassflow2d.wrapping.call_model.clean_model(self.kernel)
#
#                if os.path.exists(f"{self.bin_dir}/res/obs"):
#                    list_paths = os.listdir(f'{self.bin_dir}/res/obs')
#                    for namefile in list_paths:
#                    	a_file = open(f"{self.bin_dir}/res/obs/{namefile}", "r")
#                    	lines = a_file.readlines()
#                    	a_file.close()
#                    	new_file = open(f"{self.bin_dir}/res/obs/{namefile}", "w")
#                    	new_file.write(lines[0])
#                    	new_file.close()
#
#                self.init_all()                
#                self.kernel.dof = dassflow2d.wrapping.call_model.dof_copy(dof_init)
#                self.kernel.dof0 = dassflow2d.wrapping.call_model.dof_copy(dof_init)
#                    
#                dassflow2d.wrapping.call_model.run(self = self.kernel, arg = "direct")
                
                    
                    
        if self.run_type =="direct":
        	self.outputs = Output(bin_dir = self.bin_dir, 
                    ts = self.config["ts"], 
                    boundary_metadata = self.boundary.get_metadata() )

#=================================================================#
# Save result                                                     #
#=================================================================#
    def save_all(self):
        """
        save all elements
        """
        self.outputs.save(hdf5_path=self.hdf5_path) 
        self.config.save(hdf5_path=self.hdf5_path)
        self.meshing.save(hdf5_path=self.hdf5_path)
        self.boundary.save(hdf5_path=self.hdf5_path) 
        self.param.save(hdf5_path=self.hdf5_path) 




#================================================================================================================================#
#  PLOT result, using hdf5 file and grid object
#================================================================================================================================#




    # self : class model
    # my_hdf5_file :  <HDF5 file "simu.hdf5" (mode r)> built by build_hdf5
    # what : one of the spatial variables : "h", "u", "v", "bathy", etc..
    # what can be "zs"
    # when : id (integer) of timestep, can be also "initial", or "final"
    # title plot : the title of the plot

    def plot_var(self, path_hdf5_file = None, what = "bathy", when = 0, title_plot = " title ", save_plot= False, filename="tmp"):
        """
        Plot wished variable at wished time, sourcing hdf5 file

        Parameters
        ----------

        path_hdf5_file: str, default self.bin_dir + "/res/" + "simu.hdf5"
                  Path to the HDF5 file generated by DassFlowModel.save_res
                  The variables to plot are sourced from this file

        what: str, default="bathy"
                 which spatial variables to plot? must be one of the following :
                - h
                - u
                - v
                - bathy
                - manning_alpha
                - zs
                - vel : for velocity (sqrt(u2+v2) )
                **note that zs and vel are calculated within the function**

        when: int, default=0
                the number of the time step you want to plot
                can also be a string: **"initial"** for the initial condition (the 0th timestep) and **"final"**  for the last timestep

        title_plot: str, default = "title",
                    title of the plot

        save_plot: logical, default = False
                Do you wish to save the plot or not
        filename: str, default = "tmp"
                Path to fsave the graphic file to. (see Plotter.save_graphic)

        Examples
        --------
         no example

        See Also
        --------
        plot_grid: similar method but where you can define your own variable to plot

        """
        my_mesh = self.grid
        if path_hdf5_file is None :
            path_hdf5_file = self.bin_dir + "/res/" + "simu.hdf5"

        my_hdf5_file =  h5py.File(path_hdf5_file, "r")

        if(when == "initial") :
            when =  0
        if(when == "final") :
            when =  my_hdf5_file["spatial_var/" + str("h")].shape[1]-1

        time_corresp = my_hdf5_file["time"][:] # not used yet but could be used for title


        if what == "zs":
                my_scalar = my_hdf5_file["spatial_var/" + str("bathy")][:,when]  + my_hdf5_file["spatial_var/" + str("h")][:,when]
        elif what == "vel":
            my_scalar = np.sqrt(my_hdf5_file["spatial_var/" + str("u")][:,when] ** 2  + my_hdf5_file["spatial_var/" + str("v")][:,when] ** 2)
        else:
            my_scalar = my_hdf5_file["spatial_var/" + str(what)][:,when]


        title_scale_bar = what + "  [USI]"

        sargs = dict(vertical = True, title = title_scale_bar, position_y=0.3, position_x = 0.8  )
        # Plotter instance
        plotter = pv.Plotter()
        # mesh plot
        plotter.add_mesh(my_mesh, scalars = my_scalar, scalar_bar_args=sargs,
                 show_edges = False)
        plotter.view_xy()
        plotter.show_bounds()
        plotter.add_title(title_plot, font_size=18, color=None, font=None, shadow=False)

        if(save_plot == True):
            plotter.save_graphic(filename=f"{filename}.svg", title = title_plot)

        plotter.show()
        return(plotter)
        
        
        
        
        
        
