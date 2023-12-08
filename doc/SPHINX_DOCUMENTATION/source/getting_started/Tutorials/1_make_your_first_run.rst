.. _1_make_your_first_run:

===================================
Forward run (lake at rest)
===================================

This tutorial details how to perform a direct/forward run with `dassflow2d`. The 2D SW equations are here solved on an academic case consisting in a lake at rest. **TO_DO ref to eqs & solver used **


The goal of this case is to check whether the equilibrium state is preserved, that is preserving initial states.
More precisely the “Lake at rest” is performed to validate the well-balanced property of the numerical schemes
used in DassFlow. The equilibrium challenged here is the equilibrium between fluxes and gravity source term
:math:`S_g (U)` for a perturbed topography.


.. image:: ./images/mesh_clim.png
  :width: 600


The perturbed topography looks like:

.. image:: ./images/bathy_0.png
      :width: 600

----------------------------------------------
Run "lake at rest" case with command lines
----------------------------------------------

++++++++++++++++++++++++++++++
Set up your environment
++++++++++++++++++++++++++++++

Prepare the simulation files corresponding to "lake at rest" case with:

.. hint::

     File copy and compilation has just been performed in :ref:`getting_started` introduction, during :ref:`Installation` process. Thus, you don't need to compile the code again.

Open a terminal in (`/dassflow2d-wrap/code`) and run the following commands:

.. code-block:: bash


    # delete all files in your simulation directory before starting
    rm -r ./bin_A/*
    # Copy recursively the files provided in DassFlow case repository into your own simulation directory **code/bin_A/**.
    cp -r ../cases/tuto_case/1_lake-at-rest/bin_A/* ./bin_A
    

Note that this test case "lake at rest" can be viewed as a stability test case where the well-balancedness of the numerical scheme (see. `Math_num_doc`) is tested in terms of equilibrium preservation. It is a simple test case for which, all the boundaries are set as walls and no inflow or outflow occurs.

++++++++++++++++++++++++++++++++
Launch your first run using make
++++++++++++++++++++++++++++++++

Write the following command in your terminal:

.. .. code-block:: bash

..     make rundirect

This executes some commands coded in the Makefile and printed in the terminal at the beginning of the execution:

.. image:: ./images/tuto1_make_rundirect_initialization.png
  :width: 600


.. image:: ./images/tuto1_make_rundirect_initialization.png
  :width: 600


You should see in the terminal the successive temporal iterations of the numerical resolution of the 2D SW model.

++++++++++++++++++++++++++++++++
Investigate results
++++++++++++++++++++++++++++++++

Have a look to directory ./dassflow2d-wrap/code/bin_A/res/ and investiguate the result files produced (`result_initial.dat` and `result_final.dat`).
These dat files are in gnuplot format (cf. http://www.gnuplot.info/) ; another output format as well as plot tools are available in DassFlow and presented after.

.. NB. ON fera du vtk à partir du hdf5 car "gratos" ; à documenter plus loin

-----------------------------------
Run "lake at rest" case with Python
-----------------------------------


.. jupyter-execute::
     
     
     #####################################################################
     #####################################################################
     # PERFORM A DIRECT SIMULATION WITH  DASSFLOW2D
     # LAKE AT REST
     #
     # Introduction to basic commands of run and visualisation of results
     #####################################################################
     #####################################################################

     
     #=======================================================#
     # Source librairies
     #=======================================================#
     
     import dassflow2d as df2d
     import os
     import numpy as np
     import matplotlib
     import matplotlib.pyplot as plt
    
     #=======================================================#
     # copy of case files
     #=======================================================#
     
     os.chdir('../../')
     dassflow_dir = os.getcwd() # DassFlow directory (you can also impose your absolute path)
     
     #dassflow_dir="/home/pagarambois/Documents/Distant/dassflow2d"
     os.chdir(dassflow_dir)
     print("DassFlow directory is: ", dassflow_dir)
     
     # Define directory where case is run 
     # (its name 'bin_A' is imposed in  {dassflow_dir}/code/makefile.inc : CASEDIR='bin_A')
     run_dir = f"{dassflow_dir}/code/bin_A/" 
     
     # Define directory containing case data
     case_data_dir = f"{dassflow_dir}/cases/tuto_case/1_lake-at-rest/bin_A/"
     
     # Clean run directory
     os.system(f"rm -r {run_dir}*") 
     
     # Copy case data to runing directory
     os.system(f"cp -r {case_data_dir}* {run_dir}") # Copy of case files from existing case to bin_A
     
     # Move to code directory and clean bin directory
     os.chdir( f"{dassflow_dir}/code/")
     os.system("make cleanres cleanmin") # Clean forward run and minimization results 
     
     
.. jupyter-execute::

     #=======================================================#
     # Initialisation
     #=======================================================#

     # input file reading (simulation settings)
     df2d.wrapping.read_input(f"{run_dir}/input.txt")

     # Creation of dassflowmodel object using case data: 
     df2d.wrapping.m_mpi.init_mpi() #set the number of processes to 1
     my_model = df2d.dassflowmodel(bin_dir =  run_dir, hdf5_path = f"{dassflow_dir}/code/bin_A/res/simu.hdf5" , run_type = "direct", clean = True, custom_config = None)

     my_model.config.get()

     # Initializion of the Fortran kernel (dassflow Python library is obtained by wrapping Fortran source code)
     #initialise all fortran kernel values and source them into dassflowmodel object

     my_model.init_all()

     my_model.kernel.dof0.h[:] = 1
     my_model.kernel.dof.h[:] = 1


     
.. jupyter-execute::

     #=======================================================#
     # Run Fortran kernel and save results
     #=======================================================#
     
     my_model.run()
     my_model.save_all()

The numerical resolution is performed with variable time steps (depending on the CFL condition) and outputs are written at each writting timestep (imposed by the "dtw" parameter).
     
.. jupyter-execute::

     #=======================================================#
     # Vizualize parameters and results
     #=======================================================#

First, you can have a look at the bathymetry, friction and initial conditions (of water heigth and free surface height).     

.. jupyter-execute::
     
     # Plot of the 2D bathymetry (input parameter of the 2D shallow water model) with package plot function
     
     #plotter = my_model.plot_var(my_model.meshing.mesh_pyvista,
     #                                             what = "bathy", 
     #                                             title_plot = "Bathymetry elevation")# for a local run remove notebook option or set notebook=False 
                                        
     #plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed
     
     
.. jupyter-execute::

     # Plot the friction parameter field
     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                                  what = "manning_alpha", 
                                                  title_scale_bar ="n [m-1/3.s] ", 
                                                  title_plot = "Friction parameter (Manning coefficient)", 
                                                  notebook = True )# for a local run remove notebook option or set notebook=False 
                                        
     plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed
     
The friction is uniform as defined in this case setup. 

.. jupyter-execute::

     # Plot intial flow conditions 
        
     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                                  what = "h", 
                                                  when = 0,
                                                  title_scale_bar ="h [m] ", 
                                                  title_plot = "Initial water depth", 
                                                  notebook=True) # for a local run remove notebook option or set notebook=False 
                                        
     plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed

     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                                  what = "zs", 
                                                  when = 0,
                                                  title_scale_bar ="zs [m] ", 
                                                  title_plot = "Initial water surface elevation", 
                                                  notebook=True) # for a local run remove notebook option or set notebook=False 
                                        
     plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed
     
     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                                  what = "u", 
                                                  when = 0,
                                                  title_scale_bar ="u [m/s] ", 
                                                  title_plot = "Initial velocity u along x", 
                                                  notebook=True) # for a local run remove notebook option or set notebook=False 
                                        
     plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed

     
     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                                  what = "v", 
                                                  when = 0,
                                                  title_scale_bar ="v [m/s] ", 
                                                  title_plot = "Initial velocity v along y", 
                                                  notebook=True) # for a local run remove notebook option or set notebook=False 
                                        
     plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed


.. jupyter-execute::

     # Plot flow depth at a given time
     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                                  what = "h", 
                                                  when = 3,
                                                  title_scale_bar ="h [m] ", 
                                                  title_plot = f"Water depth at time = {my_model.outputs.result.all_time[3]}  s ", 
                                                  notebook=True) # for a local run remove notebook option or set notebook=False 
                                        
     plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed
     
     # Simulation time steps at which variables have been written 
     print(my_model.outputs.result.all_time)
     
     print("previous plot for t = ", my_model.outputs.result.all_time[3])

Let us study the flow state at the end of the simulation 

.. jupyter-execute::
     
     # Compute velocity magnitude
     
     u = my_model.outputs.result.u
     v = my_model.outputs.result.v 
     norm_vel = np.sqrt(u**2+v**2)
     
     # Print the shape of the output velocity fields
     print("shape of velocity array is : ", np.shape(norm_vel))
     print("Maximum velocity magnitude at ecah time step is : ", np.amax(norm_vel,axis=0))
     
     
     # Plot velocity magnitude at final time step
     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista, 
                                                  my_scalar = norm_vel[:,-1],
                                                  title_scale_bar ="norm(u,v) [m/s] ",
                                                  title_plot = f"Velocity magnitude at final time",
                                                  notebook=True) # for a local run remove notebook option or set notebook=False 
                                        
     plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed      
     
     # Plot water surface elevation at final time step

     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                                  what = "zs",
                                                  when = -1,
                                                  title_scale_bar ="zs [m] ", 
                                                  title_plot = f"Water surface elevation at final time", 
                                                  notebook=True) # for a local run remove notebook option or set notebook=False 
                                        
     plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed    
     
.. note::

   We can see that at the end of the simulation, the velocity magnitude can be considered as null and the water surface elevation as almost constant, hence the lake is at rest.
   
   This result, obtained on a case with non trivial bathymetry and initial state with wall lateral boundary conditions, validates the capability of the numerical scheme in preserving equilibrium.
   
You can get information on the configuration of dassflowmodel object with:

.. jupyter-execute::

   print(my_model.config)
   
   print("The numerical scheme that has been used to solve the 2D shallow water equations is:")
   print("Temporal scheme is: ", my_model.config["temp_scheme"])
   print("Spatial scheme is: ", my_model.config["spatial_scheme"])
   
.. jupyter-execute::
   
   #clean model 
   df2d.wrapping.call_model.clean_model(my_model.kernel)


.. hint::

    The above Python script is available here: :download:`1_make_your_first_run.py <../../../build/jupyter_execute/getting_started/Tutorials/1_make_your_first_run.py>`

    A Jupyter Notebook version is available here: :download:`1_make_your_first_run.ipynb <../../../build/jupyter_execute/getting_started/Tutorials/1_make_your_first_run.ipynb>`


.. warning::

       Note that the location of the dassflow directory has to be defined by setting appropriate value to **dassflow_dir**  at the begining of the above script. A relative path has been used here but you can also impose your own absolute path to run a script in terminal or from Python IDE from other directories.

.. note::
   
   The above script can be used to run any other case by simply providing case data with necessary inputs for DassFlow. 
   In next tutorial we will use this script to run another case.

