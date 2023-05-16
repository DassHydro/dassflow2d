.. _2_make_your_second_run:

===================================
Forward run (simple channel "q_in")
===================================

This tutorial details how to perform a direct/forward run with `dassflow2d` on a simple channel case provided in the package.

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
     os.chdir(dassflow_dir)
     print("DassFlow directory is: ", dassflow_dir)
     
     # Define directory where case is run 
     # (its name 'bin_A' is imposed in  {dassflow_dir}/code/makefile.inc : CASEDIR='bin_A')
     run_dir = f"{dassflow_dir}/code/bin_A/" 
     
     # Define directory containing case data
     case_data_dir = f"{dassflow_dir}/cases/tuto_case/2_qin/bin_A/"
     
     # Clean run directory
     os.system(f"rm -r {run_dir}*") 
     
     # Copy case data to runing directory
     os.system(f"cp -r {case_data_dir}* {run_dir}") # Copy of case files from existing case to bin_A
     os.system(f"cp {run_dir}/hydrograph_target.txt {run_dir}/hydrograph.txt") # Copy a "real-like" hydrograph into hydrographt.txt used as inflow
     
     # Move to 
     os.chdir( f"{dassflow_dir}/code/")
     os.system("make cleanres cleanmin") # Clean forward run and minimization results 
     
     
.. jupyter-execute::

     #=======================================================#
     # initialization
     #=======================================================#
     
     # Creation of dassflowmodel object using case data: 
     my_model = df2d.dassflowmodel(bin_dir =  f"{dassflow_dir}/code/bin_A", hdf5_path = f"{dassflow_dir}/code/bin_A/res/simu.hdf5" , run_type = "direct", clean = True)
     # Initializion of the Fortran kernel (dassflow Python library is obtained by wrapping Fortran source code)
     #initialise all fortran kernel values and source them into dassflowmodel object
     my_model.init_all()
     my_model.kernel.dof.h[:] = my_model.kernel.dof0.h[:]=1
     
First let us have a look to the case mesh and boundary conditions (BCs)
     
.. jupyter-execute:: 
     
     # Plot the mesh and BCs 
     plotter = my_model.boundary.plot(what="meshing", notebook=True) # for a local run remove notebook option or set notebook=False 
     plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed
     
     plotter = my_model.boundary.plot(what="values", notebook=True) # for a local run remove notebook option or set notebook=False 


.. jupyter-execute::

     #=======================================================#
     # Run Fortran kernel
     #=======================================================#
     
     my_model.run()

The numerical resolution is performed with variable time steps (depending on the CFL condition) and outputs are written at each writting timestep (imposed by the "dtw" parameter).
     
.. jupyter-execute::

     #=======================================================#
     # Vizualize parameters and results
     #=======================================================#

First, you can have a look at the bathymetry, friction and initial conditions (of water heigth and free surface height).     

.. jupyter-execute::
     
     # Plot of the 2D bathymetry (input parameter of the 2D shallow water model)
     
     #Plot bathymetry field
     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                                  what = "bathy", 
                                                  title_plot = "Bathymetry elevation",
                                                  notebook = True )# for a local run remove notebook option or set notebook=True 
                                        
     plotter.show(jupyter_backend='trame') # remove jupyter_backend if needed


.. jupyter-execute::

    # Plot a bathymetry longitudinal profile (at y=50m)
    allx =[]
    allz = []
    for i in range(my_model.meshing.mesh_fortran.nc):
            x =my_model.meshing.mesh_fortran.cell[i].grav.x
            y = my_model.meshing.mesh_fortran.cell[i].grav.y
            if(y==50.0):
                allx.append(x)
                allz.append(my_model.outputs.result.bathy[i-1])
                
    plt.plot(allx[1:-1],allz[1:-1])
    plt.xlabel("x [m]")
    plt.ylabel("$Z_b$ [m]")
    plt.title("Bathymetry profile at $y = 50m$")
    plt.show()
    
    # Plot a bathymetry lateral profile (at x=50m)
    ally =[]
    allz = []
    for i in range(my_model.meshing.mesh_fortran.nc):
            x =my_model.meshing.mesh_fortran.cell[i].grav.x
            y = my_model.meshing.mesh_fortran.cell[i].grav.y
            if(x==500.0):
                ally.append(y)
                allz.append(my_model.outputs.result.bathy[i-1])
                
    plt.plot(ally[1:-1],allz[1:-1])
    plt.xlabel("y [m]")
    plt.ylabel("$Z_b$ [m]")
    plt.title("Bathymetry profile at $x = 500m$")
    plt.show()
    
.. Warning::
   
   Remark that the bathymetry is more complex along x than an inclined plane as depicted with the 1D profiles. This longitudinal bathymetry variation was hardly visible on the 2D bathymetry plot only and its colorscale linearly spanned between minimum and maximum bathymetry elevations.
   
.. jupyter-execute::

     # Plot the friction parameter field     
     
     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                                  what = "manning_alpha", 
                                                  title_scale_bar ="n [m-1/3.s] ", 
                                                  title_plot = "Friction parameter (Manning coefficient)", 
                                                  notebook = True )# for a local run remove notebook option or set notebook=False 

     plotter.show(jupyter_backend='trame') # if used in Jupyter notebook, with notebook = True above
     
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

.. Note::

   Remark that this initialization corresponds to a dry channel, and with imposed hydrograph it will produce in a wet/dry front propagation at the begining of the simulation.

.. jupyter-execute::
     
     # Compute velocity magnitude 
     u = my_model.outputs.result.u
     v = my_model.outputs.result.v 
     norm_vel = np.sqrt(u**2+v**2)
     
     #compute local Froude number at each mesh cell center
     h = my_model.outputs.result.h
     g = my_model.config["g"]
     Froude = norm_vel / np.sqrt(g*h)
     
     # Print the shape of the output velocity fields
     print("The shape of the output velocity array is : \n", np.shape(norm_vel))
     print("The first number corresponds to the number of cells, \n n_cells = ", np.shape(norm_vel)[0])
     print("The second number corresponds to the number output time steps \n nt_out = ", np.shape(norm_vel)[1])
    
     print("Maximum lateral velocity v [m/s] in space at each output time step is: \n v_max_x = ", np.amax(v,axis=0))   

.. Note::

     Note that the lateral velocity v along y is close to zero and can therefore be neglected in the following analysis. This is expected since (i) the bathymetry is invariant along y, (ii) the discharge is inflowed upstream without specific lateral velocity profile, (iii) no lateral momentum exchange is modeled.
     
.. jupyter-execute::
         
     #Check the maximum Froude number
     print("Maximum Froude number in space at each output time step is: \n Fr_max_x =", np.amax(Froude,axis=0))
     
.. Note::

     Note that the maximum values of the Froude are smaller than 1 and indicate fluvial flow regimes over the whole spatial domain at output time steps.


   
.. jupyter-execute::
     
     # Plot velocity magnitude at final time step
     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista, 
                                                  my_scalar = norm_vel[:,-1],
                                                  title_scale_bar ="norm(u,v) [m/s] ",
                                                  title_plot = f"Velocity magnitude at final time",
                                                  notebook=True) # for a local run remove notebook option or set notebook=False 
     
     plotter.show(jupyter_backend='trame') # if used in Jupyter notebook, with notebook = True above    
     
     # Plot flow state at a given time
     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista,
                                                  what = "h",
                                                  when = -1,
                                                  title_scale_bar ="h [m] ", 
                                                  title_plot = f"Water depth at time = {my_model.outputs.result.all_time[3]}  s "  , 
                                                  notebook=True) # for a local run remove notebook option or set notebook=False
     
     plotter.show(jupyter_backend='trame') # if used in Jupyter notebook, with notebook = True above
     
     # Plot Froude number at final time step
     plotter = my_model.outputs.result.plot_field(my_mesh = my_model.meshing.mesh_pyvista, 
                                                  my_scalar = Froude[:,-1],
                                                  title_scale_bar ="Froude",
                                                  title_plot = f"Froude number at final time",
                                                  notebook=True) # for a local run remove notebook option or set notebook=False 
                                        
     plotter.show(jupyter_backend='trame') # if used in Jupyter notebook, with notebook = True above    
   


.. hint::

    The above Python script is available here: :download:`2_make_your_second_run.py <../../../build/jupyter_execute/getting_started/Tutorials/_2_make_your_second_run.py>`

    A Jupyter Notebook version is available here: :download:`2_make_your_second_run.ipynb <../../../build/jupyter_execute/getting_started/Tutorials/2_make_your_second_run.ipynb>`

