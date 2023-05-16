.. _2_make_your_first_4Dvar:

===================================
Inference run
===================================

This tutorial details how to perform an inverse run with `dassflow2d`. The full process is detailed for building and performing a twin experiment.  i.e. generating observations with a forward run given a true control vector, then using those observations to retrieve with the inverse VDA algorithm this true control vector starting from an a priori background value. 
Here the control vector consists in one unknown discrete parameter among inflow hydrograph or channel friction or bathymetry, and it is retrieved from water level observations in the channel.
 **TODO ref to eqs & inv algo **

The presented case correspond to the Mac Donald's type 1D  solution (see shallow water analytical solutions in SWASHES cases https://hal.archives-ouvertes.fr/hal-00628246/file/SW_analytic-complements.pdf).


----------------------------------------------
Mac Donald's type 1D solution
----------------------------------------------

In `/dassflow2d-wrap/cases/tuto_cases/2_tuto_twin-expe.py` open the script `1_tuto_twin-experiment_classic-way.py`.


.. warning::

    Define the location of your working directory by setting appropriate value to **dassflow_dir**  in the python script



.. dropdown:: 1_tuto_twin-experiment_classic-way.py

  .. code-block:: python

    ##################################################################
    ##################################################################
    # PERFORM TWIN EXPERIMENT (DIRECT RUN + INFERENCE RUN)
    #
    # Method "old way" using command system and pre-existing files
    ##################################################################
    ##################################################################


    #=======================================================#
    #  define librairies librairies
    #=======================================================#

    import dassflow2d as df2d              # dassflow2d : main package
    import os                        # os, for shell command execution (mainly for file manipulation)
    import numpy as np
    import matplotlib.pyplot as plt
    # for os.chdir() command, to open current birectory as bin_directory, defined in parameters, just below

    from subprocess import call # terminal comand

    import h5py # source hdf5 files


    # method that copy pre-existing files to bin dir, used in the script
    # source is the path of the source file and target is the path where to copy the file
    # it directly use the shell
    def cp_inference_file(source, target):
    source = f'inference_files/{source}'
    target = f'{target}'
    call(['rm', target])
    call(['cp', source, target]) # Linux


    # copy all reference files (for direct run)
    def cp_reference_file():

    cp_inference_file(source = 'hydrograph_target.txt' ,
                      target = 'hydrograph.txt')

    cp_inference_file(source = 'mesh_target.txt' ,
                      target = 'automaticaly_generated_mesh.txt')

    cp_inference_file(source = 'land_uses_target.txt' ,
                      target = 'land_uses.txt')

    cp_inference_file(source = 'input_direct.txt' ,
                      target = 'input.txt')

    # cp_dir copy files from source directory to target directory.
    # note that the directory is copied (not only the files within)
    def cp_dir(source, target):
    call(['cp', '-a', source, target]) # Linux

    # save hdf5 (results) files out of /res directory
    def cp_hdf5_file(source, target):
    source = f'res/{source}'
    target = f'hdf5/{target}'
    #    call(['mkdir', "hdf5"])
    call(['rm', target])
    call(['cp', source, target]) # Linux


    # define paths
    dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"
    bin_dir = f"{dassflow_dir}/code/bin_A"


    # delete all files in your simulation directory before starting
    os.system(f"rm -r {dassflow_dir}/code/bin_A/*")
    # Copy recursively the files provided in DassFlow case repository into your own simulation directory **code/bin_A/**.
    os.system(f"cp -r {dassflow_dir}//cases/tuto_case/7_tuto_twin-expe/bin_A/* {dassflow_dir}/code/bin_A")
    os.chdir( f"{dassflow_dir}/code/")
    os.system(f"make cleanres cleanmin ")

    #=======================================================#
    # set up true configuration
    #=======================================================#

    os.chdir(bin_dir)
    cp_reference_file()

    #=======================================================#
    # direct run of the  model
    #=======================================================#A


    # initialise fortran instance, and python corrponding data
    direct_model = df2d.DassFlowModel(bin_dir = bin_dir, run_type = "direct")
    direct_model.build_grid() # build pyvista unstructured grid from model.mesh object


    # define initial conditions
    direct_model.model.dof0.h[:] = 1
    direct_model.model.dof0.u[:] = 0
    direct_model.model.dof0.v[:] = 0
    direct_model.model.dof = direct_model.model.dof0

    direct_model.update_fortran()
    direct_model.run() # run model




    direct_model.save_res() # save simulation results in hdf5 files
    cp_hdf5_file(source = "simu.hdf5", target = "true.hdf5") #save hdf5 (results) files out of /res directory
    cp_dir('./res/obs', '.') # copy observation files


    #----------- plot
    #direct_model.plot_var(what = "bathy", when = "initial", title_plot = "bahtymetry")
    #direct_model.plot_var(what = "h", when = "initial", title_plot = "h")
    #direct_model.plot_var(what = "u", when = "initial", title_plot = "u")
    #direct_model.plot_var(what = "v", when = "initial", title_plot = "v")
    #direct_model.plot_var(what = "manning_alpha", when = "initial", title_plot = "Manning-Strickler coefficient")


    a = h5py.File(name = f"{bin_dir}/res/simu.hdf5")



    for i in range(a["spatial_var/bathy"].shape[1]):
    b=a["spatial_var/bathy"][:,i]
    zs = a["spatial_var/h"][:,i] +  a["spatial_var/bathy"][:,0]

    plt.scatter(x = np.arange(len(zs )),
                y = zs,
                linewidths = 1 )

    plt.scatter(x = np.arange(len(b )),
                y = b,
                c="red",
                linewidths = 1 )
    plt.show()





    #bc = df2d.m_model.get_bc()

    df2d.wrapping.call_model.clean_model(direct_model.model)         # deallocate correctly (necessary action)


    # clean eventual hdf5 file open
    import gc
    for obj in gc.get_objects():   # Browse through ALL objects
    if isinstance(obj, h5py.File):   # Just HDF5 files
        try:
            obj.close()
        except:
            pass # Was already closed

    ###########################################################
    #===========================================================
    # RUN INFERENCE
    #===========================================================
    ###########################################################


    #----------------------#
    #  Define Parameters
    #----------------------#

    os.chdir( f"{dassflow_dir}/code/")
    os.system(f"make cleanres cleanmin ")



    os.chdir(bin_dir)
    os.system(f"rm restart.bin ")


    # /!\ warning :: bathymetry not inferable ? --> lilian = optim not find optimum
    print("CHOOSE INFERENCE TYPE (1 hydrograph, 2 land_use, 3 = bathy)")
    inference_type = input("Enter 1,2 or 3 \n")

    if inference_type == "1":
    # -- infer hydrograph --
    cp_reference_file()
    cp_inference_file(source = 'hydrograph_prior.txt' ,
                      target = 'hydrograph.txt')

    cp_inference_file(source = 'input_hydrograph.txt' ,
                      target = 'input.txt')
    elif inference_type == "2":
    # -- same for manning -- #
    cp_reference_file()
    cp_inference_file(source = 'land_uses_prior.txt' ,
                      target = 'land_uses.txt')

    cp_inference_file(source = 'input_land_uses.txt' ,
                      target = 'input.txt')

    elif inference_type == "3":
    # -- infer hydrograph --#
    cp_reference_file()
    cp_inference_file(source = 'mesh_prior.txt' ,
                          target = 'automaticaly_generated_mesh.txt.txt')

    cp_inference_file(source = 'input_bathy.txt' ,
                          target = 'input.txt')




    #=======================================================#
    # Inference
    #=======================================================#


    my_model = df2d.DassFlowModel(bin_dir = bin_dir, run_type = "min")


    my_model.model.dof0.h[:] = 1
    my_model.model.dof0.u[:] = 0
    my_model.model.dof0.v[:] = 0

    my_model.update_fortran()
    my_model.run()

    # save simulation results in hdf5 files
    my_model.save_res()

    cp_hdf5_file(source = "simu.hdf5", target = "inference_manning.hdf5")

    inf = h5py.File(name =  f"{bin_dir}/hdf5/inference_manning.hdf5", mode = "r")


    #=======================================================#
    # LOAD HDF5 files
    #=======================================================#

    # get grid object
    my_model.build_grid()
    grid = my_model.grid
    grid.plot(show_edges =True, cpos = "xy")

    true = h5py.File(name = f"{bin_dir}/hdf5/true.hdf5", mode = "r")

    plot_title="null"
    if inference_type == "1":
    # -- infer hydrograph --
    true =  [20,20]
    inf = np.loadtxt("min/hydrograph_001.006")[:,1]
    prior = np.loadtxt("min/hydrograph_001.000")[:,1]
    plot_title="Inference of hydrograph parameter: Q(m3/s)"
    elif inference_type == "2":
    # -- same for manning -- #
    true =  true["spatial_var/manning_alpha"][:,0]
    inf = np.loadtxt("min/manning.039")[:,1]
    prior = np.loadtxt("min/manning.000")[:,1]
    plot_title="inference of manning parameter: n"

    elif inference_type == "3":
    true =  true["spatial_var/bathy"][:,0]
    inf = np.loadtxt("min/bathy.039")[:,1]
    prior = np.loadtxt("min/bathy.000")[:,1]
    plot_title="inference of bathymetry parameter (m)"



    plt.scatter( x = np.arange(len(true)), y = true, label = "target")
    plt.scatter( x = np.arange(len(inf)), y = inf, label = "infered")
    plt.scatter( x = np.arange(len(prior)), y =prior, label ="prior")
    plt.title(plot_title)
    plt.legend()

    df2d.wrapping.call_model.clean_model(my_model.model)         # deallocate correctly (necessary action)









++++++++++++++++++++++++++++++++
Script explained by block
++++++++++++++++++++++++++++++++




First, the definition necessary librairies:

.. code-block:: python

  ##################################################################
  ##################################################################
  # PERFORM TWIN EXPERIMENT (DIRECT RUN + INFERENCE RUN)
  #
  # Method "old way" using command system and pre-existing files
  ##################################################################
  ##################################################################


  #=======================================================#
  #  define librairies librairies
  #=======================================================#

  import dassflow2d as df2d              # dassflow2d : main package
  import os                        # os, for shell command execution (mainly for file manipulation)
  import numpy as np
  import matplotlib.pyplot as plt
  # for os.chdir() command, to open current birectory as bin_directory, defined in parameters, just below

  from subprocess import call # terminal comand

  import h5py # source hdf5 files


  # method that copy pre-existing files to bin dir, used in the script
  # source is the path of the source file and target is the path where to copy the file
  # it directly use the shell
  def cp_inference_file(source, target):
      source = f'inference_files/{source}'
      target = f'{target}'
      call(['rm', target])
      call(['cp', source, target]) # Linux


  # copy all reference files (for direct run)
  def cp_reference_file():

      cp_inference_file(source = 'hydrograph_target.txt' ,
                        target = 'hydrograph.txt')

      cp_inference_file(source = 'mesh_target.txt' ,
                        target = 'automaticaly_generated_mesh.txt')

      cp_inference_file(source = 'land_uses_target.txt' ,
                        target = 'land_uses.txt')

      cp_inference_file(source = 'input_direct.txt' ,
                        target = 'input.txt')

  # cp_dir copy files from source directory to target directory.
  # note that the directory is copied (not only the files within)
  def cp_dir(source, target):
    call(['cp', '-a', source, target]) # Linux

  # save hdf5 (results) files out of /res directory
  def cp_hdf5_file(source, target):
      source = f'res/{source}'
      target = f'hdf5/{target}'
  #    call(['mkdir', "hdf5"])
      call(['rm', target])
      call(['cp', source, target]) # Linux




Importing the study case and setting default input data:

.. code-block:: python


    dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"
    bin_dir = f"{dassflow_dir}/code/bin_A"


    print("Would you like to paste case files and compile DassFLow")
    args = input("Press Y or N to continue.") # ajouter exit si pas Y ou N

    if args == "Y":
    # delete all files in your simulation directory before starting
    os.system(f"rm -r {dassflow_dir}/code/bin_A/*")
    # Copy recursively the files provided in DassFlow case repository into your own simulation directory **code/bin_A/**.
    os.system(f"cp -r {dassflow_dir}//cases/tuto_case/7_tuto_twin-expe/* {dassflow_dir}/code/bin_A")
    os.chdir( f"{dassflow_dir}/code/")
    os.system(f"make cleanres cleanmin ")
    os.system(f"make install")
    #importlib.reload(df2d)

    import dassflow2d as df2d


Importing the "true" or "target" parameter configuration for the twin experiment:

.. code-block:: python

  # define paths
  dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"
  bin_dir = f"{dassflow_dir}/code/bin_A"


  # delete all files in your simulation directory before starting
  os.system(f"rm -r {dassflow_dir}/code/bin_A/*")
  # Copy recursively the files provided in DassFlow case repository into your own simulation directory **code/bin_A/**.
  os.system(f"cp -r {dassflow_dir}//cases/tuto_case/7_tuto_twin-expe/bin_A/* {dassflow_dir}/code/bin_A")
  os.chdir( f"{dassflow_dir}/code/")
  os.system(f"make cleanres cleanmin ")

  #=======================================================#
  # set up true configuration
  #=======================================================#

  os.chdir(bin_dir)
  cp_reference_file()




Performing the direct run and saving the observations:


.. code-block:: python

  #=======================================================#
  # direct run of the  model
  #=======================================================#A


  # initialise fortran instance, and python corrponding data
  direct_model = df2d.DassFlowModel(bin_dir = bin_dir, run_type = "direct")
  direct_model.build_grid() # build pyvista unstructured grid from model.mesh object


  # define initial conditions
  direct_model.model.dof0.h[:] = 1
  direct_model.model.dof0.u[:] = 0
  direct_model.model.dof0.v[:] = 0
  direct_model.model.dof = direct_model.model.dof0

  direct_model.update_fortran()
  direct_model.run() # run model




  direct_model.save_res() # save simulation results in hdf5 files
  cp_hdf5_file(source = "simu.hdf5", target = "true.hdf5") #save hdf5 (results) files out of /res directory
  cp_dir('./res/obs', '.') # copy observation files




We can have a look at water height evolution:


.. code-block:: python

  direct_model.build_grid()


    a = h5py.File(name = f"{bin_dir}/res/simu.hdf5")



    for i in range(a["spatial_var/bathy"].shape[1]):
        b=a["spatial_var/bathy"][:,i]
        zs = a["spatial_var/h"][:,i] +  a["spatial_var/bathy"][:,0]

        plt.scatter(x = np.arange(len(zs )),
                    y = zs,
                    linewidths = 1 )

        plt.scatter(x = np.arange(len(b )),
                    y = b,
                    c="red",
                    linewidths = 1 )
        plt.show()





  df2d.wrapping.call_model.clean_model(direct_model.model)         # deallocate correctly (necessary action)





Performing the inverse run:

- Importing "prior" or "first guess" data
- Running the inverse model



.. code-block:: python

  ###########################################################
  #===========================================================
  # RUN INFERENCE
  #===========================================================
  ###########################################################


  #----------------------#
  #  Define Parameters
  #----------------------#

  os.chdir( f"{dassflow_dir}/code/")
  os.system(f"make cleanres cleanmin ")



  os.chdir(bin_dir)
  os.system(f"rm restart.bin ")


  # /!\ warning :: bathymetry not inferable ? --> lilian = optim not find optimum
  print("CHOOSE INFERENCE TYPE (1 hydrograph, 2 land_use, 3 = bathy)")
  inference_type = input("Enter 1,2 or 3 \n")

  if inference_type == "1":
      # -- infer hydrograph --
      cp_reference_file()
      cp_inference_file(source = 'hydrograph_prior.txt' ,
                        target = 'hydrograph.txt')

      cp_inference_file(source = 'input_hydrograph.txt' ,
                        target = 'input.txt')
  elif inference_type == "2":
      # -- same for manning -- #
      cp_reference_file()
      cp_inference_file(source = 'land_uses_prior.txt' ,
                        target = 'land_uses.txt')

      cp_inference_file(source = 'input_land_uses.txt' ,
                        target = 'input.txt')

  elif inference_type == "3":
      # -- infer hydrograph --#
      cp_reference_file()
      cp_inference_file(source = 'mesh_prior.txt' ,
                            target = 'automaticaly_generated_mesh.txt.txt')

      cp_inference_file(source = 'input_bathy.txt' ,
                            target = 'input.txt')





Performing inference and saving the direct run and minimization results:



.. code-block:: python


  my_model = df2d.DassFlowModel(bin_dir = bin_dir, run_type = "min")

  my_model.model.dof0.h[:] = 1
  my_model.model.dof0.u[:] = 0
  my_model.model.dof0.v[:] = 0

  my_model.update_fortran()
  my_model.run()

  # save simulation results in hdf5 files
  my_model.save_res()

  cp_hdf5_file(source = "simu.hdf5", target = "inference_manning.hdf5")

  inf = h5py.File(name =  f"{bin_dir}/hdf5/inference_manning.hdf5", mode = "r")




Plotting the inferred parameters and resulting hydraulic states:


.. code-block:: python

  #=======================================================#
  # LOAD HDF5 files
  #=======================================================#

  # get grid object
  my_model.build_grid()
  grid = my_model.grid
  grid.plot(show_edges =True, cpos = "xy")

  true = h5py.File(name = f"{bin_dir}/hdf5/true.hdf5", mode = "r")

  plot_title="null"
  if inference_type == "1":
      # -- infer hydrograph --
      true =  [20,20]
      inf = np.loadtxt("min/hydrograph_001.006")[:,1]
      prior = np.loadtxt("min/hydrograph_001.000")[:,1]
      plot_title="Inference of hydrograph parameter: Q(m3/s)"
  elif inference_type == "2":
      # -- same for manning -- #
      true =  true["spatial_var/manning_alpha"][:,0]
      inf = np.loadtxt("min/manning.039")[:,1]
      prior = np.loadtxt("min/manning.000")[:,1]
      plot_title="inference of manning parameter: n"

  elif inference_type == "3":
      true =  true["spatial_var/bathy"][:,0]
      inf = np.loadtxt("min/bathy.039")[:,1]
      prior = np.loadtxt("min/bathy.000")[:,1]
      plot_title="inference of bathymetry parameter (m)"



  plt.scatter( x = np.arange(len(true)), y = true, label = "target")
  plt.scatter( x = np.arange(len(inf)), y = inf, label = "infered")
  plt.scatter( x = np.arange(len(prior)), y =prior, label ="prior")
  plt.title(plot_title)
  plt.legend()

  df2d.wrapping.call_model.clean_model(my_model.model)         # deallocate correctly (necessary action)
