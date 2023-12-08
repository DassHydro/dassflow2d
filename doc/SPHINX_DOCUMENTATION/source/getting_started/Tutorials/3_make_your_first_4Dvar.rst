.. _2_make_your_first_4Dvar:

===================================
Inference run
===================================

This tutorial details how to perform an inverse run with `dassflow2d`. The full process is detailed for building and performing a twin experiment.  i.e. generating observations with a forward run given a true control vector, then using those observations to retrieve with the inverse VDA algorithm this true control vector starting from an a priori background value.
Here the control vector consists in one unknown discrete parameter among inflow hydrograph or channel friction or bathymetry, and it is retrieved from water level observations in the channel.
**TODO ref to eqs & inv algo **

The presented case correspond to the Mac Donald's type 1D  solution (see shallow water analytical solutions in SWASHES cases https://hal.archives-ouvertes.fr/hal-00628246/file/SW_analytic-complements.pdf).

.. jupyter-execute::

    import dassflow2d as df2d              # dassflow2d : main package
    import os                        # os, for shell command execution (mainly for file manipulation)
    import numpy as np
    import matplotlib.pyplot as plt
    # for os.chdir() command, to open current birectory as bin_directory, defined in parameters, just below
    from subprocess import call # terminal comand 
    
    ##############
    # put this in package ? assim/utils ?
    ##############
   
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
    #    call(['rm', target]) 
        call(['cp', source, target]) # Linux
    
    #######################################
    #######################################
    
    os.chdir('../../')
    dassflow_dir = os.getcwd() # DassFlow directory (you can also impose your absolute path)
    os.chdir(dassflow_dir)
    print("DassFlow directory is: ", dassflow_dir)
     
    # Define directory where case is run 
    # (its name 'bin_A' is imposed in  {dassflow_dir}/code/makefile.inc : CASEDIR='bin_A')
    bin_dir = f"{dassflow_dir}/code/bin_A/" 
     
    # Define directory containing case data
    case_data_dir = f"{dassflow_dir}/cases/tuto_case/2_qin/bin_A/"

    # Delete all files in your simulation directory before starting
    os.system(f"rm -r {dassflow_dir}/code/bin_A/*")
    
    # Copy recursively the files provided in DassFlow case repository into your own simulation directory **code/bin_A/**.
    os.system(f"cp -r {dassflow_dir}//cases/tuto_case/2_tuto_twin-expe/bin_A/* {dassflow_dir}/code/bin_A")
    os.chdir( f"{dassflow_dir}/code/")
    os.system(f"make cleanres cleanmin ")

    #=======================================================#
    # set up true configuration
    #=======================================================#

    os.chdir(bin_dir)
    cp_reference_file()
        
    #=======================================================#
    # direct run of the  model
    #=======================================================#

    # initialise fortran instance, and python corrponding data
    df2d.wrapping.read_input(f"{bin_dir}/input.txt")
    df2d.wrapping.m_mpi.init_mpi()
    direct_model = df2d.dassflowmodel(bin_dir = bin_dir, hdf5_path = f"{dassflow_dir}/code/bin_A/res/simu.hdf5", run_type = "direct", clean = True)
    # then intialise meshing

    direct_model.config.get()

    direct_model.init_all()
    # define initial conditions
    direct_model.kernel.dof0.h[:] = 1
    direct_model.kernel.dof0.u[:] = 0
    direct_model.kernel.dof0.v[:] = 0
    direct_model.kernel.dof = direct_model.kernel.dof0
    direct_model.config.set({"w_obs":"1", "use_Zobs":"1"})
    direct_model.run()

    direct_model.save_all() # save simulation results in hdf5 files

    cp_hdf5_file(source = "simu.hdf5", target = "true.hdf5") #save hdf5 (results) files out of /res directory
    cp_dir('./res/obs', '.') # copy observation files
 
    df2d.wrapping.call_model.clean_model(direct_model.kernel)         # deallocate correctly (necessary action)

    #----------- some plots to add
    
    ###########################################################
    #===========================================================
    # RUN INFERENCE
    #===========================================================
    ###########################################################

    #----------------------#
    #  Define Parameters
    #----------------------#

    os.chdir( f"{dassflow_dir}/code/")
    #os.system(f"make cleanres cleanmin ")

    os.chdir(bin_dir)
    #os.system(f"rm restart.bin ")
    
    # /!\ warning :: bathymetry not inferable ? --> lilian = optim not find optimum
    #print("CHOOSE INFERENCE TYPE (1 hydrograph, 2 land_use, 3 = bathy)")
    inference_type = 1 #input("Enter 1,2 or 3 \n")

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
    df2d.wrapping.m_mpi.init_mpi()
    df2d.wrapping.read_input(f"{bin_dir}/input.txt")
    my_model = df2d.dassflowmodel(bin_dir = bin_dir, hdf5_path = f"{dassflow_dir}/code/bin_A/res/simu.hdf5", run_type = "min")
    my_model.config.get()
    my_model.config.set({"use_obs":"1", "use_Zobs":"1"})
    # then intialise meshing
    my_model.init_all()
    # define initial conditions
    my_model.kernel.dof0.h[:] = 1
    my_model.kernel.dof0.u[:] = 0
    my_model.kernel.dof0.v[:] = 0
    my_model.kernel.dof = my_model.kernel.dof0
    my_model.run() # only inference is performed




