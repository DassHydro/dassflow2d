.. _5_3_hpresc:

===============================
hpresc.txt
===============================

the hpresc.txt file enable the use of tabulated values of the relation :math:`h(t)`.



+++++++++++++++++++++++++++++
File presentation
+++++++++++++++++++++++++++++
The file is organised as follow

.. dropdown:: template of file

      ===============================================================================================
      $ comment line

      $ comment line

      $ comment line

      +---------------------------------------------------------------------------------------------+

       number-of-hpresc-curves (**must be equal to 1 presently**)

      +---------------------------------------------------------------------------------------------+

      $ comment line

      $ comment line

      $ comment line

     +---------------------------------------------------------------------------------------------+

      number-of-rows

     +---------------------------------------------------------------------------------------------+

     t11 hQ11

     t12 h12

     ...
     +---------------------------------------------------------------------------------------------+









++++++++++++++++++++++++++++++++++++++++++
Take advantage of wrapped version
++++++++++++++++++++++++++++++++++++++++++

In the directory ``/dassflow2d-wrap/Tools/2_gen_basic_channel``, execute the commands:


.. code-block:: bash

  python3 -m numpy.f2py -c -m gen_channel_case _gen_channel_case.f90

_________________________________________
generate inputs from python
_________________________________________
.. warning::

  **path to be updated**


.. warning::

  the generation of file is coded to generate a unique hpresc curve presently (case of m ulti-inflow-boundary is not coded)


.. code-block:: python

    import os
    import numpy as np

    os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/Tools/2_gen_basic_channel")
    import gen_channel_case


    hpresc = np.ndarray( (20,2) )
    hpresc[:,0] = np.arange(start=0, stop = 1000, step =1000/hpresc.shape[0])
    hpresc[:,1] = 2


    gen_channel_case.gen_bc_data(bc_typ = "hpresc",
                                 nrow=hpresc.shape[0],
                                 var1 =hpresc[:,0],
                                 var2= hpresc[:,1])


**you can paste in your bin directory the file hpresc.txt produced**, since you also change the bc_type, you **must**
replace in your bc.txt file, the type of outflow, from 'ratcurve' to 'hpresc'.

_____________________________________________________
Access and Update fortran kernel value in python
_____________________________________________________



.. code-block:: python

    import dassflow2d as df2d

    # initialise fortran instance, and python corrponding data
    my_model = df2d.DassFlowModel(bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/", run_type = "direct") # run_type can be min or direct (grad ?)
    my_model.update_fortran()

    # =========================
    # get bc values
    # ========================
    my_bc = df2d.wrapping.m_model.get_bc()
    print(my_bc.hpresc[0].h)
    print(my_bc.hpresc[0].q)

    # =========================
    # set bc values
    # ========================
    my_bc.hpresc[0].q[:] =   my_bc.hpresc[0].q[:] + 10


    # =========================
    # feed fortran kernel with new information
    # =========================

    df2d.wrapping.m_model.set_bc(my_bc)


    # =========================
    # run model
    # ========================
    my_model.run()
    my_model.save_res()
