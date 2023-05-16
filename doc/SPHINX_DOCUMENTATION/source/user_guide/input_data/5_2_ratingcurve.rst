.. _5_2_ratcurve:

===============================
ratcurve.txt
===============================

the ratcurve.txt file enable the use of tabulated values of the relation :math:`h(Q)`.



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

       number-of-rating-curves

      +---------------------------------------------------------------------------------------------+

      $ comment line

      $ comment line

      $ comment line

     +---------------------------------------------------------------------------------------------+

      number-of-tabulated-water-elevation-1 water-elevation-reference

     +---------------------------------------------------------------------------------------------+

     h11 Q11

     h12 Q12

     ...

      $ comment line

      $ comment line

      $ comment line



     +---------------------------------------------------------------------------------------------+

      number-of-tabulated-water-elevation-2 water-elevation-reference

     +---------------------------------------------------------------------------------------------+

      h21 Q21

      h22 Q22

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

  the generation of file is coded to generate a unique rating curve presently (case of m ulti-inflow-boundary is not coded)


.. code-block:: python

    import os
    import numpy as np

    os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/Tools/2_gen_basic_channel")
    import gen_channel_case


    ratcurve = np.ndarray( (20,2) )
    ratcurve[:,0] = np.arange(start=0, stop = 10, step =10/ratcurve.shape[0])
    ratcurve[:,1] = np.arange(start=0, stop = 100, step =100/ratcurve.shape[0])


    gen_channel_case.gen_bc_data(bc_typ = "ratcurve",
                                 nrow=ratcurve.shape[0],
                                 var1 =ratcurve[:,0],
                                 var2= ratcurve[:,1])


**you can paste in your bin directory the file ratcurve.txt produced**, since you also change the bc_type, you **must**
replace in your bc.txt file, the type of outflow, from 'hpresc' to 'ratcurve'.

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
    print(my_bc.rat[0].h)
    print(my_bc.rat[0].q)

    # =========================
    # set bc values
    # ========================
    my_bc.rat[0].q[:] =   my_bc.rat[0].q[:] + 10


    # =========================
    # feed fortran kernel with new information
    # =========================

    df2d.wrapping.m_model.set_bc(my_bc)


    # =========================
    # run model
    # ========================
    my_model.run()
    my_model.save_res()
