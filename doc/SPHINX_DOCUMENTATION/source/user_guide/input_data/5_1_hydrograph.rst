.. _5_1_hyd:

===============================
hydrograph.txt
===============================

the hydrograph.txt file enable the use of tabulated values of the relation :math:`Q(t)` .



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

       number-of-hydrographs

      +---------------------------------------------------------------------------------------------+

      $ comment line

      $ comment line

      $ comment line

     +---------------------------------------------------------------------------------------------+

      number-of-tabulated-times-hydrograph-1

     +---------------------------------------------------------------------------------------------+

      t11 Q11

      t12 Q12

      ...

      $ comment line

      $ comment line

      $ comment line



     +---------------------------------------------------------------------------------------------+

      number-of-tabulated-times-hydrograph-2

     +---------------------------------------------------------------------------------------------+

      t21 Q21

      t22 Q22

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

  the generation of file is coded to generate a unique hydrograph presently (case of m ulti-inflow-boundary is not coded)


.. code-block:: python

    import os
    import numpy as np

    os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/Tools/2_gen_basic_channel")
    import gen_channel_case


    qin = np.ndarray( (20,2) )
    qin[:,0] = np.arange(start=0, stop = 1000, step =1000/qin.shape[0])
    qin[:,1] = 20


    gen_channel_case.gen_bc_data(	bc_typ = "hydrograph",
    									nrow=qin.shape[0],
    									var1 =qin[:,0],
    									var2= qin[:,1])

**you can paste in your bin directory the file hydrograph.txt produced**

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
  print(my_bc.hyd[0].q)
  print(my_bc.hyd[0].t)

  # =========================
  # set bc values
  # ========================
  my_bc.hyd[0].q[:] = 500

  # =========================
  # feed fortran kernel with new information
  # =========================

  df2d.wrapping.m_model.set_bc(my_bc)


  # =========================
  # run model
  # ========================
  my_model.run()
  my_model.save_res()
