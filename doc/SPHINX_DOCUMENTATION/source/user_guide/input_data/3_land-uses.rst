.. _3_landuses:

===============================
land_uses.txt
===============================


File presentation
+++++++++++++++++++++++++++++

The friction  coefficients are defined as a land use (zoning).
The land_uses.txt file format is given in the table below :

.. dropdown:: template of file


      +---------------------------------------------------------------------------------------------+

      number-of-lands-type

      +---------------------------------------------------------------------------------------------+

      $ comment line

      land-code-1 Manning-Strickler-coefficient-1  manning_beta_coefficient-1

      land-code-2 Manning-Strickler-coefficient-2  manning_beta_coefficient-2

      +---------------------------------------------------------------------------------------------+


According to the prescribed land-code of a cell in the inhouse mesh format file (channel.geo),
the Manning-Strickler coefficient is defined with the associated value to land-code in the land_uses.txt file.

.. hint::

    macro rugosity for **Ferguson's law** to be added here ?



++++++++++++++++++++++++++++++++++++++++++
Take advantage of wrapped version
++++++++++++++++++++++++++++++++++++++++++

In the directory ``/dassflow2d-wrap/Tools/2_gen_basic_channel``, execute the commands:


.. code-block:: bash

      python3 -m numpy.f2py -c -m gen_channel_case _gen_channel_case.f90



_________________________________________
generate land_uses.txt from python
_________________________________________



As in the presentation of input.txt file go to the directory  ``/dassflow2d-wrap/Tools/2_gen_basic_channel``
You can open the script in   ``script_example/4_define_land_uses.py``, it corresponds to the following code:

.. warning::

  **path to be updated**

.. code-block:: python

  os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/Tools/2_gen_basic_channel")
  import gen_channel_case
  gen_channel_case.gen_land_uses(
                                  manning_alpha  = 0.003,
                                  manning_beta = 0)


**you can paste in your bin directory the file land_uses.txt produced**

.. note::

  gen_channel_case can only generate a uniform manning strickler coefficient. More complex generation has to be developed.


_____________________________________________________
Access and update fortran kernel value in python
_____________________________________________________
Approach 1
.........................................

.. code-block:: python

  import dassflow2d as df2d


  # initialise fortran instance, and python corrponding data
  my_model = df2d.DassFlowModel(bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/", run_type = "direct") # run_type can be min or direct (grad ?)


  print("manning value  = ", my_model.model.my_friction.manning )
  print("manning beta value  = ", my_model.model.my_friction.manning_beta)


  my_model.model.my_friction.manning[:]= 0.1
  my_model.model.my_friction.manning_beta[:]= 0.1

  my_model.update_fortran() # here, the variables my_friction.manning are updated to manning and manning_beta
  # global variables of fortran kerne that are used along the code
  my_model.run()


You can check the value of manning (not manning beta :( )  in /res/ directory in a result file.


.........................................
Approach 2
.........................................



.. code-block:: python

    # initialise fortran instance, and python corrponding data
    my_model = df2d.DassFlowModel(bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/", run_type = "direct") # run_type can be min or direct (grad ?)
    my_model.update_fortran() # here, the variables manning and manning_beta are initialised


    n=df2d.m_model.get_array_manning()    # need to be done after update_fortran() call
    nb=df2d.m_model.get_array_manning_beta()

    n[:]  = 2
    nb[:] = 5

    df2d.m_model.set_array_manning(n)
    df2d.m_model.set_array_manning_beta(nb)

    my_model.run()


You can check the value of manning (not manning beta :( )  in /res/ directory in a result file.

