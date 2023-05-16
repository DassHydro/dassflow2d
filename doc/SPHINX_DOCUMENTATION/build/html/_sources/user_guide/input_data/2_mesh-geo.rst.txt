.. _2_meshfile:

===============================
mesh.geo
===============================


++++++++++++++++++++++++++++++++++++++++++
File presentation
++++++++++++++++++++++++++++++++++++++++++


The mesh file must be created outside ``DassFlow`` by a mesh generator software [http://www.robertschneiders.de/meshgeneration/software.html].
To generate basic mesh (rectangular, with a rectangular extent) you can use the tools provided in  ``/dassflow2d-wrap/Tools/2_gen_basic_channel`` (see section **xxxx** below).




This file contains :

• A first list of the mesh nodes with (x,y) coordinates and eventually a z-coordinate defining the bed elevation (**z-coordinates for node == obsolete?**).

• A second list of the mesh elements (cells and eventually edges) defined by a suited list of the previous defined nodes.

By default, a wall type boundary is applied at the boundary :math:`\partial\Omega` of the computational domain \Omega.
Boundary edges defining a subset :math:`\Gamma` of the computational boundary `\partial\Omega` where is considered an
inflow or outflow boundary condition must be provided in the mesh file. Each subset :math:`\Gamma` is also defined by a group number and the boundary type
to apply is controlled by the **bc.txt** file.



DassFlow inhouse mesh format is described in the table below. It has been created in order to consider special issues in the Shallow Water model :
the bed elevation :math:`z_{b}` and the Manning-Strickler roughness coefficient :math:`n`. It also contains information about boundary cells correspondance with the the **bc.txt** file.



.. dropdown:: file channel.geo generated

  .. include:: input_default_files/channel.geo
     :literal:

.. dropdown:: template of file

      
      $ comment line

      number-of-nodes number-of-cells scaling

      +--------------------------------------------------------------------------------+

      $ comment line

      node-index      x-coordinate y-coordinate bed-elevation (obsolete data: must exist but is unused (can be full values of 0 ))

      ...

      +--------------------------------------------------------------------------------+

      $ comment line

      cell-index node1 node2 node3 node4 land-code bed-elevation

      ...

      +--------------------------------------------------------------------------------+

      $ comment line

      INLET number-of-cells number-of-inlet

      cell-index edge-index boundary-type ghost-cell-bed-elevation[, group-number]

      ...

      +--------------------------------------------------------------------------------+

      OUTLET number-of-cells number-of-outlet

      cell-index edge-index boundary-type ghost-cell-bed-elevation[, group-number]

      +--------------------------------------------------------------------------------+


.. note ::

    x and y coordinates of nodes must be provided in cartesian system.

.. note ::

    the bed elevation is at the cells gravity center (piecewise constant approximation over mesh cells)

    .. warning::

      the bed elevation definition using nodes is obsolete

.. note ::

    The land code gives correspondance with land_uses.txt file.

.. note ::

   Concerning INLET and OUTLET boundaries, **[, group-number]** must be defined in case in multiple INLET or OUTLET (it is the way to distinguish multiple bc).

   In case of a unique inlet and a unique outlet (and if **number-of-inlet** = 0 & **number-of-outlet** = 0 ), the value of group number is not used


 The **boundary_type** field is obsolete and must exist but is unused.




____________________________________________
Triangular or quadrangular meshing meshing
____________________________________________


The fortran kernel can process any quadrangular mesh as well as triangular meshing.

Nothing has to be updated to process quadrangular meshing, while to apply triangular meshing, the **4th** node
must have the same id as the **1st** node of the cell.





++++++++++++++++++++++++++++++++++++++++++
Take advantage of wrapped version
++++++++++++++++++++++++++++++++++++++++++

In the directory ``/dassflow2d-wrap/Tools/2_gen_basic_channel``, execute the commands:


.. code-block:: bash

  python3 -m numpy.f2py -c -m gen_channel_case _gen_channel_case.f90

_____________________________________________
Generate recatangular meshing from python
_____________________________________________

As in the presentation of input.txt file go to the directory  ``/dassflow2d-wrap/Tools/2_gen_basic_channel``
You can open the script in   ``script_example/3_define_mesh.py``, it corresponds to the following code:

.. warning::

  **path to be updated**

.. code-block:: python

  os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/Tools/2_gen_basic_channel")
  import gen_channel_case
  gen_channel_case.gen_basic_channel(
                                      nx = 10,   #number of nodes in x direction
                                      ny = 2,    #number of nodes in y direction
                                      lx = 1000, #length of channel in x direction
                                      ly = 10)   #length of channel in y direction



You can have a look at the file produced.

You can note that a bathymetry as well as a land_type code has been generated.
Their definition is made in _gen_channel_case.f90, to change bathymetry value, you can update the bathy_user() subroutine in the file.
To generate not unical land type, you have to implement yourself the needed code in  the subroutine gen_basic_channel()

.. note::

  once you modified the code, you need to recompile it. by running again the commands

  .. code-block:: bash

    python3 -m numpy.f2py -c -m gen_channel_case _gen_channel_case.f90
    python3 run.py

_____________________________________________________
Access fortran kernel value in python
_____________________________________________________


Still in the directory ``/dassflow2d-wrap/Tools/2_gen_basic_channel``, execute the commands:


.. code-block:: bash

  python3 -m numpy.f2py -c -m gen_channel_case _gen_channel_case.f90
  python3 run.py

It generates necessary files in  ``/files/`` directory (including channel.geo file). Paste them in your bin directory ( ``/dassflow2d-wrap/code/bin_A/`` ).


Then, you can execute the following code:

.. warning::

  **dassflow_dir** to be updated


.. code-block:: python3

  import dassflow2d as df2d
  dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"
  my_model = df2d.DassFlowModel(bin_dir =  f"{dassflow_dir}/code/bin_A" , arg = "direct")
  my_model.build_grid()
  my_model.plot_meshing(what = "cell")
  my_model.plot_meshing(what = "node")
  my_model.plot_meshing(what = "edge")



The plots shows the values applied and processed in the fortran kernel, you can explore the plots and compare them with the channel.geo file

To access ghost cells/nodes/edge data, you can execute the following commangs


.. code-block:: python3

    # ----------------- PRINT MESH INFO ----------------- #
    for i in range( my_model.model.mesh.neb):

        if my_model.model.mesh.edgeb[i].typlim.decode("utf-8")  != 'wall':
             print(my_model.model.mesh.edgeb[i].ind  )
             print(my_model.model.mesh.edgeb[i].typlim.decode("utf-8")  )
             print( my_model.model.mesh.edgeb[i].group)
             print(my_model.model.mesh.edgeb[i].perio)



    for i in range( my_model.model.mesh.ncb):
             print(my_model.model.mesh.cellb[i].ind,
                   my_model.model.mesh.cellb[i].typlim.decode("utf-8"),
                   my_model.model.mesh.cellb[i].group,
                   my_model.model.mesh.cellb[i].cell,
                   my_model.model.mesh.cellb[i].grav.x,
                   my_model.model.mesh.cellb[i].grav.y)

    for i in range( my_model.model.mesh.nnb):
             print(my_model.model.mesh.nodeb[i].ind,
                   my_model.model.mesh.nodeb[i].typlim.decode("utf-8"),
                   my_model.model.mesh.nodeb[i].group)

.. danger::

  note that not all the variable presented above a correctly filled



_____________________________________________________
Update fortran kernel value in python
_____________________________________________________

...........................
geometry
...........................


If you modify values of object ``my_model.model.mesh.edgeb[i].ind ``, the modification will go to the fortran kernel directly.**Nonetheless it is not recommended**. We recomand to
exclusively use the geometry defined in the mesh.geo file. The same reasonment must be applied to the boundaries correspondance definition.



...........................
Land use
...........................

Land type should be defined correctly in the mesh.geo file.


...........................
Bathymetry
...........................

you can access and modify the bathymetry with the following commands

.. warning::

  **Path to be updated**


.. danger::

  This is the correct method but it is not working (the getter is correctly working but not the setter)

.. code-block:: python3

  import dassflow2d as df2d
  dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"
  my_model = df2d.DassFlowModel(bin_dir =  f"{dassflow_dir}/code/bin_A" , arg = "direct")
  my_model.update_fortran()
  bathy= df2d.m_model.get_array_bathy_cell()
  bathy[:]=1
  df2d.m_model.set_array_bathy_cell(bathy)
  my_model.run()



The way to change bathymetry can be done this way:

.. danger::

  conflit d'interet avec léo, c'est desactivé pour le moment
  mais si réactivé, on sait que ça marche

.. warning::

  **path to be updated**

.. code-block:: python3

  import dassflow2d as df2d
  dassflow_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap"
  my_model = df2d.DassFlowModel(bin_dir =  f"{dassflow_dir}/code/bin_A" , arg = "direct")
  my_model.model.my_param_model.bathy_cell[:] = new array of bathymetry values
  my_model.update_fortran()
  my_model.run()
