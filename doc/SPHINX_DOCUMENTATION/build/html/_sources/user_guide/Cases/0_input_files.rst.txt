.. _1_input_files:

=================
Input files
=================

To perform simulations using dassflow, some input files are required.
We list the necessary files below.



+--------------------------+------------------------------------------------------+----------------------------------------+
| File                     |  Description                                         | is necessary ?                         |
+==========================+============+=========================================+========================================+
| input.txt                | Defines all variables concerning the numerical       | yes                                    |
|                          | and physical parameters  and input/output options.   |                                        |
+--------------------------+------------------------------------------------------+----------------------------------------+
| mesh.geo                 | Defines node coordinate and cells correspondance to  | yes                                    |
|                          | nodes. As well as boundaries type.                   |                                        |
|                          | Some model parameters are also defined in mesh file  |                                        |
|                          | (such as land type for friction, bathymetry)         |                                        |
+--------------------------+------------+-----------------------------------------+----------------------------------------+
| land_uses.txt            | defines correspondance between mesh's land type      | yes                                    |
|                          | and actual value                                     |                                        |
+--------------------------+------------+-----------------------------------------+----------------------------------------+
| bc.txt                   | Defines correspondance between mesh.geo file         | yes                                    |
|                          | (bc type) and   the following .txt files             |                                        |
|                          | (values applied to the bc)                           |                                        |
+--------------------------+------------+-----------------------------------------+----------------------------------------+
| hydrograph.txt           |   contains tables of values to apply on the boundary | necessary for boundaries               |
|                          |                                                      | type specified in bc. txt              |
| ratcurve.txt             |                                                      |                                        |
|                          |                                                      |                                        |
| hpresc.txt               |                                                      |                                        |
|                          |                                                      |                                        |
| zpresc.txt               |                                                      |                                        |
+--------------------------+------------+-----------------------------------------+----------------------------------------+
| obs.txt                  | Defines parameter for observed data for data         | necessary for minimization run         |
|                          | assimilation. (localisation of observed data and     |                                        |
|                          | timestep availability                                |                                        |
+--------------------------+------------+-----------------------------------------+----------------------------------------+
| rain.txt                 |to be documented                                      | no                                     |
+--------------------------+------------+-----------------------------------------+----------------------------------------+
| infil.txt                |to be documented                                      | no                                     |
+--------------------------+------------+-----------------------------------------+----------------------------------------+

.. note ::

  Most files can be generated in basic case using the script run.py after neccessary installs in
  `/dassflow2d-wrap/Tools/2_gen_basic_channel`


  Most data provided within input files can be modified a posteriori within python script
  nonetheless, only the values of variables can be modified. (Meaning that the structure and size
  of the object is fixed and can't be modified once the input files are read)


==================================
Generate input files
==================================

  Go to the directory `/dassflow2d-wrap/Tools/2_gen_basic_channel` and open the note `readme.txt`.


  run the following command in your terminal (still in directory `/dassflow2d-wrap/Tools/2_gen_basic_channel`).

  .. code-block:: bash

    python3 -m numpy.f2py -c -m gen_channel_case _gen_channel_case.f90


  This generates a wrapped package called *gen_channel_case* of subroutines written in *_gen_channel_case.f90*. Then, perform the command:

  .. code-block:: bash

    python3 run.py


  you can see that in run.py, we execute the routine defined in ``gen_dassflow.py``.
  It is in the file ``gen_dassflow.py`` that you can modify the parameters of file generation.
  The fortran subroutines  are documented within fortran script and
  gen_dassflow methods are documented within ``gen_dassflow.py`` script. (**to be done**)



==================================
Explore generated input files
==================================

You can see that multiple files have been generated in *files/* directory :


.. code-block:: bash

  ls files


You can open each files to have a look at what they look like:

- mesh data:

  - channel.geo

- Boundaries data:

  - bc.txt
  - hpresc.txt
  - hydrograph.txt

- spatial variable parameter data:

  - land_uses.txt

- inference data (not necessary to perform direct simulations)

  - obs.txt

- model parameter data:

  - input.txt (to be done)


===========================================
Detailed presentation of each input files
===========================================


- :ref:`Mesh input file <4_mesh>`
- :ref:`Boundary input files <2_bcs>`









The minimal requirements are the files ``input.txt`` which accounts for
**model and simulation parameters**. In the input.txt file you also specify the path to
the ``mesh.geo`` that define the name of the mesh file.

Additionaly, necessary files to define **boundary conditions** have to be provided.
The file ``bc.txt`` defines the boundaries type as well as their correspondance with
the ``meshing``. Depending of the type of boundry defined, you have to provide the
boundary file corresponding (such as ``hydrograph.txt``, ``ratcurve.txt``, ``hpresc.txt``,  ``zpresc.txt`` )
which correspond to tables with two entry (time and h/z/q for most files or Q/h correspondance for rating curve).

Furthermore, some parameter files can be provided, shuch as ``land_uses.txt``.

Finaly, initial conditions can be specified in ic.bin file. **to be documented by a sachant**.
