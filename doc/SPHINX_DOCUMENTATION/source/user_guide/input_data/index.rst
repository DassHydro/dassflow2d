.. _input_data_guide:


=================
Input files
=================

Input files formats and related Python scripts are listed below:

.. _table-input:

+----------------------------------+------------------------------------------------------+----------------------------------------+
| File                             |  Description                                         | is necessary ?                         |
+==================================+============+=========================================+========================================+
| :ref:`input.txt <1_inputfile>`   | Defines all variables concerning the numerical       | yes                                    |
|                                  | and physical parameters  and input/output options.   |                                        |
+----------------------------------+------------------------------------------------------+----------------------------------------+
| :ref:`mesh.geo <2_meshfile>`     | Defines node coordinate and cells correspondance to  | yes                                    |
|                                  | nodes. As well as boundaries type.                   |                                        |
|                                  | Some model parameters are also defined in mesh file  |                                        |
|                                  | (such as land type for friction, bathymetry)         |                                        |
+----------------------------------+------------------------------------------------------+----------------------------------------+
| :ref:`land_uses.txt <3_landuses>`| Defines correspondance between mesh's land type      | yes                                    |
|                                  | and actual value                                     |                                        |
+----------------------------------+------------------------------------------------------+----------------------------------------+
| :ref:`bc.txt <4_bcfile>`         | Defines correspondance between mesh.geo file         | yes                                    |
|                                  | (bc type) and   the following .txt files             |                                        |
|                                  | (values applied to the bc)                           |                                        |
+----------------------------------+------------------------------------------------------+----------------------------------------+
| :ref:`hydrograph.txt<5_1_hyd>`   | Contains tables of values to apply on the boundary   | For boundart                           |
|                                  |                                                      | types specified in bc. txt             |
| :ref:`ratcurve.txt<5_2_ratcurve>`|                                                      |                                        |
|                                  |                                                      |                                        |
| :ref:`hpresc.txt  <5_3_hpresc>`  |                                                      |                                        |
|                                  |                                                      |                                        |
| :ref:`zpresc.txt  <5_4_zpresc>`  |                                                      |                                        |
+----------------------------------+------------------------------------------------------+----------------------------------------+
| :ref:`obs.txt  <6_obs>`          | Defines parameter for observed data for data         | For minimization run                   |
|                                  | assimilation. (localisation of observed data and     |                                        |
|                                  | timestep availability                                |                                        |
+----------------------------------+------------------------------------------------------+----------------------------------------+
| :ref:`rain.txt <8_rain>`         | to be documented                                     | no                                     |
+----------------------------------+------------------------------------------------------+----------------------------------------+
| :ref:`GR4forcings_X.txt <GR4for>`| to be documented                                     | no                                     |
+----------------------------------+------------------------------------------------------+----------------------------------------+
| :ref:`GR4warmup.txt <GR4warmup>` | to be documented                                     | no                                     |
+----------------------------------+------------------------------------------------------+----------------------------------------+
| :ref:`GR4params.txt <GR4params>` | to be documented                                     | no                                     |
+----------------------------------+------------------------------------------------------+----------------------------------------+
| :ref:`infil.txt  <9_infil>`      | to be documented                                     | no                                     |
+----------------------------------+------------------------------------------------------+----------------------------------------+
| :ref:`ic.bin  <7_ic>`            | to be documented                                     | no                                     |
+----------------------------------+------------------------------------------------------+----------------------------------------+

.. note ::

  Most files can be generated in basic case using the script run.py after neccessary installs in
  `/dassflow2d-wrap/Tools/2_gen_basic_channel`


  Most data provided within input files can be modified a posteriori within python script
  nonetheless, only the values of variables can be modified. (Meaning that the structure and size
  of the object is fixed and can't be modified once the input files are read)



++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Detailed presentation of each input files
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. toctree::
   :maxdepth: 2
   :titlesonly:
   :hidden:

  1_input-txt.rst
  2_mesh-geo.rst
  3_land-uses.rst
  4_bcs-txt.rst
  5_1_hydrograph.rst
  5_2_ratingcurve.rst
  5_3_hpresc.rst
  5_4_zpresc.rst
  6_obs.rst
  7_ic.rst
  8_rain.rst
  9_infil.rst
