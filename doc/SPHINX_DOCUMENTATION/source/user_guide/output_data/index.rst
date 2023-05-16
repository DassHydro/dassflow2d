.. _output_data_guide:


=================
Output files
=================
All output result files are written in the **dassflow2d-wrap/code/bin_A/res/** directory.
The frequency and the format(s) are controlled in the user-defined :ref:`input.txt <1_inputfile>`. file respectively by the ``dtw``, ``dtp`` (and ``dta`` == obsolete) time steps
and the output switches described in the presentation of the :ref:`input.txt <1_inputfile>`.
To perform simulations using dassflow, some output files can be produced.


We list the produced files below:


**OUTPUT CONTENT THAT CAN BE SAVED**

.. _table-output:

+------------------+--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| parameter        |  file_name                     | description                                                                                                                      |
+==================+================================+==================================================================================================================================+
| model unknows and| res/                           | The Shallow Water model primitive unknows h, u, v as some model constants z_{b},                                                 |
| model constants  | result_xxxxxxE+xx.yyy'         | n can be written in output result files with the <dtw> time step in some available formats,                                      |
+------------------+--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| Post-processed   | res/post/....                  | Several post-processed variables are also written in output result files with the <dtp> ime step and in the format defined in the|
| variables        |                                | \inp{} file by the switches  n can be written in output result files with the <dtw> time step in some available formats,         |
+                  +--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
|                  | 'time_step.yyy'                | the time step dt in \left[s\right]                                                                                               |
+                  +--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
|                  | 'water_vol.yyy'                | the volume of water in the computational domain in \left[m^{3}\right]                                                            |
+                  +--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
|                  | 'water_vol_num_add.yyy '       | the volume of water in the computational domain in \left[m^{3}\right] eventually added numerically                               |
+                  +--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
|                  | 'sum_mass_flux_inflow_xxx.yyy' | **the discharges in \left[m^{3}.s^{-1}\right] at inflow**                                                                        |
+                  +--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
|                  | 'sum_mass_flux_outflow_xxx.yyy'| **the discharges in \left[m^{3}.s^{-1}\right] at outflow**                                                                       |
+                  +--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
|                  | 'sum_q_inflow_xxx.yyy'         | **the discharges in \left[m^{3}.s^{-1}\right] at inflow**                                                                        |
+                  +--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
|                  | 'sum_q_outflow_xxx.yyy'        | **the discharges in \left[m^{3}.s^{-1}\right] at outflow**                                                                       |
+------------------+--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| observations     | res/obs/....                   | An observation output result file is generated according to the \obs{} file for each defined station and                         |
| variables        |                                | two files for each defined section                                                                                               |
+                  +--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
|                  | 'obs_station_yyyy.xxx'         | model results (dof) at observed station                                                                                          |
+                  +--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
|                  | 'obs_section_yyyy.xxx'         | model results (dof) at observed section                                                                                          |
+                  +--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+
|                  | 'obs_q_h_section_yyyy.xxx'     | model results (Q) at observed section                                                                                            |
+------------------+--------------------------------+----------------------------------------------------------------------------------------------------------------------------------+


**FORMATS CAN BE SAVED**

+------------------+-----------------------------------------------------------------------------------------------------------------------------+---------------+
|format            | description                                                                                                                 | extention name|
+==================+=============================================================================================================================+===============+
|  VTK format      | ASCII encoding setting the <w_vtk> to '1' and in binary encoding setting the <w_vtk> to '2'.                                | yyy='vtk'     |
+------------------+-----------------------------------------------------------------------------------------------------------------------------+---------------+
|  tecplot format  | in ASCII encoding setting the <w_tecplot> to '1'.                                                                           | yyy='plt'     |
+------------------+-----------------------------------------------------------------------------------------------------------------------------+---------------+
|  gnuplot format  | in ASCII encoding setting the <w_gnoplot> to '1'.                                                                           | yyy='dat'     |
+------------------+-----------------------------------------------------------------------------------------------------------------------------+---------------+
|  HDF5 format     | in HDF5 encoding --> **automaticaly done setting the <w_gnoplot> to '1'**                                                   | yyy='hdf5'    |
+------------------+-----------------------------------------------------------------------------------------------------------------------------+---------------+



++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
model unknows and constants:
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

**gnuplot format**


.. dropdown:: template of file

      ==============================================

      $ comment line

      $ comment line

      id-cell x-coord y-coord h zs Manning-alpha u v


**tecplot format**

tecplot format to be documented by a user



**vtk format**

tecplot format to be documented by a user


**hdf5 format**

The results at each timesteps of  'result_xxxxxxE+xx.dat' are stored in a unical compressed file in hdf5 format.
By default, its name is simu.hdf5 and is stored in **/res/** directory.

:download:`hdf5 presentation  <files/hdf5.pdf>`



++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Post processed variables
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

All results are stored in tables (time-value).

**gnuplot format**


.. dropdown:: template of file

      ==============================================

      $ comment line

      time variable-value


The post processed variable value at specific time. **to validate**


**tecplot format**

tecplot format to be documented by a user



**vtk format**

tecplot format to be documented by a user





+++++++++++++++++++++++++++++++++
more detailed presentation
+++++++++++++++++++++++++++++++++

.. toctree::
   :maxdepth: 2
   :titlesonly:

   1_dof-and-more.rst
   2_postpro.rst
   3_assim-post.rst
