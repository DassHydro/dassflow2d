.. _3_GR4forcings_warmup:

===============================
GR4forcings.txt
===============================


+++++++++++++++++++++++++++++
File presentation
+++++++++++++++++++++++++++++

The GR4 forcings for warmup are defined per hydrological catchment (numbered X, the total number is given in GR4params.txt) for a single year.
The GR4params_X.txt file format is given in the table below:

.. dropdown:: template of file

      ===============================================================================================
      $ comment line
      
      $ comment line
            
      $ comment line
                  
      $ comment line

      +---------------------------------------------------------------------------------------------+

      number-of-time-step-to-read

      +---------------------------------------------------------------------------------------------+
      
      time(s) precipitation(mm/h) evapotransipration(mm/h)

      +---------------------------------------------------------------------------------------------+



++++++++++++++++++++++++++++++++++++++++++
Take advantage of wrapped version
++++++++++++++++++++++++++++++++++++++++++

.. code-block:: python

_____________________________________________________
Access and update fortran kernel value in python
_____________________________________________________

.........................................
Approach 1
.........................................


.........................................
Approach 2
.........................................


