.. _3_GR4params:

===============================
GR4params.txt
===============================


+++++++++++++++++++++++++++++
File presentation
+++++++++++++++++++++++++++++

The GR4 parameters are defined per hydrological catchment (4 parameters each).
The GR4params.txt file format is given in the table below :

.. dropdown:: template of file

      ===============================================================================================
      $ comment line

      $ comment line

      $ comment line

      +---------------------------------------------------------------------------------------------+

      number-of-hydrological-catchments

      +---------------------------------------------------------------------------------------------+

      $ comment line

      parameter-x1 parameter-x2 parameter-x3 parameter-x4

      $ comment line

      init-state-res1 .. init-state-res13

      $ comment line

      catchment-injection-cell (deprecated) catchment-surface(m2)

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


