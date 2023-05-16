.. _4_bcfile:

===============================
bc.txt
===============================


+++++++++++++++++++++++++++++
File presentation
+++++++++++++++++++++++++++++

The boundary condition type is prescribed in the **bc.txt** file in the format described in the table below.
A group number is assigned at each subset of the computational domain boundary :math:`\partial\Omega` as explainded
in the previous section: :ref:`mesh.geo <2_meshfile>`.
Thus, the associated type of boundary condition is read in this **bc.txt**  file that provides a relative flexibility.


.. dropdown:: template of file

      ===============================================================================================
      $ comment line

      $ comment line

      $ comment line

      +---------------------------------------------------------------------------------------------+

      number-of-tagged-bc

      +---------------------------------------------------------------------------------------------+

      $ comment line

      $ comment line

      $ comment line

      group-1 bc-type-1 [data-type-1]

      group-2 bc-type-2 [data-type-2]



bc-type take one of the following values:

+----------+------------+-----------------------------------------------------+-----------------+
| Type     | Name       | Description                                         | Data            |
+==========+============+=====================================================+=================+
| Inflow   | 'discharg1'| discharg is imposed (approximate backwater curve)   | hydrograph.txt  |
+          +------------+-----------------------------------------------------+-----------------+
|          | 'discharg2'| discharg is imposed (fixed water elevation)         | hydrograph.txt  |
+----------+------------+-----------------------------------------------------+-----------------+
| Outflow  | 'transm'   |                                                     |                 |
+          +------------+-----------------------------------------------------+-----------------+
|          | 'zspresc'  | water elevation :math:`\eta` is prescribed          | zpresc.txt      |
+          +------------+-----------------------------------------------------+-----------------+
|          | 'hpresc'   | water depth :math:`h` is prescribed                 | hpresc.txt      |
+          +------------+-----------------------------------------------------+-----------------+
|          | 'ratcurve' | rating curve based                                  | ratcurve.txt    |
+----------+------------+-----------------------------------------------------+-----------------+


data-type can take one of the following values


- 'file'
- 'internal' **tobedocumentedLeo**



+++++++++++++++++++++++++++++++++++++++
Detailed presentation of boundary types
+++++++++++++++++++++++++++++++++++++++

_________________________________________
Wall
_________________________________________
The 'wall' boundary condition type is a slip condition since there is no viscous term in the model. It can be also used as a symmetry boundary condition.
It is the default boundary condition assigned at a boundary edge not in a defined subset of the computational domain boundary :math:`\partial\Omega`.

_________________________________________
Inflow
_________________________________________


The inflow boundary condition consists of applying a discharge Q(t) to a subset \Gamma of the computational domain boundary :math:`\partial\Omega`.
The discharge relation :math:`Q(t)` is prescribed either by the hydrograph.txt file.
The data-type option must be setted to 'file' in the **bc.txt** file in order to use the tabulated values of the relation :math:`Q(t)` law in the hydrograph.txt file.

The type 'discharg1' gives a more robust method to prescribed the discharge relation :math:`Q(t)`.
It should be preferred to the 'discharg2' type which can be more precise but can generate a calculation divergence.

When the considered wet surface corresponding to :math:`\Gamma` is non trivial, practice shows that one can observe a wrong prescribed discharge  :math:`Q(t)`.
In order to overcome this drawback, a solution is to apply a feedback process on the associated ghost bed elevations
setting the variable ``feedback_inflow`` to '1' in the **input.txt** file (by default) and eventually changing the associated
feedback process coefficient ``coef_feedback`` (to '0.8' by default).


_________________________________________
Outflow
_________________________________________


The outflow boundary condition can be prescribed in different ways :

• the 'transm' type for which simple homogeneous Neumann conditions are applied to water depth and normal velocity.
  It can be used in simple cases where the ghost bed elevations are well defined (a constant slope channel for example).
  Otherwise, it can generates an improper outflow boundary condition.


• the 'zspresc' and 'hpresc' types where respectively a water surface relation :math:`z_{b}(t)`
  and a water depth relation :math:`h(t)` are imposed at the boundary


• the 'ratcurve' type for which a tabulated rating curve in the **rating_curve.txt** file is used to
  impose a relation :math:`h(Q)` at the boundary.

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

.. code-block:: python

  os.chdir("/home/livillenave/Documents/distant/dassflow2d-wrap/Tools/2_gen_basic_channel")
  import gen_channel_case
  gen_channel_case.gen_bc(
                          in_type	 = "hydrograph",
                          out_type 	 = "hpresc")


**you can replace in your bin directory the file hydrograh.txt produced**


_____________________________________________________
Access and Update fortran kernel value in python
_____________________________________________________

Should be done via the bc.txt file and not updated within python script for safety.
