.. _api_meshing:

===============
Meshing
===============
.. currentmodule:: dassflow2d.core.meshing

In Dassflow2d wrap, the meshing is stored two ways:

-  in an object: **meshing.mesh_fortran**, it corresponds to the fortran kernel mesh object.
-  in an object: **meshing.mesh_pyvista**, it corresponds toa pyvista unstructured grid object , containing mesh geometry.

++++++++++++
class
++++++++++++

.. autoclass:: Meshing


++++++++++++
Methods
++++++++++++

The following methods can be applied on objects of class ``Meshing``



  Attributes
  **********
  .. autosummary::
     :toctree: dassflow2d.core.meshing/

     Meshing.build_grid
     Meshing.plot
     Meshing.save

.. hint ::

  heavy plot method:

  .. autosummary::
    :toctree: dassflow2d.core.meshing/

      Meshing.plot_dev
