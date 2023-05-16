.. _api_boundary:

===============
Boundary
===============
.. currentmodule:: dassflow2d.core.boundary


In Dassflow2d wrap, the boundaries information are stored in **Boundary** class, it contains two main information:

-  in an object: **bc**, it corresponds to table of values of boundaries.
-  in an object: **corresp**, which gives the correspondance between meshing and table and values.


.. hint ::

  boundary inherit from **meshing.fortran_mesh**



Attributes
**********
.. autosummary::
   :toctree: dassflow2d.core.boundary/

   Boundary.__init__
   Boundary.get_mesh_corresp
   Boundary.get_metadata
   Boundary.plot
   Boundary.save
