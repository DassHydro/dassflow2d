.. _api.model:

===============
dassflowmodel
===============
.. currentmodule:: dassflow2d

**dassflowmodel** is the main class of dassflow2d packagage. It enable interfacing with fortran kernel and python librairies.

**dassflowmodel** also contains its proper data independently from fortran kernel.


++++++++++++
subclasses
++++++++++++

dassflowmodel contains subclasses:

* to manipulate the mesh:
  :class:`dassflow2d.core.meshing.Meshing`
* to manipulate the configuration file:
  :class:`dassflow2d.core.config.Config`
* to manipulate fortran simulation results:
  :class:`dassflow2d.core.output.Output`
* to manipulate boundary conditions:
  :class:`dassflow2d.core.boundary.Boundary`
* to manipulate parameters:
  :class:`dassflow2d.core.param.Param`



  ++++++++++++
  subclasses
  ++++++++++++

Attributes
**********
.. autosummary::
   :toctree: dassflow2d/

   dassflowmodel.__init__
   dassflowmodel.init_mesh
   dassflowmodel.init_dof
   dassflowmodel.init_fortran
   dassflowmodel.init_all
   dassflowmodel.run
   dassflowmodel.save_all
