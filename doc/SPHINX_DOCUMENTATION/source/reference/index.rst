.. _api:

===============
API Reference
===============

This page gives an overview of all dassflow2d objects, functions and methods.




===============
Quick Overview
===============

Dassflow2d can be seen as a package dedicated to **dassflowmodel** object.
Multiple data are made availaible in    **dassflowmodel**.
Every class (config, meshing, param boundary, output, min, obs), is built with
the same logic. We associate the following kind of functions

- get() : load fortran kernel values and fills  **dassflowmodel** by a copy
- set() : set the fortran kernel values from  **dassflowmodel**  values
- save() : save   **dassflowmodel** or kernel values to hdf5 file the data proper to the class
- load() : load from hdf5_file and fills **  **dassflowmodel** object
- source_xxx() : load data from some txt file
- plot() : represent the class, often with arguments

.. toctree::
   :maxdepth: 2

   model
   Config
   Meshing
   Boundary
   Param
   Output
   wrapping
