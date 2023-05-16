.. _api.wrapping:

.. module:: wrapping


=======================
Fortran wrapped modules
=======================

.. currentmodule:: dassflow2d.wrapping

===============
API Reference
===============

This page gives an overview of all public smash objects, functions and methods.

*****************
Derived Type
*****************

Attributes
**********
.. autosummary::
   :toctree: dassflow2d.wrapping/

   m_mesh
   m_common
   m_model
   m_adjoint
   m_mpi
   call_model










*********************
Main called routines
*********************

------------------
m_mesh
------------------


.. autosummary::
   :toctree: dassflow2d/
   
   m_model.unk
   m_model.unk.from_handle
   m_model.unk.init_array_grad_h
   m_model.unk.init_array_grad_u
   m_model.unk.init_array_grad_v
   m_model.unk.init_array_grad_z
   m_model.unk.t_display
