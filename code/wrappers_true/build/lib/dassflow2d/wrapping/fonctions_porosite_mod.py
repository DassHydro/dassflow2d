"""
Module fonctions_porosite_mod


Defined at functions_porosity.f90 lines 1-137

"""
from __future__ import print_function, absolute_import, division
import _wrapping
import f90wrap.runtime
import logging
import numpy

_arrays = {}
_objs = {}

def update_all_porosities(self, mesh):
    """
    update_all_porosities(self, mesh)
    
    
    Defined at functions_porosity.f90 lines 93-137
    
    Parameters
    ----------
    dof : Unk
    mesh : Msh
    
    =======================================================================
     Orchestre la mise \U000000e0 jour de la porosit\U000000e9 pour toutes les \
         cellules 1D-like.
    =======================================================================
     --- Arguments ---
    _COMMENT --- Variables Locales ---
    """
    _wrapping.f90wrap_fonctions_porosite_mod__update_all_porosities(dof=self._handle, \
        mesh=mesh._handle)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "fonctions_porosite_mod".')

for func in _dt_array_initialisers:
    func()
