from __future__ import print_function, absolute_import, division
import _wrapping
import f90wrap.runtime
import logging
import wrapping.call_model
import wrapping.m_adjoint
import wrapping.m_tap_vars
import wrapping.m_common
import wrapping.m_linear_algebra
import wrapping.m_model
import wrapping.m_mesh

def read_input(filename):
    """
    read_input(filename)
    
    
    Defined at input.f90 lines 1-18
    
    Parameters
    ----------
    filename : str
    
    ===================================================================================================================
    ===================================================================================================================
    """
    _wrapping.f90wrap_read_input(filename=filename)

def print_all():
    """
    print_all()
    
    
    Defined at input.f90 lines 20-26
    
    
    _COMMENT afficher toutes les valeurs lues
    """
    _wrapping.f90wrap_print_all()

def read_bc_file():
    """
    read_bc_file()
    
    
    Defined at input.f90 lines 603-658
    
    
    ===================================================================================================================
     Local Variables
    ===================================================================================================================
    """
    _wrapping.f90wrap_read_bc_file()

