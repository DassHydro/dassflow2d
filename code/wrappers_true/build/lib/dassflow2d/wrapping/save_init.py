from __future__ import print_function, absolute_import, division
import _wrapping
import f90wrap.runtime
import logging
import numpy
import wrapping.m_mesh
import wrapping.m_common
import wrapping.m_model
import wrapping.m_linear_algebra
import wrapping.m_mpi

def read_input(filename):
    """
    read_input(filename)
    
    
    Defined at input.f90 lines 1-24
    
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
    
    
    Defined at input.f90 lines 26-32
    
    
    _COMMENT afficher toutes les valeurs lues
    """
    _wrapping.f90wrap_print_all()

def read_bc_file():
    """
    line_read = read_bc_file()
    
    
    Defined at input.f90 lines 682-761
    
    
    Returns
    -------
    line_read : int
    
    ===================================================================================================================
     Local Variables
    ===================================================================================================================
    """
    line_read = _wrapping.f90wrap_read_bc_file()
    return line_read

