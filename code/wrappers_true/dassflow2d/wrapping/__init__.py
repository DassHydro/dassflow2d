"""
dassflow2d.wrapping  MODULE

enable the fortran subroutines calls

#the calls are programed in the files m_mesh.py, m_common.py, etc...
#the corresponding fortran routines are stored in _wraping.*.so file
"""


from dassflow2d.wrapping import _wrapping, m_mesh, m_common, m_linear_algebra, m_model, call_model, m_adjoint

def read_input(filename):
    """
    read_input(filename)
    
    Initialise "config" variables using "filename" file and default values in fortran kernel
    
    
    Fortran detail
    -----------------------
    Defined at input.f90 lines 1-20
    
    Parameters
    ----------
    filename : str
    """
    _wrapping.f90wrap_read_input(filename=filename)

