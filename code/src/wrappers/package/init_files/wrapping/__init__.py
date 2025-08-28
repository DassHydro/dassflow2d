"""
dassflow2d.wrapping  MODULE

enable the fortran subroutines calls

#the calls are programed in the files m_mesh.py, m_common.py, etc...
#the corresponding fortran routines are stored in _wraping.*.so file
"""


# NOUVELLE VERSION CORRIGÉE

# Importer les modules de base qui sont toujours présents
from dassflow2d.wrapping import _wrapping, m_mesh, m_common, m_linear_algebra, m_model, m_mpi, call_model

# Essayer d'importer les modules adjoints, mais ne pas planter s'ils sont absents
try:
    from dassflow2d.wrapping import m_adjoint
except ImportError:
    pass # On ignore l'erreur si m_adjoint n'a pas été compilé

try:
    from dassflow2d.wrapping import m_tap_vars
except ImportError:
    pass # On ignore aussi pour m_tap_vars


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

def print_all():
    """
    print_all()
    
    
    Defined at input.f90 lines 22-28
    
    
    """
    _wrapping.f90wrap_print_all()

def read_bc_file():
    """
    read_bc_file()
    
    
    Defined at input.f90 lines 582-633
    
    
    ===================================================================================================================
     Local Variables
    ===================================================================================================================
    """
    _wrapping.f90wrap_read_bc_file()
