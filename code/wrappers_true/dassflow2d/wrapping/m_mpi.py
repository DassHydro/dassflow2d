"""
Module m_mpi


Defined at m_mpi.f90 lines 1-416

Details: This type is used for the definition of vector of size 2( Coordinates \
    \f$(x,y)\f$ for example).
Details: This type is used for the definition of vector of size 3( Coordinates \
    \f$(x,y,z)\f$ for example).
Details: This type is used for the definition of Tensor of size 2*2( ... for \
    example).
Details: This type is used for the definition of Tensor of size 3*3( ... for \
    example).
"""
from __future__ import print_function, absolute_import, division
from dassflow2d.wrapping import _wrapping
import f90wrap.runtime
import logging
import numpy

_arrays = {}
_objs = {}

def init_mpi():
    """
    init_mpi()
    
    
    Defined at m_mpi.f90 lines 17-20
    
    
    """
    _wrapping.f90wrap_init_mpi()

def end_mpi():
    """
    end_mpi()
    
    
    Defined at m_mpi.f90 lines 22-25
    
    
    """
    _wrapping.f90wrap_end_mpi()

def fill_swap_lists(self):
    """
    fill_swap_lists(self)
    
    
    Defined at m_mpi.f90 lines 152-177
    
    Parameters
    ----------
    mesh : Msh
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_fill_swap_lists(mesh=self._handle)

def fill_swap_index(self):
    """
    fill_swap_index(self)
    
    
    Defined at m_mpi.f90 lines 179-196
    
    Parameters
    ----------
    mesh : Msh
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_fill_swap_index(mesh=self._handle)

def fill_inv_swap_index(self):
    """
    fill_inv_swap_index(self)
    
    
    Defined at m_mpi.f90 lines 198-214
    
    Parameters
    ----------
    mesh : Msh
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_fill_inv_swap_index(mesh=self._handle)

def com_var_i(var, mesh):
    """
    com_var_i(var, mesh)
    
    
    Defined at m_mpi.f90 lines 216-225
    
    Parameters
    ----------
    var : int array
    mesh : Msh
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_com_var_i(var=var, mesh=mesh._handle)

def com_var_r(var, mesh):
    """
    com_var_r(var, mesh)
    
    
    Defined at m_mpi.f90 lines 227-236
    
    Parameters
    ----------
    var : float array
    mesh : Msh
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_com_var_r(var=var, mesh=mesh._handle)

def mpi_send_recv_scal_i(to_send, proc_send, proc_recv):
    """
    to_recv = mpi_send_recv_scal_i(to_send, proc_send, proc_recv)
    
    
    Defined at m_mpi.f90 lines 238-249
    
    Parameters
    ----------
    to_send : int
    proc_send : int
    proc_recv : int
    
    Returns
    -------
    to_recv : int
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    to_recv = _wrapping.f90wrap_mpi_send_recv_scal_i(to_send=to_send, \
        proc_send=proc_send, proc_recv=proc_recv)
    return to_recv

def mpi_send_recv_scal_r(to_send, proc_send, proc_recv):
    """
    to_recv = mpi_send_recv_scal_r(to_send, proc_send, proc_recv)
    
    
    Defined at m_mpi.f90 lines 251-262
    
    Parameters
    ----------
    to_send : float
    proc_send : int
    proc_recv : int
    
    Returns
    -------
    to_recv : float
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    to_recv = _wrapping.f90wrap_mpi_send_recv_scal_r(to_send=to_send, \
        proc_send=proc_send, proc_recv=proc_recv)
    return to_recv

def mpi_send_recv_array_i(to_send, to_recv, proc_send, proc_recv):
    """
    mpi_send_recv_array_i(to_send, to_recv, proc_send, proc_recv)
    
    
    Defined at m_mpi.f90 lines 264-275
    
    Parameters
    ----------
    to_send : int array
    to_recv : int array
    proc_send : int
    proc_recv : int
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_mpi_send_recv_array_i(to_send=to_send, to_recv=to_recv, \
        proc_send=proc_send, proc_recv=proc_recv)

def mpi_send_recv_array_r(to_send, to_recv, proc_send, proc_recv):
    """
    mpi_send_recv_array_r(to_send, to_recv, proc_send, proc_recv)
    
    
    Defined at m_mpi.f90 lines 277-288
    
    Parameters
    ----------
    to_send : float array
    to_recv : float array
    proc_send : int
    proc_recv : int
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_mpi_send_recv_array_r(to_send=to_send, to_recv=to_recv, \
        proc_send=proc_send, proc_recv=proc_recv)

def mpi_sum_r(val):
    """
    mpi_sum_r(val)
    
    
    Defined at m_mpi.f90 lines 290-299
    
    Parameters
    ----------
    val : float
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_mpi_sum_r(val=val)

def mpi_sum_i(val):
    """
    mpi_sum_i(val)
    
    
    Defined at m_mpi.f90 lines 301-309
    
    Parameters
    ----------
    val : int
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_mpi_sum_i(val=val)

def mpi_max_r(val):
    """
    mpi_max_r(val)
    
    
    Defined at m_mpi.f90 lines 311-319
    
    Parameters
    ----------
    val : float
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_mpi_max_r(val=val)

def mpi_max_i(val):
    """
    mpi_max_i(val)
    
    
    Defined at m_mpi.f90 lines 321-329
    
    Parameters
    ----------
    val : int
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_mpi_max_i(val=val)

def mpi_min_r(val):
    """
    mpi_min_r(val)
    
    
    Defined at m_mpi.f90 lines 331-339
    
    Parameters
    ----------
    val : float
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_mpi_min_r(val=val)

def mpi_min_i(val):
    """
    mpi_min_i(val)
    
    
    Defined at m_mpi.f90 lines 341-349
    
    Parameters
    ----------
    val : int
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_mpi_min_i(val=val)

def mpi_bcast_r(val, pr):
    """
    mpi_bcast_r(val, pr)
    
    
    Defined at m_mpi.f90 lines 351-360
    
    Parameters
    ----------
    val : float
    pr : int
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_mpi_bcast_r(val=val, pr=pr)

def mpi_bcast_i(val, pr):
    """
    mpi_bcast_i(val, pr)
    
    
    Defined at m_mpi.f90 lines 362-371
    
    Parameters
    ----------
    val : int
    pr : int
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_mpi_bcast_i(val=val, pr=pr)

def mpi_allgather_r(val, temp):
    """
    mpi_allgather_r(val, temp)
    
    
    Defined at m_mpi.f90 lines 373-383
    
    Parameters
    ----------
    val : float
    temp : float array
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_mpi_allgather_r(val=val, temp=temp)

def mpi_allgather_i(val, temp):
    """
    mpi_allgather_i(val, temp)
    
    
    Defined at m_mpi.f90 lines 385-395
    
    Parameters
    ----------
    val : int
    temp : int array
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_mpi_allgather_i(val=val, temp=temp)

def mpi_wait_all():
    """
    mpi_wait_all()
    
    
    Defined at m_mpi.f90 lines 397-398
    
    
    """
    _wrapping.f90wrap_mpi_wait_all()

def stopping_program_sub(comment):
    """
    stopping_program_sub(comment)
    
    
    Defined at m_mpi.f90 lines 400-416
    
    Parameters
    ----------
    comment : str
    
    """
    _wrapping.f90wrap_stopping_program_sub(comment=comment)

def get_np():
    """
    Element np ftype=integer(ip) pytype=int
    
    
    Defined at m_mpi.f90 line 6
    
    """
    return _wrapping.f90wrap_m_mpi__get__np()

def set_np(np):
    _wrapping.f90wrap_m_mpi__set__np(np)

def get_proc():
    """
    Element proc ftype=integer(ip) pytype=int
    
    
    Defined at m_mpi.f90 line 7
    
    """
    return _wrapping.f90wrap_m_mpi__get__proc()

def set_proc(proc):
    _wrapping.f90wrap_m_mpi__set__proc(proc)

def get_code():
    """
    Element code ftype=integer(ip) pytype=int
    
    
    Defined at m_mpi.f90 line 8
    
    """
    return _wrapping.f90wrap_m_mpi__get__code()

def set_code(code):
    _wrapping.f90wrap_m_mpi__set__code(code)

def get_nneighb():
    """
    Element nneighb ftype=integer(ip) pytype=int
    
    
    Defined at m_mpi.f90 line 9
    
    """
    return _wrapping.f90wrap_m_mpi__get__nneighb()

def set_nneighb(nneighb):
    _wrapping.f90wrap_m_mpi__set__nneighb(nneighb)

def get_array_swap_index():
    """
    Element swap_index ftype=integer(ip) pytype=int
    
    
    Defined at m_mpi.f90 line 10
    
    """
    global swap_index
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_mpi__array__swap_index(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        swap_index = _arrays[array_handle]
    else:
        swap_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_mpi__array__swap_index)
        _arrays[array_handle] = swap_index
    return swap_index

def set_array_swap_index(swap_index):
    globals()['swap_index'][...] = swap_index

def get_array_inv_swap_index():
    """
    Element inv_swap_index ftype=integer(ip) pytype=int
    
    
    Defined at m_mpi.f90 line 10
    
    """
    global inv_swap_index
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_mpi__array__inv_swap_index(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        inv_swap_index = _arrays[array_handle]
    else:
        inv_swap_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_mpi__array__inv_swap_index)
        _arrays[array_handle] = inv_swap_index
    return inv_swap_index

def set_array_inv_swap_index(inv_swap_index):
    globals()['inv_swap_index'][...] = inv_swap_index

def get_array_part():
    """
    Element part ftype=integer(ip) pytype=int
    
    
    Defined at m_mpi.f90 line 11
    
    """
    global part
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_mpi__array__part(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        part = _arrays[array_handle]
    else:
        part = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_mpi__array__part)
        _arrays[array_handle] = part
    return part

def set_array_part(part):
    globals()['part'][...] = part

def get_array_part_size():
    """
    Element part_size ftype=integer(ip) pytype=int
    
    
    Defined at m_mpi.f90 line 12
    
    """
    global part_size
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_mpi__array__part_size(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        part_size = _arrays[array_handle]
    else:
        part_size = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_mpi__array__part_size)
        _arrays[array_handle] = part_size
    return part_size

def set_array_part_size(part_size):
    globals()['part_size'][...] = part_size

def get_array_part_neighbs():
    """
    Element part_neighbs ftype=integer(ip) pytype=int
    
    
    Defined at m_mpi.f90 line 13
    
    """
    global part_neighbs
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_mpi__array__part_neighbs(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        part_neighbs = _arrays[array_handle]
    else:
        part_neighbs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_mpi__array__part_neighbs)
        _arrays[array_handle] = part_neighbs
    return part_neighbs

def set_array_part_neighbs(part_neighbs):
    globals()['part_neighbs'][...] = part_neighbs

def get_val_tmp_r():
    """
    Element val_tmp_r ftype=real(rp) pytype=float
    
    
    Defined at m_mpi.f90 line 14
    
    """
    return _wrapping.f90wrap_m_mpi__get__val_tmp_r()

def set_val_tmp_r(val_tmp_r):
    _wrapping.f90wrap_m_mpi__set__val_tmp_r(val_tmp_r)

def get_val_tmp_i():
    """
    Element val_tmp_i ftype=integer(ip) pytype=int
    
    
    Defined at m_mpi.f90 line 15
    
    """
    return _wrapping.f90wrap_m_mpi__get__val_tmp_i()

def set_val_tmp_i(val_tmp_i):
    _wrapping.f90wrap_m_mpi__set__val_tmp_i(val_tmp_i)


_array_initialisers = [get_array_swap_index, get_array_inv_swap_index, \
    get_array_part, get_array_part_size, get_array_part_neighbs]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "m_mpi".')

for func in _dt_array_initialisers:
    func()
