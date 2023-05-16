"""
Module m_tap_vars


Defined at m_tap_vars.f90 lines 1-147

"""
from __future__ import print_function, absolute_import, division
from dassflow2d.wrapping import _wrapping
import f90wrap.runtime
import logging
from dassflow2d.wrapping.m_model import infiltration_data
from dassflow2d.wrapping.m_model import bcs

_arrays = {}
_objs = {}

def alloc_back_vars(self, dof_back, mesh):
    """
    alloc_back_vars(self, dof_back, mesh)
    
    
    Defined at m_tap_vars.f90 lines 59-127
    
    Parameters
    ----------
    dof0_back : Unk
    dof_back : Unk
    mesh : Msh
    
    ===============================================
     Classical hydraulic BCs
    ===============================================
    """
    _wrapping.f90wrap_alloc_back_vars(dof0_back=self._handle, \
        dof_back=dof_back._handle, mesh=mesh._handle)

def dealloc_back_vars():
    """
    dealloc_back_vars()
    
    
    Defined at m_tap_vars.f90 lines 129-147
    
    
    _COMMENTtype(unk) , intent(inout) :: dof0_back, dof_back
    _COMMENTif(allocated(dof0_back )) deallocate(dof0_back)
    _COMMENTif(allocated(dof_back )) deallocate(dof_back)
    """
    _wrapping.f90wrap_dealloc_back_vars()

def get_dt_diff():
    """
    Element dt_diff ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 5
    
    """
    return _wrapping.f90wrap_m_tap_vars__get__dt_diff()

def set_dt_diff(dt_diff):
    _wrapping.f90wrap_m_tap_vars__set__dt_diff(dt_diff)

def get_tc_diff():
    """
    Element tc_diff ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 5
    
    """
    return _wrapping.f90wrap_m_tap_vars__get__tc_diff()

def set_tc_diff(tc_diff):
    _wrapping.f90wrap_m_tap_vars__set__tc_diff(tc_diff)

def get_dt_back():
    """
    Element dt_back ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 6
    
    """
    return _wrapping.f90wrap_m_tap_vars__get__dt_back()

def set_dt_back(dt_back):
    _wrapping.f90wrap_m_tap_vars__set__dt_back(dt_back)

def get_tc_back():
    """
    Element tc_back ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 6
    
    """
    return _wrapping.f90wrap_m_tap_vars__get__tc_back()

def set_tc_back(tc_back):
    _wrapping.f90wrap_m_tap_vars__set__tc_back(tc_back)

def get_bc_diff():
    """
    Element bc_diff ftype=type(bcs) pytype=Bcs
    
    
    Defined at m_tap_vars.f90 line 7
    
    """
    global bc_diff
    bc_diff_handle = _wrapping.f90wrap_m_tap_vars__get__bc_diff()
    if tuple(bc_diff_handle) in _objs:
        bc_diff = _objs[tuple(bc_diff_handle)]
    else:
        bc_diff = bcs.from_handle(bc_diff_handle)
        _objs[tuple(bc_diff_handle)] = bc_diff
    return bc_diff

def set_bc_diff(bc_diff):
    bc_diff = bc_diff._handle
    _wrapping.f90wrap_m_tap_vars__set__bc_diff(bc_diff)

def get_bc_back():
    """
    Element bc_back ftype=type(bcs) pytype=Bcs
    
    
    Defined at m_tap_vars.f90 line 8
    
    """
    global bc_back
    bc_back_handle = _wrapping.f90wrap_m_tap_vars__get__bc_back()
    if tuple(bc_back_handle) in _objs:
        bc_back = _objs[tuple(bc_back_handle)]
    else:
        bc_back = bcs.from_handle(bc_back_handle)
        _objs[tuple(bc_back_handle)] = bc_back
    return bc_back

def set_bc_back(bc_back):
    bc_back = bc_back._handle
    _wrapping.f90wrap_m_tap_vars__set__bc_back(bc_back)

def get_infil_diff():
    """
    Element infil_diff ftype=type(infiltration_data) pytype=Infiltration_Data
    
    
    Defined at m_tap_vars.f90 line 9
    
    """
    global infil_diff
    infil_diff_handle = _wrapping.f90wrap_m_tap_vars__get__infil_diff()
    if tuple(infil_diff_handle) in _objs:
        infil_diff = _objs[tuple(infil_diff_handle)]
    else:
        infil_diff = infiltration_data.from_handle(infil_diff_handle)
        _objs[tuple(infil_diff_handle)] = infil_diff
    return infil_diff

def set_infil_diff(infil_diff):
    infil_diff = infil_diff._handle
    _wrapping.f90wrap_m_tap_vars__set__infil_diff(infil_diff)

def get_infil_back():
    """
    Element infil_back ftype=type(infiltration_data) pytype=Infiltration_Data
    
    
    Defined at m_tap_vars.f90 line 10
    
    """
    global infil_back
    infil_back_handle = _wrapping.f90wrap_m_tap_vars__get__infil_back()
    if tuple(infil_back_handle) in _objs:
        infil_back = _objs[tuple(infil_back_handle)]
    else:
        infil_back = infiltration_data.from_handle(infil_back_handle)
        _objs[tuple(infil_back_handle)] = infil_back
    return infil_back

def set_infil_back(infil_back):
    infil_back = infil_back._handle
    _wrapping.f90wrap_m_tap_vars__set__infil_back(infil_back)

def get_array_manning_diff():
    """
    Element manning_diff ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 13
    
    """
    global manning_diff
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_tap_vars__array__manning_diff(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        manning_diff = _arrays[array_handle]
    else:
        manning_diff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_tap_vars__array__manning_diff)
        _arrays[array_handle] = manning_diff
    return manning_diff

def set_array_manning_diff(manning_diff):
    manning_diff[...] = manning_diff

def get_array_manning_beta_diff():
    """
    Element manning_beta_diff ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 13
    
    """
    global manning_beta_diff
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_tap_vars__array__manning_beta_diff(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        manning_beta_diff = _arrays[array_handle]
    else:
        manning_beta_diff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_tap_vars__array__manning_beta_diff)
        _arrays[array_handle] = manning_beta_diff
    return manning_beta_diff

def set_array_manning_beta_diff(manning_beta_diff):
    manning_beta_diff[...] = manning_beta_diff

def get_array_bathy_cell_diff():
    """
    Element bathy_cell_diff ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 13
    
    """
    global bathy_cell_diff
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_tap_vars__array__bathy_cell_diff(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        bathy_cell_diff = _arrays[array_handle]
    else:
        bathy_cell_diff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_tap_vars__array__bathy_cell_diff)
        _arrays[array_handle] = bathy_cell_diff
    return bathy_cell_diff

def set_array_bathy_cell_diff(bathy_cell_diff):
    bathy_cell_diff[...] = bathy_cell_diff

def get_array_manning_back():
    """
    Element manning_back ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 14
    
    """
    global manning_back
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_tap_vars__array__manning_back(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        manning_back = _arrays[array_handle]
    else:
        manning_back = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_tap_vars__array__manning_back)
        _arrays[array_handle] = manning_back
    return manning_back

def set_array_manning_back(manning_back):
    manning_back[...] = manning_back

def get_array_manning_beta_back():
    """
    Element manning_beta_back ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 14
    
    """
    global manning_beta_back
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_tap_vars__array__manning_beta_back(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        manning_beta_back = _arrays[array_handle]
    else:
        manning_beta_back = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_tap_vars__array__manning_beta_back)
        _arrays[array_handle] = manning_beta_back
    return manning_beta_back

def set_array_manning_beta_back(manning_beta_back):
    manning_beta_back[...] = manning_beta_back

def get_array_bathy_cell_back():
    """
    Element bathy_cell_back ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 14
    
    """
    global bathy_cell_back
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_tap_vars__array__bathy_cell_back(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        bathy_cell_back = _arrays[array_handle]
    else:
        bathy_cell_back = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_tap_vars__array__bathy_cell_back)
        _arrays[array_handle] = bathy_cell_back
    return bathy_cell_back

def set_array_bathy_cell_back(bathy_cell_back):
    bathy_cell_back[...] = bathy_cell_back


_array_initialisers = [get_array_manning_diff, get_array_manning_beta_diff, \
    get_array_bathy_cell_diff, get_array_manning_back, \
    get_array_manning_beta_back, get_array_bathy_cell_back]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "m_tap_vars".')

for func in _dt_array_initialisers:
    func()
