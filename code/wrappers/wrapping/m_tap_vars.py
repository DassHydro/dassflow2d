"""
Module m_tap_vars


Defined at m_tap_vars.f90 lines 1-226

"""
from __future__ import print_function, absolute_import, division
import _wrapping
import f90wrap.runtime
import logging
import numpy
from wrapping.m_model import ptf_data
from wrapping.m_model import xsshp
from wrapping.m_model import infiltration_data
from wrapping.m_model import bcs

_arrays = {}
_objs = {}

def alloc_back_vars(self, dof_back, mesh):
    """
    alloc_back_vars(self, dof_back, mesh)
    
    
    Defined at m_tap_vars.f90 lines 97-198
    
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

def dealloc_back_vars(self, dof_back):
    """
    dealloc_back_vars(self, dof_back)
    
    
    Defined at m_tap_vars.f90 lines 200-226
    
    Parameters
    ----------
    dof0_back : Unk
    dof_back : Unk
    
    """
    _wrapping.f90wrap_dealloc_back_vars(dof0_back=self._handle, \
        dof_back=dof_back._handle)

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

def init_array_xsshape_back():
    global xsshape_back
    xsshape_back = \
        f90wrap.runtime.FortranDerivedTypeArray(f90wrap.runtime.empty_type,
                                    _wrapping.f90wrap_m_tap_vars__array_getitem__xsshape_back,
                                    _wrapping.f90wrap_m_tap_vars__array_setitem__xsshape_back,
                                    _wrapping.f90wrap_m_tap_vars__array_len__xsshape_back,
                                    """
    Element xsshape_back ftype=type(xsshp) pytype=Xsshp
    
    
    Defined at m_tap_vars.f90 line 9
    
    """, xsshp)
    return xsshape_back

def init_array_xsshape_diff():
    global xsshape_diff
    xsshape_diff = \
        f90wrap.runtime.FortranDerivedTypeArray(f90wrap.runtime.empty_type,
                                    _wrapping.f90wrap_m_tap_vars__array_getitem__xsshape_diff,
                                    _wrapping.f90wrap_m_tap_vars__array_setitem__xsshape_diff,
                                    _wrapping.f90wrap_m_tap_vars__array_len__xsshape_diff,
                                    """
    Element xsshape_diff ftype=type(xsshp) pytype=Xsshp
    
    
    Defined at m_tap_vars.f90 line 10
    
    """, xsshp)
    return xsshape_diff

def get_infil_diff():
    """
    Element infil_diff ftype=type(infiltration_data) pytype=Infiltration_Data
    
    
    Defined at m_tap_vars.f90 line 11
    
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
    
    
    Defined at m_tap_vars.f90 line 12
    
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

def init_array_ptf_diff():
    global ptf_diff
    ptf_diff = f90wrap.runtime.FortranDerivedTypeArray(f90wrap.runtime.empty_type,
                                    _wrapping.f90wrap_m_tap_vars__array_getitem__ptf_diff,
                                    _wrapping.f90wrap_m_tap_vars__array_setitem__ptf_diff,
                                    _wrapping.f90wrap_m_tap_vars__array_len__ptf_diff,
                                    """
    Element ptf_diff ftype=type(ptf_data) pytype=Ptf_Data
    
    
    Defined at m_tap_vars.f90 line 13
    
    """, ptf_data)
    return ptf_diff

def init_array_ptf_back():
    global ptf_back
    ptf_back = f90wrap.runtime.FortranDerivedTypeArray(f90wrap.runtime.empty_type,
                                    _wrapping.f90wrap_m_tap_vars__array_getitem__ptf_back,
                                    _wrapping.f90wrap_m_tap_vars__array_setitem__ptf_back,
                                    _wrapping.f90wrap_m_tap_vars__array_len__ptf_back,
                                    """
    Element ptf_back ftype=type(ptf_data) pytype=Ptf_Data
    
    
    Defined at m_tap_vars.f90 line 14
    
    """, ptf_data)
    return ptf_back

def get_array_manning_diff():
    """
    Element manning_diff ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 23
    
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
    globals()['manning_diff'][...] = manning_diff

def get_array_manning_beta_diff():
    """
    Element manning_beta_diff ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 23
    
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
    globals()['manning_beta_diff'][...] = manning_beta_diff

def get_array_bathy_cell_diff():
    """
    Element bathy_cell_diff ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 23
    
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
    globals()['bathy_cell_diff'][...] = bathy_cell_diff

def get_array_manning_back():
    """
    Element manning_back ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 24
    
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
    globals()['manning_back'][...] = manning_back

def get_array_manning_beta_back():
    """
    Element manning_beta_back ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 24
    
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
    globals()['manning_beta_back'][...] = manning_beta_back

def get_array_bathy_cell_back():
    """
    Element bathy_cell_back ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 24
    
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
    globals()['bathy_cell_back'][...] = bathy_cell_back

def get_array_slope_y_back():
    """
    Element slope_y_back ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 25
    
    """
    global slope_y_back
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_tap_vars__array__slope_y_back(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        slope_y_back = _arrays[array_handle]
    else:
        slope_y_back = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_tap_vars__array__slope_y_back)
        _arrays[array_handle] = slope_y_back
    return slope_y_back

def set_array_slope_y_back(slope_y_back):
    globals()['slope_y_back'][...] = slope_y_back

def get_array_slope_y_diff():
    """
    Element slope_y_diff ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 25
    
    """
    global slope_y_diff
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_tap_vars__array__slope_y_diff(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        slope_y_diff = _arrays[array_handle]
    else:
        slope_y_diff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_tap_vars__array__slope_y_diff)
        _arrays[array_handle] = slope_y_diff
    return slope_y_diff

def set_array_slope_y_diff(slope_y_diff):
    globals()['slope_y_diff'][...] = slope_y_diff

def get_array_slope_x_back():
    """
    Element slope_x_back ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 26
    
    """
    global slope_x_back
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_tap_vars__array__slope_x_back(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        slope_x_back = _arrays[array_handle]
    else:
        slope_x_back = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_tap_vars__array__slope_x_back)
        _arrays[array_handle] = slope_x_back
    return slope_x_back

def set_array_slope_x_back(slope_x_back):
    globals()['slope_x_back'][...] = slope_x_back

def get_array_slope_x_diff():
    """
    Element slope_x_diff ftype=real(rp) pytype=float
    
    
    Defined at m_tap_vars.f90 line 26
    
    """
    global slope_x_diff
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_tap_vars__array__slope_x_diff(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        slope_x_diff = _arrays[array_handle]
    else:
        slope_x_diff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_tap_vars__array__slope_x_diff)
        _arrays[array_handle] = slope_x_diff
    return slope_x_diff

def set_array_slope_x_diff(slope_x_diff):
    globals()['slope_x_diff'][...] = slope_x_diff


_array_initialisers = [get_array_manning_diff, get_array_manning_beta_diff, \
    get_array_bathy_cell_diff, get_array_manning_back, \
    get_array_manning_beta_back, get_array_bathy_cell_back, \
    get_array_slope_y_back, get_array_slope_y_diff, get_array_slope_x_back, \
    get_array_slope_x_diff]
_dt_array_initialisers = [init_array_xsshape_back, init_array_xsshape_diff, \
    init_array_ptf_diff, init_array_ptf_back]

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "m_tap_vars".')

for func in _dt_array_initialisers:
    func()
