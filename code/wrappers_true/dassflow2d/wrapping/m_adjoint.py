"""
Module m_adjoint


Defined at m_adjoint.f90 lines 1-945

"""
from __future__ import print_function, absolute_import, division
from dassflow2d.wrapping import _wrapping
import f90wrap.runtime
import logging
from dassflow2d.wrapping.m_model import unk

_arrays = {}
_objs = {}

def model_direct(self, dof0, dof):
    """
    model_direct(self, dof0, dof)
    
    
    Defined at m_adjoint.f90 lines 33-48
    
    Parameters
    ----------
    mesh : Msh
    dof0 : Unk
    dof : Unk
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_model_direct(mesh=self._handle, dof0=dof0._handle, \
        dof=dof._handle)

def model_direct_perturb(self, dof0, dof):
    """
    model_direct_perturb(self, dof0, dof)
    
    
    Defined at m_adjoint.f90 lines 50-67
    
    Parameters
    ----------
    mesh : Msh
    dof0 : Unk
    dof : Unk
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_model_direct_perturb(mesh=self._handle, dof0=dof0._handle, \
        dof=dof._handle)

def adjoint_model(self, dof0, dof):
    """
    adjoint_model(self, dof0, dof)
    
    
    Defined at m_adjoint.f90 lines 92-133
    
    Parameters
    ----------
    mesh : Msh
    dof0 : Unk
    dof : Unk
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_adjoint_model(mesh=self._handle, dof0=dof0._handle, \
        dof=dof._handle)

def write_control(self, mesh):
    """
    write_control(self, mesh)
    
    
    Defined at m_adjoint.f90 lines 288-366
    
    Parameters
    ----------
    dof0 : Unk
    mesh : Msh
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_write_control(dof0=self._handle, mesh=mesh._handle)

def write_control_back(self, mesh):
    """
    write_control_back(self, mesh)
    
    
    Defined at m_adjoint.f90 lines 441-507
    
    Parameters
    ----------
    dof0 : Unk
    mesh : Msh
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_write_control_back(dof0=self._handle, mesh=mesh._handle)

def read_control(self, mesh):
    """
    read_control(self, mesh)
    
    
    Defined at m_adjoint.f90 lines 577-651
    
    Parameters
    ----------
    dof0 : Unk
    mesh : Msh
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_read_control(dof0=self._handle, mesh=mesh._handle)

def output_control(self, mesh):
    """
    output_control(self, mesh)
    
    
    Defined at m_adjoint.f90 lines 811-901
    
    Parameters
    ----------
    dof0 : Unk
    mesh : Msh
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_output_control(dof0=self._handle, mesh=mesh._handle)

def output_control_back(self):
    """
    output_control_back(self)
    
    
    Defined at m_adjoint.f90 lines 903-945
    
    Parameters
    ----------
    mesh : Msh
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    """
    _wrapping.f90wrap_output_control_back(mesh=self._handle)

def get_dof0_diff():
    """
    Element dof0_diff ftype=type(unk) pytype=Unk
    
    
    Defined at m_adjoint.f90 line 11
    
    """
    global dof0_diff
    dof0_diff_handle = _wrapping.f90wrap_m_adjoint__get__dof0_diff()
    if tuple(dof0_diff_handle) in _objs:
        dof0_diff = _objs[tuple(dof0_diff_handle)]
    else:
        dof0_diff = unk.from_handle(dof0_diff_handle)
        _objs[tuple(dof0_diff_handle)] = dof0_diff
    return dof0_diff

def set_dof0_diff(dof0_diff):
    dof0_diff = dof0_diff._handle
    _wrapping.f90wrap_m_adjoint__set__dof0_diff(dof0_diff)

def get_dof0_back():
    """
    Element dof0_back ftype=type(unk) pytype=Unk
    
    
    Defined at m_adjoint.f90 line 12
    
    """
    global dof0_back
    dof0_back_handle = _wrapping.f90wrap_m_adjoint__get__dof0_back()
    if tuple(dof0_back_handle) in _objs:
        dof0_back = _objs[tuple(dof0_back_handle)]
    else:
        dof0_back = unk.from_handle(dof0_back_handle)
        _objs[tuple(dof0_back_handle)] = dof0_back
    return dof0_back

def set_dof0_back(dof0_back):
    dof0_back = dof0_back._handle
    _wrapping.f90wrap_m_adjoint__set__dof0_back(dof0_back)

def get_dof_diff():
    """
    Element dof_diff ftype=type(unk) pytype=Unk
    
    
    Defined at m_adjoint.f90 line 13
    
    """
    global dof_diff
    dof_diff_handle = _wrapping.f90wrap_m_adjoint__get__dof_diff()
    if tuple(dof_diff_handle) in _objs:
        dof_diff = _objs[tuple(dof_diff_handle)]
    else:
        dof_diff = unk.from_handle(dof_diff_handle)
        _objs[tuple(dof_diff_handle)] = dof_diff
    return dof_diff

def set_dof_diff(dof_diff):
    dof_diff = dof_diff._handle
    _wrapping.f90wrap_m_adjoint__set__dof_diff(dof_diff)

def get_dof_back():
    """
    Element dof_back ftype=type(unk) pytype=Unk
    
    
    Defined at m_adjoint.f90 line 14
    
    """
    global dof_back
    dof_back_handle = _wrapping.f90wrap_m_adjoint__get__dof_back()
    if tuple(dof_back_handle) in _objs:
        dof_back = _objs[tuple(dof_back_handle)]
    else:
        dof_back = unk.from_handle(dof_back_handle)
        _objs[tuple(dof_back_handle)] = dof_back
    return dof_back

def set_dof_back(dof_back):
    dof_back = dof_back._handle
    _wrapping.f90wrap_m_adjoint__set__dof_back(dof_back)

def get_array_control():
    """
    Element control ftype=real(rp) pytype=float
    
    
    Defined at m_adjoint.f90 line 15
    
    """
    global control
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_adjoint__array__control(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        control = _arrays[array_handle]
    else:
        control = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_adjoint__array__control)
        _arrays[array_handle] = control
    return control

def set_array_control(control):
    control[...] = control

def get_array_control_back():
    """
    Element control_back ftype=real(rp) pytype=float
    
    
    Defined at m_adjoint.f90 line 16
    
    """
    global control_back
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_adjoint__array__control_back(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        control_back = _arrays[array_handle]
    else:
        control_back = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_adjoint__array__control_back)
        _arrays[array_handle] = control_back
    return control_back

def set_array_control_back(control_back):
    control_back[...] = control_back

def get_array_control_diff():
    """
    Element control_diff ftype=real(rp) pytype=float
    
    
    Defined at m_adjoint.f90 line 17
    
    """
    global control_diff
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_adjoint__array__control_diff(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        control_diff = _arrays[array_handle]
    else:
        control_diff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_adjoint__array__control_diff)
        _arrays[array_handle] = control_diff
    return control_diff

def set_array_control_diff(control_diff):
    control_diff[...] = control_diff

def get_array_control_perturb():
    """
    Element control_perturb ftype=real(rp) pytype=float
    
    
    Defined at m_adjoint.f90 line 18
    
    """
    global control_perturb
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_adjoint__array__control_perturb(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        control_perturb = _arrays[array_handle]
    else:
        control_perturb = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_adjoint__array__control_perturb)
        _arrays[array_handle] = control_perturb
    return control_perturb

def set_array_control_perturb(control_perturb):
    control_perturb[...] = control_perturb

def get_array_control_lbound():
    """
    Element control_lbound ftype=real(rp) pytype=float
    
    
    Defined at m_adjoint.f90 line 19
    
    """
    global control_lbound
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_adjoint__array__control_lbound(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        control_lbound = _arrays[array_handle]
    else:
        control_lbound = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_adjoint__array__control_lbound)
        _arrays[array_handle] = control_lbound
    return control_lbound

def set_array_control_lbound(control_lbound):
    control_lbound[...] = control_lbound

def get_array_control_ubound():
    """
    Element control_ubound ftype=real(rp) pytype=float
    
    
    Defined at m_adjoint.f90 line 20
    
    """
    global control_ubound
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_adjoint__array__control_ubound(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        control_ubound = _arrays[array_handle]
    else:
        control_ubound = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_adjoint__array__control_ubound)
        _arrays[array_handle] = control_ubound
    return control_ubound

def set_array_control_ubound(control_ubound):
    control_ubound[...] = control_ubound

def get_cost():
    """
    Element cost ftype=real(rp) pytype=float
    
    
    Defined at m_adjoint.f90 line 21
    
    """
    return _wrapping.f90wrap_m_adjoint__get__cost()

def set_cost(cost):
    _wrapping.f90wrap_m_adjoint__set__cost(cost)

def get_cost_back():
    """
    Element cost_back ftype=real(rp) pytype=float
    
    
    Defined at m_adjoint.f90 line 22
    
    """
    return _wrapping.f90wrap_m_adjoint__get__cost_back()

def set_cost_back(cost_back):
    _wrapping.f90wrap_m_adjoint__set__cost_back(cost_back)

def get_cost_diff():
    """
    Element cost_diff ftype=real(rp) pytype=float
    
    
    Defined at m_adjoint.f90 line 23
    
    """
    return _wrapping.f90wrap_m_adjoint__get__cost_diff()

def set_cost_diff(cost_diff):
    _wrapping.f90wrap_m_adjoint__set__cost_diff(cost_diff)

def get_cost_perturb():
    """
    Element cost_perturb ftype=real(rp) pytype=float
    
    
    Defined at m_adjoint.f90 line 24
    
    """
    return _wrapping.f90wrap_m_adjoint__get__cost_perturb()

def set_cost_perturb(cost_perturb):
    _wrapping.f90wrap_m_adjoint__set__cost_perturb(cost_perturb)

def get_ic():
    """
    Element ic ftype=integer(ip) pytype=int
    
    
    Defined at m_adjoint.f90 line 25
    
    """
    return _wrapping.f90wrap_m_adjoint__get__ic()

def set_ic(ic):
    _wrapping.f90wrap_m_adjoint__set__ic(ic)

def get_ite_min():
    """
    Element ite_min ftype=integer(ip) pytype=int
    
    
    Defined at m_adjoint.f90 line 26
    
    """
    return _wrapping.f90wrap_m_adjoint__get__ite_min()

def set_ite_min(ite_min):
    _wrapping.f90wrap_m_adjoint__set__ite_min(ite_min)

def get_nb_vars_in_control():
    """
    Element nb_vars_in_control ftype=integer(ip) pytype=int
    
    
    Defined at m_adjoint.f90 line 27
    
    """
    return _wrapping.f90wrap_m_adjoint__get__nb_vars_in_control()

def set_nb_vars_in_control(nb_vars_in_control):
    _wrapping.f90wrap_m_adjoint__set__nb_vars_in_control(nb_vars_in_control)

def get_array_dim_vars_in_control():
    """
    Element dim_vars_in_control ftype=integer(ip) pytype=int
    
    
    Defined at m_adjoint.f90 line 28
    
    """
    global dim_vars_in_control
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_adjoint__array__dim_vars_in_control(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        dim_vars_in_control = _arrays[array_handle]
    else:
        dim_vars_in_control = \
            f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_adjoint__array__dim_vars_in_control)
        _arrays[array_handle] = dim_vars_in_control
    return dim_vars_in_control

def set_array_dim_vars_in_control(dim_vars_in_control):
    dim_vars_in_control[...] = dim_vars_in_control

def get_dim_all():
    """
    Element dim_all ftype=integer(ip) pytype=int
    
    
    Defined at m_adjoint.f90 line 29
    
    """
    return _wrapping.f90wrap_m_adjoint__get__dim_all()

def set_dim_all(dim_all):
    _wrapping.f90wrap_m_adjoint__set__dim_all(dim_all)

def get_array_bathy_cell_copy():
    """
    Element bathy_cell_copy ftype=real(rp) pytype=float
    
    
    Defined at m_adjoint.f90 line 31
    
    """
    global bathy_cell_copy
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_adjoint__array__bathy_cell_copy(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        bathy_cell_copy = _arrays[array_handle]
    else:
        bathy_cell_copy = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_adjoint__array__bathy_cell_copy)
        _arrays[array_handle] = bathy_cell_copy
    return bathy_cell_copy

def set_array_bathy_cell_copy(bathy_cell_copy):
    bathy_cell_copy[...] = bathy_cell_copy


_array_initialisers = [get_array_control, get_array_control_back, \
    get_array_control_diff, get_array_control_perturb, get_array_control_lbound, \
    get_array_control_ubound, get_array_dim_vars_in_control, \
    get_array_bathy_cell_copy]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "m_adjoint".')

for func in _dt_array_initialisers:
    func()
