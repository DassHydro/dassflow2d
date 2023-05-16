"""
Module m_common


Defined at m_common.f90 lines 1-643

"""
from __future__ import print_function, absolute_import, division
from dassflow2d.wrapping import _wrapping
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("wrapping.weights")
class weights(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=weights)
    
    
    Defined at m_common.f90 lines 80-81
    
    """
    def __init__(self, handle=None):
        """
        self = Weights()
        
        
        Defined at m_common.f90 lines 80-81
        
        
        Returns
        -------
        this : Weights
        	Object to be constructed
        
        
        Automatically generated constructor for weights
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_weights_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Weights
        
        
        Defined at m_common.f90 lines 80-81
        
        Parameters
        ----------
        this : Weights
        	Object to be destructed
        
        
        Automatically generated destructor for weights
        """
        if self._alloc:
            _wrapping.f90wrap_weights_finalise(this=self._handle)
    
    @property
    def weights(self):
        """
        Element weights ftype=real(rp) pytype=float
        
        
        Defined at m_common.f90 line 81
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_weights__array__weights(self._handle)
        if array_handle in self._arrays:
            weights = self._arrays[array_handle]
        else:
            weights = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_weights__array__weights)
            self._arrays[array_handle] = weights
        return weights
    
    @weights.setter
    def weights(self, weights):
        self.weights[...] = weights
    
    def __str__(self):
        ret = ['<weights>{\n']
        ret.append('    weights : ')
        ret.append(repr(self.weights))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def test_dt_nearest(dt_to_test):
    """
    test_dt_nearest = test_dt_nearest(dt_to_test)
    
    
    Defined at m_common.f90 lines 90-111
    
    Parameters
    ----------
    dt_to_test : float
    
    Returns
    -------
    test_dt_nearest : bool
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    _COMMENT================================================================================================================
    _COMMENT
    _COMMENT Interface Variables
    _COMMENT================================================================================================================
    _COMMENT
    _COMMENT================================================================================================================
    _COMMENT
    _COMMENT Begin
    _COMMENT================================================================================================================
    _COMMENT
    """
    test_dt_nearest = _wrapping.f90wrap_test_dt_nearest(dt_to_test=dt_to_test)
    return test_dt_nearest

def test_dt_just_after(dt_to_test):
    """
    test_dt_just_after = test_dt_just_after(dt_to_test)
    
    
    Defined at m_common.f90 lines 118-139
    
    Parameters
    ----------
    dt_to_test : float
    
    Returns
    -------
    test_dt_just_after : bool
    
    ================================================================================================================
     Interface Variables
    ================================================================================================================
    _COMMENT================================================================================================================
    _COMMENT
    _COMMENT Interface Variables
    _COMMENT================================================================================================================
    _COMMENT
    _COMMENT================================================================================================================
    _COMMENT
    _COMMENT Begin
    _COMMENT================================================================================================================
    _COMMENT
    """
    test_dt_just_after = _wrapping.f90wrap_test_dt_just_after(dt_to_test=dt_to_test)
    return test_dt_just_after

def reading_args():
    """
    reading_args()
    
    
    Defined at m_common.f90 lines 143-155
    
    
    """
    _wrapping.f90wrap_reading_args()

def machine_number_limits():
    """
    machine_number_limits()
    
    
    Defined at m_common.f90 lines 308-339
    
    
    ================================================================================================================
     Machine zero
    ================================================================================================================
    """
    _wrapping.f90wrap_machine_number_limits()

def swap_r(a, b):
    """
    swap_r(a, b)
    
    
    Defined at m_common.f90 lines 346-350
    
    Parameters
    ----------
    a : float
    b : float
    
    """
    _wrapping.f90wrap_swap_r(a=a, b=b)

def swap_i(a, b):
    """
    swap_i(a, b)
    
    
    Defined at m_common.f90 lines 357-361
    
    Parameters
    ----------
    a : int
    b : int
    
    """
    _wrapping.f90wrap_swap_i(a=a, b=b)

def swap_vec_r(vec, swap):
    """
    swap_vec_r(vec, swap)
    
    
    Defined at m_common.f90 lines 368-378
    
    Parameters
    ----------
    vec : float array
    swap : int array
    
    """
    _wrapping.f90wrap_swap_vec_r(vec=vec, swap=swap)

def swap_vec_i(vec, swap):
    """
    swap_vec_i(vec, swap)
    
    
    Defined at m_common.f90 lines 385-395
    
    Parameters
    ----------
    vec : int array
    swap : int array
    
    """
    _wrapping.f90wrap_swap_vec_i(vec=vec, swap=swap)

def div_by_except_0(a, b):
    """
    div_by_except_0 = div_by_except_0(a, b)
    
    
    Defined at m_common.f90 lines 403-410
    
    Parameters
    ----------
    a : float
    b : float
    
    Returns
    -------
    div_by_except_0 : float
    
    """
    div_by_except_0 = _wrapping.f90wrap_div_by_except_0(a=a, b=b)
    return div_by_except_0

def i4col_sort_a(m, n, a):
    """
    i4col_sort_a(m, n, a)
    
    
    Defined at m_common.f90 lines 421-449
    
    Parameters
    ----------
    m : int
    n : int
    a : int array
    
    ================================================================================================================
     Initialize.
    ================================================================================================================
    """
    _wrapping.f90wrap_i4col_sort_a(m=m, n=n, a=a)

def sort_heap_external(n, indx, i, j, isgn):
    """
    sort_heap_external(n, indx, i, j, isgn)
    
    
    Defined at m_common.f90 lines 462-556
    
    Parameters
    ----------
    n : int
    indx : int
    i : int
    j : int
    isgn : int
    
    ================================================================================================================
     INDX = 0: This is the first call.
    ================================================================================================================
    """
    _wrapping.f90wrap_sort_heap_external(n=n, indx=indx, i=i, j=j, isgn=isgn)

def i4col_swap(m, n, a, i, j):
    """
    i4col_swap(m, n, a, i, j)
    
    
    Defined at m_common.f90 lines 566-575
    
    Parameters
    ----------
    m : int
    n : int
    a : int array
    i : int
    j : int
    
    """
    _wrapping.f90wrap_i4col_swap(m=m, n=n, a=a, i=i, j=j)

def i4col_compare(m, n, a, i, j):
    """
    isgn = i4col_compare(m, n, a, i, j)
    
    
    Defined at m_common.f90 lines 586-604
    
    Parameters
    ----------
    m : int
    n : int
    a : int array
    i : int
    j : int
    
    Returns
    -------
    isgn : int
    
    """
    isgn = _wrapping.f90wrap_i4col_compare(m=m, n=n, a=a, i=i, j=j)
    return isgn

def file_name_ext(file_name, typ):
    """
    file_name_res = file_name_ext(file_name, typ)
    
    
    Defined at m_common.f90 lines 612-624
    
    Parameters
    ----------
    file_name : str
    typ : str
    
    Returns
    -------
    file_name_res : str
    
    """
    file_name_res = _wrapping.f90wrap_file_name_ext(file_name=file_name, typ=typ)
    return file_name_res

def count_lines(file_name):
    """
    nb_lines = count_lines(file_name)
    
    
    Defined at m_common.f90 lines 630-643
    
    Parameters
    ----------
    file_name : str
    
    Returns
    -------
    nb_lines : int
    
    """
    nb_lines = _wrapping.f90wrap_count_lines(file_name=file_name)
    return nb_lines

def get_ip():
    """
    Element ip ftype=integer pytype=int
    
    
    Defined at m_common.f90 line 3
    
    """
    return _wrapping.f90wrap_m_common__get__ip()

ip = get_ip()

def get_rp():
    """
    Element rp ftype=integer pytype=int
    
    
    Defined at m_common.f90 line 4
    
    """
    return _wrapping.f90wrap_m_common__get__rp()

rp = get_rp()

def get_lchar():
    """
    Element lchar ftype=integer pytype=int
    
    
    Defined at m_common.f90 line 5
    
    """
    return _wrapping.f90wrap_m_common__get__lchar()

lchar = get_lchar()

def get_i():
    """
    Element i ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 6
    
    """
    return _wrapping.f90wrap_m_common__get__i()

def set_i(i):
    _wrapping.f90wrap_m_common__set__i(i)

def get_ie():
    """
    Element ie ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 6
    
    """
    return _wrapping.f90wrap_m_common__get__ie()

def set_ie(ie):
    _wrapping.f90wrap_m_common__set__ie(ie)

def get_ik():
    """
    Element ik ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 6
    
    """
    return _wrapping.f90wrap_m_common__get__ik()

def set_ik(ik):
    _wrapping.f90wrap_m_common__set__ik(ik)

def get_ike():
    """
    Element ike ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 6
    
    """
    return _wrapping.f90wrap_m_common__get__ike()

def set_ike(ike):
    _wrapping.f90wrap_m_common__set__ike(ike)

def get_ib():
    """
    Element ib ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 6
    
    """
    return _wrapping.f90wrap_m_common__get__ib()

def set_ib(ib):
    _wrapping.f90wrap_m_common__set__ib(ib)

def get_j():
    """
    Element j ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 7
    
    """
    return _wrapping.f90wrap_m_common__get__j()

def set_j(j):
    _wrapping.f90wrap_m_common__set__j(j)

def get_je():
    """
    Element je ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 7
    
    """
    return _wrapping.f90wrap_m_common__get__je()

def set_je(je):
    _wrapping.f90wrap_m_common__set__je(je)

def get_jk():
    """
    Element jk ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 7
    
    """
    return _wrapping.f90wrap_m_common__get__jk()

def set_jk(jk):
    _wrapping.f90wrap_m_common__set__jk(jk)

def get_jke():
    """
    Element jke ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 7
    
    """
    return _wrapping.f90wrap_m_common__get__jke()

def set_jke(jke):
    _wrapping.f90wrap_m_common__set__jke(jke)

def get_jb():
    """
    Element jb ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 7
    
    """
    return _wrapping.f90wrap_m_common__get__jb()

def set_jb(jb):
    _wrapping.f90wrap_m_common__set__jb(jb)

def get_k():
    """
    Element k ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 8
    
    """
    return _wrapping.f90wrap_m_common__get__k()

def set_k(k):
    _wrapping.f90wrap_m_common__set__k(k)

def get_ke():
    """
    Element ke ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 8
    
    """
    return _wrapping.f90wrap_m_common__get__ke()

def set_ke(ke):
    _wrapping.f90wrap_m_common__set__ke(ke)

def get_kk():
    """
    Element kk ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 8
    
    """
    return _wrapping.f90wrap_m_common__get__kk()

def set_kk(kk):
    _wrapping.f90wrap_m_common__set__kk(kk)

def get_kke():
    """
    Element kke ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 8
    
    """
    return _wrapping.f90wrap_m_common__get__kke()

def set_kke(kke):
    _wrapping.f90wrap_m_common__set__kke(kke)

def get_kb():
    """
    Element kb ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 8
    
    """
    return _wrapping.f90wrap_m_common__get__kb()

def set_kb(kb):
    _wrapping.f90wrap_m_common__set__kb(kb)

def get_mesh_type():
    """
    Element mesh_type ftype=character(len=lchar) pytype=str
    
    
    Defined at m_common.f90 line 9
    
    """
    return _wrapping.f90wrap_m_common__get__mesh_type()

def set_mesh_type(mesh_type):
    _wrapping.f90wrap_m_common__set__mesh_type(mesh_type)

def get_mesh_name():
    """
    Element mesh_name ftype=character(len=lchar) pytype=str
    
    
    Defined at m_common.f90 line 10
    
    """
    return _wrapping.f90wrap_m_common__get__mesh_name()

def set_mesh_name(mesh_name):
    _wrapping.f90wrap_m_common__set__mesh_name(mesh_name)

def get_bc_n():
    """
    Element bc_n ftype=character(len=lchar) pytype=str
    
    
    Defined at m_common.f90 line 11
    
    """
    return _wrapping.f90wrap_m_common__get__bc_n()

def set_bc_n(bc_n):
    _wrapping.f90wrap_m_common__set__bc_n(bc_n)

def get_bc_s():
    """
    Element bc_s ftype=character(len=lchar) pytype=str
    
    
    Defined at m_common.f90 line 12
    
    """
    return _wrapping.f90wrap_m_common__get__bc_s()

def set_bc_s(bc_s):
    _wrapping.f90wrap_m_common__set__bc_s(bc_s)

def get_bc_w():
    """
    Element bc_w ftype=character(len=lchar) pytype=str
    
    
    Defined at m_common.f90 line 13
    
    """
    return _wrapping.f90wrap_m_common__get__bc_w()

def set_bc_w(bc_w):
    _wrapping.f90wrap_m_common__set__bc_w(bc_w)

def get_bc_e():
    """
    Element bc_e ftype=character(len=lchar) pytype=str
    
    
    Defined at m_common.f90 line 14
    
    """
    return _wrapping.f90wrap_m_common__get__bc_e()

def set_bc_e(bc_e):
    _wrapping.f90wrap_m_common__set__bc_e(bc_e)

def get_bc_rain():
    """
    Element bc_rain ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 15
    
    """
    return _wrapping.f90wrap_m_common__get__bc_rain()

def set_bc_rain(bc_rain):
    _wrapping.f90wrap_m_common__set__bc_rain(bc_rain)

def get_bc_infil():
    """
    Element bc_infil ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 16
    
    """
    return _wrapping.f90wrap_m_common__get__bc_infil()

def set_bc_infil(bc_infil):
    _wrapping.f90wrap_m_common__set__bc_infil(bc_infil)

def get_lx():
    """
    Element lx ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 17
    
    """
    return _wrapping.f90wrap_m_common__get__lx()

def set_lx(lx):
    _wrapping.f90wrap_m_common__set__lx(lx)

def get_ly():
    """
    Element ly ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 18
    
    """
    return _wrapping.f90wrap_m_common__get__ly()

def set_ly(ly):
    _wrapping.f90wrap_m_common__set__ly(ly)

def get_nx():
    """
    Element nx ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 19
    
    """
    return _wrapping.f90wrap_m_common__get__nx()

def set_nx(nx):
    _wrapping.f90wrap_m_common__set__nx(nx)

def get_ny():
    """
    Element ny ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 20
    
    """
    return _wrapping.f90wrap_m_common__get__ny()

def set_ny(ny):
    _wrapping.f90wrap_m_common__set__ny(ny)

def get_ts():
    """
    Element ts ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 21
    
    """
    return _wrapping.f90wrap_m_common__get__ts()

def set_ts(ts):
    _wrapping.f90wrap_m_common__set__ts(ts)

def get_adapt_dt():
    """
    Element adapt_dt ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 22
    
    """
    return _wrapping.f90wrap_m_common__get__adapt_dt()

def set_adapt_dt(adapt_dt):
    _wrapping.f90wrap_m_common__set__adapt_dt(adapt_dt)

def get_dt():
    """
    Element dt ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 23
    
    """
    return _wrapping.f90wrap_m_common__get__dt()

def set_dt(dt):
    _wrapping.f90wrap_m_common__set__dt(dt)

def get_cfl():
    """
    Element cfl ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 24
    
    """
    return _wrapping.f90wrap_m_common__get__cfl()

def set_cfl(cfl):
    _wrapping.f90wrap_m_common__set__cfl(cfl)

def get_dtw():
    """
    Element dtw ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 25
    
    """
    return _wrapping.f90wrap_m_common__get__dtw()

def set_dtw(dtw):
    _wrapping.f90wrap_m_common__set__dtw(dtw)

def get_dtp():
    """
    Element dtp ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 26
    
    """
    return _wrapping.f90wrap_m_common__get__dtp()

def set_dtp(dtp):
    _wrapping.f90wrap_m_common__set__dtp(dtp)

def get_dta():
    """
    Element dta ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 27
    
    """
    return _wrapping.f90wrap_m_common__get__dta()

def set_dta(dta):
    _wrapping.f90wrap_m_common__set__dta(dta)

def get_w_tecplot():
    """
    Element w_tecplot ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 28
    
    """
    return _wrapping.f90wrap_m_common__get__w_tecplot()

def set_w_tecplot(w_tecplot):
    _wrapping.f90wrap_m_common__set__w_tecplot(w_tecplot)

def get_w_vtk():
    """
    Element w_vtk ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 29
    
    """
    return _wrapping.f90wrap_m_common__get__w_vtk()

def set_w_vtk(w_vtk):
    _wrapping.f90wrap_m_common__set__w_vtk(w_vtk)

def get_w_gnuplot():
    """
    Element w_gnuplot ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 30
    
    """
    return _wrapping.f90wrap_m_common__get__w_gnuplot()

def set_w_gnuplot(w_gnuplot):
    _wrapping.f90wrap_m_common__set__w_gnuplot(w_gnuplot)

def get_w_bin():
    """
    Element w_bin ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 31
    
    """
    return _wrapping.f90wrap_m_common__get__w_bin()

def set_w_bin(w_bin):
    _wrapping.f90wrap_m_common__set__w_bin(w_bin)

def get_w_exact():
    """
    Element w_exact ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 32
    
    """
    return _wrapping.f90wrap_m_common__get__w_exact()

def set_w_exact(w_exact):
    _wrapping.f90wrap_m_common__set__w_exact(w_exact)

def get_w_norm():
    """
    Element w_norm ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 33
    
    """
    return _wrapping.f90wrap_m_common__get__w_norm()

def set_w_norm(w_norm):
    _wrapping.f90wrap_m_common__set__w_norm(w_norm)

def get_w_obs():
    """
    Element w_obs ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 34
    
    """
    return _wrapping.f90wrap_m_common__get__w_obs()

def set_w_obs(w_obs):
    _wrapping.f90wrap_m_common__set__w_obs(w_obs)

def get_use_obs():
    """
    Element use_obs ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 35
    
    """
    return _wrapping.f90wrap_m_common__get__use_obs()

def set_use_obs(use_obs):
    _wrapping.f90wrap_m_common__set__use_obs(use_obs)

def get_spatial_scheme():
    """
    Element spatial_scheme ftype=character(len=lchar) pytype=str
    
    
    Defined at m_common.f90 line 36
    
    """
    return _wrapping.f90wrap_m_common__get__spatial_scheme()

def set_spatial_scheme(spatial_scheme):
    _wrapping.f90wrap_m_common__set__spatial_scheme(spatial_scheme)

def get_temp_scheme():
    """
    Element temp_scheme ftype=character(len=lchar) pytype=str
    
    
    Defined at m_common.f90 line 37
    
    """
    return _wrapping.f90wrap_m_common__get__temp_scheme()

def set_temp_scheme(temp_scheme):
    _wrapping.f90wrap_m_common__set__temp_scheme(temp_scheme)

def get_array_args():
    """
    Element args ftype=character(len=lchar) pytype=str
    
    
    Defined at m_common.f90 line 38
    
    """
    global args
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_common__array__args(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        args = _arrays[array_handle]
    else:
        args = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_common__array__args)
        _arrays[array_handle] = args
    return args

def set_array_args(args):
    args[...] = args

def get_max_nt_for_direct():
    """
    Element max_nt_for_direct ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 39
    
    """
    return _wrapping.f90wrap_m_common__get__max_nt_for_direct()

def set_max_nt_for_direct(max_nt_for_direct):
    _wrapping.f90wrap_m_common__set__max_nt_for_direct(max_nt_for_direct)

def get_max_nt_for_adjoint():
    """
    Element max_nt_for_adjoint ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 40
    
    """
    return _wrapping.f90wrap_m_common__get__max_nt_for_adjoint()

def set_max_nt_for_adjoint(max_nt_for_adjoint):
    _wrapping.f90wrap_m_common__set__max_nt_for_adjoint(max_nt_for_adjoint)

def get_length_real():
    """
    Element length_real ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 42
    
    """
    return _wrapping.f90wrap_m_common__get__length_real()

def set_length_real(length_real):
    _wrapping.f90wrap_m_common__set__length_real(length_real)

def get_nt():
    """
    Element nt ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 43
    
    """
    return _wrapping.f90wrap_m_common__get__nt()

def set_nt(nt):
    _wrapping.f90wrap_m_common__set__nt(nt)

def get_nt0():
    """
    Element nt0 ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 43
    
    """
    return _wrapping.f90wrap_m_common__get__nt0()

def set_nt0(nt0):
    _wrapping.f90wrap_m_common__set__nt0(nt0)

def get_tc():
    """
    Element tc ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 44
    
    """
    return _wrapping.f90wrap_m_common__get__tc()

def set_tc(tc):
    _wrapping.f90wrap_m_common__set__tc(tc)

def get_tc0():
    """
    Element tc0 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 44
    
    """
    return _wrapping.f90wrap_m_common__get__tc0()

def set_tc0(tc0):
    _wrapping.f90wrap_m_common__set__tc0(tc0)

def get_end_time_loop():
    """
    Element end_time_loop ftype=logical pytype=bool
    
    
    Defined at m_common.f90 line 45
    
    """
    return _wrapping.f90wrap_m_common__get__end_time_loop()

def set_end_time_loop(end_time_loop):
    _wrapping.f90wrap_m_common__set__end_time_loop(end_time_loop)

def get_dx():
    """
    Element dx ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 46
    
    """
    return _wrapping.f90wrap_m_common__get__dx()

def set_dx(dx):
    _wrapping.f90wrap_m_common__set__dx(dx)

def get_dy():
    """
    Element dy ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 47
    
    """
    return _wrapping.f90wrap_m_common__get__dy()

def set_dy(dy):
    _wrapping.f90wrap_m_common__set__dy(dy)

def get_array_is_file_open():
    """
    Element is_file_open ftype=character(len=lchar) pytype=str
    
    
    Defined at m_common.f90 line 48
    
    """
    global is_file_open
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_common__array__is_file_open(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        is_file_open = _arrays[array_handle]
    else:
        is_file_open = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_common__array__is_file_open)
        _arrays[array_handle] = is_file_open
    return is_file_open

def set_array_is_file_open(is_file_open):
    is_file_open[...] = is_file_open

def get_file_open_counter():
    """
    Element file_open_counter ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 49
    
    """
    return _wrapping.f90wrap_m_common__get__file_open_counter()

def set_file_open_counter(file_open_counter):
    _wrapping.f90wrap_m_common__set__file_open_counter(file_open_counter)

def get_array_file_exist():
    """
    Element file_exist ftype=logical pytype=bool
    
    
    Defined at m_common.f90 line 50
    
    """
    global file_exist
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_common__array__file_exist(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        file_exist = _arrays[array_handle]
    else:
        file_exist = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_common__array__file_exist)
        _arrays[array_handle] = file_exist
    return file_exist

def set_array_file_exist(file_exist):
    file_exist[...] = file_exist

def get_buffer():
    """
    Element buffer ftype=character(len=1028) pytype=str
    
    
    Defined at m_common.f90 line 51
    
    """
    return _wrapping.f90wrap_m_common__get__buffer()

def set_buffer(buffer):
    _wrapping.f90wrap_m_common__set__buffer(buffer)

def get_logic_test():
    """
    Element logic_test ftype=logical pytype=bool
    
    
    Defined at m_common.f90 line 52
    
    """
    return _wrapping.f90wrap_m_common__get__logic_test()

def set_logic_test(logic_test):
    _wrapping.f90wrap_m_common__set__logic_test(logic_test)

def get_array_norm_inf():
    """
    Element norm_inf ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 53
    
    """
    global norm_inf
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_common__array__norm_inf(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        norm_inf = _arrays[array_handle]
    else:
        norm_inf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_common__array__norm_inf)
        _arrays[array_handle] = norm_inf
    return norm_inf

def set_array_norm_inf(norm_inf):
    norm_inf[...] = norm_inf

def get_array_norm_l1():
    """
    Element norm_l1 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 54
    
    """
    global norm_l1
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_common__array__norm_l1(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        norm_l1 = _arrays[array_handle]
    else:
        norm_l1 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_common__array__norm_l1)
        _arrays[array_handle] = norm_l1
    return norm_l1

def set_array_norm_l1(norm_l1):
    norm_l1[...] = norm_l1

def get_array_norm_l2():
    """
    Element norm_l2 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 55
    
    """
    global norm_l2
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_common__array__norm_l2(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        norm_l2 = _arrays[array_handle]
    else:
        norm_l2 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_common__array__norm_l2)
        _arrays[array_handle] = norm_l2
    return norm_l2

def set_array_norm_l2(norm_l2):
    norm_l2[...] = norm_l2

def get_verbose():
    """
    Element verbose ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 56
    
    """
    return _wrapping.f90wrap_m_common__get__verbose()

def set_verbose(verbose):
    _wrapping.f90wrap_m_common__set__verbose(verbose)

def get_restart_min():
    """
    Element restart_min ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 57
    
    """
    return _wrapping.f90wrap_m_common__get__restart_min()

def set_restart_min(restart_min):
    _wrapping.f90wrap_m_common__set__restart_min(restart_min)

def get_eps_min():
    """
    Element eps_min ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 58
    
    """
    return _wrapping.f90wrap_m_common__get__eps_min()

def set_eps_min(eps_min):
    _wrapping.f90wrap_m_common__set__eps_min(eps_min)

def get_zerom():
    """
    Element zerom ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 59
    
    """
    return _wrapping.f90wrap_m_common__get__zerom()

def set_zerom(zerom):
    _wrapping.f90wrap_m_common__set__zerom(zerom)

def get_pinfm():
    """
    Element pinfm ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 60
    
    """
    return _wrapping.f90wrap_m_common__get__pinfm()

def set_pinfm(pinfm):
    _wrapping.f90wrap_m_common__set__pinfm(pinfm)

def get_minfm():
    """
    Element minfm ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 60
    
    """
    return _wrapping.f90wrap_m_common__get__minfm()

def set_minfm(minfm):
    _wrapping.f90wrap_m_common__set__minfm(minfm)

def get_hugem():
    """
    Element hugem ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 61
    
    """
    return _wrapping.f90wrap_m_common__get__hugem()

def set_hugem(hugem):
    _wrapping.f90wrap_m_common__set__hugem(hugem)

def get_tinym():
    """
    Element tinym ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 61
    
    """
    return _wrapping.f90wrap_m_common__get__tinym()

def set_tinym(tinym):
    _wrapping.f90wrap_m_common__set__tinym(tinym)

def get_zero():
    """
    Element zero ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 62
    
    """
    return _wrapping.f90wrap_m_common__get__zero()

zero = get_zero()

def get_one():
    """
    Element one ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 63
    
    """
    return _wrapping.f90wrap_m_common__get__one()

one = get_one()

def get_two():
    """
    Element two ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 64
    
    """
    return _wrapping.f90wrap_m_common__get__two()

two = get_two()

def get_demi():
    """
    Element demi ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 65
    
    """
    return _wrapping.f90wrap_m_common__get__demi()

demi = get_demi()

def get_d1p4():
    """
    Element d1p4 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 66
    
    """
    return _wrapping.f90wrap_m_common__get__d1p4()

d1p4 = get_d1p4()

def get_d1p3():
    """
    Element d1p3 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 67
    
    """
    return _wrapping.f90wrap_m_common__get__d1p3()

d1p3 = get_d1p3()

def get_d2p3():
    """
    Element d2p3 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 68
    
    """
    return _wrapping.f90wrap_m_common__get__d2p3()

d2p3 = get_d2p3()

def get_d4p3():
    """
    Element d4p3 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 69
    
    """
    return _wrapping.f90wrap_m_common__get__d4p3()

d4p3 = get_d4p3()

def get_d5p3():
    """
    Element d5p3 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 70
    
    """
    return _wrapping.f90wrap_m_common__get__d5p3()

d5p3 = get_d5p3()

def get_d7p3():
    """
    Element d7p3 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 71
    
    """
    return _wrapping.f90wrap_m_common__get__d7p3()

d7p3 = get_d7p3()

def get_d8p3():
    """
    Element d8p3 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 72
    
    """
    return _wrapping.f90wrap_m_common__get__d8p3()

d8p3 = get_d8p3()

def get_d10p3():
    """
    Element d10p3 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 73
    
    """
    return _wrapping.f90wrap_m_common__get__d10p3()

d10p3 = get_d10p3()

def get_d3p2():
    """
    Element d3p2 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 74
    
    """
    return _wrapping.f90wrap_m_common__get__d3p2()

d3p2 = get_d3p2()

def get_d3p5():
    """
    Element d3p5 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 75
    
    """
    return _wrapping.f90wrap_m_common__get__d3p5()

d3p5 = get_d3p5()

def get_d3p8():
    """
    Element d3p8 ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 76
    
    """
    return _wrapping.f90wrap_m_common__get__d3p8()

d3p8 = get_d3p8()

def get_pi():
    """
    Element pi ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 77
    
    """
    return _wrapping.f90wrap_m_common__get__pi()

pi = get_pi()

def get_proc():
    """
    Element proc ftype=real(rp) pytype=float
    
    
    Defined at m_common.f90 line 78
    
    """
    return _wrapping.f90wrap_m_common__get__proc()

proc = get_proc()

def get_np():
    """
    Element np ftype=integer(ip) pytype=int
    
    
    Defined at m_common.f90 line 79
    
    """
    return _wrapping.f90wrap_m_common__get__np()

np = get_np()


_array_initialisers = [get_array_args, get_array_is_file_open, \
    get_array_file_exist, get_array_norm_inf, get_array_norm_l1, \
    get_array_norm_l2]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "m_common".')

for func in _dt_array_initialisers:
    func()
