"""
Module m_model


Defined at m_sw_mono.f90 lines 1-829

"""
from __future__ import print_function, absolute_import, division
from dassflow2d.wrapping import _wrapping
import f90wrap.runtime
import logging
import numpy
from dassflow2d.wrapping.m_linear_algebra import vec2d

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("wrapping.unk")
class unk(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=unk)
    
    
    Defined at m_sw_mono.f90 lines 16-26
    
    """
    def __init__(self, mesh, handle=None):
        """
        self = Unk(mesh)
        
        
        Defined at m_sw_mono.f90 lines 770-794
        
        Parameters
        ----------
        mesh : Msh
        
        Returns
        -------
        dof : Unk
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_unk_initialise(mesh=mesh._handle)
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Unk
        
        
        Defined at m_sw_mono.f90 lines 797-811
        
        Parameters
        ----------
        dof : Unk
        
        """
        if self._alloc:
            _wrapping.f90wrap_unk_finalise(dof=self._handle)
    
    @property
    def h(self):
        """
        Element h ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 18
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_unk__array__h(self._handle)
        if array_handle in self._arrays:
            h = self._arrays[array_handle]
        else:
            h = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_unk__array__h)
            self._arrays[array_handle] = h
        return h
    
    @h.setter
    def h(self, h):
        self.h[...] = h
    
    @property
    def u(self):
        """
        Element u ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 19
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_unk__array__u(self._handle)
        if array_handle in self._arrays:
            u = self._arrays[array_handle]
        else:
            u = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_unk__array__u)
            self._arrays[array_handle] = u
        return u
    
    @u.setter
    def u(self, u):
        self.u[...] = u
    
    @property
    def v(self):
        """
        Element v ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 20
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_unk__array__v(self._handle)
        if array_handle in self._arrays:
            v = self._arrays[array_handle]
        else:
            v = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_unk__array__v)
            self._arrays[array_handle] = v
        return v
    
    @v.setter
    def v(self, v):
        self.v[...] = v
    
    @property
    def infil(self):
        """
        Element infil ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 21
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_unk__array__infil(self._handle)
        if array_handle in self._arrays:
            infil = self._arrays[array_handle]
        else:
            infil = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_unk__array__infil)
            self._arrays[array_handle] = infil
        return infil
    
    @infil.setter
    def infil(self, infil):
        self.infil[...] = infil
    
    @property
    def t_display(self):
        """
        Element t_display ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 22
        
        """
        return _wrapping.f90wrap_unk__get__t_display(self._handle)
    
    @t_display.setter
    def t_display(self, t_display):
        _wrapping.f90wrap_unk__set__t_display(self._handle, t_display)
    
    def init_array_grad_h(self):
        self.grad_h = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_unk__array_getitem__grad_h,
                                        _wrapping.f90wrap_unk__array_setitem__grad_h,
                                        _wrapping.f90wrap_unk__array_len__grad_h,
                                        """
        Element grad_h ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_sw_mono.f90 line 23
        
        """, vec2d)
        return self.grad_h
    
    def init_array_grad_u(self):
        self.grad_u = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_unk__array_getitem__grad_u,
                                        _wrapping.f90wrap_unk__array_setitem__grad_u,
                                        _wrapping.f90wrap_unk__array_len__grad_u,
                                        """
        Element grad_u ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_sw_mono.f90 line 24
        
        """, vec2d)
        return self.grad_u
    
    def init_array_grad_v(self):
        self.grad_v = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_unk__array_getitem__grad_v,
                                        _wrapping.f90wrap_unk__array_setitem__grad_v,
                                        _wrapping.f90wrap_unk__array_len__grad_v,
                                        """
        Element grad_v ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_sw_mono.f90 line 25
        
        """, vec2d)
        return self.grad_v
    
    def init_array_grad_z(self):
        self.grad_z = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_unk__array_getitem__grad_z,
                                        _wrapping.f90wrap_unk__array_setitem__grad_z,
                                        _wrapping.f90wrap_unk__array_len__grad_z,
                                        """
        Element grad_z ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_sw_mono.f90 line 26
        
        """, vec2d)
        return self.grad_z
    
    def __str__(self):
        ret = ['<unk>{\n']
        ret.append('    h : ')
        ret.append(repr(self.h))
        ret.append(',\n    u : ')
        ret.append(repr(self.u))
        ret.append(',\n    v : ')
        ret.append(repr(self.v))
        ret.append(',\n    infil : ')
        ret.append(repr(self.infil))
        ret.append(',\n    t_display : ')
        ret.append(repr(self.t_display))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = [init_array_grad_h, init_array_grad_u, \
        init_array_grad_v, init_array_grad_z]
    

@f90wrap.runtime.register_class("wrapping.xsshp")
class xsshp(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=xsshp)
    
    
    Defined at m_sw_mono.f90 lines 44-50
    
    """
    def __init__(self, handle=None):
        """
        self = Xsshp()
        
        
        Defined at m_sw_mono.f90 lines 44-50
        
        
        Returns
        -------
        this : Xsshp
        	Object to be constructed
        
        
        Automatically generated constructor for xsshp
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_xsshp_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Xsshp
        
        
        Defined at m_sw_mono.f90 lines 44-50
        
        Parameters
        ----------
        this : Xsshp
        	Object to be destructed
        
        
        Automatically generated destructor for xsshp
        """
        if self._alloc:
            _wrapping.f90wrap_xsshp_finalise(this=self._handle)
    
    @property
    def xleft(self):
        """
        Element xleft ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 45
        
        """
        return _wrapping.f90wrap_xsshp__get__xleft(self._handle)
    
    @xleft.setter
    def xleft(self, xleft):
        _wrapping.f90wrap_xsshp__set__xleft(self._handle, xleft)
    
    @property
    def xcenter(self):
        """
        Element xcenter ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 46
        
        """
        return _wrapping.f90wrap_xsshp__get__xcenter(self._handle)
    
    @xcenter.setter
    def xcenter(self, xcenter):
        _wrapping.f90wrap_xsshp__set__xcenter(self._handle, xcenter)
    
    @property
    def xright(self):
        """
        Element xright ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 47
        
        """
        return _wrapping.f90wrap_xsshp__get__xright(self._handle)
    
    @xright.setter
    def xright(self, xright):
        _wrapping.f90wrap_xsshp__set__xright(self._handle, xright)
    
    @property
    def s(self):
        """
        Element s ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 48
        
        """
        return _wrapping.f90wrap_xsshp__get__s(self._handle)
    
    @s.setter
    def s(self, s):
        _wrapping.f90wrap_xsshp__set__s(self._handle, s)
    
    @property
    def hmax(self):
        """
        Element hmax ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 49
        
        """
        return _wrapping.f90wrap_xsshp__get__hmax(self._handle)
    
    @hmax.setter
    def hmax(self, hmax):
        _wrapping.f90wrap_xsshp__set__hmax(self._handle, hmax)
    
    @property
    def topz(self):
        """
        Element topz ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 50
        
        """
        return _wrapping.f90wrap_xsshp__get__topz(self._handle)
    
    @topz.setter
    def topz(self, topz):
        _wrapping.f90wrap_xsshp__set__topz(self._handle, topz)
    
    def __str__(self):
        ret = ['<xsshp>{\n']
        ret.append('    xleft : ')
        ret.append(repr(self.xleft))
        ret.append(',\n    xcenter : ')
        ret.append(repr(self.xcenter))
        ret.append(',\n    xright : ')
        ret.append(repr(self.xright))
        ret.append(',\n    s : ')
        ret.append(repr(self.s))
        ret.append(',\n    hmax : ')
        ret.append(repr(self.hmax))
        ret.append(',\n    topz : ')
        ret.append(repr(self.topz))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.greenampt")
class greenampt(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=greenampt)
    
    
    Defined at m_sw_mono.f90 lines 56-59
    
    """
    def __init__(self, handle=None):
        """
        self = Greenampt()
        
        
        Defined at m_sw_mono.f90 lines 56-59
        
        
        Returns
        -------
        this : Greenampt
        	Object to be constructed
        
        
        Automatically generated constructor for greenampt
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_greenampt_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Greenampt
        
        
        Defined at m_sw_mono.f90 lines 56-59
        
        Parameters
        ----------
        this : Greenampt
        	Object to be destructed
        
        
        Automatically generated destructor for greenampt
        """
        if self._alloc:
            _wrapping.f90wrap_greenampt_finalise(this=self._handle)
    
    @property
    def psif(self):
        """
        Element psif ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 57
        
        """
        return _wrapping.f90wrap_greenampt__get__psif(self._handle)
    
    @psif.setter
    def psif(self, psif):
        _wrapping.f90wrap_greenampt__set__psif(self._handle, psif)
    
    @property
    def ks(self):
        """
        Element ks ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 58
        
        """
        return _wrapping.f90wrap_greenampt__get__ks(self._handle)
    
    @ks.setter
    def ks(self, ks):
        _wrapping.f90wrap_greenampt__set__ks(self._handle, ks)
    
    @property
    def deltatheta(self):
        """
        Element deltatheta ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 59
        
        """
        return _wrapping.f90wrap_greenampt__get__deltatheta(self._handle)
    
    @deltatheta.setter
    def deltatheta(self, deltatheta):
        _wrapping.f90wrap_greenampt__set__deltatheta(self._handle, deltatheta)
    
    def __str__(self):
        ret = ['<greenampt>{\n']
        ret.append('    psif : ')
        ret.append(repr(self.psif))
        ret.append(',\n    ks : ')
        ret.append(repr(self.ks))
        ret.append(',\n    deltatheta : ')
        ret.append(repr(self.deltatheta))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.scs_cn")
class scs_cn(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=scs_cn)
    
    
    Defined at m_sw_mono.f90 lines 61-63
    
    """
    def __init__(self, handle=None):
        """
        self = Scs_Cn()
        
        
        Defined at m_sw_mono.f90 lines 61-63
        
        
        Returns
        -------
        this : Scs_Cn
        	Object to be constructed
        
        
        Automatically generated constructor for scs_cn
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_scs_cn_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Scs_Cn
        
        
        Defined at m_sw_mono.f90 lines 61-63
        
        Parameters
        ----------
        this : Scs_Cn
        	Object to be destructed
        
        
        Automatically generated destructor for scs_cn
        """
        if self._alloc:
            _wrapping.f90wrap_scs_cn_finalise(this=self._handle)
    
    @property
    def lambdacn(self):
        """
        Element lambdacn ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 62
        
        """
        return _wrapping.f90wrap_scs_cn__get__lambdacn(self._handle)
    
    @lambdacn.setter
    def lambdacn(self, lambdacn):
        _wrapping.f90wrap_scs_cn__set__lambdacn(self._handle, lambdacn)
    
    @property
    def cn(self):
        """
        Element cn ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 63
        
        """
        return _wrapping.f90wrap_scs_cn__get__cn(self._handle)
    
    @cn.setter
    def cn(self, cn):
        _wrapping.f90wrap_scs_cn__set__cn(self._handle, cn)
    
    def __str__(self):
        ret = ['<scs_cn>{\n']
        ret.append('    lambdacn : ')
        ret.append(repr(self.lambdacn))
        ret.append(',\n    cn : ')
        ret.append(repr(self.cn))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.infiltration_data")
class infiltration_data(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=infiltration_data)
    
    
    Defined at m_sw_mono.f90 lines 68-75
    
    """
    def __init__(self, handle=None):
        """
        self = Infiltration_Data()
        
        
        Defined at m_sw_mono.f90 lines 68-75
        
        
        Returns
        -------
        this : Infiltration_Data
        	Object to be constructed
        
        
        Automatically generated constructor for infiltration_data
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_infiltration_data_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Infiltration_Data
        
        
        Defined at m_sw_mono.f90 lines 68-75
        
        Parameters
        ----------
        this : Infiltration_Data
        	Object to be destructed
        
        
        Automatically generated destructor for infiltration_data
        """
        if self._alloc:
            _wrapping.f90wrap_infiltration_data_finalise(this=self._handle)
    
    @property
    def nland(self):
        """
        Element nland ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 69
        
        """
        return _wrapping.f90wrap_infiltration_data__get__nland(self._handle)
    
    @nland.setter
    def nland(self, nland):
        _wrapping.f90wrap_infiltration_data__set__nland(self._handle, nland)
    
    @property
    def land(self):
        """
        Element land ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 70
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_infiltration_data__array__land(self._handle)
        if array_handle in self._arrays:
            land = self._arrays[array_handle]
        else:
            land = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_infiltration_data__array__land)
            self._arrays[array_handle] = land
        return land
    
    @land.setter
    def land(self, land):
        self.land[...] = land
    
    def init_array_ga(self):
        self.ga = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_infiltration_data__array_getitem__ga,
                                        _wrapping.f90wrap_infiltration_data__array_setitem__ga,
                                        _wrapping.f90wrap_infiltration_data__array_len__ga,
                                        """
        Element ga ftype=type(greenampt) pytype=Greenampt
        
        
        Defined at m_sw_mono.f90 line 72
        
        """, greenampt)
        return self.ga
    
    def init_array_scs(self):
        self.scs = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_infiltration_data__array_getitem__scs,
                                        _wrapping.f90wrap_infiltration_data__array_setitem__scs,
                                        _wrapping.f90wrap_infiltration_data__array_len__scs,
                                        """
        Element scs ftype=type(scs_cn) pytype=Scs_Cn
        
        
        Defined at m_sw_mono.f90 line 73
        
        """, scs_cn)
        return self.scs
    
    @property
    def coord(self):
        """
        Element coord ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 74
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_infiltration_data__array__coord(self._handle)
        if array_handle in self._arrays:
            coord = self._arrays[array_handle]
        else:
            coord = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_infiltration_data__array__coord)
            self._arrays[array_handle] = coord
        return coord
    
    @coord.setter
    def coord(self, coord):
        self.coord[...] = coord
    
    @property
    def h_infil_max(self):
        """
        Element h_infil_max ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 75
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_infiltration_data__array__h_infil_max(self._handle)
        if array_handle in self._arrays:
            h_infil_max = self._arrays[array_handle]
        else:
            h_infil_max = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_infiltration_data__array__h_infil_max)
            self._arrays[array_handle] = h_infil_max
        return h_infil_max
    
    @h_infil_max.setter
    def h_infil_max(self, h_infil_max):
        self.h_infil_max[...] = h_infil_max
    
    def __str__(self):
        ret = ['<infiltration_data>{\n']
        ret.append('    nland : ')
        ret.append(repr(self.nland))
        ret.append(',\n    land : ')
        ret.append(repr(self.land))
        ret.append(',\n    coord : ')
        ret.append(repr(self.coord))
        ret.append(',\n    h_infil_max : ')
        ret.append(repr(self.h_infil_max))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = [init_array_ga, init_array_scs]
    

@f90wrap.runtime.register_class("wrapping.friction_data")
class friction_data(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=friction_data)
    
    
    Defined at m_sw_mono.f90 lines 81-86
    
    """
    def __init__(self, handle=None):
        """
        self = Friction_Data()
        
        
        Defined at m_sw_mono.f90 lines 81-86
        
        
        Returns
        -------
        this : Friction_Data
        	Object to be constructed
        
        
        Automatically generated constructor for friction_data
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_friction_data_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Friction_Data
        
        
        Defined at m_sw_mono.f90 lines 81-86
        
        Parameters
        ----------
        this : Friction_Data
        	Object to be destructed
        
        
        Automatically generated destructor for friction_data
        """
        if self._alloc:
            _wrapping.f90wrap_friction_data_finalise(this=self._handle)
    
    @property
    def nland(self):
        """
        Element nland ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 83
        
        """
        return _wrapping.f90wrap_friction_data__get__nland(self._handle)
    
    @nland.setter
    def nland(self, nland):
        _wrapping.f90wrap_friction_data__set__nland(self._handle, nland)
    
    @property
    def manning(self):
        """
        Element manning ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 84
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_friction_data__array__manning(self._handle)
        if array_handle in self._arrays:
            manning = self._arrays[array_handle]
        else:
            manning = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_friction_data__array__manning)
            self._arrays[array_handle] = manning
        return manning
    
    @manning.setter
    def manning(self, manning):
        self.manning[...] = manning
    
    @property
    def manning_beta(self):
        """
        Element manning_beta ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 85
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_friction_data__array__manning_beta(self._handle)
        if array_handle in self._arrays:
            manning_beta = self._arrays[array_handle]
        else:
            manning_beta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_friction_data__array__manning_beta)
            self._arrays[array_handle] = manning_beta
        return manning_beta
    
    @manning_beta.setter
    def manning_beta(self, manning_beta):
        self.manning_beta[...] = manning_beta
    
    @property
    def land(self):
        """
        Element land ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 86
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_friction_data__array__land(self._handle)
        if array_handle in self._arrays:
            land = self._arrays[array_handle]
        else:
            land = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_friction_data__array__land)
            self._arrays[array_handle] = land
        return land
    
    @land.setter
    def land(self, land):
        self.land[...] = land
    
    def __str__(self):
        ret = ['<friction_data>{\n']
        ret.append('    nland : ')
        ret.append(repr(self.nland))
        ret.append(',\n    manning : ')
        ret.append(repr(self.manning))
        ret.append(',\n    manning_beta : ')
        ret.append(repr(self.manning_beta))
        ret.append(',\n    land : ')
        ret.append(repr(self.land))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.porosity_data")
class porosity_data(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=porosity_data)
    
    
    Defined at m_sw_mono.f90 lines 91-95
    
    """
    def __init__(self, handle=None):
        """
        self = Porosity_Data()
        
        
        Defined at m_sw_mono.f90 lines 91-95
        
        
        Returns
        -------
        this : Porosity_Data
        	Object to be constructed
        
        
        Automatically generated constructor for porosity_data
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_porosity_data_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Porosity_Data
        
        
        Defined at m_sw_mono.f90 lines 91-95
        
        Parameters
        ----------
        this : Porosity_Data
        	Object to be destructed
        
        
        Automatically generated destructor for porosity_data
        """
        if self._alloc:
            _wrapping.f90wrap_porosity_data_finalise(this=self._handle)
    
    @property
    def nland(self):
        """
        Element nland ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 93
        
        """
        return _wrapping.f90wrap_porosity_data__get__nland(self._handle)
    
    @nland.setter
    def nland(self, nland):
        _wrapping.f90wrap_porosity_data__set__nland(self._handle, nland)
    
    @property
    def sp(self):
        """
        Element sp ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 94
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_porosity_data__array__sp(self._handle)
        if array_handle in self._arrays:
            sp = self._arrays[array_handle]
        else:
            sp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_porosity_data__array__sp)
            self._arrays[array_handle] = sp
        return sp
    
    @sp.setter
    def sp(self, sp):
        self.sp[...] = sp
    
    @property
    def land(self):
        """
        Element land ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 95
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_porosity_data__array__land(self._handle)
        if array_handle in self._arrays:
            land = self._arrays[array_handle]
        else:
            land = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_porosity_data__array__land)
            self._arrays[array_handle] = land
        return land
    
    @land.setter
    def land(self, land):
        self.land[...] = land
    
    def __str__(self):
        ret = ['<porosity_data>{\n']
        ret.append('    nland : ')
        ret.append(repr(self.nland))
        ret.append(',\n    sp : ')
        ret.append(repr(self.sp))
        ret.append(',\n    land : ')
        ret.append(repr(self.land))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.param_model")
class param_model(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=param_model)
    
    
    Defined at m_sw_mono.f90 lines 101-103
    
    """
    def __init__(self, handle=None):
        """
        self = Param_Model()
        
        
        Defined at m_sw_mono.f90 lines 101-103
        
        
        Returns
        -------
        this : Param_Model
        	Object to be constructed
        
        
        Automatically generated constructor for param_model
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_param_model_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Param_Model
        
        
        Defined at m_sw_mono.f90 lines 101-103
        
        Parameters
        ----------
        this : Param_Model
        	Object to be destructed
        
        
        Automatically generated destructor for param_model
        """
        if self._alloc:
            _wrapping.f90wrap_param_model_finalise(this=self._handle)
    
    @property
    def bathy_cell(self):
        """
        Element bathy_cell ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 103
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_param_model__array__bathy_cell(self._handle)
        if array_handle in self._arrays:
            bathy_cell = self._arrays[array_handle]
        else:
            bathy_cell = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_param_model__array__bathy_cell)
            self._arrays[array_handle] = bathy_cell
        return bathy_cell
    
    @bathy_cell.setter
    def bathy_cell(self, bathy_cell):
        self.bathy_cell[...] = bathy_cell
    
    def __str__(self):
        ret = ['<param_model>{\n']
        ret.append('    bathy_cell : ')
        ret.append(repr(self.bathy_cell))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.gr4")
class gr4(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=gr4)
    
    
    Defined at m_sw_mono.f90 lines 110-115
    
    """
    def __init__(self, handle=None):
        """
        self = Gr4()
        
        
        Defined at m_sw_mono.f90 lines 110-115
        
        
        Returns
        -------
        this : Gr4
        	Object to be constructed
        
        
        Automatically generated constructor for gr4
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_gr4_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Gr4
        
        
        Defined at m_sw_mono.f90 lines 110-115
        
        Parameters
        ----------
        this : Gr4
        	Object to be destructed
        
        
        Automatically generated destructor for gr4
        """
        if self._alloc:
            _wrapping.f90wrap_gr4_finalise(this=self._handle)
    
    @property
    def cell_id(self):
        """
        Element cell_id ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 111
        
        """
        return _wrapping.f90wrap_gr4__get__cell_id(self._handle)
    
    @cell_id.setter
    def cell_id(self, cell_id):
        _wrapping.f90wrap_gr4__set__cell_id(self._handle, cell_id)
    
    @property
    def surf(self):
        """
        Element surf ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 112
        
        """
        return _wrapping.f90wrap_gr4__get__surf(self._handle)
    
    @surf.setter
    def surf(self, surf):
        _wrapping.f90wrap_gr4__set__surf(self._handle, surf)
    
    @property
    def t(self):
        """
        Element t ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 113
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_gr4__array__t(self._handle)
        if array_handle in self._arrays:
            t = self._arrays[array_handle]
        else:
            t = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_gr4__array__t)
            self._arrays[array_handle] = t
        return t
    
    @t.setter
    def t(self, t):
        self.t[...] = t
    
    @property
    def p(self):
        """
        Element p ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 113
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_gr4__array__p(self._handle)
        if array_handle in self._arrays:
            p = self._arrays[array_handle]
        else:
            p = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_gr4__array__p)
            self._arrays[array_handle] = p
        return p
    
    @p.setter
    def p(self, p):
        self.p[...] = p
    
    @property
    def e(self):
        """
        Element e ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 113
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_gr4__array__e(self._handle)
        if array_handle in self._arrays:
            e = self._arrays[array_handle]
        else:
            e = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_gr4__array__e)
            self._arrays[array_handle] = e
        return e
    
    @e.setter
    def e(self, e):
        self.e[...] = e
    
    @property
    def q(self):
        """
        Element q ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 113
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_gr4__array__q(self._handle)
        if array_handle in self._arrays:
            q = self._arrays[array_handle]
        else:
            q = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_gr4__array__q)
            self._arrays[array_handle] = q
        return q
    
    @q.setter
    def q(self, q):
        self.q[...] = q
    
    @property
    def p0(self):
        """
        Element p0 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 113
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_gr4__array__p0(self._handle)
        if array_handle in self._arrays:
            p0 = self._arrays[array_handle]
        else:
            p0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_gr4__array__p0)
            self._arrays[array_handle] = p0
        return p0
    
    @p0.setter
    def p0(self, p0):
        self.p0[...] = p0
    
    @property
    def e0(self):
        """
        Element e0 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 113
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_gr4__array__e0(self._handle)
        if array_handle in self._arrays:
            e0 = self._arrays[array_handle]
        else:
            e0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_gr4__array__e0)
            self._arrays[array_handle] = e0
        return e0
    
    @e0.setter
    def e0(self, e0):
        self.e0[...] = e0
    
    @property
    def params(self):
        """
        Element params ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 114
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_gr4__array__params(self._handle)
        if array_handle in self._arrays:
            params = self._arrays[array_handle]
        else:
            params = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_gr4__array__params)
            self._arrays[array_handle] = params
        return params
    
    @params.setter
    def params(self, params):
        self.params[...] = params
    
    @property
    def state(self):
        """
        Element state ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 115
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_gr4__array__state(self._handle)
        if array_handle in self._arrays:
            state = self._arrays[array_handle]
        else:
            state = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_gr4__array__state)
            self._arrays[array_handle] = state
        return state
    
    @state.setter
    def state(self, state):
        self.state[...] = state
    
    def __str__(self):
        ret = ['<gr4>{\n']
        ret.append('    cell_id : ')
        ret.append(repr(self.cell_id))
        ret.append(',\n    surf : ')
        ret.append(repr(self.surf))
        ret.append(',\n    t : ')
        ret.append(repr(self.t))
        ret.append(',\n    p : ')
        ret.append(repr(self.p))
        ret.append(',\n    e : ')
        ret.append(repr(self.e))
        ret.append(',\n    q : ')
        ret.append(repr(self.q))
        ret.append(',\n    p0 : ')
        ret.append(repr(self.p0))
        ret.append(',\n    e0 : ')
        ret.append(repr(self.e0))
        ret.append(',\n    params : ')
        ret.append(repr(self.params))
        ret.append(',\n    state : ')
        ret.append(repr(self.state))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.hydrograph")
class hydrograph(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=hydrograph)
    
    
    Defined at m_sw_mono.f90 lines 118-120
    
    """
    def __init__(self, handle=None):
        """
        self = Hydrograph()
        
        
        Defined at m_sw_mono.f90 lines 118-120
        
        
        Returns
        -------
        this : Hydrograph
        	Object to be constructed
        
        
        Automatically generated constructor for hydrograph
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_hydrograph_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Hydrograph
        
        
        Defined at m_sw_mono.f90 lines 118-120
        
        Parameters
        ----------
        this : Hydrograph
        	Object to be destructed
        
        
        Automatically generated destructor for hydrograph
        """
        if self._alloc:
            _wrapping.f90wrap_hydrograph_finalise(this=self._handle)
    
    @property
    def group(self):
        """
        Element group ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 119
        
        """
        return _wrapping.f90wrap_hydrograph__get__group(self._handle)
    
    @group.setter
    def group(self, group):
        _wrapping.f90wrap_hydrograph__set__group(self._handle, group)
    
    @property
    def t(self):
        """
        Element t ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 120
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_hydrograph__array__t(self._handle)
        if array_handle in self._arrays:
            t = self._arrays[array_handle]
        else:
            t = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_hydrograph__array__t)
            self._arrays[array_handle] = t
        return t
    
    @t.setter
    def t(self, t):
        self.t[...] = t
    
    @property
    def q(self):
        """
        Element q ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 120
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_hydrograph__array__q(self._handle)
        if array_handle in self._arrays:
            q = self._arrays[array_handle]
        else:
            q = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_hydrograph__array__q)
            self._arrays[array_handle] = q
        return q
    
    @q.setter
    def q(self, q):
        self.q[...] = q
    
    def __str__(self):
        ret = ['<hydrograph>{\n']
        ret.append('    group : ')
        ret.append(repr(self.group))
        ret.append(',\n    t : ')
        ret.append(repr(self.t))
        ret.append(',\n    q : ')
        ret.append(repr(self.q))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.hpresc")
class hpresc(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=hpresc)
    
    
    Defined at m_sw_mono.f90 lines 122-124
    
    """
    def __init__(self, handle=None):
        """
        self = Hpresc()
        
        
        Defined at m_sw_mono.f90 lines 122-124
        
        
        Returns
        -------
        this : Hpresc
        	Object to be constructed
        
        
        Automatically generated constructor for hpresc
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_hpresc_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Hpresc
        
        
        Defined at m_sw_mono.f90 lines 122-124
        
        Parameters
        ----------
        this : Hpresc
        	Object to be destructed
        
        
        Automatically generated destructor for hpresc
        """
        if self._alloc:
            _wrapping.f90wrap_hpresc_finalise(this=self._handle)
    
    @property
    def group(self):
        """
        Element group ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 123
        
        """
        return _wrapping.f90wrap_hpresc__get__group(self._handle)
    
    @group.setter
    def group(self, group):
        _wrapping.f90wrap_hpresc__set__group(self._handle, group)
    
    @property
    def t(self):
        """
        Element t ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 124
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_hpresc__array__t(self._handle)
        if array_handle in self._arrays:
            t = self._arrays[array_handle]
        else:
            t = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_hpresc__array__t)
            self._arrays[array_handle] = t
        return t
    
    @t.setter
    def t(self, t):
        self.t[...] = t
    
    @property
    def h(self):
        """
        Element h ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 124
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_hpresc__array__h(self._handle)
        if array_handle in self._arrays:
            h = self._arrays[array_handle]
        else:
            h = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_hpresc__array__h)
            self._arrays[array_handle] = h
        return h
    
    @h.setter
    def h(self, h):
        self.h[...] = h
    
    def __str__(self):
        ret = ['<hpresc>{\n']
        ret.append('    group : ')
        ret.append(repr(self.group))
        ret.append(',\n    t : ')
        ret.append(repr(self.t))
        ret.append(',\n    h : ')
        ret.append(repr(self.h))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.zspresc")
class zspresc(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=zspresc)
    
    
    Defined at m_sw_mono.f90 lines 126-128
    
    """
    def __init__(self, handle=None):
        """
        self = Zspresc()
        
        
        Defined at m_sw_mono.f90 lines 126-128
        
        
        Returns
        -------
        this : Zspresc
        	Object to be constructed
        
        
        Automatically generated constructor for zspresc
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_zspresc_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Zspresc
        
        
        Defined at m_sw_mono.f90 lines 126-128
        
        Parameters
        ----------
        this : Zspresc
        	Object to be destructed
        
        
        Automatically generated destructor for zspresc
        """
        if self._alloc:
            _wrapping.f90wrap_zspresc_finalise(this=self._handle)
    
    @property
    def group(self):
        """
        Element group ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 127
        
        """
        return _wrapping.f90wrap_zspresc__get__group(self._handle)
    
    @group.setter
    def group(self, group):
        _wrapping.f90wrap_zspresc__set__group(self._handle, group)
    
    @property
    def t(self):
        """
        Element t ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 128
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_zspresc__array__t(self._handle)
        if array_handle in self._arrays:
            t = self._arrays[array_handle]
        else:
            t = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_zspresc__array__t)
            self._arrays[array_handle] = t
        return t
    
    @t.setter
    def t(self, t):
        self.t[...] = t
    
    @property
    def z(self):
        """
        Element z ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 128
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_zspresc__array__z(self._handle)
        if array_handle in self._arrays:
            z = self._arrays[array_handle]
        else:
            z = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_zspresc__array__z)
            self._arrays[array_handle] = z
        return z
    
    @z.setter
    def z(self, z):
        self.z[...] = z
    
    def __str__(self):
        ret = ['<zspresc>{\n']
        ret.append('    group : ')
        ret.append(repr(self.group))
        ret.append(',\n    t : ')
        ret.append(repr(self.t))
        ret.append(',\n    z : ')
        ret.append(repr(self.z))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.ratcurve")
class ratcurve(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=ratcurve)
    
    
    Defined at m_sw_mono.f90 lines 131-134
    
    """
    def __init__(self, handle=None):
        """
        self = Ratcurve()
        
        
        Defined at m_sw_mono.f90 lines 131-134
        
        
        Returns
        -------
        this : Ratcurve
        	Object to be constructed
        
        
        Automatically generated constructor for ratcurve
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_ratcurve_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Ratcurve
        
        
        Defined at m_sw_mono.f90 lines 131-134
        
        Parameters
        ----------
        this : Ratcurve
        	Object to be destructed
        
        
        Automatically generated destructor for ratcurve
        """
        if self._alloc:
            _wrapping.f90wrap_ratcurve_finalise(this=self._handle)
    
    @property
    def group(self):
        """
        Element group ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 132
        
        """
        return _wrapping.f90wrap_ratcurve__get__group(self._handle)
    
    @group.setter
    def group(self, group):
        _wrapping.f90wrap_ratcurve__set__group(self._handle, group)
    
    @property
    def h(self):
        """
        Element h ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 133
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_ratcurve__array__h(self._handle)
        if array_handle in self._arrays:
            h = self._arrays[array_handle]
        else:
            h = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_ratcurve__array__h)
            self._arrays[array_handle] = h
        return h
    
    @h.setter
    def h(self, h):
        self.h[...] = h
    
    @property
    def q(self):
        """
        Element q ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 133
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_ratcurve__array__q(self._handle)
        if array_handle in self._arrays:
            q = self._arrays[array_handle]
        else:
            q = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_ratcurve__array__q)
            self._arrays[array_handle] = q
        return q
    
    @q.setter
    def q(self, q):
        self.q[...] = q
    
    @property
    def z_rat_ref(self):
        """
        Element z_rat_ref ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 134
        
        """
        return _wrapping.f90wrap_ratcurve__get__z_rat_ref(self._handle)
    
    @z_rat_ref.setter
    def z_rat_ref(self, z_rat_ref):
        _wrapping.f90wrap_ratcurve__set__z_rat_ref(self._handle, z_rat_ref)
    
    @property
    def zout(self):
        """
        Element zout ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 134
        
        """
        return _wrapping.f90wrap_ratcurve__get__zout(self._handle)
    
    @zout.setter
    def zout(self, zout):
        _wrapping.f90wrap_ratcurve__set__zout(self._handle, zout)
    
    @property
    def c1(self):
        """
        Element c1 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 134
        
        """
        return _wrapping.f90wrap_ratcurve__get__c1(self._handle)
    
    @c1.setter
    def c1(self, c1):
        _wrapping.f90wrap_ratcurve__set__c1(self._handle, c1)
    
    @property
    def c2(self):
        """
        Element c2 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 134
        
        """
        return _wrapping.f90wrap_ratcurve__get__c2(self._handle)
    
    @c2.setter
    def c2(self, c2):
        _wrapping.f90wrap_ratcurve__set__c2(self._handle, c2)
    
    @property
    def pow(self):
        """
        Element pow ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 134
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_ratcurve__array__pow(self._handle)
        if array_handle in self._arrays:
            pow = self._arrays[array_handle]
        else:
            pow = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_ratcurve__array__pow)
            self._arrays[array_handle] = pow
        return pow
    
    @pow.setter
    def pow(self, pow):
        self.pow[...] = pow
    
    def __str__(self):
        ret = ['<ratcurve>{\n']
        ret.append('    group : ')
        ret.append(repr(self.group))
        ret.append(',\n    h : ')
        ret.append(repr(self.h))
        ret.append(',\n    q : ')
        ret.append(repr(self.q))
        ret.append(',\n    z_rat_ref : ')
        ret.append(repr(self.z_rat_ref))
        ret.append(',\n    zout : ')
        ret.append(repr(self.zout))
        ret.append(',\n    c1 : ')
        ret.append(repr(self.c1))
        ret.append(',\n    c2 : ')
        ret.append(repr(self.c2))
        ret.append(',\n    pow : ')
        ret.append(repr(self.pow))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.rain")
class rain(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=rain)
    
    
    Defined at m_sw_mono.f90 lines 136-144
    
    """
    def __init__(self, handle=None):
        """
        self = Rain()
        
        
        Defined at m_sw_mono.f90 lines 136-144
        
        
        Returns
        -------
        this : Rain
        	Object to be constructed
        
        
        Automatically generated constructor for rain
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_rain_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Rain
        
        
        Defined at m_sw_mono.f90 lines 136-144
        
        Parameters
        ----------
        this : Rain
        	Object to be destructed
        
        
        Automatically generated destructor for rain
        """
        if self._alloc:
            _wrapping.f90wrap_rain_finalise(this=self._handle)
    
    @property
    def x_min(self):
        """
        Element x_min ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 137
        
        """
        return _wrapping.f90wrap_rain__get__x_min(self._handle)
    
    @x_min.setter
    def x_min(self, x_min):
        _wrapping.f90wrap_rain__set__x_min(self._handle, x_min)
    
    @property
    def x_max(self):
        """
        Element x_max ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 138
        
        """
        return _wrapping.f90wrap_rain__get__x_max(self._handle)
    
    @x_max.setter
    def x_max(self, x_max):
        _wrapping.f90wrap_rain__set__x_max(self._handle, x_max)
    
    @property
    def y_min(self):
        """
        Element y_min ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 139
        
        """
        return _wrapping.f90wrap_rain__get__y_min(self._handle)
    
    @y_min.setter
    def y_min(self, y_min):
        _wrapping.f90wrap_rain__set__y_min(self._handle, y_min)
    
    @property
    def y_max(self):
        """
        Element y_max ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 140
        
        """
        return _wrapping.f90wrap_rain__get__y_max(self._handle)
    
    @y_max.setter
    def y_max(self, y_max):
        _wrapping.f90wrap_rain__set__y_max(self._handle, y_max)
    
    @property
    def tile_index(self):
        """
        Element tile_index ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 141
        
        """
        return _wrapping.f90wrap_rain__get__tile_index(self._handle)
    
    @tile_index.setter
    def tile_index(self, tile_index):
        _wrapping.f90wrap_rain__set__tile_index(self._handle, tile_index)
    
    @property
    def t(self):
        """
        Element t ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 142
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_rain__array__t(self._handle)
        if array_handle in self._arrays:
            t = self._arrays[array_handle]
        else:
            t = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_rain__array__t)
            self._arrays[array_handle] = t
        return t
    
    @t.setter
    def t(self, t):
        self.t[...] = t
    
    @property
    def q(self):
        """
        Element q ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 142
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_rain__array__q(self._handle)
        if array_handle in self._arrays:
            q = self._arrays[array_handle]
        else:
            q = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_rain__array__q)
            self._arrays[array_handle] = q
        return q
    
    @q.setter
    def q(self, q):
        self.q[...] = q
    
    @property
    def qin(self):
        """
        Element qin ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 143
        
        """
        return _wrapping.f90wrap_rain__get__qin(self._handle)
    
    @qin.setter
    def qin(self, qin):
        _wrapping.f90wrap_rain__set__qin(self._handle, qin)
    
    @property
    def cumul(self):
        """
        Element cumul ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 144
        
        """
        return _wrapping.f90wrap_rain__get__cumul(self._handle)
    
    @cumul.setter
    def cumul(self, cumul):
        _wrapping.f90wrap_rain__set__cumul(self._handle, cumul)
    
    def __str__(self):
        ret = ['<rain>{\n']
        ret.append('    x_min : ')
        ret.append(repr(self.x_min))
        ret.append(',\n    x_max : ')
        ret.append(repr(self.x_max))
        ret.append(',\n    y_min : ')
        ret.append(repr(self.y_min))
        ret.append(',\n    y_max : ')
        ret.append(repr(self.y_max))
        ret.append(',\n    tile_index : ')
        ret.append(repr(self.tile_index))
        ret.append(',\n    t : ')
        ret.append(repr(self.t))
        ret.append(',\n    q : ')
        ret.append(repr(self.q))
        ret.append(',\n    qin : ')
        ret.append(repr(self.qin))
        ret.append(',\n    cumul : ')
        ret.append(repr(self.cumul))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.bcs")
class bcs(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=bcs)
    
    
    Defined at m_sw_mono.f90 lines 147-160
    
    """
    def __init__(self, handle=None):
        """
        self = Bcs()
        
        
        Defined at m_sw_mono.f90 lines 147-160
        
        
        Returns
        -------
        this : Bcs
        	Object to be constructed
        
        
        Automatically generated constructor for bcs
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_bcs_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Bcs
        
        
        Defined at m_sw_mono.f90 lines 147-160
        
        Parameters
        ----------
        this : Bcs
        	Object to be destructed
        
        
        Automatically generated destructor for bcs
        """
        if self._alloc:
            _wrapping.f90wrap_bcs_finalise(this=self._handle)
    
    @property
    def nb(self):
        """
        Element nb ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 148
        
        """
        return _wrapping.f90wrap_bcs__get__nb(self._handle)
    
    @nb.setter
    def nb(self, nb):
        _wrapping.f90wrap_bcs__set__nb(self._handle, nb)
    
    @property
    def nb_in(self):
        """
        Element nb_in ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 148
        
        """
        return _wrapping.f90wrap_bcs__get__nb_in(self._handle)
    
    @nb_in.setter
    def nb_in(self, nb_in):
        _wrapping.f90wrap_bcs__set__nb_in(self._handle, nb_in)
    
    @property
    def nb_out(self):
        """
        Element nb_out ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 148
        
        """
        return _wrapping.f90wrap_bcs__get__nb_out(self._handle)
    
    @nb_out.setter
    def nb_out(self, nb_out):
        _wrapping.f90wrap_bcs__set__nb_out(self._handle, nb_out)
    
    @property
    def nb_rn(self):
        """
        Element nb_rn ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 148
        
        """
        return _wrapping.f90wrap_bcs__get__nb_rn(self._handle)
    
    @nb_rn.setter
    def nb_rn(self, nb_rn):
        _wrapping.f90wrap_bcs__set__nb_rn(self._handle, nb_rn)
    
    @property
    def nb_rn_t(self):
        """
        Element nb_rn_t ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 148
        
        """
        return _wrapping.f90wrap_bcs__get__nb_rn_t(self._handle)
    
    @nb_rn_t.setter
    def nb_rn_t(self, nb_rn_t):
        _wrapping.f90wrap_bcs__set__nb_rn_t(self._handle, nb_rn_t)
    
    @property
    def nb_gr4in(self):
        """
        Element nb_gr4in ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 148
        
        """
        return _wrapping.f90wrap_bcs__get__nb_gr4in(self._handle)
    
    @nb_gr4in.setter
    def nb_gr4in(self, nb_gr4in):
        _wrapping.f90wrap_bcs__set__nb_gr4in(self._handle, nb_gr4in)
    
    @property
    def typ(self):
        """
        Element typ ftype=character(len=lchar) pytype=str
        
        
        Defined at m_sw_mono.f90 line 149
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_bcs__array__typ(self._handle)
        if array_handle in self._arrays:
            typ = self._arrays[array_handle]
        else:
            typ = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_bcs__array__typ)
            self._arrays[array_handle] = typ
        return typ
    
    @typ.setter
    def typ(self, typ):
        self.typ[...] = typ
    
    @property
    def grpf(self):
        """
        Element grpf ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 150
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_bcs__array__grpf(self._handle)
        if array_handle in self._arrays:
            grpf = self._arrays[array_handle]
        else:
            grpf = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_bcs__array__grpf)
            self._arrays[array_handle] = grpf
        return grpf
    
    @grpf.setter
    def grpf(self, grpf):
        self.grpf[...] = grpf
    
    @property
    def inflow(self):
        """
        Element inflow ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 151
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_bcs__array__inflow(self._handle)
        if array_handle in self._arrays:
            inflow = self._arrays[array_handle]
        else:
            inflow = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_bcs__array__inflow)
            self._arrays[array_handle] = inflow
        return inflow
    
    @inflow.setter
    def inflow(self, inflow):
        self.inflow[...] = inflow
    
    @property
    def outflow(self):
        """
        Element outflow ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 152
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_bcs__array__outflow(self._handle)
        if array_handle in self._arrays:
            outflow = self._arrays[array_handle]
        else:
            outflow = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_bcs__array__outflow)
            self._arrays[array_handle] = outflow
        return outflow
    
    @outflow.setter
    def outflow(self, outflow):
        self.outflow[...] = outflow
    
    def init_array_hyd(self):
        self.hyd = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_bcs__array_getitem__hyd,
                                        _wrapping.f90wrap_bcs__array_setitem__hyd,
                                        _wrapping.f90wrap_bcs__array_len__hyd,
                                        """
        Element hyd ftype=type(hydrograph) pytype=Hydrograph
        
        
        Defined at m_sw_mono.f90 line 153
        
        """, hydrograph)
        return self.hyd
    
    def init_array_gr4(self):
        self.gr4 = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_bcs__array_getitem__gr4,
                                        _wrapping.f90wrap_bcs__array_setitem__gr4,
                                        _wrapping.f90wrap_bcs__array_len__gr4,
                                        """
        Element gr4 ftype=type(gr4) pytype=Gr4
        
        
        Defined at m_sw_mono.f90 line 154
        
        """, gr4)
        return self.gr4
    
    def init_array_rat(self):
        self.rat = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_bcs__array_getitem__rat,
                                        _wrapping.f90wrap_bcs__array_setitem__rat,
                                        _wrapping.f90wrap_bcs__array_len__rat,
                                        """
        Element rat ftype=type(ratcurve) pytype=Ratcurve
        
        
        Defined at m_sw_mono.f90 line 155
        
        """, ratcurve)
        return self.rat
    
    def init_array_hpresc(self):
        self.hpresc = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_bcs__array_getitem__hpresc,
                                        _wrapping.f90wrap_bcs__array_setitem__hpresc,
                                        _wrapping.f90wrap_bcs__array_len__hpresc,
                                        """
        Element hpresc ftype=type(hpresc) pytype=Hpresc
        
        
        Defined at m_sw_mono.f90 line 156
        
        """, hpresc)
        return self.hpresc
    
    def init_array_zspresc(self):
        self.zspresc = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_bcs__array_getitem__zspresc,
                                        _wrapping.f90wrap_bcs__array_setitem__zspresc,
                                        _wrapping.f90wrap_bcs__array_len__zspresc,
                                        """
        Element zspresc ftype=type(zspresc) pytype=Zspresc
        
        
        Defined at m_sw_mono.f90 line 157
        
        """, zspresc)
        return self.zspresc
    
    def init_array_rain(self):
        self.rain = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_bcs__array_getitem__rain,
                                        _wrapping.f90wrap_bcs__array_setitem__rain,
                                        _wrapping.f90wrap_bcs__array_len__rain,
                                        """
        Element rain ftype=type(rain) pytype=Rain
        
        
        Defined at m_sw_mono.f90 line 158
        
        """, rain)
        return self.rain
    
    @property
    def rain_land(self):
        """
        Element rain_land ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 159
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_bcs__array__rain_land(self._handle)
        if array_handle in self._arrays:
            rain_land = self._arrays[array_handle]
        else:
            rain_land = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_bcs__array__rain_land)
            self._arrays[array_handle] = rain_land
        return rain_land
    
    @rain_land.setter
    def rain_land(self, rain_land):
        self.rain_land[...] = rain_land
    
    @property
    def sum_mass_flux(self):
        """
        Element sum_mass_flux ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 160
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_bcs__array__sum_mass_flux(self._handle)
        if array_handle in self._arrays:
            sum_mass_flux = self._arrays[array_handle]
        else:
            sum_mass_flux = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_bcs__array__sum_mass_flux)
            self._arrays[array_handle] = sum_mass_flux
        return sum_mass_flux
    
    @sum_mass_flux.setter
    def sum_mass_flux(self, sum_mass_flux):
        self.sum_mass_flux[...] = sum_mass_flux
    
    def __str__(self):
        ret = ['<bcs>{\n']
        ret.append('    nb : ')
        ret.append(repr(self.nb))
        ret.append(',\n    nb_in : ')
        ret.append(repr(self.nb_in))
        ret.append(',\n    nb_out : ')
        ret.append(repr(self.nb_out))
        ret.append(',\n    nb_rn : ')
        ret.append(repr(self.nb_rn))
        ret.append(',\n    nb_rn_t : ')
        ret.append(repr(self.nb_rn_t))
        ret.append(',\n    nb_gr4in : ')
        ret.append(repr(self.nb_gr4in))
        ret.append(',\n    typ : ')
        ret.append(repr(self.typ))
        ret.append(',\n    grpf : ')
        ret.append(repr(self.grpf))
        ret.append(',\n    inflow : ')
        ret.append(repr(self.inflow))
        ret.append(',\n    outflow : ')
        ret.append(repr(self.outflow))
        ret.append(',\n    rain_land : ')
        ret.append(repr(self.rain_land))
        ret.append(',\n    sum_mass_flux : ')
        ret.append(repr(self.sum_mass_flux))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = [init_array_hyd, init_array_gr4, init_array_rat, \
        init_array_hpresc, init_array_zspresc, init_array_rain]
    

@f90wrap.runtime.register_class("wrapping.station_obs")
class station_obs(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=station_obs)
    
    
    Defined at m_sw_mono.f90 lines 168-177
    
    """
    def __init__(self, handle=None):
        """
        self = Station_Obs()
        
        
        Defined at m_sw_mono.f90 lines 168-177
        
        
        Returns
        -------
        this : Station_Obs
        	Object to be constructed
        
        
        Automatically generated constructor for station_obs
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_station_obs_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Station_Obs
        
        
        Defined at m_sw_mono.f90 lines 168-177
        
        Parameters
        ----------
        this : Station_Obs
        	Object to be destructed
        
        
        Automatically generated destructor for station_obs
        """
        if self._alloc:
            _wrapping.f90wrap_station_obs_finalise(this=self._handle)
    
    @property
    def weight(self):
        """
        Element weight ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 170
        
        """
        return _wrapping.f90wrap_station_obs__get__weight(self._handle)
    
    @weight.setter
    def weight(self, weight):
        _wrapping.f90wrap_station_obs__set__weight(self._handle, weight)
    
    @property
    def length(self):
        """
        Element length ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 171
        
        """
        return _wrapping.f90wrap_station_obs__get__length(self._handle)
    
    @length.setter
    def length(self, length):
        _wrapping.f90wrap_station_obs__set__length(self._handle, length)
    
    @property
    def dt_offset(self):
        """
        Element dt_offset ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 172
        
        """
        return _wrapping.f90wrap_station_obs__get__dt_offset(self._handle)
    
    @dt_offset.setter
    def dt_offset(self, dt_offset):
        _wrapping.f90wrap_station_obs__set__dt_offset(self._handle, dt_offset)
    
    @property
    def dt(self):
        """
        Element dt ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 173
        
        """
        return _wrapping.f90wrap_station_obs__get__dt(self._handle)
    
    @dt.setter
    def dt(self, dt):
        _wrapping.f90wrap_station_obs__set__dt(self._handle, dt)
    
    @property
    def dt_obs(self):
        """
        Element dt_obs ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 174
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obs__array__dt_obs(self._handle)
        if array_handle in self._arrays:
            dt_obs = self._arrays[array_handle]
        else:
            dt_obs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obs__array__dt_obs)
            self._arrays[array_handle] = dt_obs
        return dt_obs
    
    @dt_obs.setter
    def dt_obs(self, dt_obs):
        self.dt_obs[...] = dt_obs
    
    @property
    def ind_t(self):
        """
        Element ind_t ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 175
        
        """
        return _wrapping.f90wrap_station_obs__get__ind_t(self._handle)
    
    @ind_t.setter
    def ind_t(self, ind_t):
        _wrapping.f90wrap_station_obs__set__ind_t(self._handle, ind_t)
    
    @property
    def nb_dt(self):
        """
        Element nb_dt ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 176
        
        """
        return _wrapping.f90wrap_station_obs__get__nb_dt(self._handle)
    
    @nb_dt.setter
    def nb_dt(self, nb_dt):
        _wrapping.f90wrap_station_obs__set__nb_dt(self._handle, nb_dt)
    
    @property
    def t(self):
        """
        Element t ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 177
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obs__array__t(self._handle)
        if array_handle in self._arrays:
            t = self._arrays[array_handle]
        else:
            t = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obs__array__t)
            self._arrays[array_handle] = t
        return t
    
    @t.setter
    def t(self, t):
        self.t[...] = t
    
    @property
    def h(self):
        """
        Element h ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 177
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obs__array__h(self._handle)
        if array_handle in self._arrays:
            h = self._arrays[array_handle]
        else:
            h = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obs__array__h)
            self._arrays[array_handle] = h
        return h
    
    @h.setter
    def h(self, h):
        self.h[...] = h
    
    @property
    def u(self):
        """
        Element u ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 177
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obs__array__u(self._handle)
        if array_handle in self._arrays:
            u = self._arrays[array_handle]
        else:
            u = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obs__array__u)
            self._arrays[array_handle] = u
        return u
    
    @u.setter
    def u(self, u):
        self.u[...] = u
    
    @property
    def v(self):
        """
        Element v ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 177
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obs__array__v(self._handle)
        if array_handle in self._arrays:
            v = self._arrays[array_handle]
        else:
            v = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obs__array__v)
            self._arrays[array_handle] = v
        return v
    
    @v.setter
    def v(self, v):
        self.v[...] = v
    
    @property
    def q(self):
        """
        Element q ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 177
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obs__array__q(self._handle)
        if array_handle in self._arrays:
            q = self._arrays[array_handle]
        else:
            q = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obs__array__q)
            self._arrays[array_handle] = q
        return q
    
    @q.setter
    def q(self, q):
        self.q[...] = q
    
    @property
    def w(self):
        """
        Element w ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 177
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obs__array__w(self._handle)
        if array_handle in self._arrays:
            w = self._arrays[array_handle]
        else:
            w = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obs__array__w)
            self._arrays[array_handle] = w
        return w
    
    @w.setter
    def w(self, w):
        self.w[...] = w
    
    def __str__(self):
        ret = ['<station_obs>{\n']
        ret.append('    weight : ')
        ret.append(repr(self.weight))
        ret.append(',\n    length : ')
        ret.append(repr(self.length))
        ret.append(',\n    dt_offset : ')
        ret.append(repr(self.dt_offset))
        ret.append(',\n    dt : ')
        ret.append(repr(self.dt))
        ret.append(',\n    dt_obs : ')
        ret.append(repr(self.dt_obs))
        ret.append(',\n    ind_t : ')
        ret.append(repr(self.ind_t))
        ret.append(',\n    nb_dt : ')
        ret.append(repr(self.nb_dt))
        ret.append(',\n    t : ')
        ret.append(repr(self.t))
        ret.append(',\n    h : ')
        ret.append(repr(self.h))
        ret.append(',\n    u : ')
        ret.append(repr(self.u))
        ret.append(',\n    v : ')
        ret.append(repr(self.v))
        ret.append(',\n    q : ')
        ret.append(repr(self.q))
        ret.append(',\n    w : ')
        ret.append(repr(self.w))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.section_obs")
class section_obs(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=section_obs)
    
    
    Defined at m_sw_mono.f90 lines 180-184
    
    """
    def __init__(self, handle=None):
        """
        self = Section_Obs()
        
        
        Defined at m_sw_mono.f90 lines 180-184
        
        
        Returns
        -------
        this : Section_Obs
        	Object to be constructed
        
        
        Automatically generated constructor for section_obs
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_section_obs_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Section_Obs
        
        
        Defined at m_sw_mono.f90 lines 180-184
        
        Parameters
        ----------
        this : Section_Obs
        	Object to be destructed
        
        
        Automatically generated destructor for section_obs
        """
        if self._alloc:
            _wrapping.f90wrap_section_obs_finalise(this=self._handle)
    
    @property
    def dt(self):
        """
        Element dt ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 182
        
        """
        return _wrapping.f90wrap_section_obs__get__dt(self._handle)
    
    @dt.setter
    def dt(self, dt):
        _wrapping.f90wrap_section_obs__set__dt(self._handle, dt)
    
    @property
    def dx(self):
        """
        Element dx ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 182
        
        """
        return _wrapping.f90wrap_section_obs__get__dx(self._handle)
    
    @dx.setter
    def dx(self, dx):
        _wrapping.f90wrap_section_obs__set__dx(self._handle, dx)
    
    @property
    def t(self):
        """
        Element t ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 183
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_section_obs__array__t(self._handle)
        if array_handle in self._arrays:
            t = self._arrays[array_handle]
        else:
            t = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_section_obs__array__t)
            self._arrays[array_handle] = t
        return t
    
    @t.setter
    def t(self, t):
        self.t[...] = t
    
    @property
    def h(self):
        """
        Element h ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 183
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_section_obs__array__h(self._handle)
        if array_handle in self._arrays:
            h = self._arrays[array_handle]
        else:
            h = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_section_obs__array__h)
            self._arrays[array_handle] = h
        return h
    
    @h.setter
    def h(self, h):
        self.h[...] = h
    
    @property
    def u(self):
        """
        Element u ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 183
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_section_obs__array__u(self._handle)
        if array_handle in self._arrays:
            u = self._arrays[array_handle]
        else:
            u = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_section_obs__array__u)
            self._arrays[array_handle] = u
        return u
    
    @u.setter
    def u(self, u):
        self.u[...] = u
    
    @property
    def v(self):
        """
        Element v ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 183
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_section_obs__array__v(self._handle)
        if array_handle in self._arrays:
            v = self._arrays[array_handle]
        else:
            v = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_section_obs__array__v)
            self._arrays[array_handle] = v
        return v
    
    @v.setter
    def v(self, v):
        self.v[...] = v
    
    @property
    def q(self):
        """
        Element q ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 183
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_section_obs__array__q(self._handle)
        if array_handle in self._arrays:
            q = self._arrays[array_handle]
        else:
            q = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_section_obs__array__q)
            self._arrays[array_handle] = q
        return q
    
    @q.setter
    def q(self, q):
        self.q[...] = q
    
    def __str__(self):
        ret = ['<section_obs>{\n']
        ret.append('    dt : ')
        ret.append(repr(self.dt))
        ret.append(',\n    dx : ')
        ret.append(repr(self.dx))
        ret.append(',\n    t : ')
        ret.append(repr(self.t))
        ret.append(',\n    h : ')
        ret.append(repr(self.h))
        ret.append(',\n    u : ')
        ret.append(repr(self.u))
        ret.append(',\n    v : ')
        ret.append(repr(self.v))
        ret.append(',\n    q : ')
        ret.append(repr(self.q))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.station_obsQ")
class station_obsQ(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=station_obsq)
    
    
    Defined at m_sw_mono.f90 lines 186-196
    
    """
    def __init__(self, handle=None):
        """
        self = Station_Obsq()
        
        
        Defined at m_sw_mono.f90 lines 186-196
        
        
        Returns
        -------
        this : Station_Obsq
        	Object to be constructed
        
        
        Automatically generated constructor for station_obsq
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_station_obsq_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Station_Obsq
        
        
        Defined at m_sw_mono.f90 lines 186-196
        
        Parameters
        ----------
        this : Station_Obsq
        	Object to be destructed
        
        
        Automatically generated destructor for station_obsq
        """
        if self._alloc:
            _wrapping.f90wrap_station_obsq_finalise(this=self._handle)
    
    @property
    def weight(self):
        """
        Element weight ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 188
        
        """
        return _wrapping.f90wrap_station_obsq__get__weight(self._handle)
    
    @weight.setter
    def weight(self, weight):
        _wrapping.f90wrap_station_obsq__set__weight(self._handle, weight)
    
    @property
    def length(self):
        """
        Element length ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 189
        
        """
        return _wrapping.f90wrap_station_obsq__get__length(self._handle)
    
    @length.setter
    def length(self, length):
        _wrapping.f90wrap_station_obsq__set__length(self._handle, length)
    
    @property
    def dt_offset(self):
        """
        Element dt_offset ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 190
        
        """
        return _wrapping.f90wrap_station_obsq__get__dt_offset(self._handle)
    
    @dt_offset.setter
    def dt_offset(self, dt_offset):
        _wrapping.f90wrap_station_obsq__set__dt_offset(self._handle, dt_offset)
    
    @property
    def dt(self):
        """
        Element dt ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 191
        
        """
        return _wrapping.f90wrap_station_obsq__get__dt(self._handle)
    
    @dt.setter
    def dt(self, dt):
        _wrapping.f90wrap_station_obsq__set__dt(self._handle, dt)
    
    @property
    def dt_obs(self):
        """
        Element dt_obs ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 192
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obsq__array__dt_obs(self._handle)
        if array_handle in self._arrays:
            dt_obs = self._arrays[array_handle]
        else:
            dt_obs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obsq__array__dt_obs)
            self._arrays[array_handle] = dt_obs
        return dt_obs
    
    @dt_obs.setter
    def dt_obs(self, dt_obs):
        self.dt_obs[...] = dt_obs
    
    @property
    def ind_t(self):
        """
        Element ind_t ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 193
        
        """
        return _wrapping.f90wrap_station_obsq__get__ind_t(self._handle)
    
    @ind_t.setter
    def ind_t(self, ind_t):
        _wrapping.f90wrap_station_obsq__set__ind_t(self._handle, ind_t)
    
    @property
    def ind_bc(self):
        """
        Element ind_bc ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 194
        
        """
        return _wrapping.f90wrap_station_obsq__get__ind_bc(self._handle)
    
    @ind_bc.setter
    def ind_bc(self, ind_bc):
        _wrapping.f90wrap_station_obsq__set__ind_bc(self._handle, ind_bc)
    
    @property
    def nb_dt(self):
        """
        Element nb_dt ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 195
        
        """
        return _wrapping.f90wrap_station_obsq__get__nb_dt(self._handle)
    
    @nb_dt.setter
    def nb_dt(self, nb_dt):
        _wrapping.f90wrap_station_obsq__set__nb_dt(self._handle, nb_dt)
    
    @property
    def t(self):
        """
        Element t ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 196
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obsq__array__t(self._handle)
        if array_handle in self._arrays:
            t = self._arrays[array_handle]
        else:
            t = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obsq__array__t)
            self._arrays[array_handle] = t
        return t
    
    @t.setter
    def t(self, t):
        self.t[...] = t
    
    @property
    def h(self):
        """
        Element h ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 196
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obsq__array__h(self._handle)
        if array_handle in self._arrays:
            h = self._arrays[array_handle]
        else:
            h = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obsq__array__h)
            self._arrays[array_handle] = h
        return h
    
    @h.setter
    def h(self, h):
        self.h[...] = h
    
    @property
    def u(self):
        """
        Element u ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 196
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obsq__array__u(self._handle)
        if array_handle in self._arrays:
            u = self._arrays[array_handle]
        else:
            u = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obsq__array__u)
            self._arrays[array_handle] = u
        return u
    
    @u.setter
    def u(self, u):
        self.u[...] = u
    
    @property
    def v(self):
        """
        Element v ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 196
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obsq__array__v(self._handle)
        if array_handle in self._arrays:
            v = self._arrays[array_handle]
        else:
            v = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obsq__array__v)
            self._arrays[array_handle] = v
        return v
    
    @v.setter
    def v(self, v):
        self.v[...] = v
    
    @property
    def q(self):
        """
        Element q ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 196
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obsq__array__q(self._handle)
        if array_handle in self._arrays:
            q = self._arrays[array_handle]
        else:
            q = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obsq__array__q)
            self._arrays[array_handle] = q
        return q
    
    @q.setter
    def q(self, q):
        self.q[...] = q
    
    @property
    def w(self):
        """
        Element w ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 196
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_station_obsq__array__w(self._handle)
        if array_handle in self._arrays:
            w = self._arrays[array_handle]
        else:
            w = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_station_obsq__array__w)
            self._arrays[array_handle] = w
        return w
    
    @w.setter
    def w(self, w):
        self.w[...] = w
    
    def __str__(self):
        ret = ['<station_obsq>{\n']
        ret.append('    weight : ')
        ret.append(repr(self.weight))
        ret.append(',\n    length : ')
        ret.append(repr(self.length))
        ret.append(',\n    dt_offset : ')
        ret.append(repr(self.dt_offset))
        ret.append(',\n    dt : ')
        ret.append(repr(self.dt))
        ret.append(',\n    dt_obs : ')
        ret.append(repr(self.dt_obs))
        ret.append(',\n    ind_t : ')
        ret.append(repr(self.ind_t))
        ret.append(',\n    ind_bc : ')
        ret.append(repr(self.ind_bc))
        ret.append(',\n    nb_dt : ')
        ret.append(repr(self.nb_dt))
        ret.append(',\n    t : ')
        ret.append(repr(self.t))
        ret.append(',\n    h : ')
        ret.append(repr(self.h))
        ret.append(',\n    u : ')
        ret.append(repr(self.u))
        ret.append(',\n    v : ')
        ret.append(repr(self.v))
        ret.append(',\n    q : ')
        ret.append(repr(self.q))
        ret.append(',\n    w : ')
        ret.append(repr(self.w))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.soil_data")
class soil_data(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=soil_data)
    
    
    Defined at m_sw_mono.f90 lines 207-212
    
    """
    def __init__(self, handle=None):
        """
        self = Soil_Data()
        
        
        Defined at m_sw_mono.f90 lines 207-212
        
        
        Returns
        -------
        this : Soil_Data
        	Object to be constructed
        
        
        Automatically generated constructor for soil_data
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_soil_data_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Soil_Data
        
        
        Defined at m_sw_mono.f90 lines 207-212
        
        Parameters
        ----------
        this : Soil_Data
        	Object to be destructed
        
        
        Automatically generated destructor for soil_data
        """
        if self._alloc:
            _wrapping.f90wrap_soil_data_finalise(this=self._handle)
    
    @property
    def clay(self):
        """
        Element clay ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 208
        
        """
        return _wrapping.f90wrap_soil_data__get__clay(self._handle)
    
    @clay.setter
    def clay(self, clay):
        _wrapping.f90wrap_soil_data__set__clay(self._handle, clay)
    
    @property
    def silt(self):
        """
        Element silt ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 209
        
        """
        return _wrapping.f90wrap_soil_data__get__silt(self._handle)
    
    @silt.setter
    def silt(self, silt):
        _wrapping.f90wrap_soil_data__set__silt(self._handle, silt)
    
    @property
    def sand(self):
        """
        Element sand ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 210
        
        """
        return _wrapping.f90wrap_soil_data__get__sand(self._handle)
    
    @sand.setter
    def sand(self, sand):
        _wrapping.f90wrap_soil_data__set__sand(self._handle, sand)
    
    @property
    def thetas(self):
        """
        Element thetas ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 211
        
        """
        return _wrapping.f90wrap_soil_data__get__thetas(self._handle)
    
    @thetas.setter
    def thetas(self, thetas):
        _wrapping.f90wrap_soil_data__set__thetas(self._handle, thetas)
    
    @property
    def thetar(self):
        """
        Element thetar ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 212
        
        """
        return _wrapping.f90wrap_soil_data__get__thetar(self._handle)
    
    @thetar.setter
    def thetar(self, thetar):
        _wrapping.f90wrap_soil_data__set__thetar(self._handle, thetar)
    
    def __str__(self):
        ret = ['<soil_data>{\n']
        ret.append('    clay : ')
        ret.append(repr(self.clay))
        ret.append(',\n    silt : ')
        ret.append(repr(self.silt))
        ret.append(',\n    sand : ')
        ret.append(repr(self.sand))
        ret.append(',\n    thetas : ')
        ret.append(repr(self.thetas))
        ret.append(',\n    thetar : ')
        ret.append(repr(self.thetar))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.ptf_data")
class ptf_data(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=ptf_data)
    
    
    Defined at m_sw_mono.f90 lines 214-215
    
    """
    def __init__(self, handle=None):
        """
        self = Ptf_Data()
        
        
        Defined at m_sw_mono.f90 lines 214-215
        
        
        Returns
        -------
        this : Ptf_Data
        	Object to be constructed
        
        
        Automatically generated constructor for ptf_data
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_ptf_data_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Ptf_Data
        
        
        Defined at m_sw_mono.f90 lines 214-215
        
        Parameters
        ----------
        this : Ptf_Data
        	Object to be destructed
        
        
        Automatically generated destructor for ptf_data
        """
        if self._alloc:
            _wrapping.f90wrap_ptf_data_finalise(this=self._handle)
    
    @property
    def kappa(self):
        """
        Element kappa ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 215
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_ptf_data__array__kappa(self._handle)
        if array_handle in self._arrays:
            kappa = self._arrays[array_handle]
        else:
            kappa = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_ptf_data__array__kappa)
            self._arrays[array_handle] = kappa
        return kappa
    
    @kappa.setter
    def kappa(self, kappa):
        self.kappa[...] = kappa
    
    def __str__(self):
        ret = ['<ptf_data>{\n']
        ret.append('    kappa : ')
        ret.append(repr(self.kappa))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.surface_data")
class surface_data(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=surface_data)
    
    
    Defined at m_sw_mono.f90 lines 217-219
    
    """
    def __init__(self, handle=None):
        """
        self = Surface_Data()
        
        
        Defined at m_sw_mono.f90 lines 217-219
        
        
        Returns
        -------
        this : Surface_Data
        	Object to be constructed
        
        
        Automatically generated constructor for surface_data
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_surface_data_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Surface_Data
        
        
        Defined at m_sw_mono.f90 lines 217-219
        
        Parameters
        ----------
        this : Surface_Data
        	Object to be destructed
        
        
        Automatically generated destructor for surface_data
        """
        if self._alloc:
            _wrapping.f90wrap_surface_data_finalise(this=self._handle)
    
    @property
    def imperm(self):
        """
        Element imperm ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 218
        
        """
        return _wrapping.f90wrap_surface_data__get__imperm(self._handle)
    
    @imperm.setter
    def imperm(self, imperm):
        _wrapping.f90wrap_surface_data__set__imperm(self._handle, imperm)
    
    @property
    def dmax(self):
        """
        Element dmax ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 219
        
        """
        return _wrapping.f90wrap_surface_data__get__dmax(self._handle)
    
    @dmax.setter
    def dmax(self, dmax):
        _wrapping.f90wrap_surface_data__set__dmax(self._handle, dmax)
    
    def __str__(self):
        ret = ['<surface_data>{\n']
        ret.append('    imperm : ')
        ret.append(repr(self.imperm))
        ret.append(',\n    dmax : ')
        ret.append(repr(self.dmax))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.structure_data")
class structure_data(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=structure_data)
    
    
    Defined at m_sw_mono.f90 lines 221-227
    
    """
    def __init__(self, handle=None):
        """
        self = Structure_Data()
        
        
        Defined at m_sw_mono.f90 lines 221-227
        
        
        Returns
        -------
        this : Structure_Data
        	Object to be constructed
        
        
        Automatically generated constructor for structure_data
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_structure_data_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Structure_Data
        
        
        Defined at m_sw_mono.f90 lines 221-227
        
        Parameters
        ----------
        this : Structure_Data
        	Object to be destructed
        
        
        Automatically generated destructor for structure_data
        """
        if self._alloc:
            _wrapping.f90wrap_structure_data_finalise(this=self._handle)
    
    @property
    def name(self):
        """
        Element name ftype=character(len=lchar) pytype=str
        
        
        Defined at m_sw_mono.f90 line 222
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_structure_data__array__name(self._handle)
        if array_handle in self._arrays:
            name = self._arrays[array_handle]
        else:
            name = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_structure_data__array__name)
            self._arrays[array_handle] = name
        return name
    
    @name.setter
    def name(self, name):
        self.name[...] = name
    
    @property
    def c1(self):
        """
        Element c1 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 223
        
        """
        return _wrapping.f90wrap_structure_data__get__c1(self._handle)
    
    @c1.setter
    def c1(self, c1):
        _wrapping.f90wrap_structure_data__set__c1(self._handle, c1)
    
    @property
    def c2(self):
        """
        Element c2 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 224
        
        """
        return _wrapping.f90wrap_structure_data__get__c2(self._handle)
    
    @c2.setter
    def c2(self, c2):
        _wrapping.f90wrap_structure_data__set__c2(self._handle, c2)
    
    @property
    def c3(self):
        """
        Element c3 ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 225
        
        """
        return _wrapping.f90wrap_structure_data__get__c3(self._handle)
    
    @c3.setter
    def c3(self, c3):
        _wrapping.f90wrap_structure_data__set__c3(self._handle, c3)
    
    @property
    def true_x(self):
        """
        Element true_x ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 226
        
        """
        return _wrapping.f90wrap_structure_data__get__true_x(self._handle)
    
    @true_x.setter
    def true_x(self, true_x):
        _wrapping.f90wrap_structure_data__set__true_x(self._handle, true_x)
    
    @property
    def true_y(self):
        """
        Element true_y ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 227
        
        """
        return _wrapping.f90wrap_structure_data__get__true_y(self._handle)
    
    @true_y.setter
    def true_y(self, true_y):
        _wrapping.f90wrap_structure_data__set__true_y(self._handle, true_y)
    
    def __str__(self):
        ret = ['<structure_data>{\n']
        ret.append('    name : ')
        ret.append(repr(self.name))
        ret.append(',\n    c1 : ')
        ret.append(repr(self.c1))
        ret.append(',\n    c2 : ')
        ret.append(repr(self.c2))
        ret.append(',\n    c3 : ')
        ret.append(repr(self.c3))
        ret.append(',\n    true_x : ')
        ret.append(repr(self.true_x))
        ret.append(',\n    true_y : ')
        ret.append(repr(self.true_y))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.input_data")
class input_data(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=input_data)
    
    
    Defined at m_sw_mono.f90 lines 229-241
    
    """
    def __init__(self, handle=None):
        """
        self = Input_Data()
        
        
        Defined at m_sw_mono.f90 lines 229-241
        
        
        Returns
        -------
        this : Input_Data
        	Object to be constructed
        
        
        Automatically generated constructor for input_data
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_input_data_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Input_Data
        
        
        Defined at m_sw_mono.f90 lines 229-241
        
        Parameters
        ----------
        this : Input_Data
        	Object to be destructed
        
        
        Automatically generated destructor for input_data
        """
        if self._alloc:
            _wrapping.f90wrap_input_data_finalise(this=self._handle)
    
    def init_array_soil(self):
        self.soil = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_input_data__array_getitem__soil,
                                        _wrapping.f90wrap_input_data__array_setitem__soil,
                                        _wrapping.f90wrap_input_data__array_len__soil,
                                        """
        Element soil ftype=type(soil_data) pytype=Soil_Data
        
        
        Defined at m_sw_mono.f90 line 230
        
        """, soil_data)
        return self.soil
    
    @property
    def soil_nland(self):
        """
        Element soil_nland ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 231
        
        """
        return _wrapping.f90wrap_input_data__get__soil_nland(self._handle)
    
    @soil_nland.setter
    def soil_nland(self, soil_nland):
        _wrapping.f90wrap_input_data__set__soil_nland(self._handle, soil_nland)
    
    @property
    def soil_land(self):
        """
        Element soil_land ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 232
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_input_data__array__soil_land(self._handle)
        if array_handle in self._arrays:
            soil_land = self._arrays[array_handle]
        else:
            soil_land = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_input_data__array__soil_land)
            self._arrays[array_handle] = soil_land
        return soil_land
    
    @soil_land.setter
    def soil_land(self, soil_land):
        self.soil_land[...] = soil_land
    
    def init_array_ptf(self):
        self.ptf = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_input_data__array_getitem__ptf,
                                        _wrapping.f90wrap_input_data__array_setitem__ptf,
                                        _wrapping.f90wrap_input_data__array_len__ptf,
                                        """
        Element ptf ftype=type(ptf_data) pytype=Ptf_Data
        
        
        Defined at m_sw_mono.f90 line 233
        
        """, ptf_data)
        return self.ptf
    
    @property
    def ptf_nland(self):
        """
        Element ptf_nland ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 234
        
        """
        return _wrapping.f90wrap_input_data__get__ptf_nland(self._handle)
    
    @ptf_nland.setter
    def ptf_nland(self, ptf_nland):
        _wrapping.f90wrap_input_data__set__ptf_nland(self._handle, ptf_nland)
    
    @property
    def ptf_land(self):
        """
        Element ptf_land ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 235
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_input_data__array__ptf_land(self._handle)
        if array_handle in self._arrays:
            ptf_land = self._arrays[array_handle]
        else:
            ptf_land = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_input_data__array__ptf_land)
            self._arrays[array_handle] = ptf_land
        return ptf_land
    
    @ptf_land.setter
    def ptf_land(self, ptf_land):
        self.ptf_land[...] = ptf_land
    
    def init_array_surf(self):
        self.surf = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_input_data__array_getitem__surf,
                                        _wrapping.f90wrap_input_data__array_setitem__surf,
                                        _wrapping.f90wrap_input_data__array_len__surf,
                                        """
        Element surf ftype=type(surface_data) pytype=Surface_Data
        
        
        Defined at m_sw_mono.f90 line 236
        
        """, surface_data)
        return self.surf
    
    @property
    def surf_nland(self):
        """
        Element surf_nland ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 237
        
        """
        return _wrapping.f90wrap_input_data__get__surf_nland(self._handle)
    
    @surf_nland.setter
    def surf_nland(self, surf_nland):
        _wrapping.f90wrap_input_data__set__surf_nland(self._handle, surf_nland)
    
    @property
    def surf_land(self):
        """
        Element surf_land ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 238
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_input_data__array__surf_land(self._handle)
        if array_handle in self._arrays:
            surf_land = self._arrays[array_handle]
        else:
            surf_land = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_input_data__array__surf_land)
            self._arrays[array_handle] = surf_land
        return surf_land
    
    @surf_land.setter
    def surf_land(self, surf_land):
        self.surf_land[...] = surf_land
    
    def init_array_structures(self):
        self.structures = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_input_data__array_getitem__structures,
                                        _wrapping.f90wrap_input_data__array_setitem__structures,
                                        _wrapping.f90wrap_input_data__array_len__structures,
                                        """
        Element structures ftype=type(structure_data) pytype=Structure_Data
        
        
        Defined at m_sw_mono.f90 line 239
        
        """, structure_data)
        return self.structures
    
    @property
    def struct_nland(self):
        """
        Element struct_nland ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 240
        
        """
        return _wrapping.f90wrap_input_data__get__struct_nland(self._handle)
    
    @struct_nland.setter
    def struct_nland(self, struct_nland):
        _wrapping.f90wrap_input_data__set__struct_nland(self._handle, struct_nland)
    
    @property
    def struct_land(self):
        """
        Element struct_land ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 241
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_input_data__array__struct_land(self._handle)
        if array_handle in self._arrays:
            struct_land = self._arrays[array_handle]
        else:
            struct_land = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_input_data__array__struct_land)
            self._arrays[array_handle] = struct_land
        return struct_land
    
    @struct_land.setter
    def struct_land(self, struct_land):
        self.struct_land[...] = struct_land
    
    def __str__(self):
        ret = ['<input_data>{\n']
        ret.append('    soil_nland : ')
        ret.append(repr(self.soil_nland))
        ret.append(',\n    soil_land : ')
        ret.append(repr(self.soil_land))
        ret.append(',\n    ptf_nland : ')
        ret.append(repr(self.ptf_nland))
        ret.append(',\n    ptf_land : ')
        ret.append(repr(self.ptf_land))
        ret.append(',\n    surf_nland : ')
        ret.append(repr(self.surf_nland))
        ret.append(',\n    surf_land : ')
        ret.append(repr(self.surf_land))
        ret.append(',\n    struct_nland : ')
        ret.append(repr(self.struct_nland))
        ret.append(',\n    struct_land : ')
        ret.append(repr(self.struct_land))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = [init_array_soil, init_array_ptf, init_array_surf, \
        init_array_structures]
    

@f90wrap.runtime.register_class("wrapping.Input_Param")
class Input_Param(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=input_param)
    
    
    Defined at m_sw_mono.f90 lines 393-485
    
    """
    def __init__(self, handle=None):
        """
        self = Input_Param()
        
        
        Defined at m_sw_mono.f90 lines 393-485
        
        
        Returns
        -------
        this : Input_Param
        	Object to be constructed
        
        
        Automatically generated constructor for input_param
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_input_param_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Input_Param
        
        
        Defined at m_sw_mono.f90 lines 393-485
        
        Parameters
        ----------
        this : Input_Param
        	Object to be destructed
        
        
        Automatically generated destructor for input_param
        """
        if self._alloc:
            _wrapping.f90wrap_input_param_finalise(this=self._handle)
    
    @property
    def mesh_type(self):
        """
        Element mesh_type ftype=character(len=lchar) pytype=str
        
        
        Defined at m_sw_mono.f90 line 394
        
        """
        return _wrapping.f90wrap_input_param__get__mesh_type(self._handle)
    
    @mesh_type.setter
    def mesh_type(self, mesh_type):
        _wrapping.f90wrap_input_param__set__mesh_type(self._handle, mesh_type)
    
    @property
    def mesh_name(self):
        """
        Element mesh_name ftype=character(len=lchar) pytype=str
        
        
        Defined at m_sw_mono.f90 line 395
        
        """
        return _wrapping.f90wrap_input_param__get__mesh_name(self._handle)
    
    @mesh_name.setter
    def mesh_name(self, mesh_name):
        _wrapping.f90wrap_input_param__set__mesh_name(self._handle, mesh_name)
    
    @property
    def bc_n(self):
        """
        Element bc_n ftype=character(len=lchar) pytype=str
        
        
        Defined at m_sw_mono.f90 line 396
        
        """
        return _wrapping.f90wrap_input_param__get__bc_n(self._handle)
    
    @bc_n.setter
    def bc_n(self, bc_n):
        _wrapping.f90wrap_input_param__set__bc_n(self._handle, bc_n)
    
    @property
    def bc_s(self):
        """
        Element bc_s ftype=character(len=lchar) pytype=str
        
        
        Defined at m_sw_mono.f90 line 397
        
        """
        return _wrapping.f90wrap_input_param__get__bc_s(self._handle)
    
    @bc_s.setter
    def bc_s(self, bc_s):
        _wrapping.f90wrap_input_param__set__bc_s(self._handle, bc_s)
    
    @property
    def bc_w(self):
        """
        Element bc_w ftype=character(len=lchar) pytype=str
        
        
        Defined at m_sw_mono.f90 line 398
        
        """
        return _wrapping.f90wrap_input_param__get__bc_w(self._handle)
    
    @bc_w.setter
    def bc_w(self, bc_w):
        _wrapping.f90wrap_input_param__set__bc_w(self._handle, bc_w)
    
    @property
    def bc_e(self):
        """
        Element bc_e ftype=character(len=lchar) pytype=str
        
        
        Defined at m_sw_mono.f90 line 399
        
        """
        return _wrapping.f90wrap_input_param__get__bc_e(self._handle)
    
    @bc_e.setter
    def bc_e(self, bc_e):
        _wrapping.f90wrap_input_param__set__bc_e(self._handle, bc_e)
    
    @property
    def bc_rain(self):
        """
        Element bc_rain ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 400
        
        """
        return _wrapping.f90wrap_input_param__get__bc_rain(self._handle)
    
    @bc_rain.setter
    def bc_rain(self, bc_rain):
        _wrapping.f90wrap_input_param__set__bc_rain(self._handle, bc_rain)
    
    @property
    def bc_infil(self):
        """
        Element bc_infil ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 401
        
        """
        return _wrapping.f90wrap_input_param__get__bc_infil(self._handle)
    
    @bc_infil.setter
    def bc_infil(self, bc_infil):
        _wrapping.f90wrap_input_param__set__bc_infil(self._handle, bc_infil)
    
    @property
    def lx(self):
        """
        Element lx ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 402
        
        """
        return _wrapping.f90wrap_input_param__get__lx(self._handle)
    
    @lx.setter
    def lx(self, lx):
        _wrapping.f90wrap_input_param__set__lx(self._handle, lx)
    
    @property
    def ly(self):
        """
        Element ly ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 403
        
        """
        return _wrapping.f90wrap_input_param__get__ly(self._handle)
    
    @ly.setter
    def ly(self, ly):
        _wrapping.f90wrap_input_param__set__ly(self._handle, ly)
    
    @property
    def nx(self):
        """
        Element nx ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 404
        
        """
        return _wrapping.f90wrap_input_param__get__nx(self._handle)
    
    @nx.setter
    def nx(self, nx):
        _wrapping.f90wrap_input_param__set__nx(self._handle, nx)
    
    @property
    def ny(self):
        """
        Element ny ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 405
        
        """
        return _wrapping.f90wrap_input_param__get__ny(self._handle)
    
    @ny.setter
    def ny(self, ny):
        _wrapping.f90wrap_input_param__set__ny(self._handle, ny)
    
    @property
    def ts(self):
        """
        Element ts ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 406
        
        """
        return _wrapping.f90wrap_input_param__get__ts(self._handle)
    
    @ts.setter
    def ts(self, ts):
        _wrapping.f90wrap_input_param__set__ts(self._handle, ts)
    
    @property
    def adapt_dt(self):
        """
        Element adapt_dt ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 407
        
        """
        return _wrapping.f90wrap_input_param__get__adapt_dt(self._handle)
    
    @adapt_dt.setter
    def adapt_dt(self, adapt_dt):
        _wrapping.f90wrap_input_param__set__adapt_dt(self._handle, adapt_dt)
    
    @property
    def dt(self):
        """
        Element dt ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 408
        
        """
        return _wrapping.f90wrap_input_param__get__dt(self._handle)
    
    @dt.setter
    def dt(self, dt):
        _wrapping.f90wrap_input_param__set__dt(self._handle, dt)
    
    @property
    def cfl(self):
        """
        Element cfl ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 409
        
        """
        return _wrapping.f90wrap_input_param__get__cfl(self._handle)
    
    @cfl.setter
    def cfl(self, cfl):
        _wrapping.f90wrap_input_param__set__cfl(self._handle, cfl)
    
    @property
    def do_warmup(self):
        """
        Element do_warmup ftype=logical pytype=bool
        
        
        Defined at m_sw_mono.f90 line 410
        
        """
        return _wrapping.f90wrap_input_param__get__do_warmup(self._handle)
    
    @do_warmup.setter
    def do_warmup(self, do_warmup):
        _wrapping.f90wrap_input_param__set__do_warmup(self._handle, do_warmup)
    
    @property
    def dtw(self):
        """
        Element dtw ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 411
        
        """
        return _wrapping.f90wrap_input_param__get__dtw(self._handle)
    
    @dtw.setter
    def dtw(self, dtw):
        _wrapping.f90wrap_input_param__set__dtw(self._handle, dtw)
    
    @property
    def dtp(self):
        """
        Element dtp ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 412
        
        """
        return _wrapping.f90wrap_input_param__get__dtp(self._handle)
    
    @dtp.setter
    def dtp(self, dtp):
        _wrapping.f90wrap_input_param__set__dtp(self._handle, dtp)
    
    @property
    def dta(self):
        """
        Element dta ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 413
        
        """
        return _wrapping.f90wrap_input_param__get__dta(self._handle)
    
    @dta.setter
    def dta(self, dta):
        _wrapping.f90wrap_input_param__set__dta(self._handle, dta)
    
    @property
    def w_tecplot(self):
        """
        Element w_tecplot ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 415
        
        """
        return _wrapping.f90wrap_input_param__get__w_tecplot(self._handle)
    
    @w_tecplot.setter
    def w_tecplot(self, w_tecplot):
        _wrapping.f90wrap_input_param__set__w_tecplot(self._handle, w_tecplot)
    
    @property
    def w_vtk(self):
        """
        Element w_vtk ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 416
        
        """
        return _wrapping.f90wrap_input_param__get__w_vtk(self._handle)
    
    @w_vtk.setter
    def w_vtk(self, w_vtk):
        _wrapping.f90wrap_input_param__set__w_vtk(self._handle, w_vtk)
    
    @property
    def w_gnuplot(self):
        """
        Element w_gnuplot ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 417
        
        """
        return _wrapping.f90wrap_input_param__get__w_gnuplot(self._handle)
    
    @w_gnuplot.setter
    def w_gnuplot(self, w_gnuplot):
        _wrapping.f90wrap_input_param__set__w_gnuplot(self._handle, w_gnuplot)
    
    @property
    def w_bin(self):
        """
        Element w_bin ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 418
        
        """
        return _wrapping.f90wrap_input_param__get__w_bin(self._handle)
    
    @w_bin.setter
    def w_bin(self, w_bin):
        _wrapping.f90wrap_input_param__set__w_bin(self._handle, w_bin)
    
    @property
    def w_exact(self):
        """
        Element w_exact ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 419
        
        """
        return _wrapping.f90wrap_input_param__get__w_exact(self._handle)
    
    @w_exact.setter
    def w_exact(self, w_exact):
        _wrapping.f90wrap_input_param__set__w_exact(self._handle, w_exact)
    
    @property
    def w_norm(self):
        """
        Element w_norm ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 420
        
        """
        return _wrapping.f90wrap_input_param__get__w_norm(self._handle)
    
    @w_norm.setter
    def w_norm(self, w_norm):
        _wrapping.f90wrap_input_param__set__w_norm(self._handle, w_norm)
    
    @property
    def w_obs(self):
        """
        Element w_obs ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 421
        
        """
        return _wrapping.f90wrap_input_param__get__w_obs(self._handle)
    
    @w_obs.setter
    def w_obs(self, w_obs):
        _wrapping.f90wrap_input_param__set__w_obs(self._handle, w_obs)
    
    @property
    def use_obs(self):
        """
        Element use_obs ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 422
        
        """
        return _wrapping.f90wrap_input_param__get__use_obs(self._handle)
    
    @use_obs.setter
    def use_obs(self, use_obs):
        _wrapping.f90wrap_input_param__set__use_obs(self._handle, use_obs)
    
    @property
    def use_zobs(self):
        """
        Element use_zobs ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 423
        
        """
        return _wrapping.f90wrap_input_param__get__use_zobs(self._handle)
    
    @use_zobs.setter
    def use_zobs(self, use_zobs):
        _wrapping.f90wrap_input_param__set__use_zobs(self._handle, use_zobs)
    
    @property
    def use_hobs(self):
        """
        Element use_hobs ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 424
        
        """
        return _wrapping.f90wrap_input_param__get__use_hobs(self._handle)
    
    @use_hobs.setter
    def use_hobs(self, use_hobs):
        _wrapping.f90wrap_input_param__set__use_hobs(self._handle, use_hobs)
    
    @property
    def use_uvobs(self):
        """
        Element use_uvobs ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 425
        
        """
        return _wrapping.f90wrap_input_param__get__use_uvobs(self._handle)
    
    @use_uvobs.setter
    def use_uvobs(self, use_uvobs):
        _wrapping.f90wrap_input_param__set__use_uvobs(self._handle, use_uvobs)
    
    @property
    def use_huvobs(self):
        """
        Element use_huvobs ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 426
        
        """
        return _wrapping.f90wrap_input_param__get__use_huvobs(self._handle)
    
    @use_huvobs.setter
    def use_huvobs(self, use_huvobs):
        _wrapping.f90wrap_input_param__set__use_huvobs(self._handle, use_huvobs)
    
    @property
    def use_qobs(self):
        """
        Element use_qobs ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 427
        
        """
        return _wrapping.f90wrap_input_param__get__use_qobs(self._handle)
    
    @use_qobs.setter
    def use_qobs(self, use_qobs):
        _wrapping.f90wrap_input_param__set__use_qobs(self._handle, use_qobs)
    
    @property
    def use_qobs_gr4(self):
        """
        Element use_qobs_gr4 ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 428
        
        """
        return _wrapping.f90wrap_input_param__get__use_qobs_gr4(self._handle)
    
    @use_qobs_gr4.setter
    def use_qobs_gr4(self, use_qobs_gr4):
        _wrapping.f90wrap_input_param__set__use_qobs_gr4(self._handle, use_qobs_gr4)
    
    @property
    def use_nse(self):
        """
        Element use_nse ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 429
        
        """
        return _wrapping.f90wrap_input_param__get__use_nse(self._handle)
    
    @use_nse.setter
    def use_nse(self, use_nse):
        _wrapping.f90wrap_input_param__set__use_nse(self._handle, use_nse)
    
    @property
    def use_xsshp(self):
        """
        Element use_xsshp ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 430
        
        """
        return _wrapping.f90wrap_input_param__get__use_xsshp(self._handle)
    
    @use_xsshp.setter
    def use_xsshp(self, use_xsshp):
        _wrapping.f90wrap_input_param__set__use_xsshp(self._handle, use_xsshp)
    
    @property
    def xsshp_along_x(self):
        """
        Element xsshp_along_x ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 431
        
        """
        return _wrapping.f90wrap_input_param__get__xsshp_along_x(self._handle)
    
    @xsshp_along_x.setter
    def xsshp_along_x(self, xsshp_along_x):
        _wrapping.f90wrap_input_param__set__xsshp_along_x(self._handle, xsshp_along_x)
    
    @property
    def xsshp_along_y(self):
        """
        Element xsshp_along_y ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 432
        
        """
        return _wrapping.f90wrap_input_param__get__xsshp_along_y(self._handle)
    
    @xsshp_along_y.setter
    def xsshp_along_y(self, xsshp_along_y):
        _wrapping.f90wrap_input_param__set__xsshp_along_y(self._handle, xsshp_along_y)
    
    @property
    def use_ptf(self):
        """
        Element use_ptf ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 433
        
        """
        return _wrapping.f90wrap_input_param__get__use_ptf(self._handle)
    
    @use_ptf.setter
    def use_ptf(self, use_ptf):
        _wrapping.f90wrap_input_param__set__use_ptf(self._handle, use_ptf)
    
    @property
    def spatial_scheme(self):
        """
        Element spatial_scheme ftype=character(len=lchar) pytype=str
        
        
        Defined at m_sw_mono.f90 line 434
        
        """
        return _wrapping.f90wrap_input_param__get__spatial_scheme(self._handle)
    
    @spatial_scheme.setter
    def spatial_scheme(self, spatial_scheme):
        _wrapping.f90wrap_input_param__set__spatial_scheme(self._handle, spatial_scheme)
    
    @property
    def temp_scheme(self):
        """
        Element temp_scheme ftype=character(len=lchar) pytype=str
        
        
        Defined at m_sw_mono.f90 line 435
        
        """
        return _wrapping.f90wrap_input_param__get__temp_scheme(self._handle)
    
    @temp_scheme.setter
    def temp_scheme(self, temp_scheme):
        _wrapping.f90wrap_input_param__set__temp_scheme(self._handle, temp_scheme)
    
    @property
    def args(self):
        """
        Element args ftype=character(len=lchar) pytype=str
        
        
        Defined at m_sw_mono.f90 line 436
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_input_param__array__args(self._handle)
        if array_handle in self._arrays:
            args = self._arrays[array_handle]
        else:
            args = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_input_param__array__args)
            self._arrays[array_handle] = args
        return args
    
    @args.setter
    def args(self, args):
        self.args[...] = args
    
    @property
    def max_nt_for_direct(self):
        """
        Element max_nt_for_direct ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 437
        
        """
        return _wrapping.f90wrap_input_param__get__max_nt_for_direct(self._handle)
    
    @max_nt_for_direct.setter
    def max_nt_for_direct(self, max_nt_for_direct):
        _wrapping.f90wrap_input_param__set__max_nt_for_direct(self._handle, \
            max_nt_for_direct)
    
    @property
    def max_nt_for_adjoint(self):
        """
        Element max_nt_for_adjoint ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 438
        
        """
        return _wrapping.f90wrap_input_param__get__max_nt_for_adjoint(self._handle)
    
    @max_nt_for_adjoint.setter
    def max_nt_for_adjoint(self, max_nt_for_adjoint):
        _wrapping.f90wrap_input_param__set__max_nt_for_adjoint(self._handle, \
            max_nt_for_adjoint)
    
    @property
    def g(self):
        """
        Element g ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 440
        
        """
        return _wrapping.f90wrap_input_param__get__g(self._handle)
    
    @g.setter
    def g(self, g):
        _wrapping.f90wrap_input_param__set__g(self._handle, g)
    
    @property
    def heps(self):
        """
        Element heps ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 441
        
        """
        return _wrapping.f90wrap_input_param__get__heps(self._handle)
    
    @heps.setter
    def heps(self, heps):
        _wrapping.f90wrap_input_param__set__heps(self._handle, heps)
    
    @property
    def friction(self):
        """
        Element friction ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 442
        
        """
        return _wrapping.f90wrap_input_param__get__friction(self._handle)
    
    @friction.setter
    def friction(self, friction):
        _wrapping.f90wrap_input_param__set__friction(self._handle, friction)
    
    @property
    def c_shape_s(self):
        """
        Element c_shape_s ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 445
        
        """
        return _wrapping.f90wrap_input_param__get__c_shape_s(self._handle)
    
    @c_shape_s.setter
    def c_shape_s(self, c_shape_s):
        _wrapping.f90wrap_input_param__set__c_shape_s(self._handle, c_shape_s)
    
    @property
    def c_xcenter(self):
        """
        Element c_xcenter ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 446
        
        """
        return _wrapping.f90wrap_input_param__get__c_xcenter(self._handle)
    
    @c_xcenter.setter
    def c_xcenter(self, c_xcenter):
        _wrapping.f90wrap_input_param__set__c_xcenter(self._handle, c_xcenter)
    
    @property
    def c_hmax(self):
        """
        Element c_hmax ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 447
        
        """
        return _wrapping.f90wrap_input_param__get__c_hmax(self._handle)
    
    @c_hmax.setter
    def c_hmax(self, c_hmax):
        _wrapping.f90wrap_input_param__set__c_hmax(self._handle, c_hmax)
    
    @property
    def c_manning(self):
        """
        Element c_manning ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 448
        
        """
        return _wrapping.f90wrap_input_param__get__c_manning(self._handle)
    
    @c_manning.setter
    def c_manning(self, c_manning):
        _wrapping.f90wrap_input_param__set__c_manning(self._handle, c_manning)
    
    @property
    def c_manning_beta(self):
        """
        Element c_manning_beta ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 449
        
        """
        return _wrapping.f90wrap_input_param__get__c_manning_beta(self._handle)
    
    @c_manning_beta.setter
    def c_manning_beta(self, c_manning_beta):
        _wrapping.f90wrap_input_param__set__c_manning_beta(self._handle, c_manning_beta)
    
    @property
    def c_bathy(self):
        """
        Element c_bathy ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 450
        
        """
        return _wrapping.f90wrap_input_param__get__c_bathy(self._handle)
    
    @c_bathy.setter
    def c_bathy(self, c_bathy):
        _wrapping.f90wrap_input_param__set__c_bathy(self._handle, c_bathy)
    
    @property
    def c_ic(self):
        """
        Element c_ic ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 451
        
        """
        return _wrapping.f90wrap_input_param__get__c_ic(self._handle)
    
    @c_ic.setter
    def c_ic(self, c_ic):
        _wrapping.f90wrap_input_param__set__c_ic(self._handle, c_ic)
    
    @property
    def c_hydrograph(self):
        """
        Element c_hydrograph ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 452
        
        """
        return _wrapping.f90wrap_input_param__get__c_hydrograph(self._handle)
    
    @c_hydrograph.setter
    def c_hydrograph(self, c_hydrograph):
        _wrapping.f90wrap_input_param__set__c_hydrograph(self._handle, c_hydrograph)
    
    @property
    def c_ratcurve(self):
        """
        Element c_ratcurve ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 453
        
        """
        return _wrapping.f90wrap_input_param__get__c_ratcurve(self._handle)
    
    @c_ratcurve.setter
    def c_ratcurve(self, c_ratcurve):
        _wrapping.f90wrap_input_param__set__c_ratcurve(self._handle, c_ratcurve)
    
    @property
    def c_rain(self):
        """
        Element c_rain ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 454
        
        """
        return _wrapping.f90wrap_input_param__get__c_rain(self._handle)
    
    @c_rain.setter
    def c_rain(self, c_rain):
        _wrapping.f90wrap_input_param__set__c_rain(self._handle, c_rain)
    
    @property
    def c_gr4params(self):
        """
        Element c_gr4params ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 455
        
        """
        return _wrapping.f90wrap_input_param__get__c_gr4params(self._handle)
    
    @c_gr4params.setter
    def c_gr4params(self, c_gr4params):
        _wrapping.f90wrap_input_param__set__c_gr4params(self._handle, c_gr4params)
    
    @property
    def c_infil_max(self):
        """
        Element c_infil_max ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 456
        
        """
        return _wrapping.f90wrap_input_param__get__c_infil_max(self._handle)
    
    @c_infil_max.setter
    def c_infil_max(self, c_infil_max):
        _wrapping.f90wrap_input_param__set__c_infil_max(self._handle, c_infil_max)
    
    @property
    def c_ks(self):
        """
        Element c_ks ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 457
        
        """
        return _wrapping.f90wrap_input_param__get__c_ks(self._handle)
    
    @c_ks.setter
    def c_ks(self, c_ks):
        _wrapping.f90wrap_input_param__set__c_ks(self._handle, c_ks)
    
    @property
    def c_psif(self):
        """
        Element c_psif ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 458
        
        """
        return _wrapping.f90wrap_input_param__get__c_psif(self._handle)
    
    @c_psif.setter
    def c_psif(self, c_psif):
        _wrapping.f90wrap_input_param__set__c_psif(self._handle, c_psif)
    
    @property
    def c_deltatheta(self):
        """
        Element c_deltatheta ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 459
        
        """
        return _wrapping.f90wrap_input_param__get__c_deltatheta(self._handle)
    
    @c_deltatheta.setter
    def c_deltatheta(self, c_deltatheta):
        _wrapping.f90wrap_input_param__set__c_deltatheta(self._handle, c_deltatheta)
    
    @property
    def c_lambda(self):
        """
        Element c_lambda ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 460
        
        """
        return _wrapping.f90wrap_input_param__get__c_lambda(self._handle)
    
    @c_lambda.setter
    def c_lambda(self, c_lambda):
        _wrapping.f90wrap_input_param__set__c_lambda(self._handle, c_lambda)
    
    @property
    def c_cn(self):
        """
        Element c_cn ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 461
        
        """
        return _wrapping.f90wrap_input_param__get__c_cn(self._handle)
    
    @c_cn.setter
    def c_cn(self, c_cn):
        _wrapping.f90wrap_input_param__set__c_cn(self._handle, c_cn)
    
    @property
    def c_ptf(self):
        """
        Element c_ptf ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 462
        
        """
        return _wrapping.f90wrap_input_param__get__c_ptf(self._handle)
    
    @c_ptf.setter
    def c_ptf(self, c_ptf):
        _wrapping.f90wrap_input_param__set__c_ptf(self._handle, c_ptf)
    
    @property
    def eps_min(self):
        """
        Element eps_min ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 465
        
        """
        return _wrapping.f90wrap_input_param__get__eps_min(self._handle)
    
    @eps_min.setter
    def eps_min(self, eps_min):
        _wrapping.f90wrap_input_param__set__eps_min(self._handle, eps_min)
    
    @property
    def eps_manning(self):
        """
        Element eps_manning ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 466
        
        """
        return _wrapping.f90wrap_input_param__get__eps_manning(self._handle)
    
    @eps_manning.setter
    def eps_manning(self, eps_manning):
        _wrapping.f90wrap_input_param__set__eps_manning(self._handle, eps_manning)
    
    @property
    def eps_bathy(self):
        """
        Element eps_bathy ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 467
        
        """
        return _wrapping.f90wrap_input_param__get__eps_bathy(self._handle)
    
    @eps_bathy.setter
    def eps_bathy(self, eps_bathy):
        _wrapping.f90wrap_input_param__set__eps_bathy(self._handle, eps_bathy)
    
    @property
    def eps_ic(self):
        """
        Element eps_ic ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 468
        
        """
        return _wrapping.f90wrap_input_param__get__eps_ic(self._handle)
    
    @eps_ic.setter
    def eps_ic(self, eps_ic):
        _wrapping.f90wrap_input_param__set__eps_ic(self._handle, eps_ic)
    
    @property
    def eps_hydrograph(self):
        """
        Element eps_hydrograph ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 469
        
        """
        return _wrapping.f90wrap_input_param__get__eps_hydrograph(self._handle)
    
    @eps_hydrograph.setter
    def eps_hydrograph(self, eps_hydrograph):
        _wrapping.f90wrap_input_param__set__eps_hydrograph(self._handle, eps_hydrograph)
    
    @property
    def eps_ratcurve(self):
        """
        Element eps_ratcurve ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 470
        
        """
        return _wrapping.f90wrap_input_param__get__eps_ratcurve(self._handle)
    
    @eps_ratcurve.setter
    def eps_ratcurve(self, eps_ratcurve):
        _wrapping.f90wrap_input_param__set__eps_ratcurve(self._handle, eps_ratcurve)
    
    @property
    def eps_rain(self):
        """
        Element eps_rain ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 471
        
        """
        return _wrapping.f90wrap_input_param__get__eps_rain(self._handle)
    
    @eps_rain.setter
    def eps_rain(self, eps_rain):
        _wrapping.f90wrap_input_param__set__eps_rain(self._handle, eps_rain)
    
    @property
    def eps_ks(self):
        """
        Element eps_ks ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 472
        
        """
        return _wrapping.f90wrap_input_param__get__eps_ks(self._handle)
    
    @eps_ks.setter
    def eps_ks(self, eps_ks):
        _wrapping.f90wrap_input_param__set__eps_ks(self._handle, eps_ks)
    
    @property
    def eps_psif(self):
        """
        Element eps_psif ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 473
        
        """
        return _wrapping.f90wrap_input_param__get__eps_psif(self._handle)
    
    @eps_psif.setter
    def eps_psif(self, eps_psif):
        _wrapping.f90wrap_input_param__set__eps_psif(self._handle, eps_psif)
    
    @property
    def eps_deltatheta(self):
        """
        Element eps_deltatheta ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 474
        
        """
        return _wrapping.f90wrap_input_param__get__eps_deltatheta(self._handle)
    
    @eps_deltatheta.setter
    def eps_deltatheta(self, eps_deltatheta):
        _wrapping.f90wrap_input_param__set__eps_deltatheta(self._handle, eps_deltatheta)
    
    @property
    def eps_lambdacn(self):
        """
        Element eps_lambdacn ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 475
        
        """
        return _wrapping.f90wrap_input_param__get__eps_lambdacn(self._handle)
    
    @eps_lambdacn.setter
    def eps_lambdacn(self, eps_lambdacn):
        _wrapping.f90wrap_input_param__set__eps_lambdacn(self._handle, eps_lambdacn)
    
    @property
    def eps_cn(self):
        """
        Element eps_cn ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 476
        
        """
        return _wrapping.f90wrap_input_param__get__eps_cn(self._handle)
    
    @eps_cn.setter
    def eps_cn(self, eps_cn):
        _wrapping.f90wrap_input_param__set__eps_cn(self._handle, eps_cn)
    
    @property
    def eps_ptf(self):
        """
        Element eps_ptf ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 477
        
        """
        return _wrapping.f90wrap_input_param__get__eps_ptf(self._handle)
    
    @eps_ptf.setter
    def eps_ptf(self, eps_ptf):
        _wrapping.f90wrap_input_param__set__eps_ptf(self._handle, eps_ptf)
    
    @property
    def regul_manning(self):
        """
        Element regul_manning ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 479
        
        """
        return _wrapping.f90wrap_input_param__get__regul_manning(self._handle)
    
    @regul_manning.setter
    def regul_manning(self, regul_manning):
        _wrapping.f90wrap_input_param__set__regul_manning(self._handle, regul_manning)
    
    @property
    def regul_bathy(self):
        """
        Element regul_bathy ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 480
        
        """
        return _wrapping.f90wrap_input_param__get__regul_bathy(self._handle)
    
    @regul_bathy.setter
    def regul_bathy(self, regul_bathy):
        _wrapping.f90wrap_input_param__set__regul_bathy(self._handle, regul_bathy)
    
    @property
    def regul_bathy_grad(self):
        """
        Element regul_bathy_grad ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 481
        
        """
        return _wrapping.f90wrap_input_param__get__regul_bathy_grad(self._handle)
    
    @regul_bathy_grad.setter
    def regul_bathy_grad(self, regul_bathy_grad):
        _wrapping.f90wrap_input_param__set__regul_bathy_grad(self._handle, \
            regul_bathy_grad)
    
    @property
    def regul_bathy_shape(self):
        """
        Element regul_bathy_shape ftype=integer(ip) pytype=int
        
        
        Defined at m_sw_mono.f90 line 482
        
        """
        return _wrapping.f90wrap_input_param__get__regul_bathy_shape(self._handle)
    
    @regul_bathy_shape.setter
    def regul_bathy_shape(self, regul_bathy_shape):
        _wrapping.f90wrap_input_param__set__regul_bathy_shape(self._handle, \
            regul_bathy_shape)
    
    @property
    def regul_ic(self):
        """
        Element regul_ic ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 483
        
        """
        return _wrapping.f90wrap_input_param__get__regul_ic(self._handle)
    
    @regul_ic.setter
    def regul_ic(self, regul_ic):
        _wrapping.f90wrap_input_param__set__regul_ic(self._handle, regul_ic)
    
    @property
    def regul_hydrograph(self):
        """
        Element regul_hydrograph ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 484
        
        """
        return _wrapping.f90wrap_input_param__get__regul_hydrograph(self._handle)
    
    @regul_hydrograph.setter
    def regul_hydrograph(self, regul_hydrograph):
        _wrapping.f90wrap_input_param__set__regul_hydrograph(self._handle, \
            regul_hydrograph)
    
    @property
    def regul_ratcurve(self):
        """
        Element regul_ratcurve ftype=real(rp) pytype=float
        
        
        Defined at m_sw_mono.f90 line 485
        
        """
        return _wrapping.f90wrap_input_param__get__regul_ratcurve(self._handle)
    
    @regul_ratcurve.setter
    def regul_ratcurve(self, regul_ratcurve):
        _wrapping.f90wrap_input_param__set__regul_ratcurve(self._handle, regul_ratcurve)
    
    def __str__(self):
        ret = ['<input_param>{\n']
        ret.append('    mesh_type : ')
        ret.append(repr(self.mesh_type))
        ret.append(',\n    mesh_name : ')
        ret.append(repr(self.mesh_name))
        ret.append(',\n    bc_n : ')
        ret.append(repr(self.bc_n))
        ret.append(',\n    bc_s : ')
        ret.append(repr(self.bc_s))
        ret.append(',\n    bc_w : ')
        ret.append(repr(self.bc_w))
        ret.append(',\n    bc_e : ')
        ret.append(repr(self.bc_e))
        ret.append(',\n    bc_rain : ')
        ret.append(repr(self.bc_rain))
        ret.append(',\n    bc_infil : ')
        ret.append(repr(self.bc_infil))
        ret.append(',\n    lx : ')
        ret.append(repr(self.lx))
        ret.append(',\n    ly : ')
        ret.append(repr(self.ly))
        ret.append(',\n    nx : ')
        ret.append(repr(self.nx))
        ret.append(',\n    ny : ')
        ret.append(repr(self.ny))
        ret.append(',\n    ts : ')
        ret.append(repr(self.ts))
        ret.append(',\n    adapt_dt : ')
        ret.append(repr(self.adapt_dt))
        ret.append(',\n    dt : ')
        ret.append(repr(self.dt))
        ret.append(',\n    cfl : ')
        ret.append(repr(self.cfl))
        ret.append(',\n    do_warmup : ')
        ret.append(repr(self.do_warmup))
        ret.append(',\n    dtw : ')
        ret.append(repr(self.dtw))
        ret.append(',\n    dtp : ')
        ret.append(repr(self.dtp))
        ret.append(',\n    dta : ')
        ret.append(repr(self.dta))
        ret.append(',\n    w_tecplot : ')
        ret.append(repr(self.w_tecplot))
        ret.append(',\n    w_vtk : ')
        ret.append(repr(self.w_vtk))
        ret.append(',\n    w_gnuplot : ')
        ret.append(repr(self.w_gnuplot))
        ret.append(',\n    w_bin : ')
        ret.append(repr(self.w_bin))
        ret.append(',\n    w_exact : ')
        ret.append(repr(self.w_exact))
        ret.append(',\n    w_norm : ')
        ret.append(repr(self.w_norm))
        ret.append(',\n    w_obs : ')
        ret.append(repr(self.w_obs))
        ret.append(',\n    use_obs : ')
        ret.append(repr(self.use_obs))
        ret.append(',\n    use_zobs : ')
        ret.append(repr(self.use_zobs))
        ret.append(',\n    use_hobs : ')
        ret.append(repr(self.use_hobs))
        ret.append(',\n    use_uvobs : ')
        ret.append(repr(self.use_uvobs))
        ret.append(',\n    use_huvobs : ')
        ret.append(repr(self.use_huvobs))
        ret.append(',\n    use_qobs : ')
        ret.append(repr(self.use_qobs))
        ret.append(',\n    use_qobs_gr4 : ')
        ret.append(repr(self.use_qobs_gr4))
        ret.append(',\n    use_nse : ')
        ret.append(repr(self.use_nse))
        ret.append(',\n    use_xsshp : ')
        ret.append(repr(self.use_xsshp))
        ret.append(',\n    xsshp_along_x : ')
        ret.append(repr(self.xsshp_along_x))
        ret.append(',\n    xsshp_along_y : ')
        ret.append(repr(self.xsshp_along_y))
        ret.append(',\n    use_ptf : ')
        ret.append(repr(self.use_ptf))
        ret.append(',\n    spatial_scheme : ')
        ret.append(repr(self.spatial_scheme))
        ret.append(',\n    temp_scheme : ')
        ret.append(repr(self.temp_scheme))
        ret.append(',\n    args : ')
        ret.append(repr(self.args))
        ret.append(',\n    max_nt_for_direct : ')
        ret.append(repr(self.max_nt_for_direct))
        ret.append(',\n    max_nt_for_adjoint : ')
        ret.append(repr(self.max_nt_for_adjoint))
        ret.append(',\n    g : ')
        ret.append(repr(self.g))
        ret.append(',\n    heps : ')
        ret.append(repr(self.heps))
        ret.append(',\n    friction : ')
        ret.append(repr(self.friction))
        ret.append(',\n    c_shape_s : ')
        ret.append(repr(self.c_shape_s))
        ret.append(',\n    c_xcenter : ')
        ret.append(repr(self.c_xcenter))
        ret.append(',\n    c_hmax : ')
        ret.append(repr(self.c_hmax))
        ret.append(',\n    c_manning : ')
        ret.append(repr(self.c_manning))
        ret.append(',\n    c_manning_beta : ')
        ret.append(repr(self.c_manning_beta))
        ret.append(',\n    c_bathy : ')
        ret.append(repr(self.c_bathy))
        ret.append(',\n    c_ic : ')
        ret.append(repr(self.c_ic))
        ret.append(',\n    c_hydrograph : ')
        ret.append(repr(self.c_hydrograph))
        ret.append(',\n    c_ratcurve : ')
        ret.append(repr(self.c_ratcurve))
        ret.append(',\n    c_rain : ')
        ret.append(repr(self.c_rain))
        ret.append(',\n    c_gr4params : ')
        ret.append(repr(self.c_gr4params))
        ret.append(',\n    c_infil_max : ')
        ret.append(repr(self.c_infil_max))
        ret.append(',\n    c_ks : ')
        ret.append(repr(self.c_ks))
        ret.append(',\n    c_psif : ')
        ret.append(repr(self.c_psif))
        ret.append(',\n    c_deltatheta : ')
        ret.append(repr(self.c_deltatheta))
        ret.append(',\n    c_lambda : ')
        ret.append(repr(self.c_lambda))
        ret.append(',\n    c_cn : ')
        ret.append(repr(self.c_cn))
        ret.append(',\n    c_ptf : ')
        ret.append(repr(self.c_ptf))
        ret.append(',\n    eps_min : ')
        ret.append(repr(self.eps_min))
        ret.append(',\n    eps_manning : ')
        ret.append(repr(self.eps_manning))
        ret.append(',\n    eps_bathy : ')
        ret.append(repr(self.eps_bathy))
        ret.append(',\n    eps_ic : ')
        ret.append(repr(self.eps_ic))
        ret.append(',\n    eps_hydrograph : ')
        ret.append(repr(self.eps_hydrograph))
        ret.append(',\n    eps_ratcurve : ')
        ret.append(repr(self.eps_ratcurve))
        ret.append(',\n    eps_rain : ')
        ret.append(repr(self.eps_rain))
        ret.append(',\n    eps_ks : ')
        ret.append(repr(self.eps_ks))
        ret.append(',\n    eps_psif : ')
        ret.append(repr(self.eps_psif))
        ret.append(',\n    eps_deltatheta : ')
        ret.append(repr(self.eps_deltatheta))
        ret.append(',\n    eps_lambdacn : ')
        ret.append(repr(self.eps_lambdacn))
        ret.append(',\n    eps_cn : ')
        ret.append(repr(self.eps_cn))
        ret.append(',\n    eps_ptf : ')
        ret.append(repr(self.eps_ptf))
        ret.append(',\n    regul_manning : ')
        ret.append(repr(self.regul_manning))
        ret.append(',\n    regul_bathy : ')
        ret.append(repr(self.regul_bathy))
        ret.append(',\n    regul_bathy_grad : ')
        ret.append(repr(self.regul_bathy_grad))
        ret.append(',\n    regul_bathy_shape : ')
        ret.append(repr(self.regul_bathy_shape))
        ret.append(',\n    regul_ic : ')
        ret.append(repr(self.regul_ic))
        ret.append(',\n    regul_hydrograph : ')
        ret.append(repr(self.regul_hydrograph))
        ret.append(',\n    regul_ratcurve : ')
        ret.append(repr(self.regul_ratcurve))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def default_values():
    """
    default_values()
    
    
    Defined at m_sw_mono.f90 lines 489-572
    
    
    Details: h,u,v are the unknowns of the model
    """
    _wrapping.f90wrap_default_values()

def alloc_dof(self):
    """
    dof = alloc_dof(self)
    
    
    Defined at m_sw_mono.f90 lines 574-626
    
    Parameters
    ----------
    mesh : Msh
    
    Returns
    -------
    dof : Unk
    
     Allocation of Model unk
     Notes
     -----
     **alloc_dof** :
    
     - Allocate dof.
     ``common_array_allocation`` [mod_common_data.f90]
     - Run smash core.
     ``smash_core`` [mod_smash_interface.f90]
     - Write results in ``.txt`` format.
     ``write_smash_results`` [WRITE_RESULTS.f90]
    
    _COMMENT ============================= ===================================
    _COMMENT Parameters Description
    _COMMENT ============================= ===================================
    _COMMENT ``setup`` model_setup Derived Type
    _COMMENT ``domain`` mesh Derived Type
    _COMMENT ``watershed`` catchments Derived Type
    _COMMENT ``inputdata`` input_data Derived Type
    _COMMENT ``lois`` loi_ouvrage Derived Type
    _COMMENT ``param`` spatialparam Derived Type
    _COMMENT ``model_states`` spatialstates Derived Type
    _COMMENT ``model_routing_states`` spatiotemporalstates Derived Type
    _COMMENT ``outputs`` smash_outputs Derived Type
    _COMMENT ============================= ===================================
    """
    dof = _wrapping.f90wrap_alloc_dof(mesh=self._handle)
    dof = f90wrap.runtime.lookup_class("wrapping.unk").from_handle(dof, alloc=True)
    return dof

def dealloc_dof(self):
    """
    dealloc_dof(self)
    
    
    Defined at m_sw_mono.f90 lines 628-639
    
    Parameters
    ----------
    dof : Unk
    
    """
    _wrapping.f90wrap_dealloc_dof(dof=self._handle)

def dealloc_model():
    """
    dealloc_model()
    
    
    Defined at m_sw_mono.f90 lines 641-693
    
    
    ------------------------------------------------
     millascenious forgoten variables to deallocate
    ------------------------------------------------
     BC
    """
    _wrapping.f90wrap_dealloc_model()

def spatial_index_fromxy(self, xmin, xmax, ymin, ymax, spatial_index):
    """
    spatial_index_fromxy(self, xmin, xmax, ymin, ymax, spatial_index)
    
    
    Defined at m_sw_mono.f90 lines 820-828
    
    Parameters
    ----------
    mesh : Msh, [in] Mesh
    xmin : float
    xmax : float
    ymin : float
    ymax : float
    spatial_index : int, [in] Cell index
    
    Details: From a file name and a extension('tecplot', 'gnuplot', ...) complete \
        file name.
    _COMMENT \param[in] xmin, xmax, ymin, ymax Coordinates of the tile
    """
    _wrapping.f90wrap_spatial_index_fromxy(mesh=self._handle, xmin=xmin, xmax=xmax, \
        ymin=ymin, ymax=ymax, spatial_index=spatial_index)

def get_sw_nb():
    """
    Element sw_nb ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 8
    
    """
    return _wrapping.f90wrap_m_model__get__sw_nb()

def set_sw_nb(sw_nb):
    _wrapping.f90wrap_m_model__set__sw_nb(sw_nb)

def get_array_bathy_node():
    """
    Element bathy_node ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 31
    
    """
    global bathy_node
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_model__array__bathy_node(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        bathy_node = _arrays[array_handle]
    else:
        bathy_node = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_model__array__bathy_node)
        _arrays[array_handle] = bathy_node
    return bathy_node

def set_array_bathy_node(bathy_node):
    globals()['bathy_node'][...] = bathy_node

def get_array_bathy_cell():
    """
    Element bathy_cell ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 32
    
    """
    global bathy_cell
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_model__array__bathy_cell(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        bathy_cell = _arrays[array_handle]
    else:
        bathy_cell = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_model__array__bathy_cell)
        _arrays[array_handle] = bathy_cell
    return bathy_cell

def set_array_bathy_cell(bathy_cell):
    globals()['bathy_cell'][...] = bathy_cell

def get_array_manning():
    """
    Element manning ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 33
    
    """
    global manning
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_model__array__manning(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        manning = _arrays[array_handle]
    else:
        manning = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_model__array__manning)
        _arrays[array_handle] = manning
    return manning

def set_array_manning(manning):
    globals()['manning'][...] = manning

def get_array_manning_beta():
    """
    Element manning_beta ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 34
    
    """
    global manning_beta
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_model__array__manning_beta(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        manning_beta = _arrays[array_handle]
    else:
        manning_beta = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_model__array__manning_beta)
        _arrays[array_handle] = manning_beta
    return manning_beta

def set_array_manning_beta(manning_beta):
    globals()['manning_beta'][...] = manning_beta

def get_nland():
    """
    Element nland ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 35
    
    """
    return _wrapping.f90wrap_m_model__get__nland()

def set_nland(nland):
    _wrapping.f90wrap_m_model__set__nland(nland)

def get_array_land():
    """
    Element land ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 36
    
    """
    global land
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_model__array__land(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        land = _arrays[array_handle]
    else:
        land = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_model__array__land)
        _arrays[array_handle] = land
    return land

def set_array_land(land):
    globals()['land'][...] = land

def get_mass_cut():
    """
    Element mass_cut ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 40
    
    """
    return _wrapping.f90wrap_m_model__get__mass_cut()

def set_mass_cut(mass_cut):
    _wrapping.f90wrap_m_model__set__mass_cut(mass_cut)

def get_manning_data_glob():
    """
    Element manning_data_glob ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 41
    
    """
    return _wrapping.f90wrap_m_model__get__manning_data_glob()

def set_manning_data_glob(manning_data_glob):
    _wrapping.f90wrap_m_model__set__manning_data_glob(manning_data_glob)

def get_array_slope_y():
    """
    Element slope_y ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 42
    
    """
    global slope_y
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_model__array__slope_y(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        slope_y = _arrays[array_handle]
    else:
        slope_y = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_model__array__slope_y)
        _arrays[array_handle] = slope_y
    return slope_y

def set_array_slope_y(slope_y):
    globals()['slope_y'][...] = slope_y

def get_array_slope_x():
    """
    Element slope_x ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 43
    
    """
    global slope_x
    array_ndim, array_type, array_shape, array_handle = \
        _wrapping.f90wrap_m_model__array__slope_x(f90wrap.runtime.empty_handle)
    if array_handle in _arrays:
        slope_x = _arrays[array_handle]
    else:
        slope_x = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                f90wrap.runtime.empty_handle,
                                _wrapping.f90wrap_m_model__array__slope_x)
        _arrays[array_handle] = slope_x
    return slope_x

def set_array_slope_x(slope_x):
    globals()['slope_x'][...] = slope_x

def get_single_porosity():
    """
    Element single_porosity ftype=type(porosity_data) pytype=Porosity_Data
    
    
    Defined at m_sw_mono.f90 line 97
    
    """
    global single_porosity
    single_porosity_handle = _wrapping.f90wrap_m_model__get__single_porosity()
    if tuple(single_porosity_handle) in _objs:
        single_porosity = _objs[tuple(single_porosity_handle)]
    else:
        single_porosity = porosity_data.from_handle(single_porosity_handle)
        _objs[tuple(single_porosity_handle)] = single_porosity
    return single_porosity

def set_single_porosity(single_porosity):
    single_porosity = single_porosity._handle
    _wrapping.f90wrap_m_model__set__single_porosity(single_porosity)

def get_feedback_inflow():
    """
    Element feedback_inflow ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 108
    
    """
    return _wrapping.f90wrap_m_model__get__feedback_inflow()

def set_feedback_inflow(feedback_inflow):
    _wrapping.f90wrap_m_model__set__feedback_inflow(feedback_inflow)

def get_coef_feedback():
    """
    Element coef_feedback ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 109
    
    """
    return _wrapping.f90wrap_m_model__get__coef_feedback()

def set_coef_feedback(coef_feedback):
    _wrapping.f90wrap_m_model__set__coef_feedback(coef_feedback)

def get_bc():
    """
    Element bc ftype=type(bcs) pytype=Bcs
    
    
    Defined at m_sw_mono.f90 line 163
    
    """
    global bc
    bc_handle = _wrapping.f90wrap_m_model__get__bc()
    if tuple(bc_handle) in _objs:
        bc = _objs[tuple(bc_handle)]
    else:
        bc = bcs.from_handle(bc_handle)
        _objs[tuple(bc_handle)] = bc
    return bc

def set_bc(bc):
    bc = bc._handle
    _wrapping.f90wrap_m_model__set__bc(bc)

def get_phys_desc():
    """
    Element phys_desc ftype=type(input_data) pytype=Input_Data
    
    
    Defined at m_sw_mono.f90 line 243
    
    """
    global phys_desc
    phys_desc_handle = _wrapping.f90wrap_m_model__get__phys_desc()
    if tuple(phys_desc_handle) in _objs:
        phys_desc = _objs[tuple(phys_desc_handle)]
    else:
        phys_desc = input_data.from_handle(phys_desc_handle)
        _objs[tuple(phys_desc_handle)] = phys_desc
    return phys_desc

def set_phys_desc(phys_desc):
    phys_desc = phys_desc._handle
    _wrapping.f90wrap_m_model__set__phys_desc(phys_desc)

def get_g():
    """
    Element g ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 248
    
    """
    return _wrapping.f90wrap_m_model__get__g()

def set_g(g):
    _wrapping.f90wrap_m_model__set__g(g)

def get_heps():
    """
    Element heps ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 249
    
    """
    return _wrapping.f90wrap_m_model__get__heps()

def set_heps(heps):
    _wrapping.f90wrap_m_model__set__heps(heps)

def get_friction():
    """
    Element friction ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 250
    
    """
    return _wrapping.f90wrap_m_model__get__friction()

def set_friction(friction):
    _wrapping.f90wrap_m_model__set__friction(friction)

def get_c_shape_s():
    """
    Element c_shape_s ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 253
    
    """
    return _wrapping.f90wrap_m_model__get__c_shape_s()

def set_c_shape_s(c_shape_s):
    _wrapping.f90wrap_m_model__set__c_shape_s(c_shape_s)

def get_c_xcenter():
    """
    Element c_xcenter ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 254
    
    """
    return _wrapping.f90wrap_m_model__get__c_xcenter()

def set_c_xcenter(c_xcenter):
    _wrapping.f90wrap_m_model__set__c_xcenter(c_xcenter)

def get_c_hmax():
    """
    Element c_hmax ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 255
    
    """
    return _wrapping.f90wrap_m_model__get__c_hmax()

def set_c_hmax(c_hmax):
    _wrapping.f90wrap_m_model__set__c_hmax(c_hmax)

def get_c_manning():
    """
    Element c_manning ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 256
    
    """
    return _wrapping.f90wrap_m_model__get__c_manning()

def set_c_manning(c_manning):
    _wrapping.f90wrap_m_model__set__c_manning(c_manning)

def get_c_manning_beta():
    """
    Element c_manning_beta ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 257
    
    """
    return _wrapping.f90wrap_m_model__get__c_manning_beta()

def set_c_manning_beta(c_manning_beta):
    _wrapping.f90wrap_m_model__set__c_manning_beta(c_manning_beta)

def get_c_bathy():
    """
    Element c_bathy ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 258
    
    """
    return _wrapping.f90wrap_m_model__get__c_bathy()

def set_c_bathy(c_bathy):
    _wrapping.f90wrap_m_model__set__c_bathy(c_bathy)

def get_c_slope_y():
    """
    Element c_slope_y ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 259
    
    """
    return _wrapping.f90wrap_m_model__get__c_slope_y()

def set_c_slope_y(c_slope_y):
    _wrapping.f90wrap_m_model__set__c_slope_y(c_slope_y)

def get_c_slope_x():
    """
    Element c_slope_x ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 260
    
    """
    return _wrapping.f90wrap_m_model__get__c_slope_x()

def set_c_slope_x(c_slope_x):
    _wrapping.f90wrap_m_model__set__c_slope_x(c_slope_x)

def get_c_ic():
    """
    Element c_ic ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 261
    
    """
    return _wrapping.f90wrap_m_model__get__c_ic()

def set_c_ic(c_ic):
    _wrapping.f90wrap_m_model__set__c_ic(c_ic)

def get_c_hydrograph():
    """
    Element c_hydrograph ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 262
    
    """
    return _wrapping.f90wrap_m_model__get__c_hydrograph()

def set_c_hydrograph(c_hydrograph):
    _wrapping.f90wrap_m_model__set__c_hydrograph(c_hydrograph)

def get_c_ratcurve():
    """
    Element c_ratcurve ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 263
    
    """
    return _wrapping.f90wrap_m_model__get__c_ratcurve()

def set_c_ratcurve(c_ratcurve):
    _wrapping.f90wrap_m_model__set__c_ratcurve(c_ratcurve)

def get_c_gr4params():
    """
    Element c_gr4params ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 264
    
    """
    return _wrapping.f90wrap_m_model__get__c_gr4params()

def set_c_gr4params(c_gr4params):
    _wrapping.f90wrap_m_model__set__c_gr4params(c_gr4params)

def get_c_rain():
    """
    Element c_rain ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 265
    
    """
    return _wrapping.f90wrap_m_model__get__c_rain()

def set_c_rain(c_rain):
    _wrapping.f90wrap_m_model__set__c_rain(c_rain)

def get_c_infil_max():
    """
    Element c_infil_max ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 266
    
    """
    return _wrapping.f90wrap_m_model__get__c_infil_max()

def set_c_infil_max(c_infil_max):
    _wrapping.f90wrap_m_model__set__c_infil_max(c_infil_max)

def get_c_ks():
    """
    Element c_ks ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 267
    
    """
    return _wrapping.f90wrap_m_model__get__c_ks()

def set_c_ks(c_ks):
    _wrapping.f90wrap_m_model__set__c_ks(c_ks)

def get_c_psif():
    """
    Element c_psif ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 268
    
    """
    return _wrapping.f90wrap_m_model__get__c_psif()

def set_c_psif(c_psif):
    _wrapping.f90wrap_m_model__set__c_psif(c_psif)

def get_c_deltatheta():
    """
    Element c_deltatheta ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 269
    
    """
    return _wrapping.f90wrap_m_model__get__c_deltatheta()

def set_c_deltatheta(c_deltatheta):
    _wrapping.f90wrap_m_model__set__c_deltatheta(c_deltatheta)

def get_c_lambda():
    """
    Element c_lambda ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 270
    
    """
    return _wrapping.f90wrap_m_model__get__c_lambda()

def set_c_lambda(c_lambda):
    _wrapping.f90wrap_m_model__set__c_lambda(c_lambda)

def get_c_cn():
    """
    Element c_cn ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 271
    
    """
    return _wrapping.f90wrap_m_model__get__c_cn()

def set_c_cn(c_cn):
    _wrapping.f90wrap_m_model__set__c_cn(c_cn)

def get_c_ptf():
    """
    Element c_ptf ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 272
    
    """
    return _wrapping.f90wrap_m_model__get__c_ptf()

def set_c_ptf(c_ptf):
    _wrapping.f90wrap_m_model__set__c_ptf(c_ptf)

def get_eps_manning():
    """
    Element eps_manning ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 275
    
    """
    return _wrapping.f90wrap_m_model__get__eps_manning()

def set_eps_manning(eps_manning):
    _wrapping.f90wrap_m_model__set__eps_manning(eps_manning)

def get_eps_bathy():
    """
    Element eps_bathy ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 276
    
    """
    return _wrapping.f90wrap_m_model__get__eps_bathy()

def set_eps_bathy(eps_bathy):
    _wrapping.f90wrap_m_model__set__eps_bathy(eps_bathy)

def get_eps_ic():
    """
    Element eps_ic ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 277
    
    """
    return _wrapping.f90wrap_m_model__get__eps_ic()

def set_eps_ic(eps_ic):
    _wrapping.f90wrap_m_model__set__eps_ic(eps_ic)

def get_eps_hydrograph():
    """
    Element eps_hydrograph ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 278
    
    """
    return _wrapping.f90wrap_m_model__get__eps_hydrograph()

def set_eps_hydrograph(eps_hydrograph):
    _wrapping.f90wrap_m_model__set__eps_hydrograph(eps_hydrograph)

def get_eps_ratcurve():
    """
    Element eps_ratcurve ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 279
    
    """
    return _wrapping.f90wrap_m_model__get__eps_ratcurve()

def set_eps_ratcurve(eps_ratcurve):
    _wrapping.f90wrap_m_model__set__eps_ratcurve(eps_ratcurve)

def get_eps_gr4params():
    """
    Element eps_gr4params ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 280
    
    """
    return _wrapping.f90wrap_m_model__get__eps_gr4params()

def set_eps_gr4params(eps_gr4params):
    _wrapping.f90wrap_m_model__set__eps_gr4params(eps_gr4params)

def get_eps_rain():
    """
    Element eps_rain ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 281
    
    """
    return _wrapping.f90wrap_m_model__get__eps_rain()

def set_eps_rain(eps_rain):
    _wrapping.f90wrap_m_model__set__eps_rain(eps_rain)

def get_eps_ks():
    """
    Element eps_ks ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 282
    
    """
    return _wrapping.f90wrap_m_model__get__eps_ks()

def set_eps_ks(eps_ks):
    _wrapping.f90wrap_m_model__set__eps_ks(eps_ks)

def get_eps_psif():
    """
    Element eps_psif ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 283
    
    """
    return _wrapping.f90wrap_m_model__get__eps_psif()

def set_eps_psif(eps_psif):
    _wrapping.f90wrap_m_model__set__eps_psif(eps_psif)

def get_eps_deltatheta():
    """
    Element eps_deltatheta ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 284
    
    """
    return _wrapping.f90wrap_m_model__get__eps_deltatheta()

def set_eps_deltatheta(eps_deltatheta):
    _wrapping.f90wrap_m_model__set__eps_deltatheta(eps_deltatheta)

def get_eps_lambdacn():
    """
    Element eps_lambdacn ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 285
    
    """
    return _wrapping.f90wrap_m_model__get__eps_lambdacn()

def set_eps_lambdacn(eps_lambdacn):
    _wrapping.f90wrap_m_model__set__eps_lambdacn(eps_lambdacn)

def get_eps_cn():
    """
    Element eps_cn ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 286
    
    """
    return _wrapping.f90wrap_m_model__get__eps_cn()

def set_eps_cn(eps_cn):
    _wrapping.f90wrap_m_model__set__eps_cn(eps_cn)

def get_eps_ptf():
    """
    Element eps_ptf ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 287
    
    """
    return _wrapping.f90wrap_m_model__get__eps_ptf()

def set_eps_ptf(eps_ptf):
    _wrapping.f90wrap_m_model__set__eps_ptf(eps_ptf)

def get_regul_manning():
    """
    Element regul_manning ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 288
    
    """
    return _wrapping.f90wrap_m_model__get__regul_manning()

def set_regul_manning(regul_manning):
    _wrapping.f90wrap_m_model__set__regul_manning(regul_manning)

def get_regul_bathy():
    """
    Element regul_bathy ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 289
    
    """
    return _wrapping.f90wrap_m_model__get__regul_bathy()

def set_regul_bathy(regul_bathy):
    _wrapping.f90wrap_m_model__set__regul_bathy(regul_bathy)

def get_regul_bathy_grad():
    """
    Element regul_bathy_grad ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 290
    
    """
    return _wrapping.f90wrap_m_model__get__regul_bathy_grad()

def set_regul_bathy_grad(regul_bathy_grad):
    _wrapping.f90wrap_m_model__set__regul_bathy_grad(regul_bathy_grad)

def get_regul_bathy_shape():
    """
    Element regul_bathy_shape ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 291
    
    """
    return _wrapping.f90wrap_m_model__get__regul_bathy_shape()

def set_regul_bathy_shape(regul_bathy_shape):
    _wrapping.f90wrap_m_model__set__regul_bathy_shape(regul_bathy_shape)

def get_regul_ic():
    """
    Element regul_ic ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 292
    
    """
    return _wrapping.f90wrap_m_model__get__regul_ic()

def set_regul_ic(regul_ic):
    _wrapping.f90wrap_m_model__set__regul_ic(regul_ic)

def get_regul_hydrograph():
    """
    Element regul_hydrograph ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 293
    
    """
    return _wrapping.f90wrap_m_model__get__regul_hydrograph()

def set_regul_hydrograph(regul_hydrograph):
    _wrapping.f90wrap_m_model__set__regul_hydrograph(regul_hydrograph)

def get_regul_ratcurve():
    """
    Element regul_ratcurve ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 294
    
    """
    return _wrapping.f90wrap_m_model__get__regul_ratcurve()

def set_regul_ratcurve(regul_ratcurve):
    _wrapping.f90wrap_m_model__set__regul_ratcurve(regul_ratcurve)

def get_regul_gr4params():
    """
    Element regul_gr4params ftype=real(rp) pytype=float
    
    
    Defined at m_sw_mono.f90 line 295
    
    """
    return _wrapping.f90wrap_m_model__get__regul_gr4params()

def set_regul_gr4params(regul_gr4params):
    _wrapping.f90wrap_m_model__set__regul_gr4params(regul_gr4params)

def get_fix_time_step_serie():
    """
    Element fix_time_step_serie ftype=integer(ip) pytype=int
    
    
    Defined at m_sw_mono.f90 line 296
    
    """
    return _wrapping.f90wrap_m_model__get__fix_time_step_serie()

def set_fix_time_step_serie(fix_time_step_serie):
    _wrapping.f90wrap_m_model__set__fix_time_step_serie(fix_time_step_serie)


_array_initialisers = [get_array_bathy_node, get_array_bathy_cell, \
    get_array_manning, get_array_manning_beta, get_array_land, \
    get_array_slope_y, get_array_slope_x]
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "m_model".')

for func in _dt_array_initialisers:
    func()
