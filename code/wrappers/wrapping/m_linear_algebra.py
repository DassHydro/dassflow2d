"""
Module m_linear_algebra


Defined at m_linear_algebra.f90 lines 1-47

"""
from __future__ import print_function, absolute_import, division
import _wrapping
import f90wrap.runtime
import logging

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("wrapping.vec2d")
class vec2d(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=vec2d)
    
    
    Defined at m_linear_algebra.f90 lines 7-9
    
    """
    def __init__(self, handle=None):
        """
        self = Vec2D()
        
        
        Defined at m_linear_algebra.f90 lines 7-9
        
        
        Returns
        -------
        this : Vec2D
        	Object to be constructed
        
        
        Automatically generated constructor for vec2d
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_vec2d_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Vec2D
        
        
        Defined at m_linear_algebra.f90 lines 7-9
        
        Parameters
        ----------
        this : Vec2D
        	Object to be destructed
        
        
        Automatically generated destructor for vec2d
        """
        if self._alloc:
            _wrapping.f90wrap_vec2d_finalise(this=self._handle)
    
    @property
    def x(self):
        """
        Element x ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 8
        
        """
        return _wrapping.f90wrap_vec2d__get__x(self._handle)
    
    @x.setter
    def x(self, x):
        _wrapping.f90wrap_vec2d__set__x(self._handle, x)
    
    @property
    def y(self):
        """
        Element y ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 9
        
        """
        return _wrapping.f90wrap_vec2d__get__y(self._handle)
    
    @y.setter
    def y(self, y):
        _wrapping.f90wrap_vec2d__set__y(self._handle, y)
    
    def __str__(self):
        ret = ['<vec2d>{\n']
        ret.append('    x : ')
        ret.append(repr(self.x))
        ret.append(',\n    y : ')
        ret.append(repr(self.y))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.vec3d")
class vec3d(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=vec3d)
    
    
    Defined at m_linear_algebra.f90 lines 14-15
    
    """
    def __init__(self, handle=None):
        """
        self = Vec3D()
        
        
        Defined at m_linear_algebra.f90 lines 14-15
        
        
        Returns
        -------
        this : Vec3D
        	Object to be constructed
        
        
        Automatically generated constructor for vec3d
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_vec3d_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Vec3D
        
        
        Defined at m_linear_algebra.f90 lines 14-15
        
        Parameters
        ----------
        this : Vec3D
        	Object to be destructed
        
        
        Automatically generated destructor for vec3d
        """
        if self._alloc:
            _wrapping.f90wrap_vec3d_finalise(this=self._handle)
    
    @property
    def x(self):
        """
        Element x ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 15
        
        """
        return _wrapping.f90wrap_vec3d__get__x(self._handle)
    
    @x.setter
    def x(self, x):
        _wrapping.f90wrap_vec3d__set__x(self._handle, x)
    
    @property
    def y(self):
        """
        Element y ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 15
        
        """
        return _wrapping.f90wrap_vec3d__get__y(self._handle)
    
    @y.setter
    def y(self, y):
        _wrapping.f90wrap_vec3d__set__y(self._handle, y)
    
    @property
    def z(self):
        """
        Element z ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 15
        
        """
        return _wrapping.f90wrap_vec3d__get__z(self._handle)
    
    @z.setter
    def z(self, z):
        _wrapping.f90wrap_vec3d__set__z(self._handle, z)
    
    def __str__(self):
        ret = ['<vec3d>{\n']
        ret.append('    x : ')
        ret.append(repr(self.x))
        ret.append(',\n    y : ')
        ret.append(repr(self.y))
        ret.append(',\n    z : ')
        ret.append(repr(self.z))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.tens2d")
class tens2d(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=tens2d)
    
    
    Defined at m_linear_algebra.f90 lines 20-21
    
    """
    def __init__(self, handle=None):
        """
        self = Tens2D()
        
        
        Defined at m_linear_algebra.f90 lines 20-21
        
        
        Returns
        -------
        this : Tens2D
        	Object to be constructed
        
        
        Automatically generated constructor for tens2d
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_tens2d_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Tens2D
        
        
        Defined at m_linear_algebra.f90 lines 20-21
        
        Parameters
        ----------
        this : Tens2D
        	Object to be destructed
        
        
        Automatically generated destructor for tens2d
        """
        if self._alloc:
            _wrapping.f90wrap_tens2d_finalise(this=self._handle)
    
    @property
    def xx(self):
        """
        Element xx ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 21
        
        """
        return _wrapping.f90wrap_tens2d__get__xx(self._handle)
    
    @xx.setter
    def xx(self, xx):
        _wrapping.f90wrap_tens2d__set__xx(self._handle, xx)
    
    @property
    def xy(self):
        """
        Element xy ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 21
        
        """
        return _wrapping.f90wrap_tens2d__get__xy(self._handle)
    
    @xy.setter
    def xy(self, xy):
        _wrapping.f90wrap_tens2d__set__xy(self._handle, xy)
    
    @property
    def yx(self):
        """
        Element yx ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 21
        
        """
        return _wrapping.f90wrap_tens2d__get__yx(self._handle)
    
    @yx.setter
    def yx(self, yx):
        _wrapping.f90wrap_tens2d__set__yx(self._handle, yx)
    
    @property
    def yy(self):
        """
        Element yy ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 21
        
        """
        return _wrapping.f90wrap_tens2d__get__yy(self._handle)
    
    @yy.setter
    def yy(self, yy):
        _wrapping.f90wrap_tens2d__set__yy(self._handle, yy)
    
    def __str__(self):
        ret = ['<tens2d>{\n']
        ret.append('    xx : ')
        ret.append(repr(self.xx))
        ret.append(',\n    xy : ')
        ret.append(repr(self.xy))
        ret.append(',\n    yx : ')
        ret.append(repr(self.yx))
        ret.append(',\n    yy : ')
        ret.append(repr(self.yy))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.tens3d")
class tens3d(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=tens3d)
    
    
    Defined at m_linear_algebra.f90 lines 27-28
    
    """
    def __init__(self, handle=None):
        """
        self = Tens3D()
        
        
        Defined at m_linear_algebra.f90 lines 27-28
        
        
        Returns
        -------
        this : Tens3D
        	Object to be constructed
        
        
        Automatically generated constructor for tens3d
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_tens3d_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Tens3D
        
        
        Defined at m_linear_algebra.f90 lines 27-28
        
        Parameters
        ----------
        this : Tens3D
        	Object to be destructed
        
        
        Automatically generated destructor for tens3d
        """
        if self._alloc:
            _wrapping.f90wrap_tens3d_finalise(this=self._handle)
    
    @property
    def xx(self):
        """
        Element xx ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 28
        
        """
        return _wrapping.f90wrap_tens3d__get__xx(self._handle)
    
    @xx.setter
    def xx(self, xx):
        _wrapping.f90wrap_tens3d__set__xx(self._handle, xx)
    
    @property
    def xy(self):
        """
        Element xy ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 28
        
        """
        return _wrapping.f90wrap_tens3d__get__xy(self._handle)
    
    @xy.setter
    def xy(self, xy):
        _wrapping.f90wrap_tens3d__set__xy(self._handle, xy)
    
    @property
    def xz(self):
        """
        Element xz ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 28
        
        """
        return _wrapping.f90wrap_tens3d__get__xz(self._handle)
    
    @xz.setter
    def xz(self, xz):
        _wrapping.f90wrap_tens3d__set__xz(self._handle, xz)
    
    @property
    def yx(self):
        """
        Element yx ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 28
        
        """
        return _wrapping.f90wrap_tens3d__get__yx(self._handle)
    
    @yx.setter
    def yx(self, yx):
        _wrapping.f90wrap_tens3d__set__yx(self._handle, yx)
    
    @property
    def yy(self):
        """
        Element yy ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 28
        
        """
        return _wrapping.f90wrap_tens3d__get__yy(self._handle)
    
    @yy.setter
    def yy(self, yy):
        _wrapping.f90wrap_tens3d__set__yy(self._handle, yy)
    
    @property
    def yz(self):
        """
        Element yz ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 28
        
        """
        return _wrapping.f90wrap_tens3d__get__yz(self._handle)
    
    @yz.setter
    def yz(self, yz):
        _wrapping.f90wrap_tens3d__set__yz(self._handle, yz)
    
    @property
    def zx(self):
        """
        Element zx ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 28
        
        """
        return _wrapping.f90wrap_tens3d__get__zx(self._handle)
    
    @zx.setter
    def zx(self, zx):
        _wrapping.f90wrap_tens3d__set__zx(self._handle, zx)
    
    @property
    def zy(self):
        """
        Element zy ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 28
        
        """
        return _wrapping.f90wrap_tens3d__get__zy(self._handle)
    
    @zy.setter
    def zy(self, zy):
        _wrapping.f90wrap_tens3d__set__zy(self._handle, zy)
    
    @property
    def zz(self):
        """
        Element zz ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 28
        
        """
        return _wrapping.f90wrap_tens3d__get__zz(self._handle)
    
    @zz.setter
    def zz(self, zz):
        _wrapping.f90wrap_tens3d__set__zz(self._handle, zz)
    
    def __str__(self):
        ret = ['<tens3d>{\n']
        ret.append('    xx : ')
        ret.append(repr(self.xx))
        ret.append(',\n    xy : ')
        ret.append(repr(self.xy))
        ret.append(',\n    xz : ')
        ret.append(repr(self.xz))
        ret.append(',\n    yx : ')
        ret.append(repr(self.yx))
        ret.append(',\n    yy : ')
        ret.append(repr(self.yy))
        ret.append(',\n    yz : ')
        ret.append(repr(self.yz))
        ret.append(',\n    zx : ')
        ret.append(repr(self.zx))
        ret.append(',\n    zy : ')
        ret.append(repr(self.zy))
        ret.append(',\n    zz : ')
        ret.append(repr(self.zz))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.sys_lin_full")
class sys_lin_full(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=sys_lin_full)
    
    
    Defined at m_linear_algebra.f90 lines 32-40
    
    """
    def __init__(self, handle=None):
        """
        self = Sys_Lin_Full()
        
        
        Defined at m_linear_algebra.f90 lines 32-40
        
        
        Returns
        -------
        this : Sys_Lin_Full
        	Object to be constructed
        
        
        Automatically generated constructor for sys_lin_full
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_sys_lin_full_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Sys_Lin_Full
        
        
        Defined at m_linear_algebra.f90 lines 32-40
        
        Parameters
        ----------
        this : Sys_Lin_Full
        	Object to be destructed
        
        
        Automatically generated destructor for sys_lin_full
        """
        if self._alloc:
            _wrapping.f90wrap_sys_lin_full_finalise(this=self._handle)
    
    @property
    def n(self):
        """
        Element n ftype=integer(ip) pytype=int
        
        
        Defined at m_linear_algebra.f90 line 34
        
        """
        return _wrapping.f90wrap_sys_lin_full__get__n(self._handle)
    
    @n.setter
    def n(self, n):
        _wrapping.f90wrap_sys_lin_full__set__n(self._handle, n)
    
    @property
    def a(self):
        """
        Element a ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 36
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_sys_lin_full__array__a(self._handle)
        if array_handle in self._arrays:
            a = self._arrays[array_handle]
        else:
            a = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_sys_lin_full__array__a)
            self._arrays[array_handle] = a
        return a
    
    @a.setter
    def a(self, a):
        self.a[...] = a
    
    @property
    def rhs(self):
        """
        Element rhs ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 38
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_sys_lin_full__array__rhs(self._handle)
        if array_handle in self._arrays:
            rhs = self._arrays[array_handle]
        else:
            rhs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_sys_lin_full__array__rhs)
            self._arrays[array_handle] = rhs
        return rhs
    
    @rhs.setter
    def rhs(self, rhs):
        self.rhs[...] = rhs
    
    @property
    def det(self):
        """
        Element det ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 40
        
        """
        return _wrapping.f90wrap_sys_lin_full__get__det(self._handle)
    
    @det.setter
    def det(self, det):
        _wrapping.f90wrap_sys_lin_full__set__det(self._handle, det)
    
    def __str__(self):
        ret = ['<sys_lin_full>{\n']
        ret.append('    n : ')
        ret.append(repr(self.n))
        ret.append(',\n    a : ')
        ret.append(repr(self.a))
        ret.append(',\n    rhs : ')
        ret.append(repr(self.rhs))
        ret.append(',\n    det : ')
        ret.append(repr(self.det))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.sys_lin_sparse")
class sys_lin_sparse(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=sys_lin_sparse)
    
    
    Defined at m_linear_algebra.f90 lines 43-46
    
    """
    def __init__(self, handle=None):
        """
        self = Sys_Lin_Sparse()
        
        
        Defined at m_linear_algebra.f90 lines 43-46
        
        
        Returns
        -------
        this : Sys_Lin_Sparse
        	Object to be constructed
        
        
        Automatically generated constructor for sys_lin_sparse
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_sys_lin_sparse_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Sys_Lin_Sparse
        
        
        Defined at m_linear_algebra.f90 lines 43-46
        
        Parameters
        ----------
        this : Sys_Lin_Sparse
        	Object to be destructed
        
        
        Automatically generated destructor for sys_lin_sparse
        """
        if self._alloc:
            _wrapping.f90wrap_sys_lin_sparse_finalise(this=self._handle)
    
    @property
    def n(self):
        """
        Element n ftype=integer(ip) pytype=int
        
        
        Defined at m_linear_algebra.f90 line 44
        
        """
        return _wrapping.f90wrap_sys_lin_sparse__get__n(self._handle)
    
    @n.setter
    def n(self, n):
        _wrapping.f90wrap_sys_lin_sparse__set__n(self._handle, n)
    
    @property
    def nz(self):
        """
        Element nz ftype=integer(ip) pytype=int
        
        
        Defined at m_linear_algebra.f90 line 44
        
        """
        return _wrapping.f90wrap_sys_lin_sparse__get__nz(self._handle)
    
    @nz.setter
    def nz(self, nz):
        _wrapping.f90wrap_sys_lin_sparse__set__nz(self._handle, nz)
    
    @property
    def ia(self):
        """
        Element ia ftype=integer(ip) pytype=int
        
        
        Defined at m_linear_algebra.f90 line 45
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_sys_lin_sparse__array__ia(self._handle)
        if array_handle in self._arrays:
            ia = self._arrays[array_handle]
        else:
            ia = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_sys_lin_sparse__array__ia)
            self._arrays[array_handle] = ia
        return ia
    
    @ia.setter
    def ia(self, ia):
        self.ia[...] = ia
    
    @property
    def ja(self):
        """
        Element ja ftype=integer(ip) pytype=int
        
        
        Defined at m_linear_algebra.f90 line 45
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_sys_lin_sparse__array__ja(self._handle)
        if array_handle in self._arrays:
            ja = self._arrays[array_handle]
        else:
            ja = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_sys_lin_sparse__array__ja)
            self._arrays[array_handle] = ja
        return ja
    
    @ja.setter
    def ja(self, ja):
        self.ja[...] = ja
    
    @property
    def swap(self):
        """
        Element swap ftype=integer(ip) pytype=int
        
        
        Defined at m_linear_algebra.f90 line 45
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_sys_lin_sparse__array__swap(self._handle)
        if array_handle in self._arrays:
            swap = self._arrays[array_handle]
        else:
            swap = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_sys_lin_sparse__array__swap)
            self._arrays[array_handle] = swap
        return swap
    
    @swap.setter
    def swap(self, swap):
        self.swap[...] = swap
    
    @property
    def a(self):
        """
        Element a ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 46
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_sys_lin_sparse__array__a(self._handle)
        if array_handle in self._arrays:
            a = self._arrays[array_handle]
        else:
            a = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_sys_lin_sparse__array__a)
            self._arrays[array_handle] = a
        return a
    
    @a.setter
    def a(self, a):
        self.a[...] = a
    
    @property
    def x(self):
        """
        Element x ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 46
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_sys_lin_sparse__array__x(self._handle)
        if array_handle in self._arrays:
            x = self._arrays[array_handle]
        else:
            x = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_sys_lin_sparse__array__x)
            self._arrays[array_handle] = x
        return x
    
    @x.setter
    def x(self, x):
        self.x[...] = x
    
    @property
    def x0(self):
        """
        Element x0 ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 46
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_sys_lin_sparse__array__x0(self._handle)
        if array_handle in self._arrays:
            x0 = self._arrays[array_handle]
        else:
            x0 = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_sys_lin_sparse__array__x0)
            self._arrays[array_handle] = x0
        return x0
    
    @x0.setter
    def x0(self, x0):
        self.x0[...] = x0
    
    @property
    def rhs(self):
        """
        Element rhs ftype=real(rp) pytype=float
        
        
        Defined at m_linear_algebra.f90 line 46
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_sys_lin_sparse__array__rhs(self._handle)
        if array_handle in self._arrays:
            rhs = self._arrays[array_handle]
        else:
            rhs = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_sys_lin_sparse__array__rhs)
            self._arrays[array_handle] = rhs
        return rhs
    
    @rhs.setter
    def rhs(self, rhs):
        self.rhs[...] = rhs
    
    def __str__(self):
        ret = ['<sys_lin_sparse>{\n']
        ret.append('    n : ')
        ret.append(repr(self.n))
        ret.append(',\n    nz : ')
        ret.append(repr(self.nz))
        ret.append(',\n    ia : ')
        ret.append(repr(self.ia))
        ret.append(',\n    ja : ')
        ret.append(repr(self.ja))
        ret.append(',\n    swap : ')
        ret.append(repr(self.swap))
        ret.append(',\n    a : ')
        ret.append(repr(self.a))
        ret.append(',\n    x : ')
        ret.append(repr(self.x))
        ret.append(',\n    x0 : ')
        ret.append(repr(self.x0))
        ret.append(',\n    rhs : ')
        ret.append(repr(self.rhs))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module \
        "m_linear_algebra".')

for func in _dt_array_initialisers:
    func()
