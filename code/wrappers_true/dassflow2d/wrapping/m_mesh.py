"""
Module m_mesh


Defined at m_mesh.f90 lines 1-389

"""
from __future__ import print_function, absolute_import, division
from dassflow2d.wrapping import _wrapping
import f90wrap.runtime
import logging
from dassflow2d.wrapping.m_linear_algebra import vec2d

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("wrapping.NodeType")
class NodeType(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=nodetype)
    
    
    Defined at m_mesh.f90 lines 13-18
    
    """
    def __init__(self, dim1, dim2, handle=None):
        """
        self = Nodetype(dim1, dim2)
        
        
        Defined at m_mesh.f90 lines 352-362
        
        Parameters
        ----------
        dim1 : int
        dim2 : int
        
        Returns
        -------
        node : Nodetype
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_nodetype_initialise(dim1=dim1, dim2=dim2)
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Nodetype
        
        
        Defined at m_mesh.f90 lines 13-18
        
        Parameters
        ----------
        this : Nodetype
        	Object to be destructed
        
        
        Automatically generated destructor for nodetype
        """
        if self._alloc:
            _wrapping.f90wrap_nodetype_finalise(this=self._handle)
    
    @property
    def cell(self):
        """
        Element cell ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_nodetype__array__cell(self._handle)
        if array_handle in self._arrays:
            cell = self._arrays[array_handle]
        else:
            cell = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_nodetype__array__cell)
            self._arrays[array_handle] = cell
        return cell
    
    @cell.setter
    def cell(self, cell):
        self.cell[...] = cell
    
    @property
    def edge(self):
        """
        Element edge ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_nodetype__array__edge(self._handle)
        if array_handle in self._arrays:
            edge = self._arrays[array_handle]
        else:
            edge = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_nodetype__array__edge)
            self._arrays[array_handle] = edge
        return edge
    
    @edge.setter
    def edge(self, edge):
        self.edge[...] = edge
    
    @property
    def coord(self):
        """
        Element coord ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_mesh.f90 line 16
        
        """
        coord_handle = _wrapping.f90wrap_nodetype__get__coord(self._handle)
        if tuple(coord_handle) in self._objs:
            coord = self._objs[tuple(coord_handle)]
        else:
            coord = vec2d.from_handle(coord_handle)
            self._objs[tuple(coord_handle)] = coord
        return coord
    
    @coord.setter
    def coord(self, coord):
        coord = coord._handle
        _wrapping.f90wrap_nodetype__set__coord(self._handle, coord)
    
    @property
    def boundary(self):
        """
        Element boundary ftype=logical pytype=bool
        
        
        Defined at m_mesh.f90 line 17
        
        """
        return _wrapping.f90wrap_nodetype__get__boundary(self._handle)
    
    @boundary.setter
    def boundary(self, boundary):
        _wrapping.f90wrap_nodetype__set__boundary(self._handle, boundary)
    
    @property
    def lim(self):
        """
        Element lim ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 18
        
        """
        return _wrapping.f90wrap_nodetype__get__lim(self._handle)
    
    @lim.setter
    def lim(self, lim):
        _wrapping.f90wrap_nodetype__set__lim(self._handle, lim)
    
    def __str__(self):
        ret = ['<nodetype>{\n']
        ret.append('    cell : ')
        ret.append(repr(self.cell))
        ret.append(',\n    edge : ')
        ret.append(repr(self.edge))
        ret.append(',\n    coord : ')
        ret.append(repr(self.coord))
        ret.append(',\n    boundary : ')
        ret.append(repr(self.boundary))
        ret.append(',\n    lim : ')
        ret.append(repr(self.lim))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.NodeTypeLim")
class NodeTypeLim(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=nodetypelim)
    
    
    Defined at m_mesh.f90 lines 24-27
    
    """
    def __init__(self, ind, group, handle=None):
        """
        self = Nodetypelim(ind, group)
        
        
        Defined at m_mesh.f90 lines 365-372
        
        Parameters
        ----------
        ind : int
        group : int
        
        Returns
        -------
        nodelim : Nodetypelim
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_nodetypelim_initialise(ind=ind, group=group)
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Nodetypelim
        
        
        Defined at m_mesh.f90 lines 24-27
        
        Parameters
        ----------
        this : Nodetypelim
        	Object to be destructed
        
        
        Automatically generated destructor for nodetypelim
        """
        if self._alloc:
            _wrapping.f90wrap_nodetypelim_finalise(this=self._handle)
    
    @property
    def ind(self):
        """
        Element ind ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 25
        
        """
        return _wrapping.f90wrap_nodetypelim__get__ind(self._handle)
    
    @ind.setter
    def ind(self, ind):
        _wrapping.f90wrap_nodetypelim__set__ind(self._handle, ind)
    
    @property
    def typlim(self):
        """
        Element typlim ftype=character(len=lchar) pytype=str
        
        
        Defined at m_mesh.f90 line 26
        
        """
        return _wrapping.f90wrap_nodetypelim__get__typlim(self._handle)
    
    @typlim.setter
    def typlim(self, typlim):
        _wrapping.f90wrap_nodetypelim__set__typlim(self._handle, typlim)
    
    @property
    def group(self):
        """
        Element group ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 27
        
        """
        return _wrapping.f90wrap_nodetypelim__get__group(self._handle)
    
    @group.setter
    def group(self, group):
        _wrapping.f90wrap_nodetypelim__set__group(self._handle, group)
    
    def __str__(self):
        ret = ['<nodetypelim>{\n']
        ret.append('    ind : ')
        ret.append(repr(self.ind))
        ret.append(',\n    typlim : ')
        ret.append(repr(self.typlim))
        ret.append(',\n    group : ')
        ret.append(repr(self.group))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.CellType")
class CellType(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=celltype)
    
    
    Defined at m_mesh.f90 lines 33-43
    
    """
    def __init__(self, handle=None):
        """
        self = Celltype()
        
        
        Defined at m_mesh.f90 lines 33-43
        
        
        Returns
        -------
        this : Celltype
        	Object to be constructed
        
        
        Automatically generated constructor for celltype
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_celltype_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Celltype
        
        
        Defined at m_mesh.f90 lines 33-43
        
        Parameters
        ----------
        this : Celltype
        	Object to be destructed
        
        
        Automatically generated destructor for celltype
        """
        if self._alloc:
            _wrapping.f90wrap_celltype_finalise(this=self._handle)
    
    @property
    def node(self):
        """
        Element node ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 34
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_celltype__array__node(self._handle)
        if array_handle in self._arrays:
            node = self._arrays[array_handle]
        else:
            node = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_celltype__array__node)
            self._arrays[array_handle] = node
        return node
    
    @node.setter
    def node(self, node):
        self.node[...] = node
    
    @property
    def cell(self):
        """
        Element cell ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 35
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_celltype__array__cell(self._handle)
        if array_handle in self._arrays:
            cell = self._arrays[array_handle]
        else:
            cell = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_celltype__array__cell)
            self._arrays[array_handle] = cell
        return cell
    
    @cell.setter
    def cell(self, cell):
        self.cell[...] = cell
    
    @property
    def edge(self):
        """
        Element edge ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 36
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_celltype__array__edge(self._handle)
        if array_handle in self._arrays:
            edge = self._arrays[array_handle]
        else:
            edge = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_celltype__array__edge)
            self._arrays[array_handle] = edge
        return edge
    
    @edge.setter
    def edge(self, edge):
        self.edge[...] = edge
    
    @property
    def nbed(self):
        """
        Element nbed ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 37
        
        """
        return _wrapping.f90wrap_celltype__get__nbed(self._handle)
    
    @nbed.setter
    def nbed(self, nbed):
        _wrapping.f90wrap_celltype__set__nbed(self._handle, nbed)
    
    @property
    def boundary(self):
        """
        Element boundary ftype=logical pytype=bool
        
        
        Defined at m_mesh.f90 line 38
        
        """
        return _wrapping.f90wrap_celltype__get__boundary(self._handle)
    
    @boundary.setter
    def boundary(self, boundary):
        _wrapping.f90wrap_celltype__set__boundary(self._handle, boundary)
    
    @property
    def surf(self):
        """
        Element surf ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 39
        
        """
        return _wrapping.f90wrap_celltype__get__surf(self._handle)
    
    @surf.setter
    def surf(self, surf):
        _wrapping.f90wrap_celltype__set__surf(self._handle, surf)
    
    @property
    def invsurf(self):
        """
        Element invsurf ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 40
        
        """
        return _wrapping.f90wrap_celltype__get__invsurf(self._handle)
    
    @invsurf.setter
    def invsurf(self, invsurf):
        _wrapping.f90wrap_celltype__set__invsurf(self._handle, invsurf)
    
    @property
    def peri(self):
        """
        Element peri ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 41
        
        """
        return _wrapping.f90wrap_celltype__get__peri(self._handle)
    
    @peri.setter
    def peri(self, peri):
        _wrapping.f90wrap_celltype__set__peri(self._handle, peri)
    
    @property
    def grav(self):
        """
        Element grav ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_mesh.f90 line 42
        
        """
        grav_handle = _wrapping.f90wrap_celltype__get__grav(self._handle)
        if tuple(grav_handle) in self._objs:
            grav = self._objs[tuple(grav_handle)]
        else:
            grav = vec2d.from_handle(grav_handle)
            self._objs[tuple(grav_handle)] = grav
        return grav
    
    @grav.setter
    def grav(self, grav):
        grav = grav._handle
        _wrapping.f90wrap_celltype__set__grav(self._handle, grav)
    
    @property
    def rain(self):
        """
        Element rain ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 43
        
        """
        return _wrapping.f90wrap_celltype__get__rain(self._handle)
    
    @rain.setter
    def rain(self, rain):
        _wrapping.f90wrap_celltype__set__rain(self._handle, rain)
    
    def __str__(self):
        ret = ['<celltype>{\n']
        ret.append('    node : ')
        ret.append(repr(self.node))
        ret.append(',\n    cell : ')
        ret.append(repr(self.cell))
        ret.append(',\n    edge : ')
        ret.append(repr(self.edge))
        ret.append(',\n    nbed : ')
        ret.append(repr(self.nbed))
        ret.append(',\n    boundary : ')
        ret.append(repr(self.boundary))
        ret.append(',\n    surf : ')
        ret.append(repr(self.surf))
        ret.append(',\n    invsurf : ')
        ret.append(repr(self.invsurf))
        ret.append(',\n    peri : ')
        ret.append(repr(self.peri))
        ret.append(',\n    grav : ')
        ret.append(repr(self.grav))
        ret.append(',\n    rain : ')
        ret.append(repr(self.rain))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.CellTypeLim")
class CellTypeLim(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=celltypelim)
    
    
    Defined at m_mesh.f90 lines 49-54
    
    """
    def __init__(self, handle=None):
        """
        self = Celltypelim()
        
        
        Defined at m_mesh.f90 lines 49-54
        
        
        Returns
        -------
        this : Celltypelim
        	Object to be constructed
        
        
        Automatically generated constructor for celltypelim
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_celltypelim_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Celltypelim
        
        
        Defined at m_mesh.f90 lines 49-54
        
        Parameters
        ----------
        this : Celltypelim
        	Object to be destructed
        
        
        Automatically generated destructor for celltypelim
        """
        if self._alloc:
            _wrapping.f90wrap_celltypelim_finalise(this=self._handle)
    
    @property
    def ind(self):
        """
        Element ind ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 50
        
        """
        return _wrapping.f90wrap_celltypelim__get__ind(self._handle)
    
    @ind.setter
    def ind(self, ind):
        _wrapping.f90wrap_celltypelim__set__ind(self._handle, ind)
    
    @property
    def typlim(self):
        """
        Element typlim ftype=character(len=lchar) pytype=str
        
        
        Defined at m_mesh.f90 line 51
        
        """
        return _wrapping.f90wrap_celltypelim__get__typlim(self._handle)
    
    @typlim.setter
    def typlim(self, typlim):
        _wrapping.f90wrap_celltypelim__set__typlim(self._handle, typlim)
    
    @property
    def group(self):
        """
        Element group ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 52
        
        """
        return _wrapping.f90wrap_celltypelim__get__group(self._handle)
    
    @group.setter
    def group(self, group):
        _wrapping.f90wrap_celltypelim__set__group(self._handle, group)
    
    @property
    def cell(self):
        """
        Element cell ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 53
        
        """
        return _wrapping.f90wrap_celltypelim__get__cell(self._handle)
    
    @cell.setter
    def cell(self, cell):
        _wrapping.f90wrap_celltypelim__set__cell(self._handle, cell)
    
    @property
    def grav(self):
        """
        Element grav ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_mesh.f90 line 54
        
        """
        grav_handle = _wrapping.f90wrap_celltypelim__get__grav(self._handle)
        if tuple(grav_handle) in self._objs:
            grav = self._objs[tuple(grav_handle)]
        else:
            grav = vec2d.from_handle(grav_handle)
            self._objs[tuple(grav_handle)] = grav
        return grav
    
    @grav.setter
    def grav(self, grav):
        grav = grav._handle
        _wrapping.f90wrap_celltypelim__set__grav(self._handle, grav)
    
    def __str__(self):
        ret = ['<celltypelim>{\n']
        ret.append('    ind : ')
        ret.append(repr(self.ind))
        ret.append(',\n    typlim : ')
        ret.append(repr(self.typlim))
        ret.append(',\n    group : ')
        ret.append(repr(self.group))
        ret.append(',\n    cell : ')
        ret.append(repr(self.cell))
        ret.append(',\n    grav : ')
        ret.append(repr(self.grav))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.EdgeType")
class EdgeType(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=edgetype)
    
    
    Defined at m_mesh.f90 lines 60-72
    
    """
    def __init__(self, handle=None):
        """
        self = Edgetype()
        
        
        Defined at m_mesh.f90 lines 60-72
        
        
        Returns
        -------
        this : Edgetype
        	Object to be constructed
        
        
        Automatically generated constructor for edgetype
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_edgetype_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Edgetype
        
        
        Defined at m_mesh.f90 lines 60-72
        
        Parameters
        ----------
        this : Edgetype
        	Object to be destructed
        
        
        Automatically generated destructor for edgetype
        """
        if self._alloc:
            _wrapping.f90wrap_edgetype_finalise(this=self._handle)
    
    @property
    def node(self):
        """
        Element node ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 61
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_edgetype__array__node(self._handle)
        if array_handle in self._arrays:
            node = self._arrays[array_handle]
        else:
            node = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_edgetype__array__node)
            self._arrays[array_handle] = node
        return node
    
    @node.setter
    def node(self, node):
        self.node[...] = node
    
    @property
    def cell(self):
        """
        Element cell ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 62
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _wrapping.f90wrap_edgetype__array__cell(self._handle)
        if array_handle in self._arrays:
            cell = self._arrays[array_handle]
        else:
            cell = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    self._handle,
                                    _wrapping.f90wrap_edgetype__array__cell)
            self._arrays[array_handle] = cell
        return cell
    
    @cell.setter
    def cell(self, cell):
        self.cell[...] = cell
    
    @property
    def cell1d2d(self):
        """
        Element cell1d2d ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 63
        
        """
        return _wrapping.f90wrap_edgetype__get__cell1d2d(self._handle)
    
    @cell1d2d.setter
    def cell1d2d(self, cell1d2d):
        _wrapping.f90wrap_edgetype__set__cell1d2d(self._handle, cell1d2d)
    
    @property
    def boundary(self):
        """
        Element boundary ftype=logical pytype=bool
        
        
        Defined at m_mesh.f90 line 64
        
        """
        return _wrapping.f90wrap_edgetype__get__boundary(self._handle)
    
    @boundary.setter
    def boundary(self, boundary):
        _wrapping.f90wrap_edgetype__set__boundary(self._handle, boundary)
    
    @property
    def subdomain(self):
        """
        Element subdomain ftype=logical pytype=bool
        
        
        Defined at m_mesh.f90 line 65
        
        """
        return _wrapping.f90wrap_edgetype__get__subdomain(self._handle)
    
    @subdomain.setter
    def subdomain(self, subdomain):
        _wrapping.f90wrap_edgetype__set__subdomain(self._handle, subdomain)
    
    @property
    def lim(self):
        """
        Element lim ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 66
        
        """
        return _wrapping.f90wrap_edgetype__get__lim(self._handle)
    
    @lim.setter
    def lim(self, lim):
        _wrapping.f90wrap_edgetype__set__lim(self._handle, lim)
    
    @property
    def length(self):
        """
        Element length ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 67
        
        """
        return _wrapping.f90wrap_edgetype__get__length(self._handle)
    
    @length.setter
    def length(self, length):
        _wrapping.f90wrap_edgetype__set__length(self._handle, length)
    
    @property
    def center(self):
        """
        Element center ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_mesh.f90 line 68
        
        """
        center_handle = _wrapping.f90wrap_edgetype__get__center(self._handle)
        if tuple(center_handle) in self._objs:
            center = self._objs[tuple(center_handle)]
        else:
            center = vec2d.from_handle(center_handle)
            self._objs[tuple(center_handle)] = center
        return center
    
    @center.setter
    def center(self, center):
        center = center._handle
        _wrapping.f90wrap_edgetype__set__center(self._handle, center)
    
    @property
    def normal(self):
        """
        Element normal ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_mesh.f90 line 69
        
        """
        normal_handle = _wrapping.f90wrap_edgetype__get__normal(self._handle)
        if tuple(normal_handle) in self._objs:
            normal = self._objs[tuple(normal_handle)]
        else:
            normal = vec2d.from_handle(normal_handle)
            self._objs[tuple(normal_handle)] = normal
        return normal
    
    @normal.setter
    def normal(self, normal):
        normal = normal._handle
        _wrapping.f90wrap_edgetype__set__normal(self._handle, normal)
    
    @property
    def tangent(self):
        """
        Element tangent ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_mesh.f90 line 70
        
        """
        tangent_handle = _wrapping.f90wrap_edgetype__get__tangent(self._handle)
        if tuple(tangent_handle) in self._objs:
            tangent = self._objs[tuple(tangent_handle)]
        else:
            tangent = vec2d.from_handle(tangent_handle)
            self._objs[tuple(tangent_handle)] = tangent
        return tangent
    
    @tangent.setter
    def tangent(self, tangent):
        tangent = tangent._handle
        _wrapping.f90wrap_edgetype__set__tangent(self._handle, tangent)
    
    @property
    def vcell(self):
        """
        Element vcell ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_mesh.f90 line 71
        
        """
        vcell_handle = _wrapping.f90wrap_edgetype__get__vcell(self._handle)
        if tuple(vcell_handle) in self._objs:
            vcell = self._objs[tuple(vcell_handle)]
        else:
            vcell = vec2d.from_handle(vcell_handle)
            self._objs[tuple(vcell_handle)] = vcell
        return vcell
    
    @vcell.setter
    def vcell(self, vcell):
        vcell = vcell._handle
        _wrapping.f90wrap_edgetype__set__vcell(self._handle, vcell)
    
    def init_array_v_edge_cell(self):
        self.v_edge_cell = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_edgetype__array_getitem__v_edge_cell,
                                        _wrapping.f90wrap_edgetype__array_setitem__v_edge_cell,
                                        _wrapping.f90wrap_edgetype__array_len__v_edge_cell,
                                        """
        Element v_edge_cell ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_mesh.f90 line 72
        
        """, vec2d)
        return self.v_edge_cell
    
    def __str__(self):
        ret = ['<edgetype>{\n']
        ret.append('    node : ')
        ret.append(repr(self.node))
        ret.append(',\n    cell : ')
        ret.append(repr(self.cell))
        ret.append(',\n    cell1d2d : ')
        ret.append(repr(self.cell1d2d))
        ret.append(',\n    boundary : ')
        ret.append(repr(self.boundary))
        ret.append(',\n    subdomain : ')
        ret.append(repr(self.subdomain))
        ret.append(',\n    lim : ')
        ret.append(repr(self.lim))
        ret.append(',\n    length : ')
        ret.append(repr(self.length))
        ret.append(',\n    center : ')
        ret.append(repr(self.center))
        ret.append(',\n    normal : ')
        ret.append(repr(self.normal))
        ret.append(',\n    tangent : ')
        ret.append(repr(self.tangent))
        ret.append(',\n    vcell : ')
        ret.append(repr(self.vcell))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = [init_array_v_edge_cell]
    

@f90wrap.runtime.register_class("wrapping.EdgeTypeLim")
class EdgeTypeLim(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=edgetypelim)
    
    
    Defined at m_mesh.f90 lines 78-82
    
    """
    def __init__(self, handle=None):
        """
        self = Edgetypelim()
        
        
        Defined at m_mesh.f90 lines 78-82
        
        
        Returns
        -------
        this : Edgetypelim
        	Object to be constructed
        
        
        Automatically generated constructor for edgetypelim
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_edgetypelim_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Edgetypelim
        
        
        Defined at m_mesh.f90 lines 78-82
        
        Parameters
        ----------
        this : Edgetypelim
        	Object to be destructed
        
        
        Automatically generated destructor for edgetypelim
        """
        if self._alloc:
            _wrapping.f90wrap_edgetypelim_finalise(this=self._handle)
    
    @property
    def ind(self):
        """
        Element ind ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 79
        
        """
        return _wrapping.f90wrap_edgetypelim__get__ind(self._handle)
    
    @ind.setter
    def ind(self, ind):
        _wrapping.f90wrap_edgetypelim__set__ind(self._handle, ind)
    
    @property
    def typlim(self):
        """
        Element typlim ftype=character(len=lchar) pytype=str
        
        
        Defined at m_mesh.f90 line 80
        
        """
        return _wrapping.f90wrap_edgetypelim__get__typlim(self._handle)
    
    @typlim.setter
    def typlim(self, typlim):
        _wrapping.f90wrap_edgetypelim__set__typlim(self._handle, typlim)
    
    @property
    def group(self):
        """
        Element group ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 81
        
        """
        return _wrapping.f90wrap_edgetypelim__get__group(self._handle)
    
    @group.setter
    def group(self, group):
        _wrapping.f90wrap_edgetypelim__set__group(self._handle, group)
    
    @property
    def perio(self):
        """
        Element perio ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 82
        
        """
        return _wrapping.f90wrap_edgetypelim__get__perio(self._handle)
    
    @perio.setter
    def perio(self, perio):
        _wrapping.f90wrap_edgetypelim__set__perio(self._handle, perio)
    
    def __str__(self):
        ret = ['<edgetypelim>{\n']
        ret.append('    ind : ')
        ret.append(repr(self.ind))
        ret.append(',\n    typlim : ')
        ret.append(repr(self.typlim))
        ret.append(',\n    group : ')
        ret.append(repr(self.group))
        ret.append(',\n    perio : ')
        ret.append(repr(self.perio))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

@f90wrap.runtime.register_class("wrapping.msh")
class msh(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=msh)
    
    
    Defined at m_mesh.f90 lines 88-103
    
    """
    def __init__(self, handle=None):
        """
        self = Msh()
        
        
        Defined at m_mesh.f90 lines 375-378
        
        
        Returns
        -------
        mesh : Msh
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_msh_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Msh
        
        
        Defined at m_mesh.f90 lines 381-389
        
        Parameters
        ----------
        mesh : Msh
        
        """
        if self._alloc:
            _wrapping.f90wrap_msh_finalise(mesh=self._handle)
    
    @property
    def nn(self):
        """
        Element nn ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 89
        
        """
        return _wrapping.f90wrap_msh__get__nn(self._handle)
    
    @nn.setter
    def nn(self, nn):
        _wrapping.f90wrap_msh__set__nn(self._handle, nn)
    
    @property
    def nnb(self):
        """
        Element nnb ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 90
        
        """
        return _wrapping.f90wrap_msh__get__nnb(self._handle)
    
    @nnb.setter
    def nnb(self, nnb):
        _wrapping.f90wrap_msh__set__nnb(self._handle, nnb)
    
    @property
    def nc(self):
        """
        Element nc ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 91
        
        """
        return _wrapping.f90wrap_msh__get__nc(self._handle)
    
    @nc.setter
    def nc(self, nc):
        _wrapping.f90wrap_msh__set__nc(self._handle, nc)
    
    @property
    def ncb(self):
        """
        Element ncb ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 92
        
        """
        return _wrapping.f90wrap_msh__get__ncb(self._handle)
    
    @ncb.setter
    def ncb(self, ncb):
        _wrapping.f90wrap_msh__set__ncb(self._handle, ncb)
    
    @property
    def ne(self):
        """
        Element ne ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 93
        
        """
        return _wrapping.f90wrap_msh__get__ne(self._handle)
    
    @ne.setter
    def ne(self, ne):
        _wrapping.f90wrap_msh__set__ne(self._handle, ne)
    
    @property
    def neb(self):
        """
        Element neb ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 94
        
        """
        return _wrapping.f90wrap_msh__get__neb(self._handle)
    
    @neb.setter
    def neb(self, neb):
        _wrapping.f90wrap_msh__set__neb(self._handle, neb)
    
    def init_array_node(self):
        self.node = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_msh__array_getitem__node,
                                        _wrapping.f90wrap_msh__array_setitem__node,
                                        _wrapping.f90wrap_msh__array_len__node,
                                        """
        Element node ftype=type(nodetype) pytype=Nodetype
        
        
        Defined at m_mesh.f90 line 95
        
        """, NodeType)
        return self.node
    
    def init_array_nodeb(self):
        self.nodeb = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_msh__array_getitem__nodeb,
                                        _wrapping.f90wrap_msh__array_setitem__nodeb,
                                        _wrapping.f90wrap_msh__array_len__nodeb,
                                        """
        Element nodeb ftype=type(nodetypelim) pytype=Nodetypelim
        
        
        Defined at m_mesh.f90 line 96
        
        """, NodeTypeLim)
        return self.nodeb
    
    def init_array_cell(self):
        self.cell = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_msh__array_getitem__cell,
                                        _wrapping.f90wrap_msh__array_setitem__cell,
                                        _wrapping.f90wrap_msh__array_len__cell,
                                        """
        Element cell ftype=type(celltype) pytype=Celltype
        
        
        Defined at m_mesh.f90 line 97
        
        """, CellType)
        return self.cell
    
    def init_array_cellb(self):
        self.cellb = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_msh__array_getitem__cellb,
                                        _wrapping.f90wrap_msh__array_setitem__cellb,
                                        _wrapping.f90wrap_msh__array_len__cellb,
                                        """
        Element cellb ftype=type(celltypelim) pytype=Celltypelim
        
        
        Defined at m_mesh.f90 line 98
        
        """, CellTypeLim)
        return self.cellb
    
    def init_array_edge(self):
        self.edge = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_msh__array_getitem__edge,
                                        _wrapping.f90wrap_msh__array_setitem__edge,
                                        _wrapping.f90wrap_msh__array_len__edge,
                                        """
        Element edge ftype=type(edgetype) pytype=Edgetype
        
        
        Defined at m_mesh.f90 line 99
        
        """, EdgeType)
        return self.edge
    
    def init_array_edgeb(self):
        self.edgeb = f90wrap.runtime.FortranDerivedTypeArray(self,
                                        _wrapping.f90wrap_msh__array_getitem__edgeb,
                                        _wrapping.f90wrap_msh__array_setitem__edgeb,
                                        _wrapping.f90wrap_msh__array_len__edgeb,
                                        """
        Element edgeb ftype=type(edgetypelim) pytype=Edgetypelim
        
        
        Defined at m_mesh.f90 line 100
        
        """, EdgeTypeLim)
        return self.edgeb
    
    @property
    def file_name(self):
        """
        Element file_name ftype=character(len=lchar) pytype=str
        
        
        Defined at m_mesh.f90 line 101
        
        """
        return _wrapping.f90wrap_msh__get__file_name(self._handle)
    
    @file_name.setter
    def file_name(self, file_name):
        _wrapping.f90wrap_msh__set__file_name(self._handle, file_name)
    
    @property
    def scal(self):
        """
        Element scal ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 102
        
        """
        return _wrapping.f90wrap_msh__get__scal(self._handle)
    
    @scal.setter
    def scal(self, scal):
        _wrapping.f90wrap_msh__set__scal(self._handle, scal)
    
    @property
    def surf(self):
        """
        Element surf ftype=real(rp) pytype=float
        
        
        Defined at m_mesh.f90 line 103
        
        """
        return _wrapping.f90wrap_msh__get__surf(self._handle)
    
    @surf.setter
    def surf(self, surf):
        _wrapping.f90wrap_msh__set__surf(self._handle, surf)
    
    def __str__(self):
        ret = ['<msh>{\n']
        ret.append('    nn : ')
        ret.append(repr(self.nn))
        ret.append(',\n    nnb : ')
        ret.append(repr(self.nnb))
        ret.append(',\n    nc : ')
        ret.append(repr(self.nc))
        ret.append(',\n    ncb : ')
        ret.append(repr(self.ncb))
        ret.append(',\n    ne : ')
        ret.append(repr(self.ne))
        ret.append(',\n    neb : ')
        ret.append(repr(self.neb))
        ret.append(',\n    file_name : ')
        ret.append(repr(self.file_name))
        ret.append(',\n    scal : ')
        ret.append(repr(self.scal))
        ret.append(',\n    surf : ')
        ret.append(repr(self.surf))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = [init_array_node, init_array_nodeb, init_array_cell, \
        init_array_cellb, init_array_edge, init_array_edgeb]
    

@f90wrap.runtime.register_class("wrapping.point_in_mesh")
class point_in_mesh(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=point_in_mesh)
    
    
    Defined at m_mesh.f90 lines 110-112
    
    """
    def __init__(self, handle=None):
        """
        self = Point_In_Mesh()
        
        
        Defined at m_mesh.f90 lines 110-112
        
        
        Returns
        -------
        this : Point_In_Mesh
        	Object to be constructed
        
        
        Automatically generated constructor for point_in_mesh
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_point_in_mesh_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Point_In_Mesh
        
        
        Defined at m_mesh.f90 lines 110-112
        
        Parameters
        ----------
        this : Point_In_Mesh
        	Object to be destructed
        
        
        Automatically generated destructor for point_in_mesh
        """
        if self._alloc:
            _wrapping.f90wrap_point_in_mesh_finalise(this=self._handle)
    
    @property
    def coord(self):
        """
        Element coord ftype=type(vec2d) pytype=Vec2D
        
        
        Defined at m_mesh.f90 line 111
        
        """
        coord_handle = _wrapping.f90wrap_point_in_mesh__get__coord(self._handle)
        if tuple(coord_handle) in self._objs:
            coord = self._objs[tuple(coord_handle)]
        else:
            coord = vec2d.from_handle(coord_handle)
            self._objs[tuple(coord_handle)] = coord
        return coord
    
    @coord.setter
    def coord(self, coord):
        coord = coord._handle
        _wrapping.f90wrap_point_in_mesh__set__coord(self._handle, coord)
    
    @property
    def cell(self):
        """
        Element cell ftype=integer(ip) pytype=int
        
        
        Defined at m_mesh.f90 line 112
        
        """
        return _wrapping.f90wrap_point_in_mesh__get__cell(self._handle)
    
    @cell.setter
    def cell(self, cell):
        _wrapping.f90wrap_point_in_mesh__set__cell(self._handle, cell)
    
    def __str__(self):
        ret = ['<point_in_mesh>{\n']
        ret.append('    coord : ')
        ret.append(repr(self.coord))
        ret.append(',\n    cell : ')
        ret.append(repr(self.cell))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def display_mesh_cell(self):
    """
    display_mesh_cell(self)
    
    
    Defined at m_mesh.f90 lines 126-134
    
    Parameters
    ----------
    mesh : Msh
    
    """
    _wrapping.f90wrap_display_mesh_cell(mesh=self._handle)

def calc_cells_connectivity(self):
    """
    calc_cells_connectivity(self)
    
    
    Defined at m_mesh.f90 lines 290-349
    
    Parameters
    ----------
    mesh : Msh
    
    """
    _wrapping.f90wrap_calc_cells_connectivity(mesh=self._handle)

def get_maxed():
    """
    Element maxed ftype=integer pytype=int
    
    
    Defined at m_mesh.f90 line 8
    
    """
    return _wrapping.f90wrap_m_mesh__get__maxed()

maxed = get_maxed()


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "m_mesh".')

for func in _dt_array_initialisers:
    func()
