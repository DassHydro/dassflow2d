"""
Module call_model


Defined at call_run_model.f90 lines 1-794

"""
from __future__ import print_function, absolute_import, division
import _wrapping
import f90wrap.runtime
import logging
from wrapping.m_model import infiltration_data
from wrapping.m_mesh import msh
from wrapping.m_model import input_data
from wrapping.m_model import friction_data
from wrapping.m_model import param_model
from wrapping.m_model import Input_Param
from wrapping.m_model import unk

_arrays = {}
_objs = {}

@f90wrap.runtime.register_class("wrapping.Model")
class Model(f90wrap.runtime.FortranDerivedType):
    """
    Type(name=model)
    
    
    Defined at call_run_model.f90 lines 13-23
    
    """
    def __init__(self, handle=None):
        """
        self = Model()
        
        
        Defined at call_run_model.f90 lines 28-31
        
        
        Returns
        -------
        mdl : Model
        
        """
        f90wrap.runtime.FortranDerivedType.__init__(self)
        result = _wrapping.f90wrap_model_initialise()
        self._handle = result[0] if isinstance(result, tuple) else result
    
    def __del__(self):
        """
        Destructor for class Model
        
        
        Defined at call_run_model.f90 lines 35-45
        
        Parameters
        ----------
        mdl : Model
        
        """
        if self._alloc:
            _wrapping.f90wrap_model_finalise(mdl=self._handle)
    
    @property
    def mesh(self):
        """
        Element mesh ftype=type(msh) pytype=Msh
        
        
        Defined at call_run_model.f90 line 14
        
        """
        mesh_handle = _wrapping.f90wrap_model__get__mesh(self._handle)
        if tuple(mesh_handle) in self._objs:
            mesh = self._objs[tuple(mesh_handle)]
        else:
            mesh = msh.from_handle(mesh_handle)
            self._objs[tuple(mesh_handle)] = mesh
        return mesh
    
    @mesh.setter
    def mesh(self, mesh):
        mesh = mesh._handle
        _wrapping.f90wrap_model__set__mesh(self._handle, mesh)
    
    @property
    def dof0(self):
        """
        Element dof0 ftype=type(unk) pytype=Unk
        
        
        Defined at call_run_model.f90 line 15
        
        """
        dof0_handle = _wrapping.f90wrap_model__get__dof0(self._handle)
        if tuple(dof0_handle) in self._objs:
            dof0 = self._objs[tuple(dof0_handle)]
        else:
            dof0 = unk.from_handle(dof0_handle)
            self._objs[tuple(dof0_handle)] = dof0
        return dof0
    
    @dof0.setter
    def dof0(self, dof0):
        dof0 = dof0._handle
        _wrapping.f90wrap_model__set__dof0(self._handle, dof0)
    
    @property
    def dof(self):
        """
        Element dof ftype=type(unk) pytype=Unk
        
        
        Defined at call_run_model.f90 line 16
        
        """
        dof_handle = _wrapping.f90wrap_model__get__dof(self._handle)
        if tuple(dof_handle) in self._objs:
            dof = self._objs[tuple(dof_handle)]
        else:
            dof = unk.from_handle(dof_handle)
            self._objs[tuple(dof_handle)] = dof
        return dof
    
    @dof.setter
    def dof(self, dof):
        dof = dof._handle
        _wrapping.f90wrap_model__set__dof(self._handle, dof)
    
    @property
    def cost(self):
        """
        Element cost ftype=real(rp) pytype=float
        
        
        Defined at call_run_model.f90 line 17
        
        """
        return _wrapping.f90wrap_model__get__cost(self._handle)
    
    @cost.setter
    def cost(self, cost):
        _wrapping.f90wrap_model__set__cost(self._handle, cost)
    
    @property
    def param(self):
        """
        Element param ftype=type(input_param) pytype=Input_Param
        
        
        Defined at call_run_model.f90 line 18
        
        """
        param_handle = _wrapping.f90wrap_model__get__param(self._handle)
        if tuple(param_handle) in self._objs:
            param = self._objs[tuple(param_handle)]
        else:
            param = Input_Param.from_handle(param_handle)
            self._objs[tuple(param_handle)] = param
        return param
    
    @param.setter
    def param(self, param):
        param = param._handle
        _wrapping.f90wrap_model__set__param(self._handle, param)
    
    @property
    def my_friction(self):
        """
        Element my_friction ftype=type(friction_data) pytype=Friction_Data
        
        
        Defined at call_run_model.f90 line 19
        
        """
        my_friction_handle = _wrapping.f90wrap_model__get__my_friction(self._handle)
        if tuple(my_friction_handle) in self._objs:
            my_friction = self._objs[tuple(my_friction_handle)]
        else:
            my_friction = friction_data.from_handle(my_friction_handle)
            self._objs[tuple(my_friction_handle)] = my_friction
        return my_friction
    
    @my_friction.setter
    def my_friction(self, my_friction):
        my_friction = my_friction._handle
        _wrapping.f90wrap_model__set__my_friction(self._handle, my_friction)
    
    @property
    def my_infiltration(self):
        """
        Element my_infiltration ftype=type(infiltration_data) pytype=Infiltration_Data
        
        
        Defined at call_run_model.f90 line 20
        
        """
        my_infiltration_handle = \
            _wrapping.f90wrap_model__get__my_infiltration(self._handle)
        if tuple(my_infiltration_handle) in self._objs:
            my_infiltration = self._objs[tuple(my_infiltration_handle)]
        else:
            my_infiltration = infiltration_data.from_handle(my_infiltration_handle)
            self._objs[tuple(my_infiltration_handle)] = my_infiltration
        return my_infiltration
    
    @my_infiltration.setter
    def my_infiltration(self, my_infiltration):
        my_infiltration = my_infiltration._handle
        _wrapping.f90wrap_model__set__my_infiltration(self._handle, my_infiltration)
    
    @property
    def my_phys_desc(self):
        """
        Element my_phys_desc ftype=type(input_data) pytype=Input_Data
        
        
        Defined at call_run_model.f90 line 21
        
        """
        my_phys_desc_handle = _wrapping.f90wrap_model__get__my_phys_desc(self._handle)
        if tuple(my_phys_desc_handle) in self._objs:
            my_phys_desc = self._objs[tuple(my_phys_desc_handle)]
        else:
            my_phys_desc = input_data.from_handle(my_phys_desc_handle)
            self._objs[tuple(my_phys_desc_handle)] = my_phys_desc
        return my_phys_desc
    
    @my_phys_desc.setter
    def my_phys_desc(self, my_phys_desc):
        my_phys_desc = my_phys_desc._handle
        _wrapping.f90wrap_model__set__my_phys_desc(self._handle, my_phys_desc)
    
    @property
    def my_param_model(self):
        """
        Element my_param_model ftype=type(param_model) pytype=Param_Model
        
        
        Defined at call_run_model.f90 line 22
        
        """
        my_param_model_handle = \
            _wrapping.f90wrap_model__get__my_param_model(self._handle)
        if tuple(my_param_model_handle) in self._objs:
            my_param_model = self._objs[tuple(my_param_model_handle)]
        else:
            my_param_model = param_model.from_handle(my_param_model_handle)
            self._objs[tuple(my_param_model_handle)] = my_param_model
        return my_param_model
    
    @my_param_model.setter
    def my_param_model(self, my_param_model):
        my_param_model = my_param_model._handle
        _wrapping.f90wrap_model__set__my_param_model(self._handle, my_param_model)
    
    def __str__(self):
        ret = ['<model>{\n']
        ret.append('    mesh : ')
        ret.append(repr(self.mesh))
        ret.append(',\n    dof0 : ')
        ret.append(repr(self.dof0))
        ret.append(',\n    dof : ')
        ret.append(repr(self.dof))
        ret.append(',\n    cost : ')
        ret.append(repr(self.cost))
        ret.append(',\n    param : ')
        ret.append(repr(self.param))
        ret.append(',\n    my_friction : ')
        ret.append(repr(self.my_friction))
        ret.append(',\n    my_infiltration : ')
        ret.append(repr(self.my_infiltration))
        ret.append(',\n    my_phys_desc : ')
        ret.append(repr(self.my_phys_desc))
        ret.append(',\n    my_param_model : ')
        ret.append(repr(self.my_param_model))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

def init_solver(self):
    """
    init_solver(self)
    
    
    Defined at call_run_model.f90 lines 47-75
    
    Parameters
    ----------
    mdl : Model
    
    ++++++++++++++++++++++++++++++++++++++++++++++++
    FORTRAN DOCUMENTATION
    ++++++++++++++++++++++++++++++++++++++++++++++++
     init_solver
    
     initialize solver
    
     -------------------
     details
     -------------------
    
     calls Init_Linear_Solver(model) from m_numeric.f90
    """
    _wrapping.f90wrap_init_solver(mdl=self._handle)

def init_friction(self):
    """
    init_friction(self)
    
    
    Defined at call_run_model.f90 lines 77-88
    
    Parameters
    ----------
    mdl : Model
    
    > init_friction solver
    call f90 routine friction_initialise, infiltration_initialise,
    blablas
    blablasS
    blablassss
    """
    _wrapping.f90wrap_init_friction(mdl=self._handle)

def init_fortran(self):
    """
    init_fortran(self)
    
    
    Defined at call_run_model.f90 lines 90-117
    
    Parameters
    ----------
    mdl : Model
    
    ++++++++++++++++++++++++++++++++++++++++++++++++
    FORTRAN DOCUMENTATION
    ++++++++++++++++++++++++++++++++++++++++++++++++
    init_fortran
    
    Initialize fortran derived types
    - bc
    - land use
    - eventual bc value(bc%hyd, bc%rat, bc%zpresc, bc%hpresc, etc...)
     (manning, manning_beta, bathy_cell are ALREADY defined with mesh)
    
    -------------------
    details
    -------------------
    
    call not interfaced f90 subroutines
    - Initial() from initialization.f90
    - Init_Schemes() from m_numeric.f90
    """
    _wrapping.f90wrap_init_fortran(mdl=self._handle)

def init_all(self):
    """
    init_all(self)
    
    
    Defined at call_run_model.f90 lines 119-126
    
    Parameters
    ----------
    mdl : Model
    
    > Init_all solver
     call f90 routine friction_initialise, infiltration_initialise,
    """
    _wrapping.f90wrap_init_all(mdl=self._handle)

def init_back(self):
    """
    init_back(self)
    
    
    Defined at call_run_model.f90 lines 132-136
    
    Parameters
    ----------
    mdl : Model
    
    """
    _wrapping.f90wrap_init_back(mdl=self._handle)

def run(self, arg=None):
    """
    run(self[, arg])
    
    
    Defined at call_run_model.f90 lines 141-193
    
    Parameters
    ----------
    mdl : Model
    arg : str
    
    ++++++++++++++++++++++++++++++++++++++++++++++++
    FORTRAN DOCUMENTATION
    ++++++++++++++++++++++++++++++++++++++++++++++++
    run
    -------------------
    details
    -------------------
    _COMMENT
    _COMMENTrun the model depending on arg value, the different cases are treated
    _COMMENT
    _COMMENT- 'direct' : direct simulation
    _COMMENT
    _COMMENT- 'testadj' : for debug ??
    _COMMENT
    _COMMENT- 'grad' : for returning grad_cost for sensibility analisys
    _COMMENT
    _COMMENT- 'min' : perform a minimization
    _COMMENT
    _COMMENT- else : f90wrap_abort
    """
    _wrapping.f90wrap_run(mdl=self._handle, arg=arg)

def clean_model(self):
    """
    clean_model(self)
    
    
    Defined at call_run_model.f90 lines 196-210
    
    Parameters
    ----------
    mdl : Model
    
    ++++++++++++++++++++++++++++++++++++++++++++++++
    FORTRAN DOCUMENTATION
    ++++++++++++++++++++++++++++++++++++++++++++++++
    run
    -------------------
    details
    -------------------
     Dealocate the model variable
     so that it will be safe to run another simulation
     WARNING: if not done, allocation fault will happen
    """
    _wrapping.f90wrap_clean_model(mdl=self._handle)

def output_gnu(self, dof, filename):
    """
    output_gnu(self, dof, filename)
    
    
    Defined at call_run_model.f90 lines 216-221
    
    Parameters
    ----------
    mdl : Model
    dof : Unk
    filename : str
    
    """
    _wrapping.f90wrap_output_gnu(mdl=self._handle, dof=dof._handle, \
        filename=filename)

def output_vtk(self, dof, filename):
    """
    output_vtk(self, dof, filename)
    
    
    Defined at call_run_model.f90 lines 224-229
    
    Parameters
    ----------
    mdl : Model
    dof : Unk
    filename : str
    
    """
    _wrapping.f90wrap_output_vtk(mdl=self._handle, dof=dof._handle, \
        filename=filename)

def output_tec(self, dof, filename):
    """
    output_tec(self, dof, filename)
    
    
    Defined at call_run_model.f90 lines 232-237
    
    Parameters
    ----------
    mdl : Model
    dof : Unk
    filename : str
    
    """
    _wrapping.f90wrap_output_tec(mdl=self._handle, dof=dof._handle, \
        filename=filename)

def __init__(self, mesh, handle=None):
    """
    self = Friction_Data(mesh)
    
    
    Defined at call_run_model.f90 lines 239-295
    
    Parameters
    ----------
    mesh : Msh
    
    Returns
    -------
    my_friction : Friction_Data
    
    ++++++++++++++++++++++++++++++++++++++++++++++++
    FORTRAN DOCUMENTATION
    ++++++++++++++++++++++++++++++++++++++++++++++++
     define my_friction fortran type(linked to mdl%my_friction)
     allocate my_friction%xxx and set values from text file
    """
    f90wrap.runtime.FortranDerivedType.__init__(self)
    result = _wrapping.f90wrap_friction_initialise(mesh=mesh._handle)
    self._handle = result[0] if isinstance(result, tuple) else result

def __del__(self):
    """
    Destructor for class Friction_Data
    
    
    Defined at call_run_model.f90 lines 297-305
    
    Parameters
    ----------
    my_friction : Friction_Data
    
    ++++++++++++++++++++++++++++++++++++++++++++++++
    FORTRAN DOCUMENTATION
    ++++++++++++++++++++++++++++++++++++++++++++++++
     deallocate my_friction%xxx derived type
    """
    if self._alloc:
        _wrapping.f90wrap_friction_finalise(my_friction=self._handle)

@classmethod
def infiltration_initialise(cls, mesh):
    """
    self = Infiltration_Data(mesh)
    
    
    Defined at call_run_model.f90 lines 307-410
    
    Parameters
    ----------
    mesh : Msh
    
    Returns
    -------
    my_infiltration : Infiltration_Data
    
    ===================================================================================================================
     set the correspondance between cell and land value
    ===================================================================================================================
    ++++++++++++++++++++++++++++++++++++++++++++++++
    FORTRAN DOCUMENTATION
    ++++++++++++++++++++++++++++++++++++++++++++++++
    _COMMENT initialise friction(allocate and set values from file)
    """
    bare_class = cls.__new__(cls)
    f90wrap.runtime.FortranDerivedType.__init__(bare_class)
    result = _wrapping.f90wrap_infiltration_initialise(mesh=mesh._handle)
    bare_class._handle = result[0] if isinstance(result, tuple) else result
    return bare_class

def __del__(self):
    """
    Destructor for class Infiltration_Data
    
    
    Defined at call_run_model.f90 lines 412-421
    
    Parameters
    ----------
    my_infiltration : Infiltration_Data
    
    """
    if self._alloc:
        _wrapping.f90wrap_infiltration_finalise(my_infiltration=self._handle)

@classmethod
def phys_desc_initialise(cls, mesh):
    """
    self = Input_Data(mesh)
    
    
    Defined at call_run_model.f90 lines 423-429
    
    Parameters
    ----------
    mesh : Msh
    
    Returns
    -------
    my_phys_desc : Input_Data
    
    """
    bare_class = cls.__new__(cls)
    f90wrap.runtime.FortranDerivedType.__init__(bare_class)
    result = _wrapping.f90wrap_phys_desc_initialise(mesh=mesh._handle)
    bare_class._handle = result[0] if isinstance(result, tuple) else result
    return bare_class

def __del__(self):
    """
    Destructor for class Input_Data
    
    
    Defined at call_run_model.f90 lines 431-436
    
    Parameters
    ----------
    my_phys_desc : Input_Data
    
    """
    if self._alloc:
        _wrapping.f90wrap_phys_desc_finalise(my_phys_desc=self._handle)

def func(self, ctrl_in, grad_func):
    """
    cost_func = func(self, ctrl_in, grad_func)
    
    
    Defined at call_run_model.f90 lines 438-466
    
    Parameters
    ----------
    mdl : Model
    ctrl_in : float array
    grad_func : float array
    
    Returns
    -------
    cost_func : float
    
    ++++++++++++++++++++++++++++++++++++++++++++++++
    FORTRAN DOCUMENTATION
    ++++++++++++++++++++++++++++++++++++++++++++++++
     func
    compute the cost function and the gradient of cost function
    -------------------
    details
    -------------------
     set control = ctrl_in, and compute the cost function and the gradient of cost \
         function
     dim_all is defined in m_adjoint, is calculated once the control vector is \
         allocated(write_control)
     \param[inout] mdl : Model, the model instance
     \param[in] ctrl_in : real(rp), , the control vector to define as new control \
         vector
     \param[out] cost_func : real(rp),dimension(1), the curent value of cost \
         funtion(cost in fortran, cost_func in python)
     \param[out] grad_func : real(rp),dimension(control),the curent value of gradient \
         of cost funtion(control_back in fortran, grad_func in python)
    """
    cost_func = _wrapping.f90wrap_func(mdl=self._handle, ctrl_in=ctrl_in, \
        grad_func=grad_func)
    return cost_func

def write_land_uses():
    """
    write_land_uses()
    
    
    Defined at call_run_model.f90 lines 609-625
    
    
    =================================================================================================
    "
    """
    _wrapping.f90wrap_write_land_uses()

def write_bc():
    """
    write_bc()
    
    
    Defined at call_run_model.f90 lines 632-707
    
    
    ===================================================================================================================
     Interface Variables
    ===================================================================================================================
     all information needed is in object bc.
     the type bcs is defined in module m_model
     bcs is initialized in either
     - read_dass_mesh() --> call read_bc_file()
     - Initial()
     most subroutines modifying and using bc information are located in in \
         src/sw_mono/boundary.f90 file
    ===================================================================================================================
     Local variables
    ===================================================================================================================
    """
    _wrapping.f90wrap_write_bc()

def write_hydrograph():
    """
    write_hydrograph()
    
    
    Defined at call_run_model.f90 lines 709-752
    
    
    ===================================================================================================================
     Interface Variables
    ===================================================================================================================
     all information needed is in object bc.
     the type bcs is defined in module m_model
     bcs is initialized in either
     - read_dass_mesh() --> call read_bc_file()
     - Initial()
     most subroutines modifying and using bc information are located in in \
         src/sw_mono/boundary.f90 file
    ===================================================================================================================
     Local variables
    ===================================================================================================================
    """
    _wrapping.f90wrap_write_hydrograph()

def boundaries_copy(self):
    """
    o = boundaries_copy(self)
    
    
    Defined at call_run_model.f90 lines 754-758
    
    Parameters
    ----------
    i : Bcs
    
    Returns
    -------
    o : Bcs
    
    """
    o = _wrapping.f90wrap_boundaries_copy(i=self._handle)
    o = f90wrap.runtime.lookup_class("wrapping.bcs").from_handle(o, alloc=True)
    return o

def dof_copy(self):
    """
    o = dof_copy(self)
    
    
    Defined at call_run_model.f90 lines 763-767
    
    Parameters
    ----------
    i : Unk
    
    Returns
    -------
    o : Unk
    
    """
    o = _wrapping.f90wrap_dof_copy(i=self._handle)
    o = f90wrap.runtime.lookup_class("wrapping.unk").from_handle(o, alloc=True)
    return o

def mesh_copy(self):
    """
    o = mesh_copy(self)
    
    
    Defined at call_run_model.f90 lines 772-776
    
    Parameters
    ----------
    i : Msh
    
    Returns
    -------
    o : Msh
    
    """
    o = _wrapping.f90wrap_mesh_copy(i=self._handle)
    o = f90wrap.runtime.lookup_class("wrapping.msh").from_handle(o, alloc=True)
    return o

def reallocate_manning(new_size):
    """
    reallocate_manning(new_size)
    
    
    Defined at call_run_model.f90 lines 782-794
    
    Parameters
    ----------
    new_size : int
    
    _COMMENT, manning)
    _COMMENTsize_in, size_out, manning_out )
    """
    _wrapping.f90wrap_reallocate_manning(new_size=new_size)


_array_initialisers = []
_dt_array_initialisers = []

try:
    for func in _array_initialisers:
        func()
except ValueError:
    logging.debug('unallocated array(s) detected on import of module "call_model".')

for func in _dt_array_initialisers:
    func()
