! Module call_model defined in file call_run_model.f90

subroutine f90wrap_model__get__mesh(this, f90wrap_mesh)
    use call_model, only: model
    use m_mesh, only: msh
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_mesh(2)
    type(msh_ptr_type) :: mesh_ptr
    
    this_ptr = transfer(this, this_ptr)
    mesh_ptr%p => this_ptr%p%mesh
    f90wrap_mesh = transfer(mesh_ptr,f90wrap_mesh)
end subroutine f90wrap_model__get__mesh

subroutine f90wrap_model__set__mesh(this, f90wrap_mesh)
    use call_model, only: model
    use m_mesh, only: msh
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_mesh(2)
    type(msh_ptr_type) :: mesh_ptr
    
    this_ptr = transfer(this, this_ptr)
    mesh_ptr = transfer(f90wrap_mesh,mesh_ptr)
    this_ptr%p%mesh = mesh_ptr%p
end subroutine f90wrap_model__set__mesh

subroutine f90wrap_model__get__dof0(this, f90wrap_dof0)
    use call_model, only: model
    use m_model, only: unk
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_dof0(2)
    type(unk_ptr_type) :: dof0_ptr
    
    this_ptr = transfer(this, this_ptr)
    dof0_ptr%p => this_ptr%p%dof0
    f90wrap_dof0 = transfer(dof0_ptr,f90wrap_dof0)
end subroutine f90wrap_model__get__dof0

subroutine f90wrap_model__set__dof0(this, f90wrap_dof0)
    use call_model, only: model
    use m_model, only: unk
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_dof0(2)
    type(unk_ptr_type) :: dof0_ptr
    
    this_ptr = transfer(this, this_ptr)
    dof0_ptr = transfer(f90wrap_dof0,dof0_ptr)
    this_ptr%p%dof0 = dof0_ptr%p
end subroutine f90wrap_model__set__dof0

subroutine f90wrap_model__get__dof(this, f90wrap_dof)
    use call_model, only: model
    use m_model, only: unk
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_dof(2)
    type(unk_ptr_type) :: dof_ptr
    
    this_ptr = transfer(this, this_ptr)
    dof_ptr%p => this_ptr%p%dof
    f90wrap_dof = transfer(dof_ptr,f90wrap_dof)
end subroutine f90wrap_model__get__dof

subroutine f90wrap_model__set__dof(this, f90wrap_dof)
    use call_model, only: model
    use m_model, only: unk
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_dof(2)
    type(unk_ptr_type) :: dof_ptr
    
    this_ptr = transfer(this, this_ptr)
    dof_ptr = transfer(f90wrap_dof,dof_ptr)
    this_ptr%p%dof = dof_ptr%p
end subroutine f90wrap_model__set__dof

subroutine f90wrap_model__get__cost(this, f90wrap_cost)
    use call_model, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_cost
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_cost = this_ptr%p%cost
end subroutine f90wrap_model__get__cost

subroutine f90wrap_model__set__cost(this, f90wrap_cost)
    use call_model, only: model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_cost
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%cost = f90wrap_cost
end subroutine f90wrap_model__set__cost

subroutine f90wrap_model__get__param(this, f90wrap_param)
    use call_model, only: model
    use m_model, only: input_param
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_param(2)
    type(input_param_ptr_type) :: param_ptr
    
    this_ptr = transfer(this, this_ptr)
    param_ptr%p => this_ptr%p%param
    f90wrap_param = transfer(param_ptr,f90wrap_param)
end subroutine f90wrap_model__get__param

subroutine f90wrap_model__set__param(this, f90wrap_param)
    use call_model, only: model
    use m_model, only: input_param
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_param(2)
    type(input_param_ptr_type) :: param_ptr
    
    this_ptr = transfer(this, this_ptr)
    param_ptr = transfer(f90wrap_param,param_ptr)
    this_ptr%p%param = param_ptr%p
end subroutine f90wrap_model__set__param

subroutine f90wrap_model__get__my_friction(this, f90wrap_my_friction)
    use call_model, only: model
    use m_model, only: friction_data
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type friction_data_ptr_type
        type(friction_data), pointer :: p => NULL()
    end type friction_data_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_my_friction(2)
    type(friction_data_ptr_type) :: my_friction_ptr
    
    this_ptr = transfer(this, this_ptr)
    my_friction_ptr%p => this_ptr%p%my_friction
    f90wrap_my_friction = transfer(my_friction_ptr,f90wrap_my_friction)
end subroutine f90wrap_model__get__my_friction

subroutine f90wrap_model__set__my_friction(this, f90wrap_my_friction)
    use call_model, only: model
    use m_model, only: friction_data
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type friction_data_ptr_type
        type(friction_data), pointer :: p => NULL()
    end type friction_data_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_my_friction(2)
    type(friction_data_ptr_type) :: my_friction_ptr
    
    this_ptr = transfer(this, this_ptr)
    my_friction_ptr = transfer(f90wrap_my_friction,my_friction_ptr)
    this_ptr%p%my_friction = my_friction_ptr%p
end subroutine f90wrap_model__set__my_friction

subroutine f90wrap_model__get__my_infiltration(this, f90wrap_my_infiltration)
    use call_model, only: model
    use m_model, only: infiltration_data
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_my_infiltration(2)
    type(infiltration_data_ptr_type) :: my_infiltration_ptr
    
    this_ptr = transfer(this, this_ptr)
    my_infiltration_ptr%p => this_ptr%p%my_infiltration
    f90wrap_my_infiltration = transfer(my_infiltration_ptr,f90wrap_my_infiltration)
end subroutine f90wrap_model__get__my_infiltration

subroutine f90wrap_model__set__my_infiltration(this, f90wrap_my_infiltration)
    use call_model, only: model
    use m_model, only: infiltration_data
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_my_infiltration(2)
    type(infiltration_data_ptr_type) :: my_infiltration_ptr
    
    this_ptr = transfer(this, this_ptr)
    my_infiltration_ptr = transfer(f90wrap_my_infiltration,my_infiltration_ptr)
    this_ptr%p%my_infiltration = my_infiltration_ptr%p
end subroutine f90wrap_model__set__my_infiltration

subroutine f90wrap_model__get__my_phys_desc(this, f90wrap_my_phys_desc)
    use call_model, only: model
    use m_model, only: input_data
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_my_phys_desc(2)
    type(input_data_ptr_type) :: my_phys_desc_ptr
    
    this_ptr = transfer(this, this_ptr)
    my_phys_desc_ptr%p => this_ptr%p%my_phys_desc
    f90wrap_my_phys_desc = transfer(my_phys_desc_ptr,f90wrap_my_phys_desc)
end subroutine f90wrap_model__get__my_phys_desc

subroutine f90wrap_model__set__my_phys_desc(this, f90wrap_my_phys_desc)
    use call_model, only: model
    use m_model, only: input_data
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_my_phys_desc(2)
    type(input_data_ptr_type) :: my_phys_desc_ptr
    
    this_ptr = transfer(this, this_ptr)
    my_phys_desc_ptr = transfer(f90wrap_my_phys_desc,my_phys_desc_ptr)
    this_ptr%p%my_phys_desc = my_phys_desc_ptr%p
end subroutine f90wrap_model__set__my_phys_desc

subroutine f90wrap_model__get__my_param_model(this, f90wrap_my_param_model)
    use call_model, only: model
    use m_model, only: param_model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type param_model_ptr_type
        type(param_model), pointer :: p => NULL()
    end type param_model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_my_param_model(2)
    type(param_model_ptr_type) :: my_param_model_ptr
    
    this_ptr = transfer(this, this_ptr)
    my_param_model_ptr%p => this_ptr%p%my_param_model
    f90wrap_my_param_model = transfer(my_param_model_ptr,f90wrap_my_param_model)
end subroutine f90wrap_model__get__my_param_model

subroutine f90wrap_model__set__my_param_model(this, f90wrap_my_param_model)
    use call_model, only: model
    use m_model, only: param_model
    implicit none
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type param_model_ptr_type
        type(param_model), pointer :: p => NULL()
    end type param_model_ptr_type
    integer, intent(in)   :: this(2)
    type(model_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_my_param_model(2)
    type(param_model_ptr_type) :: my_param_model_ptr
    
    this_ptr = transfer(this, this_ptr)
    my_param_model_ptr = transfer(f90wrap_my_param_model,my_param_model_ptr)
    this_ptr%p%my_param_model = my_param_model_ptr%p
end subroutine f90wrap_model__set__my_param_model

subroutine f90wrap_model_initialise(mdl)
    use call_model, only: model, model_initialise
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(out), dimension(2) :: mdl
    allocate(mdl_ptr%p)
    call model_initialise(mdl=mdl_ptr%p)
    mdl = transfer(mdl_ptr, mdl)
end subroutine f90wrap_model_initialise

subroutine f90wrap_model_finalise(mdl)
    use call_model, only: model, model_finalise
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    mdl_ptr = transfer(mdl, mdl_ptr)
    call model_finalise(mdl=mdl_ptr%p)
    deallocate(mdl_ptr%p)
end subroutine f90wrap_model_finalise

subroutine f90wrap_init_solver(mdl)
    use call_model, only: model, init_solver
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    mdl_ptr = transfer(mdl, mdl_ptr)
    call init_solver(mdl=mdl_ptr%p)
end subroutine f90wrap_init_solver

subroutine f90wrap_init_friction(mdl)
    use call_model, only: init_friction, model
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    mdl_ptr = transfer(mdl, mdl_ptr)
    call init_friction(mdl=mdl_ptr%p)
end subroutine f90wrap_init_friction

subroutine f90wrap_init_fortran(mdl)
    use call_model, only: model, init_fortran
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    mdl_ptr = transfer(mdl, mdl_ptr)
    call init_fortran(mdl=mdl_ptr%p)
end subroutine f90wrap_init_fortran

subroutine f90wrap_init_all(mdl)
    use call_model, only: model, init_all
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    mdl_ptr = transfer(mdl, mdl_ptr)
    call init_all(mdl=mdl_ptr%p)
end subroutine f90wrap_init_all

subroutine f90wrap_init_back(mdl)
    use call_model, only: model, init_back
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    mdl_ptr = transfer(mdl, mdl_ptr)
    call init_back(mdl=mdl_ptr%p)
end subroutine f90wrap_init_back

subroutine f90wrap_run(mdl, arg)
    use call_model, only: model, run
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    character*(*), optional, intent(in) :: arg
    mdl_ptr = transfer(mdl, mdl_ptr)
    call run(mdl=mdl_ptr%p, arg=arg)
end subroutine f90wrap_run

subroutine f90wrap_clean_model(mdl)
    use call_model, only: model, clean_model
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    mdl_ptr = transfer(mdl, mdl_ptr)
    call clean_model(mdl=mdl_ptr%p)
end subroutine f90wrap_clean_model

subroutine f90wrap_output_gnu(mdl, dof, filename)
    use call_model, only: model, output_gnu
    use m_model, only: unk
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    type(unk_ptr_type) :: dof_ptr
    integer, intent(in), dimension(2) :: dof
    character(128), intent(in) :: filename
    mdl_ptr = transfer(mdl, mdl_ptr)
    dof_ptr = transfer(dof, dof_ptr)
    call output_gnu(mdl=mdl_ptr%p, dof=dof_ptr%p, filename=filename)
end subroutine f90wrap_output_gnu

subroutine f90wrap_output_vtk(mdl, dof, filename)
    use call_model, only: model, output_vtk
    use m_model, only: unk
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    type(unk_ptr_type) :: dof_ptr
    integer, intent(in), dimension(2) :: dof
    character(128), intent(in) :: filename
    mdl_ptr = transfer(mdl, mdl_ptr)
    dof_ptr = transfer(dof, dof_ptr)
    call output_vtk(mdl=mdl_ptr%p, dof=dof_ptr%p, filename=filename)
end subroutine f90wrap_output_vtk

subroutine f90wrap_output_tec(mdl, dof, filename)
    use call_model, only: model, output_tec
    use m_model, only: unk
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    type(unk_ptr_type) :: dof_ptr
    integer, intent(in), dimension(2) :: dof
    character(128), intent(in) :: filename
    mdl_ptr = transfer(mdl, mdl_ptr)
    dof_ptr = transfer(dof, dof_ptr)
    call output_tec(mdl=mdl_ptr%p, dof=dof_ptr%p, filename=filename)
end subroutine f90wrap_output_tec

subroutine f90wrap_friction_initialise(my_friction, mesh)
    use call_model, only: friction_initialise
    use m_model, only: friction_data
    use m_mesh, only: msh
    implicit none
    
    type friction_data_ptr_type
        type(friction_data), pointer :: p => NULL()
    end type friction_data_ptr_type
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type(friction_data_ptr_type) :: my_friction_ptr
    integer, intent(out), dimension(2) :: my_friction
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    mesh_ptr = transfer(mesh, mesh_ptr)
    allocate(my_friction_ptr%p)
    call friction_initialise(my_friction=my_friction_ptr%p, mesh=mesh_ptr%p)
    my_friction = transfer(my_friction_ptr, my_friction)
end subroutine f90wrap_friction_initialise

subroutine f90wrap_friction_finalise(my_friction)
    use m_model, only: friction_data
    use call_model, only: friction_finalise
    implicit none
    
    type friction_data_ptr_type
        type(friction_data), pointer :: p => NULL()
    end type friction_data_ptr_type
    type(friction_data_ptr_type) :: my_friction_ptr
    integer, intent(in), dimension(2) :: my_friction
    my_friction_ptr = transfer(my_friction, my_friction_ptr)
    call friction_finalise(my_friction=my_friction_ptr%p)
    deallocate(my_friction_ptr%p)
end subroutine f90wrap_friction_finalise

subroutine f90wrap_infiltration_initialise(my_infiltration, mesh)
    use m_mesh, only: msh
    use call_model, only: infiltration_initialise
    use m_model, only: infiltration_data
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    type(infiltration_data_ptr_type) :: my_infiltration_ptr
    integer, intent(out), dimension(2) :: my_infiltration
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    mesh_ptr = transfer(mesh, mesh_ptr)
    allocate(my_infiltration_ptr%p)
    call infiltration_initialise(my_infiltration=my_infiltration_ptr%p, mesh=mesh_ptr%p)
    my_infiltration = transfer(my_infiltration_ptr, my_infiltration)
end subroutine f90wrap_infiltration_initialise

subroutine f90wrap_infiltration_finalise(my_infiltration)
    use call_model, only: infiltration_finalise
    use m_model, only: infiltration_data
    implicit none
    
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    type(infiltration_data_ptr_type) :: my_infiltration_ptr
    integer, intent(in), dimension(2) :: my_infiltration
    my_infiltration_ptr = transfer(my_infiltration, my_infiltration_ptr)
    call infiltration_finalise(my_infiltration=my_infiltration_ptr%p)
    deallocate(my_infiltration_ptr%p)
end subroutine f90wrap_infiltration_finalise

subroutine f90wrap_phys_desc_initialise(my_phys_desc, mesh)
    use m_mesh, only: msh
    use m_model, only: input_data
    use call_model, only: phys_desc_initialise
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type(input_data_ptr_type) :: my_phys_desc_ptr
    integer, intent(out), dimension(2) :: my_phys_desc
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    mesh_ptr = transfer(mesh, mesh_ptr)
    allocate(my_phys_desc_ptr%p)
    call phys_desc_initialise(my_phys_desc=my_phys_desc_ptr%p, mesh=mesh_ptr%p)
    my_phys_desc = transfer(my_phys_desc_ptr, my_phys_desc)
end subroutine f90wrap_phys_desc_initialise

subroutine f90wrap_phys_desc_finalise(my_phys_desc)
    use m_model, only: input_data
    use call_model, only: phys_desc_finalise
    implicit none
    
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type(input_data_ptr_type) :: my_phys_desc_ptr
    integer, intent(in), dimension(2) :: my_phys_desc
    my_phys_desc_ptr = transfer(my_phys_desc, my_phys_desc_ptr)
    call phys_desc_finalise(my_phys_desc=my_phys_desc_ptr%p)
    deallocate(my_phys_desc_ptr%p)
end subroutine f90wrap_phys_desc_finalise

subroutine f90wrap_func(mdl, ctrl_in, cost_func, grad_func, n0, n1)
    use call_model, only: model, func
    implicit none
    
    type model_ptr_type
        type(model), pointer :: p => NULL()
    end type model_ptr_type
    type(model_ptr_type) :: mdl_ptr
    integer, intent(in), dimension(2) :: mdl
    real(8), intent(in), dimension(n0) :: ctrl_in
    real(8), intent(out) :: cost_func
    real(8), intent(inout), dimension(n1) :: grad_func
    integer :: n0
    !f2py intent(hide), depend(ctrl_in) :: n0 = shape(ctrl_in,0)
    integer :: n1
    !f2py intent(hide), depend(grad_func) :: n1 = shape(grad_func,0)
    mdl_ptr = transfer(mdl, mdl_ptr)
    call func(mdl=mdl_ptr%p, ctrl_in=ctrl_in, cost_func=cost_func, grad_func=grad_func)
end subroutine f90wrap_func

subroutine f90wrap_write_land_uses
    use call_model, only: write_land_uses
    implicit none
    
    call write_land_uses()
end subroutine f90wrap_write_land_uses

subroutine f90wrap_write_bc
    use call_model, only: write_bc
    implicit none
    
    call write_bc()
end subroutine f90wrap_write_bc

subroutine f90wrap_write_hydrograph
    use call_model, only: write_hydrograph
    implicit none
    
    call write_hydrograph()
end subroutine f90wrap_write_hydrograph

subroutine f90wrap_boundaries_copy(i, o)
    use call_model, only: boundaries_copy
    use m_model, only: bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type(bcs_ptr_type) :: i_ptr
    integer, intent(in), dimension(2) :: i
    type(bcs_ptr_type) :: o_ptr
    integer, intent(out), dimension(2) :: o
    i_ptr = transfer(i, i_ptr)
    allocate(o_ptr%p)
    call boundaries_copy(i=i_ptr%p, o=o_ptr%p)
    o = transfer(o_ptr, o)
end subroutine f90wrap_boundaries_copy

subroutine f90wrap_dof_copy(i, o)
    use call_model, only: dof_copy
    use m_model, only: unk
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(unk_ptr_type) :: i_ptr
    integer, intent(in), dimension(2) :: i
    type(unk_ptr_type) :: o_ptr
    integer, intent(out), dimension(2) :: o
    i_ptr = transfer(i, i_ptr)
    allocate(o_ptr%p)
    call dof_copy(i=i_ptr%p, o=o_ptr%p)
    o = transfer(o_ptr, o)
end subroutine f90wrap_dof_copy

subroutine f90wrap_mesh_copy(i, o)
    use m_mesh, only: msh
    use call_model, only: mesh_copy
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type(msh_ptr_type) :: i_ptr
    integer, intent(in), dimension(2) :: i
    type(msh_ptr_type) :: o_ptr
    integer, intent(out), dimension(2) :: o
    i_ptr = transfer(i, i_ptr)
    allocate(o_ptr%p)
    call mesh_copy(i=i_ptr%p, o=o_ptr%p)
    o = transfer(o_ptr, o)
end subroutine f90wrap_mesh_copy

subroutine f90wrap_reallocate_manning(new_size)
    use call_model, only: reallocate_manning
    implicit none
    
    integer, intent(in) :: new_size
    call reallocate_manning(new_size=new_size)
end subroutine f90wrap_reallocate_manning

! End of module call_model defined in file call_run_model.f90

