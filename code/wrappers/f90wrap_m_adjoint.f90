! Module m_adjoint defined in file m_adjoint.f90

subroutine f90wrap_model_direct(mesh, dof0, dof)
    use m_mesh, only: msh
    use m_adjoint, only: model_direct
    use m_model, only: unk
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    type(unk_ptr_type) :: dof0_ptr
    integer, intent(in), dimension(2) :: dof0
    type(unk_ptr_type) :: dof_ptr
    integer, intent(in), dimension(2) :: dof
    mesh_ptr = transfer(mesh, mesh_ptr)
    dof0_ptr = transfer(dof0, dof0_ptr)
    dof_ptr = transfer(dof, dof_ptr)
    call model_direct(mesh=mesh_ptr%p, dof0=dof0_ptr%p, dof=dof_ptr%p)
end subroutine f90wrap_model_direct

subroutine f90wrap_model_direct_perturb(mesh, dof0, dof)
    use m_mesh, only: msh
    use m_adjoint, only: model_direct_perturb
    use m_model, only: unk
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    type(unk_ptr_type) :: dof0_ptr
    integer, intent(in), dimension(2) :: dof0
    type(unk_ptr_type) :: dof_ptr
    integer, intent(in), dimension(2) :: dof
    mesh_ptr = transfer(mesh, mesh_ptr)
    dof0_ptr = transfer(dof0, dof0_ptr)
    dof_ptr = transfer(dof, dof_ptr)
    call model_direct_perturb(mesh=mesh_ptr%p, dof0=dof0_ptr%p, dof=dof_ptr%p)
end subroutine f90wrap_model_direct_perturb

subroutine f90wrap_adjoint_model(mesh, dof0, dof)
    use m_mesh, only: msh
    use m_adjoint, only: adjoint_model
    use m_model, only: unk
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    type(unk_ptr_type) :: dof0_ptr
    integer, intent(in), dimension(2) :: dof0
    type(unk_ptr_type) :: dof_ptr
    integer, intent(in), dimension(2) :: dof
    mesh_ptr = transfer(mesh, mesh_ptr)
    dof0_ptr = transfer(dof0, dof0_ptr)
    dof_ptr = transfer(dof, dof_ptr)
    call adjoint_model(mesh=mesh_ptr%p, dof0=dof0_ptr%p, dof=dof_ptr%p)
end subroutine f90wrap_adjoint_model

subroutine f90wrap_write_control(dof0, mesh)
    use m_mesh, only: msh
    use m_adjoint, only: write_control
    use m_model, only: unk
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(unk_ptr_type) :: dof0_ptr
    integer, intent(in), dimension(2) :: dof0
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    dof0_ptr = transfer(dof0, dof0_ptr)
    mesh_ptr = transfer(mesh, mesh_ptr)
    call write_control(dof0=dof0_ptr%p, mesh=mesh_ptr%p)
end subroutine f90wrap_write_control

subroutine f90wrap_write_control_back(dof0, mesh)
    use m_adjoint, only: write_control_back
    use m_model, only: unk
    use m_mesh, only: msh
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(unk_ptr_type) :: dof0_ptr
    integer, intent(in), dimension(2) :: dof0
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    dof0_ptr = transfer(dof0, dof0_ptr)
    mesh_ptr = transfer(mesh, mesh_ptr)
    call write_control_back(dof0=dof0_ptr%p, mesh=mesh_ptr%p)
end subroutine f90wrap_write_control_back

subroutine f90wrap_read_control(dof0, mesh)
    use m_mesh, only: msh
    use m_adjoint, only: read_control
    use m_model, only: unk
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(unk_ptr_type) :: dof0_ptr
    integer, intent(in), dimension(2) :: dof0
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    dof0_ptr = transfer(dof0, dof0_ptr)
    mesh_ptr = transfer(mesh, mesh_ptr)
    call read_control(dof0=dof0_ptr%p, mesh=mesh_ptr%p)
end subroutine f90wrap_read_control

subroutine f90wrap_output_control(dof0, mesh)
    use m_mesh, only: msh
    use m_model, only: unk
    use m_adjoint, only: output_control
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(unk_ptr_type) :: dof0_ptr
    integer, intent(in), dimension(2) :: dof0
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    dof0_ptr = transfer(dof0, dof0_ptr)
    mesh_ptr = transfer(mesh, mesh_ptr)
    call output_control(dof0=dof0_ptr%p, mesh=mesh_ptr%p)
end subroutine f90wrap_output_control

subroutine f90wrap_output_control_back(mesh)
    use m_mesh, only: msh
    use m_adjoint, only: output_control_back
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    mesh_ptr = transfer(mesh, mesh_ptr)
    call output_control_back(mesh=mesh_ptr%p)
end subroutine f90wrap_output_control_back

subroutine f90wrap_m_adjoint__get__dof0_diff(f90wrap_dof0_diff)
    use m_model, only: unk
    use m_adjoint, only: m_adjoint_dof0_diff => dof0_diff
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(out) :: f90wrap_dof0_diff(2)
    type(unk_ptr_type) :: dof0_diff_ptr
    
    dof0_diff_ptr%p => m_adjoint_dof0_diff
    f90wrap_dof0_diff = transfer(dof0_diff_ptr,f90wrap_dof0_diff)
end subroutine f90wrap_m_adjoint__get__dof0_diff

subroutine f90wrap_m_adjoint__set__dof0_diff(f90wrap_dof0_diff)
    use m_model, only: unk
    use m_adjoint, only: m_adjoint_dof0_diff => dof0_diff
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(in) :: f90wrap_dof0_diff(2)
    type(unk_ptr_type) :: dof0_diff_ptr
    
    dof0_diff_ptr = transfer(f90wrap_dof0_diff,dof0_diff_ptr)
    m_adjoint_dof0_diff = dof0_diff_ptr%p
end subroutine f90wrap_m_adjoint__set__dof0_diff

subroutine f90wrap_m_adjoint__get__dof0_back(f90wrap_dof0_back)
    use m_model, only: unk
    use m_adjoint, only: m_adjoint_dof0_back => dof0_back
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(out) :: f90wrap_dof0_back(2)
    type(unk_ptr_type) :: dof0_back_ptr
    
    dof0_back_ptr%p => m_adjoint_dof0_back
    f90wrap_dof0_back = transfer(dof0_back_ptr,f90wrap_dof0_back)
end subroutine f90wrap_m_adjoint__get__dof0_back

subroutine f90wrap_m_adjoint__set__dof0_back(f90wrap_dof0_back)
    use m_model, only: unk
    use m_adjoint, only: m_adjoint_dof0_back => dof0_back
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(in) :: f90wrap_dof0_back(2)
    type(unk_ptr_type) :: dof0_back_ptr
    
    dof0_back_ptr = transfer(f90wrap_dof0_back,dof0_back_ptr)
    m_adjoint_dof0_back = dof0_back_ptr%p
end subroutine f90wrap_m_adjoint__set__dof0_back

subroutine f90wrap_m_adjoint__get__dof_diff(f90wrap_dof_diff)
    use m_model, only: unk
    use m_adjoint, only: m_adjoint_dof_diff => dof_diff
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(out) :: f90wrap_dof_diff(2)
    type(unk_ptr_type) :: dof_diff_ptr
    
    dof_diff_ptr%p => m_adjoint_dof_diff
    f90wrap_dof_diff = transfer(dof_diff_ptr,f90wrap_dof_diff)
end subroutine f90wrap_m_adjoint__get__dof_diff

subroutine f90wrap_m_adjoint__set__dof_diff(f90wrap_dof_diff)
    use m_model, only: unk
    use m_adjoint, only: m_adjoint_dof_diff => dof_diff
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(in) :: f90wrap_dof_diff(2)
    type(unk_ptr_type) :: dof_diff_ptr
    
    dof_diff_ptr = transfer(f90wrap_dof_diff,dof_diff_ptr)
    m_adjoint_dof_diff = dof_diff_ptr%p
end subroutine f90wrap_m_adjoint__set__dof_diff

subroutine f90wrap_m_adjoint__get__dof_back(f90wrap_dof_back)
    use m_model, only: unk
    use m_adjoint, only: m_adjoint_dof_back => dof_back
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(out) :: f90wrap_dof_back(2)
    type(unk_ptr_type) :: dof_back_ptr
    
    dof_back_ptr%p => m_adjoint_dof_back
    f90wrap_dof_back = transfer(dof_back_ptr,f90wrap_dof_back)
end subroutine f90wrap_m_adjoint__get__dof_back

subroutine f90wrap_m_adjoint__set__dof_back(f90wrap_dof_back)
    use m_model, only: unk
    use m_adjoint, only: m_adjoint_dof_back => dof_back
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(in) :: f90wrap_dof_back(2)
    type(unk_ptr_type) :: dof_back_ptr
    
    dof_back_ptr = transfer(f90wrap_dof_back,dof_back_ptr)
    m_adjoint_dof_back = dof_back_ptr%p
end subroutine f90wrap_m_adjoint__set__dof_back

subroutine f90wrap_m_adjoint__array__control(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_random
    use m_time_screen
    use m_tap_vars
    use m_obs
    use m_model
    use m_adjoint, only: m_adjoint_control => control
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_adjoint_control)) then
        dshape(1:1) = shape(m_adjoint_control)
        dloc = loc(m_adjoint_control)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_adjoint__array__control

subroutine f90wrap_m_adjoint__array__control_back(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_random
    use m_time_screen
    use m_tap_vars
    use m_obs
    use m_model
    use m_adjoint, only: m_adjoint_control_back => control_back
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_adjoint_control_back)) then
        dshape(1:1) = shape(m_adjoint_control_back)
        dloc = loc(m_adjoint_control_back)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_adjoint__array__control_back

subroutine f90wrap_m_adjoint__array__control_diff(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_random
    use m_time_screen
    use m_tap_vars
    use m_obs
    use m_model
    use m_adjoint, only: m_adjoint_control_diff => control_diff
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_adjoint_control_diff)) then
        dshape(1:1) = shape(m_adjoint_control_diff)
        dloc = loc(m_adjoint_control_diff)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_adjoint__array__control_diff

subroutine f90wrap_m_adjoint__array__control_perturb(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_random
    use m_time_screen
    use m_tap_vars
    use m_obs
    use m_model
    use m_adjoint, only: m_adjoint_control_perturb => control_perturb
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_adjoint_control_perturb)) then
        dshape(1:1) = shape(m_adjoint_control_perturb)
        dloc = loc(m_adjoint_control_perturb)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_adjoint__array__control_perturb

subroutine f90wrap_m_adjoint__array__control_lbound(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_random
    use m_time_screen
    use m_tap_vars
    use m_obs
    use m_model
    use m_adjoint, only: m_adjoint_control_lbound => control_lbound
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_adjoint_control_lbound)) then
        dshape(1:1) = shape(m_adjoint_control_lbound)
        dloc = loc(m_adjoint_control_lbound)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_adjoint__array__control_lbound

subroutine f90wrap_m_adjoint__array__control_ubound(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_random
    use m_time_screen
    use m_tap_vars
    use m_obs
    use m_model
    use m_adjoint, only: m_adjoint_control_ubound => control_ubound
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_adjoint_control_ubound)) then
        dshape(1:1) = shape(m_adjoint_control_ubound)
        dloc = loc(m_adjoint_control_ubound)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_adjoint__array__control_ubound

subroutine f90wrap_m_adjoint__get__cost(f90wrap_cost)
    use m_adjoint, only: m_adjoint_cost => cost
    implicit none
    real(8), intent(out) :: f90wrap_cost
    
    f90wrap_cost = m_adjoint_cost
end subroutine f90wrap_m_adjoint__get__cost

subroutine f90wrap_m_adjoint__set__cost(f90wrap_cost)
    use m_adjoint, only: m_adjoint_cost => cost
    implicit none
    real(8), intent(in) :: f90wrap_cost
    
    m_adjoint_cost = f90wrap_cost
end subroutine f90wrap_m_adjoint__set__cost

subroutine f90wrap_m_adjoint__get__cost_back(f90wrap_cost_back)
    use m_adjoint, only: m_adjoint_cost_back => cost_back
    implicit none
    real(8), intent(out) :: f90wrap_cost_back
    
    f90wrap_cost_back = m_adjoint_cost_back
end subroutine f90wrap_m_adjoint__get__cost_back

subroutine f90wrap_m_adjoint__set__cost_back(f90wrap_cost_back)
    use m_adjoint, only: m_adjoint_cost_back => cost_back
    implicit none
    real(8), intent(in) :: f90wrap_cost_back
    
    m_adjoint_cost_back = f90wrap_cost_back
end subroutine f90wrap_m_adjoint__set__cost_back

subroutine f90wrap_m_adjoint__get__cost_diff(f90wrap_cost_diff)
    use m_adjoint, only: m_adjoint_cost_diff => cost_diff
    implicit none
    real(8), intent(out) :: f90wrap_cost_diff
    
    f90wrap_cost_diff = m_adjoint_cost_diff
end subroutine f90wrap_m_adjoint__get__cost_diff

subroutine f90wrap_m_adjoint__set__cost_diff(f90wrap_cost_diff)
    use m_adjoint, only: m_adjoint_cost_diff => cost_diff
    implicit none
    real(8), intent(in) :: f90wrap_cost_diff
    
    m_adjoint_cost_diff = f90wrap_cost_diff
end subroutine f90wrap_m_adjoint__set__cost_diff

subroutine f90wrap_m_adjoint__get__cost_perturb(f90wrap_cost_perturb)
    use m_adjoint, only: m_adjoint_cost_perturb => cost_perturb
    implicit none
    real(8), intent(out) :: f90wrap_cost_perturb
    
    f90wrap_cost_perturb = m_adjoint_cost_perturb
end subroutine f90wrap_m_adjoint__get__cost_perturb

subroutine f90wrap_m_adjoint__set__cost_perturb(f90wrap_cost_perturb)
    use m_adjoint, only: m_adjoint_cost_perturb => cost_perturb
    implicit none
    real(8), intent(in) :: f90wrap_cost_perturb
    
    m_adjoint_cost_perturb = f90wrap_cost_perturb
end subroutine f90wrap_m_adjoint__set__cost_perturb

subroutine f90wrap_m_adjoint__get__ic(f90wrap_ic)
    use m_adjoint, only: m_adjoint_ic => ic
    implicit none
    integer(4), intent(out) :: f90wrap_ic
    
    f90wrap_ic = m_adjoint_ic
end subroutine f90wrap_m_adjoint__get__ic

subroutine f90wrap_m_adjoint__set__ic(f90wrap_ic)
    use m_adjoint, only: m_adjoint_ic => ic
    implicit none
    integer(4), intent(in) :: f90wrap_ic
    
    m_adjoint_ic = f90wrap_ic
end subroutine f90wrap_m_adjoint__set__ic

subroutine f90wrap_m_adjoint__get__ite_min(f90wrap_ite_min)
    use m_adjoint, only: m_adjoint_ite_min => ite_min
    implicit none
    integer(4), intent(out) :: f90wrap_ite_min
    
    f90wrap_ite_min = m_adjoint_ite_min
end subroutine f90wrap_m_adjoint__get__ite_min

subroutine f90wrap_m_adjoint__set__ite_min(f90wrap_ite_min)
    use m_adjoint, only: m_adjoint_ite_min => ite_min
    implicit none
    integer(4), intent(in) :: f90wrap_ite_min
    
    m_adjoint_ite_min = f90wrap_ite_min
end subroutine f90wrap_m_adjoint__set__ite_min

subroutine f90wrap_m_adjoint__get__nb_vars_in_control(f90wrap_nb_vars_in_control)
    use m_adjoint, only: m_adjoint_nb_vars_in_control => nb_vars_in_control
    implicit none
    integer(4), intent(out) :: f90wrap_nb_vars_in_control
    
    f90wrap_nb_vars_in_control = m_adjoint_nb_vars_in_control
end subroutine f90wrap_m_adjoint__get__nb_vars_in_control

subroutine f90wrap_m_adjoint__set__nb_vars_in_control(f90wrap_nb_vars_in_control)
    use m_adjoint, only: m_adjoint_nb_vars_in_control => nb_vars_in_control
    implicit none
    integer(4), intent(in) :: f90wrap_nb_vars_in_control
    
    m_adjoint_nb_vars_in_control = f90wrap_nb_vars_in_control
end subroutine f90wrap_m_adjoint__set__nb_vars_in_control

subroutine f90wrap_m_adjoint__array__dim_vars_in_control(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_random
    use m_time_screen
    use m_tap_vars
    use m_obs
    use m_model
    use m_adjoint, only: m_adjoint_dim_vars_in_control => dim_vars_in_control
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    dshape(1:1) = shape(m_adjoint_dim_vars_in_control)
    dloc = loc(m_adjoint_dim_vars_in_control)
end subroutine f90wrap_m_adjoint__array__dim_vars_in_control

subroutine f90wrap_m_adjoint__get__dim_all(f90wrap_dim_all)
    use m_adjoint, only: m_adjoint_dim_all => dim_all
    implicit none
    integer(4), intent(out) :: f90wrap_dim_all
    
    f90wrap_dim_all = m_adjoint_dim_all
end subroutine f90wrap_m_adjoint__get__dim_all

subroutine f90wrap_m_adjoint__set__dim_all(f90wrap_dim_all)
    use m_adjoint, only: m_adjoint_dim_all => dim_all
    implicit none
    integer(4), intent(in) :: f90wrap_dim_all
    
    m_adjoint_dim_all = f90wrap_dim_all
end subroutine f90wrap_m_adjoint__set__dim_all

subroutine f90wrap_m_adjoint__array__bathy_cell_copy(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_random
    use m_time_screen
    use m_tap_vars
    use m_obs
    use m_model
    use m_adjoint, only: m_adjoint_bathy_cell_copy => bathy_cell_copy
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_adjoint_bathy_cell_copy)) then
        dshape(1:1) = shape(m_adjoint_bathy_cell_copy)
        dloc = loc(m_adjoint_bathy_cell_copy)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_adjoint__array__bathy_cell_copy

! End of module m_adjoint defined in file m_adjoint.f90

