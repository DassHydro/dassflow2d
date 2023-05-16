! Module m_tap_vars defined in file m_tap_vars.f90

subroutine f90wrap_alloc_back_vars(dof0_back, dof_back, mesh)
    use m_tap_vars, only: alloc_back_vars
    use m_model, only: unk
    use m_mesh, only: msh
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(unk_ptr_type) :: dof0_back_ptr
    integer, intent(in), dimension(2) :: dof0_back
    type(unk_ptr_type) :: dof_back_ptr
    integer, intent(in), dimension(2) :: dof_back
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    dof0_back_ptr = transfer(dof0_back, dof0_back_ptr)
    dof_back_ptr = transfer(dof_back, dof_back_ptr)
    mesh_ptr = transfer(mesh, mesh_ptr)
    call alloc_back_vars(dof0_back=dof0_back_ptr%p, dof_back=dof_back_ptr%p, mesh=mesh_ptr%p)
end subroutine f90wrap_alloc_back_vars

subroutine f90wrap_dealloc_back_vars
    use m_tap_vars, only: dealloc_back_vars
    implicit none
    
    call dealloc_back_vars()
end subroutine f90wrap_dealloc_back_vars

subroutine f90wrap_m_tap_vars__get__dt_diff(f90wrap_dt_diff)
    use m_tap_vars, only: m_tap_vars_dt_diff => dt_diff
    implicit none
    real(8), intent(out) :: f90wrap_dt_diff
    
    f90wrap_dt_diff = m_tap_vars_dt_diff
end subroutine f90wrap_m_tap_vars__get__dt_diff

subroutine f90wrap_m_tap_vars__set__dt_diff(f90wrap_dt_diff)
    use m_tap_vars, only: m_tap_vars_dt_diff => dt_diff
    implicit none
    real(8), intent(in) :: f90wrap_dt_diff
    
    m_tap_vars_dt_diff = f90wrap_dt_diff
end subroutine f90wrap_m_tap_vars__set__dt_diff

subroutine f90wrap_m_tap_vars__get__tc_diff(f90wrap_tc_diff)
    use m_tap_vars, only: m_tap_vars_tc_diff => tc_diff
    implicit none
    real(8), intent(out) :: f90wrap_tc_diff
    
    f90wrap_tc_diff = m_tap_vars_tc_diff
end subroutine f90wrap_m_tap_vars__get__tc_diff

subroutine f90wrap_m_tap_vars__set__tc_diff(f90wrap_tc_diff)
    use m_tap_vars, only: m_tap_vars_tc_diff => tc_diff
    implicit none
    real(8), intent(in) :: f90wrap_tc_diff
    
    m_tap_vars_tc_diff = f90wrap_tc_diff
end subroutine f90wrap_m_tap_vars__set__tc_diff

subroutine f90wrap_m_tap_vars__get__dt_back(f90wrap_dt_back)
    use m_tap_vars, only: m_tap_vars_dt_back => dt_back
    implicit none
    real(8), intent(out) :: f90wrap_dt_back
    
    f90wrap_dt_back = m_tap_vars_dt_back
end subroutine f90wrap_m_tap_vars__get__dt_back

subroutine f90wrap_m_tap_vars__set__dt_back(f90wrap_dt_back)
    use m_tap_vars, only: m_tap_vars_dt_back => dt_back
    implicit none
    real(8), intent(in) :: f90wrap_dt_back
    
    m_tap_vars_dt_back = f90wrap_dt_back
end subroutine f90wrap_m_tap_vars__set__dt_back

subroutine f90wrap_m_tap_vars__get__tc_back(f90wrap_tc_back)
    use m_tap_vars, only: m_tap_vars_tc_back => tc_back
    implicit none
    real(8), intent(out) :: f90wrap_tc_back
    
    f90wrap_tc_back = m_tap_vars_tc_back
end subroutine f90wrap_m_tap_vars__get__tc_back

subroutine f90wrap_m_tap_vars__set__tc_back(f90wrap_tc_back)
    use m_tap_vars, only: m_tap_vars_tc_back => tc_back
    implicit none
    real(8), intent(in) :: f90wrap_tc_back
    
    m_tap_vars_tc_back = f90wrap_tc_back
end subroutine f90wrap_m_tap_vars__set__tc_back

subroutine f90wrap_m_tap_vars__get__bc_diff(f90wrap_bc_diff)
    use m_model, only: bcs
    use m_tap_vars, only: m_tap_vars_bc_diff => bc_diff
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(out) :: f90wrap_bc_diff(2)
    type(bcs_ptr_type) :: bc_diff_ptr
    
    bc_diff_ptr%p => m_tap_vars_bc_diff
    f90wrap_bc_diff = transfer(bc_diff_ptr,f90wrap_bc_diff)
end subroutine f90wrap_m_tap_vars__get__bc_diff

subroutine f90wrap_m_tap_vars__set__bc_diff(f90wrap_bc_diff)
    use m_model, only: bcs
    use m_tap_vars, only: m_tap_vars_bc_diff => bc_diff
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(in) :: f90wrap_bc_diff(2)
    type(bcs_ptr_type) :: bc_diff_ptr
    
    bc_diff_ptr = transfer(f90wrap_bc_diff,bc_diff_ptr)
    m_tap_vars_bc_diff = bc_diff_ptr%p
end subroutine f90wrap_m_tap_vars__set__bc_diff

subroutine f90wrap_m_tap_vars__get__bc_back(f90wrap_bc_back)
    use m_model, only: bcs
    use m_tap_vars, only: m_tap_vars_bc_back => bc_back
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(out) :: f90wrap_bc_back(2)
    type(bcs_ptr_type) :: bc_back_ptr
    
    bc_back_ptr%p => m_tap_vars_bc_back
    f90wrap_bc_back = transfer(bc_back_ptr,f90wrap_bc_back)
end subroutine f90wrap_m_tap_vars__get__bc_back

subroutine f90wrap_m_tap_vars__set__bc_back(f90wrap_bc_back)
    use m_model, only: bcs
    use m_tap_vars, only: m_tap_vars_bc_back => bc_back
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(in) :: f90wrap_bc_back(2)
    type(bcs_ptr_type) :: bc_back_ptr
    
    bc_back_ptr = transfer(f90wrap_bc_back,bc_back_ptr)
    m_tap_vars_bc_back = bc_back_ptr%p
end subroutine f90wrap_m_tap_vars__set__bc_back

subroutine f90wrap_m_tap_vars__get__infil_diff(f90wrap_infil_diff)
    use m_model, only: infiltration_data
    use m_tap_vars, only: m_tap_vars_infil_diff => infil_diff
    implicit none
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer, intent(out) :: f90wrap_infil_diff(2)
    type(infiltration_data_ptr_type) :: infil_diff_ptr
    
    infil_diff_ptr%p => m_tap_vars_infil_diff
    f90wrap_infil_diff = transfer(infil_diff_ptr,f90wrap_infil_diff)
end subroutine f90wrap_m_tap_vars__get__infil_diff

subroutine f90wrap_m_tap_vars__set__infil_diff(f90wrap_infil_diff)
    use m_model, only: infiltration_data
    use m_tap_vars, only: m_tap_vars_infil_diff => infil_diff
    implicit none
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer, intent(in) :: f90wrap_infil_diff(2)
    type(infiltration_data_ptr_type) :: infil_diff_ptr
    
    infil_diff_ptr = transfer(f90wrap_infil_diff,infil_diff_ptr)
    m_tap_vars_infil_diff = infil_diff_ptr%p
end subroutine f90wrap_m_tap_vars__set__infil_diff

subroutine f90wrap_m_tap_vars__get__infil_back(f90wrap_infil_back)
    use m_model, only: infiltration_data
    use m_tap_vars, only: m_tap_vars_infil_back => infil_back
    implicit none
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer, intent(out) :: f90wrap_infil_back(2)
    type(infiltration_data_ptr_type) :: infil_back_ptr
    
    infil_back_ptr%p => m_tap_vars_infil_back
    f90wrap_infil_back = transfer(infil_back_ptr,f90wrap_infil_back)
end subroutine f90wrap_m_tap_vars__get__infil_back

subroutine f90wrap_m_tap_vars__set__infil_back(f90wrap_infil_back)
    use m_model, only: infiltration_data
    use m_tap_vars, only: m_tap_vars_infil_back => infil_back
    implicit none
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer, intent(in) :: f90wrap_infil_back(2)
    type(infiltration_data_ptr_type) :: infil_back_ptr
    
    infil_back_ptr = transfer(f90wrap_infil_back,infil_back_ptr)
    m_tap_vars_infil_back = infil_back_ptr%p
end subroutine f90wrap_m_tap_vars__set__infil_back

subroutine f90wrap_m_tap_vars__array__manning_diff(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_obs
    use m_model
    use m_tap_vars, only: m_tap_vars_manning_diff => manning_diff
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_tap_vars_manning_diff)) then
        dshape(1:1) = shape(m_tap_vars_manning_diff)
        dloc = loc(m_tap_vars_manning_diff)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_tap_vars__array__manning_diff

subroutine f90wrap_m_tap_vars__array__manning_beta_diff(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_obs
    use m_model
    use m_tap_vars, only: m_tap_vars_manning_beta_diff => manning_beta_diff
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_tap_vars_manning_beta_diff)) then
        dshape(1:1) = shape(m_tap_vars_manning_beta_diff)
        dloc = loc(m_tap_vars_manning_beta_diff)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_tap_vars__array__manning_beta_diff

subroutine f90wrap_m_tap_vars__array__bathy_cell_diff(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_obs
    use m_model
    use m_tap_vars, only: m_tap_vars_bathy_cell_diff => bathy_cell_diff
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_tap_vars_bathy_cell_diff)) then
        dshape(1:1) = shape(m_tap_vars_bathy_cell_diff)
        dloc = loc(m_tap_vars_bathy_cell_diff)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_tap_vars__array__bathy_cell_diff

subroutine f90wrap_m_tap_vars__array__manning_back(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_obs
    use m_model
    use m_tap_vars, only: m_tap_vars_manning_back => manning_back
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_tap_vars_manning_back)) then
        dshape(1:1) = shape(m_tap_vars_manning_back)
        dloc = loc(m_tap_vars_manning_back)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_tap_vars__array__manning_back

subroutine f90wrap_m_tap_vars__array__manning_beta_back(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_obs
    use m_model
    use m_tap_vars, only: m_tap_vars_manning_beta_back => manning_beta_back
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_tap_vars_manning_beta_back)) then
        dshape(1:1) = shape(m_tap_vars_manning_beta_back)
        dloc = loc(m_tap_vars_manning_beta_back)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_tap_vars__array__manning_beta_back

subroutine f90wrap_m_tap_vars__array__bathy_cell_back(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_obs
    use m_model
    use m_tap_vars, only: m_tap_vars_bathy_cell_back => bathy_cell_back
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_tap_vars_bathy_cell_back)) then
        dshape(1:1) = shape(m_tap_vars_bathy_cell_back)
        dloc = loc(m_tap_vars_bathy_cell_back)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_tap_vars__array__bathy_cell_back

! End of module m_tap_vars defined in file m_tap_vars.f90

