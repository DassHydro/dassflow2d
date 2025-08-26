! Module m_mpi defined in file m_mpi.f90

subroutine f90wrap_m_mpi__init_mpi
    use m_mpi, only: init_mpi
    implicit none
    
    call Init_MPI()
end subroutine f90wrap_m_mpi__init_mpi

subroutine f90wrap_m_mpi__end_mpi
    use m_mpi, only: end_mpi
    implicit none
    
    call End_MPI()
end subroutine f90wrap_m_mpi__end_mpi

subroutine f90wrap_m_mpi__fill_swap_lists(mesh)
    use m_mesh, only: msh
    use m_mpi, only: fill_swap_lists
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    mesh_ptr = transfer(mesh, mesh_ptr)
    call fill_swap_lists(mesh=mesh_ptr%p)
end subroutine f90wrap_m_mpi__fill_swap_lists

subroutine f90wrap_m_mpi__fill_swap_index(mesh)
    use m_mesh, only: msh
    use m_mpi, only: fill_swap_index
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    mesh_ptr = transfer(mesh, mesh_ptr)
    call fill_swap_index(mesh=mesh_ptr%p)
end subroutine f90wrap_m_mpi__fill_swap_index

subroutine f90wrap_m_mpi__fill_inv_swap_index(mesh)
    use m_mesh, only: msh
    use m_mpi, only: fill_inv_swap_index
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    mesh_ptr = transfer(mesh, mesh_ptr)
    call fill_inv_swap_index(mesh=mesh_ptr%p)
end subroutine f90wrap_m_mpi__fill_inv_swap_index

subroutine f90wrap_m_mpi__com_var_i(var, mesh, n0)
    use m_mpi, only: com_var_i
    use m_mesh, only: msh
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer(4), intent(inout), dimension(n0) :: var
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    integer :: n0
    !f2py intent(hide), depend(var) :: n0 = shape(var,0)
    mesh_ptr = transfer(mesh, mesh_ptr)
    call com_var_i(var=var, mesh=mesh_ptr%p)
end subroutine f90wrap_m_mpi__com_var_i

subroutine f90wrap_m_mpi__com_var_r(var, mesh, n0)
    use m_mpi, only: com_var_r
    use m_mesh, only: msh
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    real(8), intent(inout), dimension(n0) :: var
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    integer :: n0
    !f2py intent(hide), depend(var) :: n0 = shape(var,0)
    mesh_ptr = transfer(mesh, mesh_ptr)
    call com_var_r(var=var, mesh=mesh_ptr%p)
end subroutine f90wrap_m_mpi__com_var_r

subroutine f90wrap_m_mpi__mpi_send_recv_scal_i(to_send, to_recv, proc_send, proc_recv)
    use m_mpi, only: mpi_send_recv_scal_i
    implicit none
    
    integer(4) :: to_send
    integer(4), intent(out) :: to_recv
    integer(4), intent(in) :: proc_send
    integer(4), intent(in) :: proc_recv
    call mpi_send_recv_scal_i(to_send=to_send, to_recv=to_recv, proc_send=proc_send, proc_recv=proc_recv)
end subroutine f90wrap_m_mpi__mpi_send_recv_scal_i

subroutine f90wrap_m_mpi__mpi_send_recv_scal_r(to_send, to_recv, proc_send, proc_recv)
    use m_mpi, only: mpi_send_recv_scal_r
    implicit none
    
    real(8) :: to_send
    real(8), intent(out) :: to_recv
    integer(4), intent(in) :: proc_send
    integer(4), intent(in) :: proc_recv
    call mpi_send_recv_scal_r(to_send=to_send, to_recv=to_recv, proc_send=proc_send, proc_recv=proc_recv)
end subroutine f90wrap_m_mpi__mpi_send_recv_scal_r

subroutine f90wrap_m_mpi__mpi_send_recv_array_i(to_send, to_recv, proc_send, proc_recv, n0, n1)
    use m_mpi, only: mpi_send_recv_array_i
    implicit none
    
    integer(4), dimension(n0) :: to_send
    integer(4), intent(inout), dimension(n1) :: to_recv
    integer(4), intent(in) :: proc_send
    integer(4), intent(in) :: proc_recv
    integer :: n0
    !f2py intent(hide), depend(to_send) :: n0 = shape(to_send,0)
    integer :: n1
    !f2py intent(hide), depend(to_recv) :: n1 = shape(to_recv,0)
    call mpi_send_recv_array_i(to_send=to_send, to_recv=to_recv, proc_send=proc_send, proc_recv=proc_recv)
end subroutine f90wrap_m_mpi__mpi_send_recv_array_i

subroutine f90wrap_m_mpi__mpi_send_recv_array_r(to_send, to_recv, proc_send, proc_recv, n0, n1)
    use m_mpi, only: mpi_send_recv_array_r
    implicit none
    
    real(8), dimension(n0) :: to_send
    real(8), intent(inout), dimension(n1) :: to_recv
    integer(4), intent(in) :: proc_send
    integer(4), intent(in) :: proc_recv
    integer :: n0
    !f2py intent(hide), depend(to_send) :: n0 = shape(to_send,0)
    integer :: n1
    !f2py intent(hide), depend(to_recv) :: n1 = shape(to_recv,0)
    call mpi_send_recv_array_r(to_send=to_send, to_recv=to_recv, proc_send=proc_send, proc_recv=proc_recv)
end subroutine f90wrap_m_mpi__mpi_send_recv_array_r

subroutine f90wrap_m_mpi__mpi_sum_r(val)
    use m_mpi, only: mpi_sum_r
    implicit none
    
    real(8), intent(inout) :: val
    call mpi_sum_r(val=val)
end subroutine f90wrap_m_mpi__mpi_sum_r

subroutine f90wrap_m_mpi__mpi_sum_i(val)
    use m_mpi, only: mpi_sum_i
    implicit none
    
    integer(4), intent(inout) :: val
    call mpi_sum_i(val=val)
end subroutine f90wrap_m_mpi__mpi_sum_i

subroutine f90wrap_m_mpi__mpi_max_r(val)
    use m_mpi, only: mpi_max_r
    implicit none
    
    real(8), intent(inout) :: val
    call mpi_max_r(val=val)
end subroutine f90wrap_m_mpi__mpi_max_r

subroutine f90wrap_m_mpi__mpi_max_i(val)
    use m_mpi, only: mpi_max_i
    implicit none
    
    integer(4), intent(inout) :: val
    call mpi_max_i(val=val)
end subroutine f90wrap_m_mpi__mpi_max_i

subroutine f90wrap_m_mpi__mpi_min_r(val)
    use m_mpi, only: mpi_min_r
    implicit none
    
    real(8), intent(inout) :: val
    call mpi_min_r(val=val)
end subroutine f90wrap_m_mpi__mpi_min_r

subroutine f90wrap_m_mpi__mpi_min_i(val)
    use m_mpi, only: mpi_min_i
    implicit none
    
    integer(4), intent(inout) :: val
    call mpi_min_i(val=val)
end subroutine f90wrap_m_mpi__mpi_min_i

subroutine f90wrap_m_mpi__mpi_bcast_r(val, pr)
    use m_mpi, only: mpi_bcast_r
    implicit none
    
    real(8), intent(inout) :: val
    integer(4), intent(in) :: pr
    call mpi_bcast_r(val=val, pr=pr)
end subroutine f90wrap_m_mpi__mpi_bcast_r

subroutine f90wrap_m_mpi__mpi_bcast_i(val, pr)
    use m_mpi, only: mpi_bcast_i
    implicit none
    
    integer(4), intent(inout) :: val
    integer(4), intent(in) :: pr
    call mpi_bcast_i(val=val, pr=pr)
end subroutine f90wrap_m_mpi__mpi_bcast_i

subroutine f90wrap_m_mpi__mpi_allgather_r(val, temp, n0)
    use m_mpi, only: mpi_allgather_r
    implicit none
    
    real(8), intent(inout) :: val
    real(8), intent(inout), dimension(n0) :: temp
    integer :: n0
    !f2py intent(hide), depend(temp) :: n0 = shape(temp,0)
    call mpi_allgather_r(val=val, temp=temp)
end subroutine f90wrap_m_mpi__mpi_allgather_r

subroutine f90wrap_m_mpi__mpi_allgather_i(val, temp, n0)
    use m_mpi, only: mpi_allgather_i
    implicit none
    
    integer(4), intent(inout) :: val
    integer(4), intent(inout), dimension(n0) :: temp
    integer :: n0
    !f2py intent(hide), depend(temp) :: n0 = shape(temp,0)
    call mpi_allgather_i(val=val, temp=temp)
end subroutine f90wrap_m_mpi__mpi_allgather_i

subroutine f90wrap_m_mpi__mpi_wait_all
    use m_mpi, only: mpi_wait_all
    implicit none
    
    call mpi_wait_all()
end subroutine f90wrap_m_mpi__mpi_wait_all

subroutine f90wrap_m_mpi__stopping_program_sub(comment)
    use m_mpi, only: stopping_program_sub
    implicit none
    
    character*(*), intent(in) :: comment
    call Stopping_Program_Sub(comment=comment)
end subroutine f90wrap_m_mpi__stopping_program_sub

subroutine f90wrap_m_mpi__get__np(f90wrap_np)
    use m_mpi, only: m_mpi_np => np
    implicit none
    integer(4), intent(out) :: f90wrap_np
    
    f90wrap_np = m_mpi_np
end subroutine f90wrap_m_mpi__get__np

subroutine f90wrap_m_mpi__set__np(f90wrap_np)
    use m_mpi, only: m_mpi_np => np
    implicit none
    integer(4), intent(in) :: f90wrap_np
    
    m_mpi_np = f90wrap_np
end subroutine f90wrap_m_mpi__set__np

subroutine f90wrap_m_mpi__get__proc(f90wrap_proc)
    use m_mpi, only: m_mpi_proc => proc
    implicit none
    integer(4), intent(out) :: f90wrap_proc
    
    f90wrap_proc = m_mpi_proc
end subroutine f90wrap_m_mpi__get__proc

subroutine f90wrap_m_mpi__set__proc(f90wrap_proc)
    use m_mpi, only: m_mpi_proc => proc
    implicit none
    integer(4), intent(in) :: f90wrap_proc
    
    m_mpi_proc = f90wrap_proc
end subroutine f90wrap_m_mpi__set__proc

subroutine f90wrap_m_mpi__get__code(f90wrap_code)
    use m_mpi, only: m_mpi_code => code
    implicit none
    integer(4), intent(out) :: f90wrap_code
    
    f90wrap_code = m_mpi_code
end subroutine f90wrap_m_mpi__get__code

subroutine f90wrap_m_mpi__set__code(f90wrap_code)
    use m_mpi, only: m_mpi_code => code
    implicit none
    integer(4), intent(in) :: f90wrap_code
    
    m_mpi_code = f90wrap_code
end subroutine f90wrap_m_mpi__set__code

subroutine f90wrap_m_mpi__get__nneighb(f90wrap_nneighb)
    use m_mpi, only: m_mpi_nneighb => nneighb
    implicit none
    integer(4), intent(out) :: f90wrap_nneighb
    
    f90wrap_nneighb = m_mpi_nneighb
end subroutine f90wrap_m_mpi__get__nneighb

subroutine f90wrap_m_mpi__set__nneighb(f90wrap_nneighb)
    use m_mpi, only: m_mpi_nneighb => nneighb
    implicit none
    integer(4), intent(in) :: f90wrap_nneighb
    
    m_mpi_nneighb = f90wrap_nneighb
end subroutine f90wrap_m_mpi__set__nneighb

subroutine f90wrap_m_mpi__array__swap_index(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_mpi, only: m_mpi_swap_index => swap_index
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(m_mpi_swap_index)) then
        dshape(1:1) = shape(m_mpi_swap_index)
        dloc = loc(m_mpi_swap_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_mpi__array__swap_index

subroutine f90wrap_m_mpi__array__inv_swap_index(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_mpi, only: m_mpi_inv_swap_index => inv_swap_index
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(m_mpi_inv_swap_index)) then
        dshape(1:1) = shape(m_mpi_inv_swap_index)
        dloc = loc(m_mpi_inv_swap_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_mpi__array__inv_swap_index

subroutine f90wrap_m_mpi__array__part(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_mpi, only: m_mpi_part => part
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(m_mpi_part)) then
        dshape(1:1) = shape(m_mpi_part)
        dloc = loc(m_mpi_part)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_mpi__array__part

subroutine f90wrap_m_mpi__array__part_size(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_mpi, only: m_mpi_part_size => part_size
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(m_mpi_part_size)) then
        dshape(1:1) = shape(m_mpi_part_size)
        dloc = loc(m_mpi_part_size)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_mpi__array__part_size

subroutine f90wrap_m_mpi__array__part_neighbs(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_mpi, only: m_mpi_part_neighbs => part_neighbs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    if (allocated(m_mpi_part_neighbs)) then
        dshape(1:2) = shape(m_mpi_part_neighbs)
        dloc = loc(m_mpi_part_neighbs)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_mpi__array__part_neighbs

subroutine f90wrap_m_mpi__get__val_tmp_r(f90wrap_val_tmp_r)
    use m_mpi, only: m_mpi_val_tmp_r => val_tmp_r
    implicit none
    real(8), intent(out) :: f90wrap_val_tmp_r
    
    f90wrap_val_tmp_r = m_mpi_val_tmp_r
end subroutine f90wrap_m_mpi__get__val_tmp_r

subroutine f90wrap_m_mpi__set__val_tmp_r(f90wrap_val_tmp_r)
    use m_mpi, only: m_mpi_val_tmp_r => val_tmp_r
    implicit none
    real(8), intent(in) :: f90wrap_val_tmp_r
    
    m_mpi_val_tmp_r = f90wrap_val_tmp_r
end subroutine f90wrap_m_mpi__set__val_tmp_r

subroutine f90wrap_m_mpi__get__val_tmp_i(f90wrap_val_tmp_i)
    use m_mpi, only: m_mpi_val_tmp_i => val_tmp_i
    implicit none
    integer(4), intent(out) :: f90wrap_val_tmp_i
    
    f90wrap_val_tmp_i = m_mpi_val_tmp_i
end subroutine f90wrap_m_mpi__get__val_tmp_i

subroutine f90wrap_m_mpi__set__val_tmp_i(f90wrap_val_tmp_i)
    use m_mpi, only: m_mpi_val_tmp_i => val_tmp_i
    implicit none
    integer(4), intent(in) :: f90wrap_val_tmp_i
    
    m_mpi_val_tmp_i = f90wrap_val_tmp_i
end subroutine f90wrap_m_mpi__set__val_tmp_i

! End of module m_mpi defined in file m_mpi.f90

