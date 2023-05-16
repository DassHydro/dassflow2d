! Module m_common defined in file m_common.f90

subroutine f90wrap_weights__array__weights(this, nd, dtype, dshape, dloc)
    use m_common, only: weights
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type weights_ptr_type
        type(weights), pointer :: p => NULL()
    end type weights_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(weights_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%weights)) then
        dshape(1:1) = shape(this_ptr%p%weights)
        dloc = loc(this_ptr%p%weights)
    else
        dloc = 0
    end if
end subroutine f90wrap_weights__array__weights

subroutine f90wrap_weights_initialise(this)
    use m_common, only: weights
    implicit none
    
    type weights_ptr_type
        type(weights), pointer :: p => NULL()
    end type weights_ptr_type
    type(weights_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_weights_initialise

subroutine f90wrap_weights_finalise(this)
    use m_common, only: weights
    implicit none
    
    type weights_ptr_type
        type(weights), pointer :: p => NULL()
    end type weights_ptr_type
    type(weights_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_weights_finalise

subroutine f90wrap_test_dt_nearest(ret_test_dt_nearest, dt_to_test)
    use m_common, only: test_dt_nearest
    implicit none
    
    logical, intent(out) :: ret_test_dt_nearest
    real(8), intent(in) :: dt_to_test
    ret_test_dt_nearest = test_dt_nearest(dt_to_test=dt_to_test)
end subroutine f90wrap_test_dt_nearest

subroutine f90wrap_test_dt_just_after(ret_test_dt_just_after, dt_to_test)
    use m_common, only: test_dt_just_after
    implicit none
    
    logical, intent(out) :: ret_test_dt_just_after
    real(8), intent(in) :: dt_to_test
    ret_test_dt_just_after = test_dt_just_after(dt_to_test=dt_to_test)
end subroutine f90wrap_test_dt_just_after

subroutine f90wrap_reading_args
    use m_common, only: reading_args
    implicit none
    
    call reading_args()
end subroutine f90wrap_reading_args

subroutine f90wrap_machine_number_limits
    use m_common, only: machine_number_limits
    implicit none
    
    call machine_number_limits()
end subroutine f90wrap_machine_number_limits

subroutine f90wrap_swap_r(a, b)
    use m_common, only: swap_r
    implicit none
    
    real(8), intent(inout) :: a
    real(8), intent(inout) :: b
    call swap_r(a=a, b=b)
end subroutine f90wrap_swap_r

subroutine f90wrap_swap_i(a, b)
    use m_common, only: swap_i
    implicit none
    
    integer(4), intent(inout) :: a
    integer(4), intent(inout) :: b
    call swap_i(a=a, b=b)
end subroutine f90wrap_swap_i

subroutine f90wrap_swap_vec_r(vec, swap, n0, n1)
    use m_common, only: swap_vec_r
    implicit none
    
    real(8), intent(inout), dimension(n0) :: vec
    integer(4), dimension(n1) :: swap
    integer :: n0
    !f2py intent(hide), depend(vec) :: n0 = shape(vec,0)
    integer :: n1
    !f2py intent(hide), depend(swap) :: n1 = shape(swap,0)
    call swap_vec_r(vec=vec, swap=swap)
end subroutine f90wrap_swap_vec_r

subroutine f90wrap_swap_vec_i(vec, swap, n0, n1)
    use m_common, only: swap_vec_i
    implicit none
    
    integer(4), intent(inout), dimension(n0) :: vec
    integer(4), dimension(n1) :: swap
    integer :: n0
    !f2py intent(hide), depend(vec) :: n0 = shape(vec,0)
    integer :: n1
    !f2py intent(hide), depend(swap) :: n1 = shape(swap,0)
    call swap_vec_i(vec=vec, swap=swap)
end subroutine f90wrap_swap_vec_i

subroutine f90wrap_div_by_except_0(a, ret_div_by_except_0, b)
    use m_common, only: div_by_except_0
    implicit none
    
    real(8), intent(in) :: a
    real(8), intent(out) :: ret_div_by_except_0
    real(8), intent(in) :: b
    ret_div_by_except_0 = div_by_except_0(a=a, b=b)
end subroutine f90wrap_div_by_except_0

subroutine f90wrap_i4col_sort_a(m, n, a, n0, n1)
    use m_common, only: i4col_sort_a
    implicit none
    
    integer(4) :: m
    integer(4) :: n
    integer(4), intent(inout), dimension(n0,n1) :: a
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,0)
    integer :: n1
    !f2py intent(hide), depend(a) :: n1 = shape(a,1)
    call i4col_sort_a(m=m, n=n, a=a)
end subroutine f90wrap_i4col_sort_a

subroutine f90wrap_sort_heap_external(n, indx, i, j, isgn)
    use m_common, only: sort_heap_external
    implicit none
    
    integer(4) :: n
    integer(4), intent(inout) :: indx
    integer(4) :: i
    integer(4) :: j
    integer(4) :: isgn
    call sort_heap_external(n=n, indx=indx, i=i, j=j, isgn=isgn)
end subroutine f90wrap_sort_heap_external

subroutine f90wrap_i4col_swap(m, n, a, i, j, n0, n1)
    use m_common, only: i4col_swap
    implicit none
    
    integer(4) :: m
    integer(4) :: n
    integer(4), intent(inout), dimension(n0,n1) :: a
    integer(4) :: i
    integer(4) :: j
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,0)
    integer :: n1
    !f2py intent(hide), depend(a) :: n1 = shape(a,1)
    call i4col_swap(m=m, n=n, a=a, i=i, j=j)
end subroutine f90wrap_i4col_swap

subroutine f90wrap_i4col_compare(m, n, a, i, j, isgn, n0, n1)
    use m_common, only: i4col_compare
    implicit none
    
    integer(4) :: m
    integer(4) :: n
    integer(4), dimension(n0,n1) :: a
    integer(4) :: i
    integer(4) :: j
    integer(4), intent(out) :: isgn
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,0)
    integer :: n1
    !f2py intent(hide), depend(a) :: n1 = shape(a,1)
    call i4col_compare(m=m, n=n, a=a, i=i, j=j, isgn=isgn)
end subroutine f90wrap_i4col_compare

subroutine f90wrap_file_name_ext(file_name, ret_file_name_res, typ)
    use m_common, only: file_name_ext
    implicit none
    
    character*(*), intent(in) :: file_name
    character(1024), intent(out) :: ret_file_name_res
    character*(*), intent(in) :: typ
    ret_file_name_res = file_name_ext(file_name=file_name, typ=typ)
end subroutine f90wrap_file_name_ext

subroutine f90wrap_count_lines(ret_nb_lines, file_name)
    use m_common, only: count_lines
    implicit none
    
    integer(4), intent(out) :: ret_nb_lines
    character*(*), intent(in) :: file_name
    ret_nb_lines = count_lines(file_name=file_name)
end subroutine f90wrap_count_lines

subroutine f90wrap_m_common__get__ip(f90wrap_ip)
    use m_common, only: m_common_ip => ip
    implicit none
    integer, intent(out) :: f90wrap_ip
    
    f90wrap_ip = m_common_ip
end subroutine f90wrap_m_common__get__ip

subroutine f90wrap_m_common__get__rp(f90wrap_rp)
    use m_common, only: m_common_rp => rp
    implicit none
    integer, intent(out) :: f90wrap_rp
    
    f90wrap_rp = m_common_rp
end subroutine f90wrap_m_common__get__rp

subroutine f90wrap_m_common__get__lchar(f90wrap_lchar)
    use m_common, only: m_common_lchar => lchar
    implicit none
    integer, intent(out) :: f90wrap_lchar
    
    f90wrap_lchar = m_common_lchar
end subroutine f90wrap_m_common__get__lchar

subroutine f90wrap_m_common__get__i(f90wrap_i)
    use m_common, only: m_common_i => i
    implicit none
    integer(4), intent(out) :: f90wrap_i
    
    f90wrap_i = m_common_i
end subroutine f90wrap_m_common__get__i

subroutine f90wrap_m_common__set__i(f90wrap_i)
    use m_common, only: m_common_i => i
    implicit none
    integer(4), intent(in) :: f90wrap_i
    
    m_common_i = f90wrap_i
end subroutine f90wrap_m_common__set__i

subroutine f90wrap_m_common__get__ie(f90wrap_ie)
    use m_common, only: m_common_ie => ie
    implicit none
    integer(4), intent(out) :: f90wrap_ie
    
    f90wrap_ie = m_common_ie
end subroutine f90wrap_m_common__get__ie

subroutine f90wrap_m_common__set__ie(f90wrap_ie)
    use m_common, only: m_common_ie => ie
    implicit none
    integer(4), intent(in) :: f90wrap_ie
    
    m_common_ie = f90wrap_ie
end subroutine f90wrap_m_common__set__ie

subroutine f90wrap_m_common__get__iK(f90wrap_iK)
    use m_common, only: m_common_iK => iK
    implicit none
    integer(4), intent(out) :: f90wrap_iK
    
    f90wrap_iK = m_common_iK
end subroutine f90wrap_m_common__get__iK

subroutine f90wrap_m_common__set__iK(f90wrap_iK)
    use m_common, only: m_common_iK => iK
    implicit none
    integer(4), intent(in) :: f90wrap_iK
    
    m_common_iK = f90wrap_iK
end subroutine f90wrap_m_common__set__iK

subroutine f90wrap_m_common__get__iKe(f90wrap_iKe)
    use m_common, only: m_common_iKe => iKe
    implicit none
    integer(4), intent(out) :: f90wrap_iKe
    
    f90wrap_iKe = m_common_iKe
end subroutine f90wrap_m_common__get__iKe

subroutine f90wrap_m_common__set__iKe(f90wrap_iKe)
    use m_common, only: m_common_iKe => iKe
    implicit none
    integer(4), intent(in) :: f90wrap_iKe
    
    m_common_iKe = f90wrap_iKe
end subroutine f90wrap_m_common__set__iKe

subroutine f90wrap_m_common__get__ib(f90wrap_ib)
    use m_common, only: m_common_ib => ib
    implicit none
    integer(4), intent(out) :: f90wrap_ib
    
    f90wrap_ib = m_common_ib
end subroutine f90wrap_m_common__get__ib

subroutine f90wrap_m_common__set__ib(f90wrap_ib)
    use m_common, only: m_common_ib => ib
    implicit none
    integer(4), intent(in) :: f90wrap_ib
    
    m_common_ib = f90wrap_ib
end subroutine f90wrap_m_common__set__ib

subroutine f90wrap_m_common__get__j(f90wrap_j)
    use m_common, only: m_common_j => j
    implicit none
    integer(4), intent(out) :: f90wrap_j
    
    f90wrap_j = m_common_j
end subroutine f90wrap_m_common__get__j

subroutine f90wrap_m_common__set__j(f90wrap_j)
    use m_common, only: m_common_j => j
    implicit none
    integer(4), intent(in) :: f90wrap_j
    
    m_common_j = f90wrap_j
end subroutine f90wrap_m_common__set__j

subroutine f90wrap_m_common__get__je(f90wrap_je)
    use m_common, only: m_common_je => je
    implicit none
    integer(4), intent(out) :: f90wrap_je
    
    f90wrap_je = m_common_je
end subroutine f90wrap_m_common__get__je

subroutine f90wrap_m_common__set__je(f90wrap_je)
    use m_common, only: m_common_je => je
    implicit none
    integer(4), intent(in) :: f90wrap_je
    
    m_common_je = f90wrap_je
end subroutine f90wrap_m_common__set__je

subroutine f90wrap_m_common__get__jK(f90wrap_jK)
    use m_common, only: m_common_jK => jK
    implicit none
    integer(4), intent(out) :: f90wrap_jK
    
    f90wrap_jK = m_common_jK
end subroutine f90wrap_m_common__get__jK

subroutine f90wrap_m_common__set__jK(f90wrap_jK)
    use m_common, only: m_common_jK => jK
    implicit none
    integer(4), intent(in) :: f90wrap_jK
    
    m_common_jK = f90wrap_jK
end subroutine f90wrap_m_common__set__jK

subroutine f90wrap_m_common__get__jKe(f90wrap_jKe)
    use m_common, only: m_common_jKe => jKe
    implicit none
    integer(4), intent(out) :: f90wrap_jKe
    
    f90wrap_jKe = m_common_jKe
end subroutine f90wrap_m_common__get__jKe

subroutine f90wrap_m_common__set__jKe(f90wrap_jKe)
    use m_common, only: m_common_jKe => jKe
    implicit none
    integer(4), intent(in) :: f90wrap_jKe
    
    m_common_jKe = f90wrap_jKe
end subroutine f90wrap_m_common__set__jKe

subroutine f90wrap_m_common__get__jb(f90wrap_jb)
    use m_common, only: m_common_jb => jb
    implicit none
    integer(4), intent(out) :: f90wrap_jb
    
    f90wrap_jb = m_common_jb
end subroutine f90wrap_m_common__get__jb

subroutine f90wrap_m_common__set__jb(f90wrap_jb)
    use m_common, only: m_common_jb => jb
    implicit none
    integer(4), intent(in) :: f90wrap_jb
    
    m_common_jb = f90wrap_jb
end subroutine f90wrap_m_common__set__jb

subroutine f90wrap_m_common__get__k(f90wrap_k)
    use m_common, only: m_common_k => k
    implicit none
    integer(4), intent(out) :: f90wrap_k
    
    f90wrap_k = m_common_k
end subroutine f90wrap_m_common__get__k

subroutine f90wrap_m_common__set__k(f90wrap_k)
    use m_common, only: m_common_k => k
    implicit none
    integer(4), intent(in) :: f90wrap_k
    
    m_common_k = f90wrap_k
end subroutine f90wrap_m_common__set__k

subroutine f90wrap_m_common__get__ke(f90wrap_ke)
    use m_common, only: m_common_ke => ke
    implicit none
    integer(4), intent(out) :: f90wrap_ke
    
    f90wrap_ke = m_common_ke
end subroutine f90wrap_m_common__get__ke

subroutine f90wrap_m_common__set__ke(f90wrap_ke)
    use m_common, only: m_common_ke => ke
    implicit none
    integer(4), intent(in) :: f90wrap_ke
    
    m_common_ke = f90wrap_ke
end subroutine f90wrap_m_common__set__ke

subroutine f90wrap_m_common__get__kK(f90wrap_kK)
    use m_common, only: m_common_kK => kK
    implicit none
    integer(4), intent(out) :: f90wrap_kK
    
    f90wrap_kK = m_common_kK
end subroutine f90wrap_m_common__get__kK

subroutine f90wrap_m_common__set__kK(f90wrap_kK)
    use m_common, only: m_common_kK => kK
    implicit none
    integer(4), intent(in) :: f90wrap_kK
    
    m_common_kK = f90wrap_kK
end subroutine f90wrap_m_common__set__kK

subroutine f90wrap_m_common__get__kKe(f90wrap_kKe)
    use m_common, only: m_common_kKe => kKe
    implicit none
    integer(4), intent(out) :: f90wrap_kKe
    
    f90wrap_kKe = m_common_kKe
end subroutine f90wrap_m_common__get__kKe

subroutine f90wrap_m_common__set__kKe(f90wrap_kKe)
    use m_common, only: m_common_kKe => kKe
    implicit none
    integer(4), intent(in) :: f90wrap_kKe
    
    m_common_kKe = f90wrap_kKe
end subroutine f90wrap_m_common__set__kKe

subroutine f90wrap_m_common__get__kb(f90wrap_kb)
    use m_common, only: m_common_kb => kb
    implicit none
    integer(4), intent(out) :: f90wrap_kb
    
    f90wrap_kb = m_common_kb
end subroutine f90wrap_m_common__get__kb

subroutine f90wrap_m_common__set__kb(f90wrap_kb)
    use m_common, only: m_common_kb => kb
    implicit none
    integer(4), intent(in) :: f90wrap_kb
    
    m_common_kb = f90wrap_kb
end subroutine f90wrap_m_common__set__kb

subroutine f90wrap_m_common__get__mesh_type(f90wrap_mesh_type)
    use m_common, only: m_common_mesh_type => mesh_type
    implicit none
    character(1024), intent(out) :: f90wrap_mesh_type
    
    f90wrap_mesh_type = m_common_mesh_type
end subroutine f90wrap_m_common__get__mesh_type

subroutine f90wrap_m_common__set__mesh_type(f90wrap_mesh_type)
    use m_common, only: m_common_mesh_type => mesh_type
    implicit none
    character(1024), intent(in) :: f90wrap_mesh_type
    
    m_common_mesh_type = f90wrap_mesh_type
end subroutine f90wrap_m_common__set__mesh_type

subroutine f90wrap_m_common__get__mesh_name(f90wrap_mesh_name)
    use m_common, only: m_common_mesh_name => mesh_name
    implicit none
    character(1024), intent(out) :: f90wrap_mesh_name
    
    f90wrap_mesh_name = m_common_mesh_name
end subroutine f90wrap_m_common__get__mesh_name

subroutine f90wrap_m_common__set__mesh_name(f90wrap_mesh_name)
    use m_common, only: m_common_mesh_name => mesh_name
    implicit none
    character(1024), intent(in) :: f90wrap_mesh_name
    
    m_common_mesh_name = f90wrap_mesh_name
end subroutine f90wrap_m_common__set__mesh_name

subroutine f90wrap_m_common__get__bc_N(f90wrap_bc_N)
    use m_common, only: m_common_bc_N => bc_N
    implicit none
    character(1024), intent(out) :: f90wrap_bc_N
    
    f90wrap_bc_N = m_common_bc_N
end subroutine f90wrap_m_common__get__bc_N

subroutine f90wrap_m_common__set__bc_N(f90wrap_bc_N)
    use m_common, only: m_common_bc_N => bc_N
    implicit none
    character(1024), intent(in) :: f90wrap_bc_N
    
    m_common_bc_N = f90wrap_bc_N
end subroutine f90wrap_m_common__set__bc_N

subroutine f90wrap_m_common__get__bc_S(f90wrap_bc_S)
    use m_common, only: m_common_bc_S => bc_S
    implicit none
    character(1024), intent(out) :: f90wrap_bc_S
    
    f90wrap_bc_S = m_common_bc_S
end subroutine f90wrap_m_common__get__bc_S

subroutine f90wrap_m_common__set__bc_S(f90wrap_bc_S)
    use m_common, only: m_common_bc_S => bc_S
    implicit none
    character(1024), intent(in) :: f90wrap_bc_S
    
    m_common_bc_S = f90wrap_bc_S
end subroutine f90wrap_m_common__set__bc_S

subroutine f90wrap_m_common__get__bc_W(f90wrap_bc_W)
    use m_common, only: m_common_bc_W => bc_W
    implicit none
    character(1024), intent(out) :: f90wrap_bc_W
    
    f90wrap_bc_W = m_common_bc_W
end subroutine f90wrap_m_common__get__bc_W

subroutine f90wrap_m_common__set__bc_W(f90wrap_bc_W)
    use m_common, only: m_common_bc_W => bc_W
    implicit none
    character(1024), intent(in) :: f90wrap_bc_W
    
    m_common_bc_W = f90wrap_bc_W
end subroutine f90wrap_m_common__set__bc_W

subroutine f90wrap_m_common__get__bc_E(f90wrap_bc_E)
    use m_common, only: m_common_bc_E => bc_E
    implicit none
    character(1024), intent(out) :: f90wrap_bc_E
    
    f90wrap_bc_E = m_common_bc_E
end subroutine f90wrap_m_common__get__bc_E

subroutine f90wrap_m_common__set__bc_E(f90wrap_bc_E)
    use m_common, only: m_common_bc_E => bc_E
    implicit none
    character(1024), intent(in) :: f90wrap_bc_E
    
    m_common_bc_E = f90wrap_bc_E
end subroutine f90wrap_m_common__set__bc_E

subroutine f90wrap_m_common__get__bc_rain(f90wrap_bc_rain)
    use m_common, only: m_common_bc_rain => bc_rain
    implicit none
    integer(4), intent(out) :: f90wrap_bc_rain
    
    f90wrap_bc_rain = m_common_bc_rain
end subroutine f90wrap_m_common__get__bc_rain

subroutine f90wrap_m_common__set__bc_rain(f90wrap_bc_rain)
    use m_common, only: m_common_bc_rain => bc_rain
    implicit none
    integer(4), intent(in) :: f90wrap_bc_rain
    
    m_common_bc_rain = f90wrap_bc_rain
end subroutine f90wrap_m_common__set__bc_rain

subroutine f90wrap_m_common__get__bc_infil(f90wrap_bc_infil)
    use m_common, only: m_common_bc_infil => bc_infil
    implicit none
    integer(4), intent(out) :: f90wrap_bc_infil
    
    f90wrap_bc_infil = m_common_bc_infil
end subroutine f90wrap_m_common__get__bc_infil

subroutine f90wrap_m_common__set__bc_infil(f90wrap_bc_infil)
    use m_common, only: m_common_bc_infil => bc_infil
    implicit none
    integer(4), intent(in) :: f90wrap_bc_infil
    
    m_common_bc_infil = f90wrap_bc_infil
end subroutine f90wrap_m_common__set__bc_infil

subroutine f90wrap_m_common__get__lx(f90wrap_lx)
    use m_common, only: m_common_lx => lx
    implicit none
    real(8), intent(out) :: f90wrap_lx
    
    f90wrap_lx = m_common_lx
end subroutine f90wrap_m_common__get__lx

subroutine f90wrap_m_common__set__lx(f90wrap_lx)
    use m_common, only: m_common_lx => lx
    implicit none
    real(8), intent(in) :: f90wrap_lx
    
    m_common_lx = f90wrap_lx
end subroutine f90wrap_m_common__set__lx

subroutine f90wrap_m_common__get__ly(f90wrap_ly)
    use m_common, only: m_common_ly => ly
    implicit none
    real(8), intent(out) :: f90wrap_ly
    
    f90wrap_ly = m_common_ly
end subroutine f90wrap_m_common__get__ly

subroutine f90wrap_m_common__set__ly(f90wrap_ly)
    use m_common, only: m_common_ly => ly
    implicit none
    real(8), intent(in) :: f90wrap_ly
    
    m_common_ly = f90wrap_ly
end subroutine f90wrap_m_common__set__ly

subroutine f90wrap_m_common__get__nx(f90wrap_nx)
    use m_common, only: m_common_nx => nx
    implicit none
    integer(4), intent(out) :: f90wrap_nx
    
    f90wrap_nx = m_common_nx
end subroutine f90wrap_m_common__get__nx

subroutine f90wrap_m_common__set__nx(f90wrap_nx)
    use m_common, only: m_common_nx => nx
    implicit none
    integer(4), intent(in) :: f90wrap_nx
    
    m_common_nx = f90wrap_nx
end subroutine f90wrap_m_common__set__nx

subroutine f90wrap_m_common__get__ny(f90wrap_ny)
    use m_common, only: m_common_ny => ny
    implicit none
    integer(4), intent(out) :: f90wrap_ny
    
    f90wrap_ny = m_common_ny
end subroutine f90wrap_m_common__get__ny

subroutine f90wrap_m_common__set__ny(f90wrap_ny)
    use m_common, only: m_common_ny => ny
    implicit none
    integer(4), intent(in) :: f90wrap_ny
    
    m_common_ny = f90wrap_ny
end subroutine f90wrap_m_common__set__ny

subroutine f90wrap_m_common__get__ts(f90wrap_ts)
    use m_common, only: m_common_ts => ts
    implicit none
    real(8), intent(out) :: f90wrap_ts
    
    f90wrap_ts = m_common_ts
end subroutine f90wrap_m_common__get__ts

subroutine f90wrap_m_common__set__ts(f90wrap_ts)
    use m_common, only: m_common_ts => ts
    implicit none
    real(8), intent(in) :: f90wrap_ts
    
    m_common_ts = f90wrap_ts
end subroutine f90wrap_m_common__set__ts

subroutine f90wrap_m_common__get__adapt_dt(f90wrap_adapt_dt)
    use m_common, only: m_common_adapt_dt => adapt_dt
    implicit none
    integer(4), intent(out) :: f90wrap_adapt_dt
    
    f90wrap_adapt_dt = m_common_adapt_dt
end subroutine f90wrap_m_common__get__adapt_dt

subroutine f90wrap_m_common__set__adapt_dt(f90wrap_adapt_dt)
    use m_common, only: m_common_adapt_dt => adapt_dt
    implicit none
    integer(4), intent(in) :: f90wrap_adapt_dt
    
    m_common_adapt_dt = f90wrap_adapt_dt
end subroutine f90wrap_m_common__set__adapt_dt

subroutine f90wrap_m_common__get__dt(f90wrap_dt)
    use m_common, only: m_common_dt => dt
    implicit none
    real(8), intent(out) :: f90wrap_dt
    
    f90wrap_dt = m_common_dt
end subroutine f90wrap_m_common__get__dt

subroutine f90wrap_m_common__set__dt(f90wrap_dt)
    use m_common, only: m_common_dt => dt
    implicit none
    real(8), intent(in) :: f90wrap_dt
    
    m_common_dt = f90wrap_dt
end subroutine f90wrap_m_common__set__dt

subroutine f90wrap_m_common__get__cfl(f90wrap_cfl)
    use m_common, only: m_common_cfl => cfl
    implicit none
    real(8), intent(out) :: f90wrap_cfl
    
    f90wrap_cfl = m_common_cfl
end subroutine f90wrap_m_common__get__cfl

subroutine f90wrap_m_common__set__cfl(f90wrap_cfl)
    use m_common, only: m_common_cfl => cfl
    implicit none
    real(8), intent(in) :: f90wrap_cfl
    
    m_common_cfl = f90wrap_cfl
end subroutine f90wrap_m_common__set__cfl

subroutine f90wrap_m_common__get__dtw(f90wrap_dtw)
    use m_common, only: m_common_dtw => dtw
    implicit none
    real(8), intent(out) :: f90wrap_dtw
    
    f90wrap_dtw = m_common_dtw
end subroutine f90wrap_m_common__get__dtw

subroutine f90wrap_m_common__set__dtw(f90wrap_dtw)
    use m_common, only: m_common_dtw => dtw
    implicit none
    real(8), intent(in) :: f90wrap_dtw
    
    m_common_dtw = f90wrap_dtw
end subroutine f90wrap_m_common__set__dtw

subroutine f90wrap_m_common__get__dtp(f90wrap_dtp)
    use m_common, only: m_common_dtp => dtp
    implicit none
    real(8), intent(out) :: f90wrap_dtp
    
    f90wrap_dtp = m_common_dtp
end subroutine f90wrap_m_common__get__dtp

subroutine f90wrap_m_common__set__dtp(f90wrap_dtp)
    use m_common, only: m_common_dtp => dtp
    implicit none
    real(8), intent(in) :: f90wrap_dtp
    
    m_common_dtp = f90wrap_dtp
end subroutine f90wrap_m_common__set__dtp

subroutine f90wrap_m_common__get__dta(f90wrap_dta)
    use m_common, only: m_common_dta => dta
    implicit none
    real(8), intent(out) :: f90wrap_dta
    
    f90wrap_dta = m_common_dta
end subroutine f90wrap_m_common__get__dta

subroutine f90wrap_m_common__set__dta(f90wrap_dta)
    use m_common, only: m_common_dta => dta
    implicit none
    real(8), intent(in) :: f90wrap_dta
    
    m_common_dta = f90wrap_dta
end subroutine f90wrap_m_common__set__dta

subroutine f90wrap_m_common__get__w_tecplot(f90wrap_w_tecplot)
    use m_common, only: m_common_w_tecplot => w_tecplot
    implicit none
    integer(4), intent(out) :: f90wrap_w_tecplot
    
    f90wrap_w_tecplot = m_common_w_tecplot
end subroutine f90wrap_m_common__get__w_tecplot

subroutine f90wrap_m_common__set__w_tecplot(f90wrap_w_tecplot)
    use m_common, only: m_common_w_tecplot => w_tecplot
    implicit none
    integer(4), intent(in) :: f90wrap_w_tecplot
    
    m_common_w_tecplot = f90wrap_w_tecplot
end subroutine f90wrap_m_common__set__w_tecplot

subroutine f90wrap_m_common__get__w_vtk(f90wrap_w_vtk)
    use m_common, only: m_common_w_vtk => w_vtk
    implicit none
    integer(4), intent(out) :: f90wrap_w_vtk
    
    f90wrap_w_vtk = m_common_w_vtk
end subroutine f90wrap_m_common__get__w_vtk

subroutine f90wrap_m_common__set__w_vtk(f90wrap_w_vtk)
    use m_common, only: m_common_w_vtk => w_vtk
    implicit none
    integer(4), intent(in) :: f90wrap_w_vtk
    
    m_common_w_vtk = f90wrap_w_vtk
end subroutine f90wrap_m_common__set__w_vtk

subroutine f90wrap_m_common__get__w_gnuplot(f90wrap_w_gnuplot)
    use m_common, only: m_common_w_gnuplot => w_gnuplot
    implicit none
    integer(4), intent(out) :: f90wrap_w_gnuplot
    
    f90wrap_w_gnuplot = m_common_w_gnuplot
end subroutine f90wrap_m_common__get__w_gnuplot

subroutine f90wrap_m_common__set__w_gnuplot(f90wrap_w_gnuplot)
    use m_common, only: m_common_w_gnuplot => w_gnuplot
    implicit none
    integer(4), intent(in) :: f90wrap_w_gnuplot
    
    m_common_w_gnuplot = f90wrap_w_gnuplot
end subroutine f90wrap_m_common__set__w_gnuplot

subroutine f90wrap_m_common__get__w_bin(f90wrap_w_bin)
    use m_common, only: m_common_w_bin => w_bin
    implicit none
    integer(4), intent(out) :: f90wrap_w_bin
    
    f90wrap_w_bin = m_common_w_bin
end subroutine f90wrap_m_common__get__w_bin

subroutine f90wrap_m_common__set__w_bin(f90wrap_w_bin)
    use m_common, only: m_common_w_bin => w_bin
    implicit none
    integer(4), intent(in) :: f90wrap_w_bin
    
    m_common_w_bin = f90wrap_w_bin
end subroutine f90wrap_m_common__set__w_bin

subroutine f90wrap_m_common__get__w_exact(f90wrap_w_exact)
    use m_common, only: m_common_w_exact => w_exact
    implicit none
    integer(4), intent(out) :: f90wrap_w_exact
    
    f90wrap_w_exact = m_common_w_exact
end subroutine f90wrap_m_common__get__w_exact

subroutine f90wrap_m_common__set__w_exact(f90wrap_w_exact)
    use m_common, only: m_common_w_exact => w_exact
    implicit none
    integer(4), intent(in) :: f90wrap_w_exact
    
    m_common_w_exact = f90wrap_w_exact
end subroutine f90wrap_m_common__set__w_exact

subroutine f90wrap_m_common__get__w_norm(f90wrap_w_norm)
    use m_common, only: m_common_w_norm => w_norm
    implicit none
    integer(4), intent(out) :: f90wrap_w_norm
    
    f90wrap_w_norm = m_common_w_norm
end subroutine f90wrap_m_common__get__w_norm

subroutine f90wrap_m_common__set__w_norm(f90wrap_w_norm)
    use m_common, only: m_common_w_norm => w_norm
    implicit none
    integer(4), intent(in) :: f90wrap_w_norm
    
    m_common_w_norm = f90wrap_w_norm
end subroutine f90wrap_m_common__set__w_norm

subroutine f90wrap_m_common__get__w_obs(f90wrap_w_obs)
    use m_common, only: m_common_w_obs => w_obs
    implicit none
    integer(4), intent(out) :: f90wrap_w_obs
    
    f90wrap_w_obs = m_common_w_obs
end subroutine f90wrap_m_common__get__w_obs

subroutine f90wrap_m_common__set__w_obs(f90wrap_w_obs)
    use m_common, only: m_common_w_obs => w_obs
    implicit none
    integer(4), intent(in) :: f90wrap_w_obs
    
    m_common_w_obs = f90wrap_w_obs
end subroutine f90wrap_m_common__set__w_obs

subroutine f90wrap_m_common__get__use_obs(f90wrap_use_obs)
    use m_common, only: m_common_use_obs => use_obs
    implicit none
    integer(4), intent(out) :: f90wrap_use_obs
    
    f90wrap_use_obs = m_common_use_obs
end subroutine f90wrap_m_common__get__use_obs

subroutine f90wrap_m_common__set__use_obs(f90wrap_use_obs)
    use m_common, only: m_common_use_obs => use_obs
    implicit none
    integer(4), intent(in) :: f90wrap_use_obs
    
    m_common_use_obs = f90wrap_use_obs
end subroutine f90wrap_m_common__set__use_obs

subroutine f90wrap_m_common__get__spatial_scheme(f90wrap_spatial_scheme)
    use m_common, only: m_common_spatial_scheme => spatial_scheme
    implicit none
    character(1024), intent(out) :: f90wrap_spatial_scheme
    
    f90wrap_spatial_scheme = m_common_spatial_scheme
end subroutine f90wrap_m_common__get__spatial_scheme

subroutine f90wrap_m_common__set__spatial_scheme(f90wrap_spatial_scheme)
    use m_common, only: m_common_spatial_scheme => spatial_scheme
    implicit none
    character(1024), intent(in) :: f90wrap_spatial_scheme
    
    m_common_spatial_scheme = f90wrap_spatial_scheme
end subroutine f90wrap_m_common__set__spatial_scheme

subroutine f90wrap_m_common__get__temp_scheme(f90wrap_temp_scheme)
    use m_common, only: m_common_temp_scheme => temp_scheme
    implicit none
    character(1024), intent(out) :: f90wrap_temp_scheme
    
    f90wrap_temp_scheme = m_common_temp_scheme
end subroutine f90wrap_m_common__get__temp_scheme

subroutine f90wrap_m_common__set__temp_scheme(f90wrap_temp_scheme)
    use m_common, only: m_common_temp_scheme => temp_scheme
    implicit none
    character(1024), intent(in) :: f90wrap_temp_scheme
    
    m_common_temp_scheme = f90wrap_temp_scheme
end subroutine f90wrap_m_common__set__temp_scheme

subroutine f90wrap_m_common__array__args(dummy_this, nd, dtype, dshape, dloc)
    use m_common, only: m_common_args => args
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    if (allocated(m_common_args)) then
        dshape(1:2) = (/len(m_common_args(1)), shape(m_common_args)/)
        dloc = loc(m_common_args)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_common__array__args

subroutine f90wrap_m_common__get__max_nt_for_direct(f90wrap_max_nt_for_direct)
    use m_common, only: m_common_max_nt_for_direct => max_nt_for_direct
    implicit none
    integer(4), intent(out) :: f90wrap_max_nt_for_direct
    
    f90wrap_max_nt_for_direct = m_common_max_nt_for_direct
end subroutine f90wrap_m_common__get__max_nt_for_direct

subroutine f90wrap_m_common__set__max_nt_for_direct(f90wrap_max_nt_for_direct)
    use m_common, only: m_common_max_nt_for_direct => max_nt_for_direct
    implicit none
    integer(4), intent(in) :: f90wrap_max_nt_for_direct
    
    m_common_max_nt_for_direct = f90wrap_max_nt_for_direct
end subroutine f90wrap_m_common__set__max_nt_for_direct

subroutine f90wrap_m_common__get__max_nt_for_adjoint(f90wrap_max_nt_for_adjoint)
    use m_common, only: m_common_max_nt_for_adjoint => max_nt_for_adjoint
    implicit none
    integer(4), intent(out) :: f90wrap_max_nt_for_adjoint
    
    f90wrap_max_nt_for_adjoint = m_common_max_nt_for_adjoint
end subroutine f90wrap_m_common__get__max_nt_for_adjoint

subroutine f90wrap_m_common__set__max_nt_for_adjoint(f90wrap_max_nt_for_adjoint)
    use m_common, only: m_common_max_nt_for_adjoint => max_nt_for_adjoint
    implicit none
    integer(4), intent(in) :: f90wrap_max_nt_for_adjoint
    
    m_common_max_nt_for_adjoint = f90wrap_max_nt_for_adjoint
end subroutine f90wrap_m_common__set__max_nt_for_adjoint

subroutine f90wrap_m_common__get__length_real(f90wrap_length_real)
    use m_common, only: m_common_length_real => length_real
    implicit none
    integer(4), intent(out) :: f90wrap_length_real
    
    f90wrap_length_real = m_common_length_real
end subroutine f90wrap_m_common__get__length_real

subroutine f90wrap_m_common__set__length_real(f90wrap_length_real)
    use m_common, only: m_common_length_real => length_real
    implicit none
    integer(4), intent(in) :: f90wrap_length_real
    
    m_common_length_real = f90wrap_length_real
end subroutine f90wrap_m_common__set__length_real

subroutine f90wrap_m_common__get__nt(f90wrap_nt)
    use m_common, only: m_common_nt => nt
    implicit none
    integer(4), intent(out) :: f90wrap_nt
    
    f90wrap_nt = m_common_nt
end subroutine f90wrap_m_common__get__nt

subroutine f90wrap_m_common__set__nt(f90wrap_nt)
    use m_common, only: m_common_nt => nt
    implicit none
    integer(4), intent(in) :: f90wrap_nt
    
    m_common_nt = f90wrap_nt
end subroutine f90wrap_m_common__set__nt

subroutine f90wrap_m_common__get__nt0(f90wrap_nt0)
    use m_common, only: m_common_nt0 => nt0
    implicit none
    integer(4), intent(out) :: f90wrap_nt0
    
    f90wrap_nt0 = m_common_nt0
end subroutine f90wrap_m_common__get__nt0

subroutine f90wrap_m_common__set__nt0(f90wrap_nt0)
    use m_common, only: m_common_nt0 => nt0
    implicit none
    integer(4), intent(in) :: f90wrap_nt0
    
    m_common_nt0 = f90wrap_nt0
end subroutine f90wrap_m_common__set__nt0

subroutine f90wrap_m_common__get__tc(f90wrap_tc)
    use m_common, only: m_common_tc => tc
    implicit none
    real(8), intent(out) :: f90wrap_tc
    
    f90wrap_tc = m_common_tc
end subroutine f90wrap_m_common__get__tc

subroutine f90wrap_m_common__set__tc(f90wrap_tc)
    use m_common, only: m_common_tc => tc
    implicit none
    real(8), intent(in) :: f90wrap_tc
    
    m_common_tc = f90wrap_tc
end subroutine f90wrap_m_common__set__tc

subroutine f90wrap_m_common__get__tc0(f90wrap_tc0)
    use m_common, only: m_common_tc0 => tc0
    implicit none
    real(8), intent(out) :: f90wrap_tc0
    
    f90wrap_tc0 = m_common_tc0
end subroutine f90wrap_m_common__get__tc0

subroutine f90wrap_m_common__set__tc0(f90wrap_tc0)
    use m_common, only: m_common_tc0 => tc0
    implicit none
    real(8), intent(in) :: f90wrap_tc0
    
    m_common_tc0 = f90wrap_tc0
end subroutine f90wrap_m_common__set__tc0

subroutine f90wrap_m_common__get__end_time_loop(f90wrap_end_time_loop)
    use m_common, only: m_common_end_time_loop => end_time_loop
    implicit none
    logical, intent(out) :: f90wrap_end_time_loop
    
    f90wrap_end_time_loop = m_common_end_time_loop
end subroutine f90wrap_m_common__get__end_time_loop

subroutine f90wrap_m_common__set__end_time_loop(f90wrap_end_time_loop)
    use m_common, only: m_common_end_time_loop => end_time_loop
    implicit none
    logical, intent(in) :: f90wrap_end_time_loop
    
    m_common_end_time_loop = f90wrap_end_time_loop
end subroutine f90wrap_m_common__set__end_time_loop

subroutine f90wrap_m_common__get__dx(f90wrap_dx)
    use m_common, only: m_common_dx => dx
    implicit none
    real(8), intent(out) :: f90wrap_dx
    
    f90wrap_dx = m_common_dx
end subroutine f90wrap_m_common__get__dx

subroutine f90wrap_m_common__set__dx(f90wrap_dx)
    use m_common, only: m_common_dx => dx
    implicit none
    real(8), intent(in) :: f90wrap_dx
    
    m_common_dx = f90wrap_dx
end subroutine f90wrap_m_common__set__dx

subroutine f90wrap_m_common__get__dy(f90wrap_dy)
    use m_common, only: m_common_dy => dy
    implicit none
    real(8), intent(out) :: f90wrap_dy
    
    f90wrap_dy = m_common_dy
end subroutine f90wrap_m_common__get__dy

subroutine f90wrap_m_common__set__dy(f90wrap_dy)
    use m_common, only: m_common_dy => dy
    implicit none
    real(8), intent(in) :: f90wrap_dy
    
    m_common_dy = f90wrap_dy
end subroutine f90wrap_m_common__set__dy

subroutine f90wrap_m_common__array__is_file_open(dummy_this, nd, dtype, dshape, dloc)
    use m_common, only: m_common_is_file_open => is_file_open
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    dshape(1:2) = (/len(m_common_is_file_open(1)), shape(m_common_is_file_open)/)
    dloc = loc(m_common_is_file_open)
end subroutine f90wrap_m_common__array__is_file_open

subroutine f90wrap_m_common__get__file_open_counter(f90wrap_file_open_counter)
    use m_common, only: m_common_file_open_counter => file_open_counter
    implicit none
    integer(4), intent(out) :: f90wrap_file_open_counter
    
    f90wrap_file_open_counter = m_common_file_open_counter
end subroutine f90wrap_m_common__get__file_open_counter

subroutine f90wrap_m_common__set__file_open_counter(f90wrap_file_open_counter)
    use m_common, only: m_common_file_open_counter => file_open_counter
    implicit none
    integer(4), intent(in) :: f90wrap_file_open_counter
    
    m_common_file_open_counter = f90wrap_file_open_counter
end subroutine f90wrap_m_common__set__file_open_counter

subroutine f90wrap_m_common__array__file_exist(dummy_this, nd, dtype, dshape, dloc)
    use m_common, only: m_common_file_exist => file_exist
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    dshape(1:1) = shape(m_common_file_exist)
    dloc = loc(m_common_file_exist)
end subroutine f90wrap_m_common__array__file_exist

subroutine f90wrap_m_common__get__buffer(f90wrap_buffer)
    use m_common, only: m_common_buffer => buffer
    implicit none
    character(1028), intent(out) :: f90wrap_buffer
    
    f90wrap_buffer = m_common_buffer
end subroutine f90wrap_m_common__get__buffer

subroutine f90wrap_m_common__set__buffer(f90wrap_buffer)
    use m_common, only: m_common_buffer => buffer
    implicit none
    character(1028), intent(in) :: f90wrap_buffer
    
    m_common_buffer = f90wrap_buffer
end subroutine f90wrap_m_common__set__buffer

subroutine f90wrap_m_common__get__logic_test(f90wrap_logic_test)
    use m_common, only: m_common_logic_test => logic_test
    implicit none
    logical, intent(out) :: f90wrap_logic_test
    
    f90wrap_logic_test = m_common_logic_test
end subroutine f90wrap_m_common__get__logic_test

subroutine f90wrap_m_common__set__logic_test(f90wrap_logic_test)
    use m_common, only: m_common_logic_test => logic_test
    implicit none
    logical, intent(in) :: f90wrap_logic_test
    
    m_common_logic_test = f90wrap_logic_test
end subroutine f90wrap_m_common__set__logic_test

subroutine f90wrap_m_common__array__norm_inf(dummy_this, nd, dtype, dshape, dloc)
    use m_common, only: m_common_norm_inf => norm_inf
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(m_common_norm_inf)
    dloc = loc(m_common_norm_inf)
end subroutine f90wrap_m_common__array__norm_inf

subroutine f90wrap_m_common__array__norm_L1(dummy_this, nd, dtype, dshape, dloc)
    use m_common, only: m_common_norm_l1 => norm_l1
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(m_common_norm_L1)
    dloc = loc(m_common_norm_L1)
end subroutine f90wrap_m_common__array__norm_L1

subroutine f90wrap_m_common__array__norm_L2(dummy_this, nd, dtype, dshape, dloc)
    use m_common, only: m_common_norm_l2 => norm_l2
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(m_common_norm_L2)
    dloc = loc(m_common_norm_L2)
end subroutine f90wrap_m_common__array__norm_L2

subroutine f90wrap_m_common__get__verbose(f90wrap_verbose)
    use m_common, only: m_common_verbose => verbose
    implicit none
    integer(4), intent(out) :: f90wrap_verbose
    
    f90wrap_verbose = m_common_verbose
end subroutine f90wrap_m_common__get__verbose

subroutine f90wrap_m_common__set__verbose(f90wrap_verbose)
    use m_common, only: m_common_verbose => verbose
    implicit none
    integer(4), intent(in) :: f90wrap_verbose
    
    m_common_verbose = f90wrap_verbose
end subroutine f90wrap_m_common__set__verbose

subroutine f90wrap_m_common__get__restart_min(f90wrap_restart_min)
    use m_common, only: m_common_restart_min => restart_min
    implicit none
    integer(4), intent(out) :: f90wrap_restart_min
    
    f90wrap_restart_min = m_common_restart_min
end subroutine f90wrap_m_common__get__restart_min

subroutine f90wrap_m_common__set__restart_min(f90wrap_restart_min)
    use m_common, only: m_common_restart_min => restart_min
    implicit none
    integer(4), intent(in) :: f90wrap_restart_min
    
    m_common_restart_min = f90wrap_restart_min
end subroutine f90wrap_m_common__set__restart_min

subroutine f90wrap_m_common__get__eps_min(f90wrap_eps_min)
    use m_common, only: m_common_eps_min => eps_min
    implicit none
    real(8), intent(out) :: f90wrap_eps_min
    
    f90wrap_eps_min = m_common_eps_min
end subroutine f90wrap_m_common__get__eps_min

subroutine f90wrap_m_common__set__eps_min(f90wrap_eps_min)
    use m_common, only: m_common_eps_min => eps_min
    implicit none
    real(8), intent(in) :: f90wrap_eps_min
    
    m_common_eps_min = f90wrap_eps_min
end subroutine f90wrap_m_common__set__eps_min

subroutine f90wrap_m_common__get__zerom(f90wrap_zerom)
    use m_common, only: m_common_zerom => zerom
    implicit none
    real(8), intent(out) :: f90wrap_zerom
    
    f90wrap_zerom = m_common_zerom
end subroutine f90wrap_m_common__get__zerom

subroutine f90wrap_m_common__set__zerom(f90wrap_zerom)
    use m_common, only: m_common_zerom => zerom
    implicit none
    real(8), intent(in) :: f90wrap_zerom
    
    m_common_zerom = f90wrap_zerom
end subroutine f90wrap_m_common__set__zerom

subroutine f90wrap_m_common__get__pinfm(f90wrap_pinfm)
    use m_common, only: m_common_pinfm => pinfm
    implicit none
    real(8), intent(out) :: f90wrap_pinfm
    
    f90wrap_pinfm = m_common_pinfm
end subroutine f90wrap_m_common__get__pinfm

subroutine f90wrap_m_common__set__pinfm(f90wrap_pinfm)
    use m_common, only: m_common_pinfm => pinfm
    implicit none
    real(8), intent(in) :: f90wrap_pinfm
    
    m_common_pinfm = f90wrap_pinfm
end subroutine f90wrap_m_common__set__pinfm

subroutine f90wrap_m_common__get__minfm(f90wrap_minfm)
    use m_common, only: m_common_minfm => minfm
    implicit none
    real(8), intent(out) :: f90wrap_minfm
    
    f90wrap_minfm = m_common_minfm
end subroutine f90wrap_m_common__get__minfm

subroutine f90wrap_m_common__set__minfm(f90wrap_minfm)
    use m_common, only: m_common_minfm => minfm
    implicit none
    real(8), intent(in) :: f90wrap_minfm
    
    m_common_minfm = f90wrap_minfm
end subroutine f90wrap_m_common__set__minfm

subroutine f90wrap_m_common__get__hugem(f90wrap_hugem)
    use m_common, only: m_common_hugem => hugem
    implicit none
    real(8), intent(out) :: f90wrap_hugem
    
    f90wrap_hugem = m_common_hugem
end subroutine f90wrap_m_common__get__hugem

subroutine f90wrap_m_common__set__hugem(f90wrap_hugem)
    use m_common, only: m_common_hugem => hugem
    implicit none
    real(8), intent(in) :: f90wrap_hugem
    
    m_common_hugem = f90wrap_hugem
end subroutine f90wrap_m_common__set__hugem

subroutine f90wrap_m_common__get__tinym(f90wrap_tinym)
    use m_common, only: m_common_tinym => tinym
    implicit none
    real(8), intent(out) :: f90wrap_tinym
    
    f90wrap_tinym = m_common_tinym
end subroutine f90wrap_m_common__get__tinym

subroutine f90wrap_m_common__set__tinym(f90wrap_tinym)
    use m_common, only: m_common_tinym => tinym
    implicit none
    real(8), intent(in) :: f90wrap_tinym
    
    m_common_tinym = f90wrap_tinym
end subroutine f90wrap_m_common__set__tinym

subroutine f90wrap_m_common__get__zero(f90wrap_zero)
    use m_common, only: m_common_zero => zero
    implicit none
    real(8), intent(out) :: f90wrap_zero
    
    f90wrap_zero = m_common_zero
end subroutine f90wrap_m_common__get__zero

subroutine f90wrap_m_common__get__one(f90wrap_one)
    use m_common, only: m_common_one => one
    implicit none
    real(8), intent(out) :: f90wrap_one
    
    f90wrap_one = m_common_one
end subroutine f90wrap_m_common__get__one

subroutine f90wrap_m_common__get__two(f90wrap_two)
    use m_common, only: m_common_two => two
    implicit none
    real(8), intent(out) :: f90wrap_two
    
    f90wrap_two = m_common_two
end subroutine f90wrap_m_common__get__two

subroutine f90wrap_m_common__get__demi(f90wrap_demi)
    use m_common, only: m_common_demi => demi
    implicit none
    real(8), intent(out) :: f90wrap_demi
    
    f90wrap_demi = m_common_demi
end subroutine f90wrap_m_common__get__demi

subroutine f90wrap_m_common__get__d1p4(f90wrap_d1p4)
    use m_common, only: m_common_d1p4 => d1p4
    implicit none
    real(8), intent(out) :: f90wrap_d1p4
    
    f90wrap_d1p4 = m_common_d1p4
end subroutine f90wrap_m_common__get__d1p4

subroutine f90wrap_m_common__get__d1p3(f90wrap_d1p3)
    use m_common, only: m_common_d1p3 => d1p3
    implicit none
    real(8), intent(out) :: f90wrap_d1p3
    
    f90wrap_d1p3 = m_common_d1p3
end subroutine f90wrap_m_common__get__d1p3

subroutine f90wrap_m_common__get__d2p3(f90wrap_d2p3)
    use m_common, only: m_common_d2p3 => d2p3
    implicit none
    real(8), intent(out) :: f90wrap_d2p3
    
    f90wrap_d2p3 = m_common_d2p3
end subroutine f90wrap_m_common__get__d2p3

subroutine f90wrap_m_common__get__d4p3(f90wrap_d4p3)
    use m_common, only: m_common_d4p3 => d4p3
    implicit none
    real(8), intent(out) :: f90wrap_d4p3
    
    f90wrap_d4p3 = m_common_d4p3
end subroutine f90wrap_m_common__get__d4p3

subroutine f90wrap_m_common__get__d5p3(f90wrap_d5p3)
    use m_common, only: m_common_d5p3 => d5p3
    implicit none
    real(8), intent(out) :: f90wrap_d5p3
    
    f90wrap_d5p3 = m_common_d5p3
end subroutine f90wrap_m_common__get__d5p3

subroutine f90wrap_m_common__get__d7p3(f90wrap_d7p3)
    use m_common, only: m_common_d7p3 => d7p3
    implicit none
    real(8), intent(out) :: f90wrap_d7p3
    
    f90wrap_d7p3 = m_common_d7p3
end subroutine f90wrap_m_common__get__d7p3

subroutine f90wrap_m_common__get__d8p3(f90wrap_d8p3)
    use m_common, only: m_common_d8p3 => d8p3
    implicit none
    real(8), intent(out) :: f90wrap_d8p3
    
    f90wrap_d8p3 = m_common_d8p3
end subroutine f90wrap_m_common__get__d8p3

subroutine f90wrap_m_common__get__d10p3(f90wrap_d10p3)
    use m_common, only: m_common_d10p3 => d10p3
    implicit none
    real(8), intent(out) :: f90wrap_d10p3
    
    f90wrap_d10p3 = m_common_d10p3
end subroutine f90wrap_m_common__get__d10p3

subroutine f90wrap_m_common__get__d3p2(f90wrap_d3p2)
    use m_common, only: m_common_d3p2 => d3p2
    implicit none
    real(8), intent(out) :: f90wrap_d3p2
    
    f90wrap_d3p2 = m_common_d3p2
end subroutine f90wrap_m_common__get__d3p2

subroutine f90wrap_m_common__get__d3p5(f90wrap_d3p5)
    use m_common, only: m_common_d3p5 => d3p5
    implicit none
    real(8), intent(out) :: f90wrap_d3p5
    
    f90wrap_d3p5 = m_common_d3p5
end subroutine f90wrap_m_common__get__d3p5

subroutine f90wrap_m_common__get__d3p8(f90wrap_d3p8)
    use m_common, only: m_common_d3p8 => d3p8
    implicit none
    real(8), intent(out) :: f90wrap_d3p8
    
    f90wrap_d3p8 = m_common_d3p8
end subroutine f90wrap_m_common__get__d3p8

subroutine f90wrap_m_common__get__pi(f90wrap_pi)
    use m_common, only: m_common_pi => pi
    implicit none
    real(8), intent(out) :: f90wrap_pi
    
    f90wrap_pi = m_common_pi
end subroutine f90wrap_m_common__get__pi

subroutine f90wrap_m_common__get__proc(f90wrap_proc)
    use m_common, only: m_common_proc => proc
    implicit none
    real(8), intent(out) :: f90wrap_proc
    
    f90wrap_proc = m_common_proc
end subroutine f90wrap_m_common__get__proc

subroutine f90wrap_m_common__get__np(f90wrap_np)
    use m_common, only: m_common_np => np
    implicit none
    integer(4), intent(out) :: f90wrap_np
    
    f90wrap_np = m_common_np
end subroutine f90wrap_m_common__get__np

! End of module m_common defined in file m_common.f90

