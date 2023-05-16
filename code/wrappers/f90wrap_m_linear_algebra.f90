! Module m_linear_algebra defined in file m_linear_algebra.f90

subroutine f90wrap_vec2d__get__x(this, f90wrap_x)
    use m_linear_algebra, only: vec2d
    implicit none
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec2d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_x
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_x = this_ptr%p%x
end subroutine f90wrap_vec2d__get__x

subroutine f90wrap_vec2d__set__x(this, f90wrap_x)
    use m_linear_algebra, only: vec2d
    implicit none
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec2d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_x
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%x = f90wrap_x
end subroutine f90wrap_vec2d__set__x

subroutine f90wrap_vec2d__get__y(this, f90wrap_y)
    use m_linear_algebra, only: vec2d
    implicit none
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec2d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_y
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_y = this_ptr%p%y
end subroutine f90wrap_vec2d__get__y

subroutine f90wrap_vec2d__set__y(this, f90wrap_y)
    use m_linear_algebra, only: vec2d
    implicit none
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec2d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_y
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%y = f90wrap_y
end subroutine f90wrap_vec2d__set__y

subroutine f90wrap_vec2d_initialise(this)
    use m_linear_algebra, only: vec2d
    implicit none
    
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    type(vec2d_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_vec2d_initialise

subroutine f90wrap_vec2d_finalise(this)
    use m_linear_algebra, only: vec2d
    implicit none
    
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    type(vec2d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_vec2d_finalise

subroutine f90wrap_vec3d__get__x(this, f90wrap_x)
    use m_linear_algebra, only: vec3d
    implicit none
    type vec3d_ptr_type
        type(vec3d), pointer :: p => NULL()
    end type vec3d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec3d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_x
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_x = this_ptr%p%x
end subroutine f90wrap_vec3d__get__x

subroutine f90wrap_vec3d__set__x(this, f90wrap_x)
    use m_linear_algebra, only: vec3d
    implicit none
    type vec3d_ptr_type
        type(vec3d), pointer :: p => NULL()
    end type vec3d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec3d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_x
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%x = f90wrap_x
end subroutine f90wrap_vec3d__set__x

subroutine f90wrap_vec3d__get__y(this, f90wrap_y)
    use m_linear_algebra, only: vec3d
    implicit none
    type vec3d_ptr_type
        type(vec3d), pointer :: p => NULL()
    end type vec3d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec3d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_y
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_y = this_ptr%p%y
end subroutine f90wrap_vec3d__get__y

subroutine f90wrap_vec3d__set__y(this, f90wrap_y)
    use m_linear_algebra, only: vec3d
    implicit none
    type vec3d_ptr_type
        type(vec3d), pointer :: p => NULL()
    end type vec3d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec3d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_y
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%y = f90wrap_y
end subroutine f90wrap_vec3d__set__y

subroutine f90wrap_vec3d__get__z(this, f90wrap_z)
    use m_linear_algebra, only: vec3d
    implicit none
    type vec3d_ptr_type
        type(vec3d), pointer :: p => NULL()
    end type vec3d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec3d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_z
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_z = this_ptr%p%z
end subroutine f90wrap_vec3d__get__z

subroutine f90wrap_vec3d__set__z(this, f90wrap_z)
    use m_linear_algebra, only: vec3d
    implicit none
    type vec3d_ptr_type
        type(vec3d), pointer :: p => NULL()
    end type vec3d_ptr_type
    integer, intent(in)   :: this(2)
    type(vec3d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_z
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%z = f90wrap_z
end subroutine f90wrap_vec3d__set__z

subroutine f90wrap_vec3d_initialise(this)
    use m_linear_algebra, only: vec3d
    implicit none
    
    type vec3d_ptr_type
        type(vec3d), pointer :: p => NULL()
    end type vec3d_ptr_type
    type(vec3d_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_vec3d_initialise

subroutine f90wrap_vec3d_finalise(this)
    use m_linear_algebra, only: vec3d
    implicit none
    
    type vec3d_ptr_type
        type(vec3d), pointer :: p => NULL()
    end type vec3d_ptr_type
    type(vec3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_vec3d_finalise

subroutine f90wrap_tens2d__get__xx(this, f90wrap_xx)
    use m_linear_algebra, only: tens2d
    implicit none
    type tens2d_ptr_type
        type(tens2d), pointer :: p => NULL()
    end type tens2d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens2d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_xx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_xx = this_ptr%p%xx
end subroutine f90wrap_tens2d__get__xx

subroutine f90wrap_tens2d__set__xx(this, f90wrap_xx)
    use m_linear_algebra, only: tens2d
    implicit none
    type tens2d_ptr_type
        type(tens2d), pointer :: p => NULL()
    end type tens2d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens2d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_xx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%xx = f90wrap_xx
end subroutine f90wrap_tens2d__set__xx

subroutine f90wrap_tens2d__get__xy(this, f90wrap_xy)
    use m_linear_algebra, only: tens2d
    implicit none
    type tens2d_ptr_type
        type(tens2d), pointer :: p => NULL()
    end type tens2d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens2d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_xy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_xy = this_ptr%p%xy
end subroutine f90wrap_tens2d__get__xy

subroutine f90wrap_tens2d__set__xy(this, f90wrap_xy)
    use m_linear_algebra, only: tens2d
    implicit none
    type tens2d_ptr_type
        type(tens2d), pointer :: p => NULL()
    end type tens2d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens2d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_xy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%xy = f90wrap_xy
end subroutine f90wrap_tens2d__set__xy

subroutine f90wrap_tens2d__get__yx(this, f90wrap_yx)
    use m_linear_algebra, only: tens2d
    implicit none
    type tens2d_ptr_type
        type(tens2d), pointer :: p => NULL()
    end type tens2d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens2d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_yx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_yx = this_ptr%p%yx
end subroutine f90wrap_tens2d__get__yx

subroutine f90wrap_tens2d__set__yx(this, f90wrap_yx)
    use m_linear_algebra, only: tens2d
    implicit none
    type tens2d_ptr_type
        type(tens2d), pointer :: p => NULL()
    end type tens2d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens2d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_yx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%yx = f90wrap_yx
end subroutine f90wrap_tens2d__set__yx

subroutine f90wrap_tens2d__get__yy(this, f90wrap_yy)
    use m_linear_algebra, only: tens2d
    implicit none
    type tens2d_ptr_type
        type(tens2d), pointer :: p => NULL()
    end type tens2d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens2d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_yy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_yy = this_ptr%p%yy
end subroutine f90wrap_tens2d__get__yy

subroutine f90wrap_tens2d__set__yy(this, f90wrap_yy)
    use m_linear_algebra, only: tens2d
    implicit none
    type tens2d_ptr_type
        type(tens2d), pointer :: p => NULL()
    end type tens2d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens2d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_yy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%yy = f90wrap_yy
end subroutine f90wrap_tens2d__set__yy

subroutine f90wrap_tens2d_initialise(this)
    use m_linear_algebra, only: tens2d
    implicit none
    
    type tens2d_ptr_type
        type(tens2d), pointer :: p => NULL()
    end type tens2d_ptr_type
    type(tens2d_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_tens2d_initialise

subroutine f90wrap_tens2d_finalise(this)
    use m_linear_algebra, only: tens2d
    implicit none
    
    type tens2d_ptr_type
        type(tens2d), pointer :: p => NULL()
    end type tens2d_ptr_type
    type(tens2d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_tens2d_finalise

subroutine f90wrap_tens3d__get__xx(this, f90wrap_xx)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_xx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_xx = this_ptr%p%xx
end subroutine f90wrap_tens3d__get__xx

subroutine f90wrap_tens3d__set__xx(this, f90wrap_xx)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_xx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%xx = f90wrap_xx
end subroutine f90wrap_tens3d__set__xx

subroutine f90wrap_tens3d__get__xy(this, f90wrap_xy)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_xy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_xy = this_ptr%p%xy
end subroutine f90wrap_tens3d__get__xy

subroutine f90wrap_tens3d__set__xy(this, f90wrap_xy)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_xy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%xy = f90wrap_xy
end subroutine f90wrap_tens3d__set__xy

subroutine f90wrap_tens3d__get__xz(this, f90wrap_xz)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_xz
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_xz = this_ptr%p%xz
end subroutine f90wrap_tens3d__get__xz

subroutine f90wrap_tens3d__set__xz(this, f90wrap_xz)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_xz
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%xz = f90wrap_xz
end subroutine f90wrap_tens3d__set__xz

subroutine f90wrap_tens3d__get__yx(this, f90wrap_yx)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_yx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_yx = this_ptr%p%yx
end subroutine f90wrap_tens3d__get__yx

subroutine f90wrap_tens3d__set__yx(this, f90wrap_yx)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_yx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%yx = f90wrap_yx
end subroutine f90wrap_tens3d__set__yx

subroutine f90wrap_tens3d__get__yy(this, f90wrap_yy)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_yy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_yy = this_ptr%p%yy
end subroutine f90wrap_tens3d__get__yy

subroutine f90wrap_tens3d__set__yy(this, f90wrap_yy)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_yy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%yy = f90wrap_yy
end subroutine f90wrap_tens3d__set__yy

subroutine f90wrap_tens3d__get__yz(this, f90wrap_yz)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_yz
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_yz = this_ptr%p%yz
end subroutine f90wrap_tens3d__get__yz

subroutine f90wrap_tens3d__set__yz(this, f90wrap_yz)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_yz
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%yz = f90wrap_yz
end subroutine f90wrap_tens3d__set__yz

subroutine f90wrap_tens3d__get__zx(this, f90wrap_zx)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_zx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_zx = this_ptr%p%zx
end subroutine f90wrap_tens3d__get__zx

subroutine f90wrap_tens3d__set__zx(this, f90wrap_zx)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_zx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%zx = f90wrap_zx
end subroutine f90wrap_tens3d__set__zx

subroutine f90wrap_tens3d__get__zy(this, f90wrap_zy)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_zy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_zy = this_ptr%p%zy
end subroutine f90wrap_tens3d__get__zy

subroutine f90wrap_tens3d__set__zy(this, f90wrap_zy)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_zy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%zy = f90wrap_zy
end subroutine f90wrap_tens3d__set__zy

subroutine f90wrap_tens3d__get__zz(this, f90wrap_zz)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_zz
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_zz = this_ptr%p%zz
end subroutine f90wrap_tens3d__get__zz

subroutine f90wrap_tens3d__set__zz(this, f90wrap_zz)
    use m_linear_algebra, only: tens3d
    implicit none
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    integer, intent(in)   :: this(2)
    type(tens3d_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_zz
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%zz = f90wrap_zz
end subroutine f90wrap_tens3d__set__zz

subroutine f90wrap_tens3d_initialise(this)
    use m_linear_algebra, only: tens3d
    implicit none
    
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    type(tens3d_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_tens3d_initialise

subroutine f90wrap_tens3d_finalise(this)
    use m_linear_algebra, only: tens3d
    implicit none
    
    type tens3d_ptr_type
        type(tens3d), pointer :: p => NULL()
    end type tens3d_ptr_type
    type(tens3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_tens3d_finalise

subroutine f90wrap_sys_lin_full__get__n(this, f90wrap_n)
    use m_linear_algebra, only: sys_lin_full
    implicit none
    type sys_lin_full_ptr_type
        type(sys_lin_full), pointer :: p => NULL()
    end type sys_lin_full_ptr_type
    integer, intent(in)   :: this(2)
    type(sys_lin_full_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_n
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_n = this_ptr%p%n
end subroutine f90wrap_sys_lin_full__get__n

subroutine f90wrap_sys_lin_full__set__n(this, f90wrap_n)
    use m_linear_algebra, only: sys_lin_full
    implicit none
    type sys_lin_full_ptr_type
        type(sys_lin_full), pointer :: p => NULL()
    end type sys_lin_full_ptr_type
    integer, intent(in)   :: this(2)
    type(sys_lin_full_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_n
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%n = f90wrap_n
end subroutine f90wrap_sys_lin_full__set__n

subroutine f90wrap_sys_lin_full__array__a(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: sys_lin_full
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type sys_lin_full_ptr_type
        type(sys_lin_full), pointer :: p => NULL()
    end type sys_lin_full_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(sys_lin_full_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%a)) then
        dshape(1:2) = shape(this_ptr%p%a)
        dloc = loc(this_ptr%p%a)
    else
        dloc = 0
    end if
end subroutine f90wrap_sys_lin_full__array__a

subroutine f90wrap_sys_lin_full__array__rhs(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: sys_lin_full
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type sys_lin_full_ptr_type
        type(sys_lin_full), pointer :: p => NULL()
    end type sys_lin_full_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(sys_lin_full_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%rhs)) then
        dshape(1:1) = shape(this_ptr%p%rhs)
        dloc = loc(this_ptr%p%rhs)
    else
        dloc = 0
    end if
end subroutine f90wrap_sys_lin_full__array__rhs

subroutine f90wrap_sys_lin_full__get__det(this, f90wrap_det)
    use m_linear_algebra, only: sys_lin_full
    implicit none
    type sys_lin_full_ptr_type
        type(sys_lin_full), pointer :: p => NULL()
    end type sys_lin_full_ptr_type
    integer, intent(in)   :: this(2)
    type(sys_lin_full_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_det
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_det = this_ptr%p%det
end subroutine f90wrap_sys_lin_full__get__det

subroutine f90wrap_sys_lin_full__set__det(this, f90wrap_det)
    use m_linear_algebra, only: sys_lin_full
    implicit none
    type sys_lin_full_ptr_type
        type(sys_lin_full), pointer :: p => NULL()
    end type sys_lin_full_ptr_type
    integer, intent(in)   :: this(2)
    type(sys_lin_full_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_det
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%det = f90wrap_det
end subroutine f90wrap_sys_lin_full__set__det

subroutine f90wrap_sys_lin_full_initialise(this)
    use m_linear_algebra, only: sys_lin_full
    implicit none
    
    type sys_lin_full_ptr_type
        type(sys_lin_full), pointer :: p => NULL()
    end type sys_lin_full_ptr_type
    type(sys_lin_full_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_sys_lin_full_initialise

subroutine f90wrap_sys_lin_full_finalise(this)
    use m_linear_algebra, only: sys_lin_full
    implicit none
    
    type sys_lin_full_ptr_type
        type(sys_lin_full), pointer :: p => NULL()
    end type sys_lin_full_ptr_type
    type(sys_lin_full_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_sys_lin_full_finalise

subroutine f90wrap_sys_lin_sparse__get__n(this, f90wrap_n)
    use m_linear_algebra, only: sys_lin_sparse
    implicit none
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    integer, intent(in)   :: this(2)
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_n
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_n = this_ptr%p%n
end subroutine f90wrap_sys_lin_sparse__get__n

subroutine f90wrap_sys_lin_sparse__set__n(this, f90wrap_n)
    use m_linear_algebra, only: sys_lin_sparse
    implicit none
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    integer, intent(in)   :: this(2)
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_n
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%n = f90wrap_n
end subroutine f90wrap_sys_lin_sparse__set__n

subroutine f90wrap_sys_lin_sparse__get__nz(this, f90wrap_nz)
    use m_linear_algebra, only: sys_lin_sparse
    implicit none
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    integer, intent(in)   :: this(2)
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nz
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nz = this_ptr%p%nz
end subroutine f90wrap_sys_lin_sparse__get__nz

subroutine f90wrap_sys_lin_sparse__set__nz(this, f90wrap_nz)
    use m_linear_algebra, only: sys_lin_sparse
    implicit none
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    integer, intent(in)   :: this(2)
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nz
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nz = f90wrap_nz
end subroutine f90wrap_sys_lin_sparse__set__nz

subroutine f90wrap_sys_lin_sparse__array__ia(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: sys_lin_sparse
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ia)) then
        dshape(1:1) = shape(this_ptr%p%ia)
        dloc = loc(this_ptr%p%ia)
    else
        dloc = 0
    end if
end subroutine f90wrap_sys_lin_sparse__array__ia

subroutine f90wrap_sys_lin_sparse__array__ja(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: sys_lin_sparse
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%ja)) then
        dshape(1:1) = shape(this_ptr%p%ja)
        dloc = loc(this_ptr%p%ja)
    else
        dloc = 0
    end if
end subroutine f90wrap_sys_lin_sparse__array__ja

subroutine f90wrap_sys_lin_sparse__array__swap(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: sys_lin_sparse
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%swap)) then
        dshape(1:1) = shape(this_ptr%p%swap)
        dloc = loc(this_ptr%p%swap)
    else
        dloc = 0
    end if
end subroutine f90wrap_sys_lin_sparse__array__swap

subroutine f90wrap_sys_lin_sparse__array__a(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: sys_lin_sparse
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%a)) then
        dshape(1:1) = shape(this_ptr%p%a)
        dloc = loc(this_ptr%p%a)
    else
        dloc = 0
    end if
end subroutine f90wrap_sys_lin_sparse__array__a

subroutine f90wrap_sys_lin_sparse__array__x(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: sys_lin_sparse
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%x)) then
        dshape(1:1) = shape(this_ptr%p%x)
        dloc = loc(this_ptr%p%x)
    else
        dloc = 0
    end if
end subroutine f90wrap_sys_lin_sparse__array__x

subroutine f90wrap_sys_lin_sparse__array__x0(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: sys_lin_sparse
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%x0)) then
        dshape(1:1) = shape(this_ptr%p%x0)
        dloc = loc(this_ptr%p%x0)
    else
        dloc = 0
    end if
end subroutine f90wrap_sys_lin_sparse__array__x0

subroutine f90wrap_sys_lin_sparse__array__rhs(this, nd, dtype, dshape, dloc)
    use m_linear_algebra, only: sys_lin_sparse
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%rhs)) then
        dshape(1:1) = shape(this_ptr%p%rhs)
        dloc = loc(this_ptr%p%rhs)
    else
        dloc = 0
    end if
end subroutine f90wrap_sys_lin_sparse__array__rhs

subroutine f90wrap_sys_lin_sparse_initialise(this)
    use m_linear_algebra, only: sys_lin_sparse
    implicit none
    
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_sys_lin_sparse_initialise

subroutine f90wrap_sys_lin_sparse_finalise(this)
    use m_linear_algebra, only: sys_lin_sparse
    implicit none
    
    type sys_lin_sparse_ptr_type
        type(sys_lin_sparse), pointer :: p => NULL()
    end type sys_lin_sparse_ptr_type
    type(sys_lin_sparse_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_sys_lin_sparse_finalise

! End of module m_linear_algebra defined in file m_linear_algebra.f90

