! Module m_model defined in file m_sw_mono.f90

subroutine f90wrap_unk__array__h(this, nd, dtype, dshape, dloc)
    use m_model, only: unk
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(unk_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%h)) then
        dshape(1:1) = shape(this_ptr%p%h)
        dloc = loc(this_ptr%p%h)
    else
        dloc = 0
    end if
end subroutine f90wrap_unk__array__h

subroutine f90wrap_unk__array__u(this, nd, dtype, dshape, dloc)
    use m_model, only: unk
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(unk_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%u)) then
        dshape(1:1) = shape(this_ptr%p%u)
        dloc = loc(this_ptr%p%u)
    else
        dloc = 0
    end if
end subroutine f90wrap_unk__array__u

subroutine f90wrap_unk__array__v(this, nd, dtype, dshape, dloc)
    use m_model, only: unk
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(unk_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%v)) then
        dshape(1:1) = shape(this_ptr%p%v)
        dloc = loc(this_ptr%p%v)
    else
        dloc = 0
    end if
end subroutine f90wrap_unk__array__v

subroutine f90wrap_unk__array__infil(this, nd, dtype, dshape, dloc)
    use m_model, only: unk
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(unk_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%infil)) then
        dshape(1:1) = shape(this_ptr%p%infil)
        dloc = loc(this_ptr%p%infil)
    else
        dloc = 0
    end if
end subroutine f90wrap_unk__array__infil

subroutine f90wrap_unk__get__t_display(this, f90wrap_t_display)
    use m_model, only: unk
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(in)   :: this(2)
    type(unk_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_t_display
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_t_display = this_ptr%p%t_display
end subroutine f90wrap_unk__get__t_display

subroutine f90wrap_unk__set__t_display(this, f90wrap_t_display)
    use m_model, only: unk
    implicit none
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    integer, intent(in)   :: this(2)
    type(unk_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_t_display
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%t_display = f90wrap_t_display
end subroutine f90wrap_unk__set__t_display

subroutine f90wrap_unk__array_getitem__grad_h(f90wrap_this, f90wrap_i, grad_hitem)
    
    use m_model, only: unk
    use m_linear_algebra, only: vec2d
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(unk_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: grad_hitem(2)
    type(vec2d_ptr_type) :: grad_h_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%grad_h)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%grad_h)) then
            call f90wrap_abort("array index out of range")
        else
            grad_h_ptr%p => this_ptr%p%grad_h(f90wrap_i)
            grad_hitem = transfer(grad_h_ptr,grad_hitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_unk__array_getitem__grad_h

subroutine f90wrap_unk__array_setitem__grad_h(f90wrap_this, f90wrap_i, grad_hitem)
    
    use m_model, only: unk
    use m_linear_algebra, only: vec2d
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(unk_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: grad_hitem(2)
    type(vec2d_ptr_type) :: grad_h_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%grad_h)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%grad_h)) then
            call f90wrap_abort("array index out of range")
        else
            grad_h_ptr = transfer(grad_hitem,grad_h_ptr)
            this_ptr%p%grad_h(f90wrap_i) = grad_h_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_unk__array_setitem__grad_h

subroutine f90wrap_unk__array_len__grad_h(f90wrap_this, f90wrap_n)
    
    use m_model, only: unk
    use m_linear_algebra, only: vec2d
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(unk_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%grad_h)) then
        f90wrap_n = size(this_ptr%p%grad_h)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_unk__array_len__grad_h

subroutine f90wrap_unk__array_getitem__grad_u(f90wrap_this, f90wrap_i, grad_uitem)
    
    use m_model, only: unk
    use m_linear_algebra, only: vec2d
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(unk_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: grad_uitem(2)
    type(vec2d_ptr_type) :: grad_u_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%grad_u)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%grad_u)) then
            call f90wrap_abort("array index out of range")
        else
            grad_u_ptr%p => this_ptr%p%grad_u(f90wrap_i)
            grad_uitem = transfer(grad_u_ptr,grad_uitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_unk__array_getitem__grad_u

subroutine f90wrap_unk__array_setitem__grad_u(f90wrap_this, f90wrap_i, grad_uitem)
    
    use m_model, only: unk
    use m_linear_algebra, only: vec2d
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(unk_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: grad_uitem(2)
    type(vec2d_ptr_type) :: grad_u_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%grad_u)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%grad_u)) then
            call f90wrap_abort("array index out of range")
        else
            grad_u_ptr = transfer(grad_uitem,grad_u_ptr)
            this_ptr%p%grad_u(f90wrap_i) = grad_u_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_unk__array_setitem__grad_u

subroutine f90wrap_unk__array_len__grad_u(f90wrap_this, f90wrap_n)
    
    use m_model, only: unk
    use m_linear_algebra, only: vec2d
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(unk_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%grad_u)) then
        f90wrap_n = size(this_ptr%p%grad_u)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_unk__array_len__grad_u

subroutine f90wrap_unk__array_getitem__grad_v(f90wrap_this, f90wrap_i, grad_vitem)
    
    use m_model, only: unk
    use m_linear_algebra, only: vec2d
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(unk_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: grad_vitem(2)
    type(vec2d_ptr_type) :: grad_v_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%grad_v)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%grad_v)) then
            call f90wrap_abort("array index out of range")
        else
            grad_v_ptr%p => this_ptr%p%grad_v(f90wrap_i)
            grad_vitem = transfer(grad_v_ptr,grad_vitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_unk__array_getitem__grad_v

subroutine f90wrap_unk__array_setitem__grad_v(f90wrap_this, f90wrap_i, grad_vitem)
    
    use m_model, only: unk
    use m_linear_algebra, only: vec2d
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(unk_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: grad_vitem(2)
    type(vec2d_ptr_type) :: grad_v_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%grad_v)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%grad_v)) then
            call f90wrap_abort("array index out of range")
        else
            grad_v_ptr = transfer(grad_vitem,grad_v_ptr)
            this_ptr%p%grad_v(f90wrap_i) = grad_v_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_unk__array_setitem__grad_v

subroutine f90wrap_unk__array_len__grad_v(f90wrap_this, f90wrap_n)
    
    use m_model, only: unk
    use m_linear_algebra, only: vec2d
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(unk_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%grad_v)) then
        f90wrap_n = size(this_ptr%p%grad_v)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_unk__array_len__grad_v

subroutine f90wrap_unk__array_getitem__grad_z(f90wrap_this, f90wrap_i, grad_zitem)
    
    use m_model, only: unk
    use m_linear_algebra, only: vec2d
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(unk_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: grad_zitem(2)
    type(vec2d_ptr_type) :: grad_z_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%grad_z)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%grad_z)) then
            call f90wrap_abort("array index out of range")
        else
            grad_z_ptr%p => this_ptr%p%grad_z(f90wrap_i)
            grad_zitem = transfer(grad_z_ptr,grad_zitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_unk__array_getitem__grad_z

subroutine f90wrap_unk__array_setitem__grad_z(f90wrap_this, f90wrap_i, grad_zitem)
    
    use m_model, only: unk
    use m_linear_algebra, only: vec2d
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(unk_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: grad_zitem(2)
    type(vec2d_ptr_type) :: grad_z_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%grad_z)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%grad_z)) then
            call f90wrap_abort("array index out of range")
        else
            grad_z_ptr = transfer(grad_zitem,grad_z_ptr)
            this_ptr%p%grad_z(f90wrap_i) = grad_z_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_unk__array_setitem__grad_z

subroutine f90wrap_unk__array_len__grad_z(f90wrap_this, f90wrap_n)
    
    use m_model, only: unk
    use m_linear_algebra, only: vec2d
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(unk_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%grad_z)) then
        f90wrap_n = size(this_ptr%p%grad_z)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_unk__array_len__grad_z

subroutine f90wrap_unk_initialise(dof, mesh)
    use m_mesh, only: msh
    use m_model, only: unk, unk_initialise
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(unk_ptr_type) :: dof_ptr
    integer, intent(out), dimension(2) :: dof
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    mesh_ptr = transfer(mesh, mesh_ptr)
    allocate(dof_ptr%p)
    call unk_initialise(dof=dof_ptr%p, mesh=mesh_ptr%p)
    dof = transfer(dof_ptr, dof)
end subroutine f90wrap_unk_initialise

subroutine f90wrap_unk_finalise(dof)
    use m_model, only: unk, unk_finalise
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(unk_ptr_type) :: dof_ptr
    integer, intent(in), dimension(2) :: dof
    dof_ptr = transfer(dof, dof_ptr)
    call unk_finalise(dof=dof_ptr%p)
    deallocate(dof_ptr%p)
end subroutine f90wrap_unk_finalise

subroutine f90wrap_greenampt__get__PsiF(this, f90wrap_PsiF)
    use m_model, only: greenampt
    implicit none
    type greenampt_ptr_type
        type(greenampt), pointer :: p => NULL()
    end type greenampt_ptr_type
    integer, intent(in)   :: this(2)
    type(greenampt_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_PsiF
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_PsiF = this_ptr%p%PsiF
end subroutine f90wrap_greenampt__get__PsiF

subroutine f90wrap_greenampt__set__PsiF(this, f90wrap_PsiF)
    use m_model, only: greenampt
    implicit none
    type greenampt_ptr_type
        type(greenampt), pointer :: p => NULL()
    end type greenampt_ptr_type
    integer, intent(in)   :: this(2)
    type(greenampt_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_PsiF
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%PsiF = f90wrap_PsiF
end subroutine f90wrap_greenampt__set__PsiF

subroutine f90wrap_greenampt__get__Ks(this, f90wrap_Ks)
    use m_model, only: greenampt
    implicit none
    type greenampt_ptr_type
        type(greenampt), pointer :: p => NULL()
    end type greenampt_ptr_type
    integer, intent(in)   :: this(2)
    type(greenampt_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_Ks
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Ks = this_ptr%p%Ks
end subroutine f90wrap_greenampt__get__Ks

subroutine f90wrap_greenampt__set__Ks(this, f90wrap_Ks)
    use m_model, only: greenampt
    implicit none
    type greenampt_ptr_type
        type(greenampt), pointer :: p => NULL()
    end type greenampt_ptr_type
    integer, intent(in)   :: this(2)
    type(greenampt_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_Ks
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Ks = f90wrap_Ks
end subroutine f90wrap_greenampt__set__Ks

subroutine f90wrap_greenampt__get__DeltaTheta(this, f90wrap_DeltaTheta)
    use m_model, only: greenampt
    implicit none
    type greenampt_ptr_type
        type(greenampt), pointer :: p => NULL()
    end type greenampt_ptr_type
    integer, intent(in)   :: this(2)
    type(greenampt_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_DeltaTheta
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_DeltaTheta = this_ptr%p%DeltaTheta
end subroutine f90wrap_greenampt__get__DeltaTheta

subroutine f90wrap_greenampt__set__DeltaTheta(this, f90wrap_DeltaTheta)
    use m_model, only: greenampt
    implicit none
    type greenampt_ptr_type
        type(greenampt), pointer :: p => NULL()
    end type greenampt_ptr_type
    integer, intent(in)   :: this(2)
    type(greenampt_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_DeltaTheta
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%DeltaTheta = f90wrap_DeltaTheta
end subroutine f90wrap_greenampt__set__DeltaTheta

subroutine f90wrap_greenampt_initialise(this)
    use m_model, only: greenampt
    implicit none
    
    type greenampt_ptr_type
        type(greenampt), pointer :: p => NULL()
    end type greenampt_ptr_type
    type(greenampt_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_greenampt_initialise

subroutine f90wrap_greenampt_finalise(this)
    use m_model, only: greenampt
    implicit none
    
    type greenampt_ptr_type
        type(greenampt), pointer :: p => NULL()
    end type greenampt_ptr_type
    type(greenampt_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_greenampt_finalise

subroutine f90wrap_scs_cn__get__lambda(this, f90wrap_lambda)
    use m_model, only: scs_cn
    implicit none
    type scs_cn_ptr_type
        type(scs_cn), pointer :: p => NULL()
    end type scs_cn_ptr_type
    integer, intent(in)   :: this(2)
    type(scs_cn_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_lambda
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_lambda = this_ptr%p%lambda
end subroutine f90wrap_scs_cn__get__lambda

subroutine f90wrap_scs_cn__set__lambda(this, f90wrap_lambda)
    use m_model, only: scs_cn
    implicit none
    type scs_cn_ptr_type
        type(scs_cn), pointer :: p => NULL()
    end type scs_cn_ptr_type
    integer, intent(in)   :: this(2)
    type(scs_cn_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_lambda
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lambda = f90wrap_lambda
end subroutine f90wrap_scs_cn__set__lambda

subroutine f90wrap_scs_cn__get__CN(this, f90wrap_CN)
    use m_model, only: scs_cn
    implicit none
    type scs_cn_ptr_type
        type(scs_cn), pointer :: p => NULL()
    end type scs_cn_ptr_type
    integer, intent(in)   :: this(2)
    type(scs_cn_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_CN
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_CN = this_ptr%p%CN
end subroutine f90wrap_scs_cn__get__CN

subroutine f90wrap_scs_cn__set__CN(this, f90wrap_CN)
    use m_model, only: scs_cn
    implicit none
    type scs_cn_ptr_type
        type(scs_cn), pointer :: p => NULL()
    end type scs_cn_ptr_type
    integer, intent(in)   :: this(2)
    type(scs_cn_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_CN
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%CN = f90wrap_CN
end subroutine f90wrap_scs_cn__set__CN

subroutine f90wrap_scs_cn_initialise(this)
    use m_model, only: scs_cn
    implicit none
    
    type scs_cn_ptr_type
        type(scs_cn), pointer :: p => NULL()
    end type scs_cn_ptr_type
    type(scs_cn_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_scs_cn_initialise

subroutine f90wrap_scs_cn_finalise(this)
    use m_model, only: scs_cn
    implicit none
    
    type scs_cn_ptr_type
        type(scs_cn), pointer :: p => NULL()
    end type scs_cn_ptr_type
    type(scs_cn_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_scs_cn_finalise

subroutine f90wrap_infiltration_data__get__nland(this, f90wrap_nland)
    use m_model, only: infiltration_data
    implicit none
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer, intent(in)   :: this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nland
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nland = this_ptr%p%nland
end subroutine f90wrap_infiltration_data__get__nland

subroutine f90wrap_infiltration_data__set__nland(this, f90wrap_nland)
    use m_model, only: infiltration_data
    implicit none
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer, intent(in)   :: this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nland
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nland = f90wrap_nland
end subroutine f90wrap_infiltration_data__set__nland

subroutine f90wrap_infiltration_data__array__land(this, nd, dtype, dshape, dloc)
    use m_model, only: infiltration_data
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%land)) then
        dshape(1:1) = shape(this_ptr%p%land)
        dloc = loc(this_ptr%p%land)
    else
        dloc = 0
    end if
end subroutine f90wrap_infiltration_data__array__land

subroutine f90wrap_infiltration_data__array__infil_qty(this, nd, dtype, dshape, dloc)
    use m_model, only: infiltration_data
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%infil_qty)) then
        dshape(1:1) = shape(this_ptr%p%infil_qty)
        dloc = loc(this_ptr%p%infil_qty)
    else
        dloc = 0
    end if
end subroutine f90wrap_infiltration_data__array__infil_qty

subroutine f90wrap_infiltration_data__array_getitem__ga(f90wrap_this, f90wrap_i, gaitem)
    
    use m_model, only: greenampt, infiltration_data
    implicit none
    
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    type greenampt_ptr_type
        type(greenampt), pointer :: p => NULL()
    end type greenampt_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: gaitem(2)
    type(greenampt_ptr_type) :: ga_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%ga)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%ga)) then
            call f90wrap_abort("array index out of range")
        else
            ga_ptr%p => this_ptr%p%ga(f90wrap_i)
            gaitem = transfer(ga_ptr,gaitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_infiltration_data__array_getitem__ga

subroutine f90wrap_infiltration_data__array_setitem__ga(f90wrap_this, f90wrap_i, gaitem)
    
    use m_model, only: greenampt, infiltration_data
    implicit none
    
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    type greenampt_ptr_type
        type(greenampt), pointer :: p => NULL()
    end type greenampt_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: gaitem(2)
    type(greenampt_ptr_type) :: ga_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%ga)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%ga)) then
            call f90wrap_abort("array index out of range")
        else
            ga_ptr = transfer(gaitem,ga_ptr)
            this_ptr%p%ga(f90wrap_i) = ga_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_infiltration_data__array_setitem__ga

subroutine f90wrap_infiltration_data__array_len__ga(f90wrap_this, f90wrap_n)
    
    use m_model, only: greenampt, infiltration_data
    implicit none
    
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    type greenampt_ptr_type
        type(greenampt), pointer :: p => NULL()
    end type greenampt_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%ga)) then
        f90wrap_n = size(this_ptr%p%ga)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_infiltration_data__array_len__ga

subroutine f90wrap_infiltration_data__array_getitem__scs(f90wrap_this, f90wrap_i, scsitem)
    
    use m_model, only: infiltration_data, scs_cn
    implicit none
    
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    type scs_cn_ptr_type
        type(scs_cn), pointer :: p => NULL()
    end type scs_cn_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: scsitem(2)
    type(scs_cn_ptr_type) :: scs_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%scs)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%scs)) then
            call f90wrap_abort("array index out of range")
        else
            scs_ptr%p => this_ptr%p%scs(f90wrap_i)
            scsitem = transfer(scs_ptr,scsitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_infiltration_data__array_getitem__scs

subroutine f90wrap_infiltration_data__array_setitem__scs(f90wrap_this, f90wrap_i, scsitem)
    
    use m_model, only: infiltration_data, scs_cn
    implicit none
    
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    type scs_cn_ptr_type
        type(scs_cn), pointer :: p => NULL()
    end type scs_cn_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: scsitem(2)
    type(scs_cn_ptr_type) :: scs_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%scs)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%scs)) then
            call f90wrap_abort("array index out of range")
        else
            scs_ptr = transfer(scsitem,scs_ptr)
            this_ptr%p%scs(f90wrap_i) = scs_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_infiltration_data__array_setitem__scs

subroutine f90wrap_infiltration_data__array_len__scs(f90wrap_this, f90wrap_n)
    
    use m_model, only: infiltration_data, scs_cn
    implicit none
    
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    type scs_cn_ptr_type
        type(scs_cn), pointer :: p => NULL()
    end type scs_cn_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%scs)) then
        f90wrap_n = size(this_ptr%p%scs)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_infiltration_data__array_len__scs

subroutine f90wrap_infiltration_data__array__x_min(this, nd, dtype, dshape, dloc)
    use m_model, only: infiltration_data
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%x_min)) then
        dshape(1:1) = shape(this_ptr%p%x_min)
        dloc = loc(this_ptr%p%x_min)
    else
        dloc = 0
    end if
end subroutine f90wrap_infiltration_data__array__x_min

subroutine f90wrap_infiltration_data__array__x_max(this, nd, dtype, dshape, dloc)
    use m_model, only: infiltration_data
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%x_max)) then
        dshape(1:1) = shape(this_ptr%p%x_max)
        dloc = loc(this_ptr%p%x_max)
    else
        dloc = 0
    end if
end subroutine f90wrap_infiltration_data__array__x_max

subroutine f90wrap_infiltration_data__array__y_min(this, nd, dtype, dshape, dloc)
    use m_model, only: infiltration_data
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%y_min)) then
        dshape(1:1) = shape(this_ptr%p%y_min)
        dloc = loc(this_ptr%p%y_min)
    else
        dloc = 0
    end if
end subroutine f90wrap_infiltration_data__array__y_min

subroutine f90wrap_infiltration_data__array__y_max(this, nd, dtype, dshape, dloc)
    use m_model, only: infiltration_data
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(infiltration_data_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%y_max)) then
        dshape(1:1) = shape(this_ptr%p%y_max)
        dloc = loc(this_ptr%p%y_max)
    else
        dloc = 0
    end if
end subroutine f90wrap_infiltration_data__array__y_max

subroutine f90wrap_infiltration_data_initialise(this)
    use m_model, only: infiltration_data
    implicit none
    
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    type(infiltration_data_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_infiltration_data_initialise

subroutine f90wrap_infiltration_data_finalise(this)
    use m_model, only: infiltration_data
    implicit none
    
    type infiltration_data_ptr_type
        type(infiltration_data), pointer :: p => NULL()
    end type infiltration_data_ptr_type
    type(infiltration_data_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_infiltration_data_finalise

subroutine f90wrap_friction_data__get__nland(this, f90wrap_nland)
    use m_model, only: friction_data
    implicit none
    type friction_data_ptr_type
        type(friction_data), pointer :: p => NULL()
    end type friction_data_ptr_type
    integer, intent(in)   :: this(2)
    type(friction_data_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nland
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nland = this_ptr%p%nland
end subroutine f90wrap_friction_data__get__nland

subroutine f90wrap_friction_data__set__nland(this, f90wrap_nland)
    use m_model, only: friction_data
    implicit none
    type friction_data_ptr_type
        type(friction_data), pointer :: p => NULL()
    end type friction_data_ptr_type
    integer, intent(in)   :: this(2)
    type(friction_data_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nland
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nland = f90wrap_nland
end subroutine f90wrap_friction_data__set__nland

subroutine f90wrap_friction_data__array__manning(this, nd, dtype, dshape, dloc)
    use m_model, only: friction_data
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type friction_data_ptr_type
        type(friction_data), pointer :: p => NULL()
    end type friction_data_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(friction_data_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%manning)) then
        dshape(1:1) = shape(this_ptr%p%manning)
        dloc = loc(this_ptr%p%manning)
    else
        dloc = 0
    end if
end subroutine f90wrap_friction_data__array__manning

subroutine f90wrap_friction_data__array__manning_beta(this, nd, dtype, dshape, dloc)
    use m_model, only: friction_data
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type friction_data_ptr_type
        type(friction_data), pointer :: p => NULL()
    end type friction_data_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(friction_data_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%manning_beta)) then
        dshape(1:1) = shape(this_ptr%p%manning_beta)
        dloc = loc(this_ptr%p%manning_beta)
    else
        dloc = 0
    end if
end subroutine f90wrap_friction_data__array__manning_beta

subroutine f90wrap_friction_data__array__land(this, nd, dtype, dshape, dloc)
    use m_model, only: friction_data
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type friction_data_ptr_type
        type(friction_data), pointer :: p => NULL()
    end type friction_data_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(friction_data_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%land)) then
        dshape(1:1) = shape(this_ptr%p%land)
        dloc = loc(this_ptr%p%land)
    else
        dloc = 0
    end if
end subroutine f90wrap_friction_data__array__land

subroutine f90wrap_friction_data_initialise(this)
    use m_model, only: friction_data
    implicit none
    
    type friction_data_ptr_type
        type(friction_data), pointer :: p => NULL()
    end type friction_data_ptr_type
    type(friction_data_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_friction_data_initialise

subroutine f90wrap_friction_data_finalise(this)
    use m_model, only: friction_data
    implicit none
    
    type friction_data_ptr_type
        type(friction_data), pointer :: p => NULL()
    end type friction_data_ptr_type
    type(friction_data_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_friction_data_finalise

subroutine f90wrap_param_model__array__bathy_cell(this, nd, dtype, dshape, dloc)
    use m_model, only: param_model
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type param_model_ptr_type
        type(param_model), pointer :: p => NULL()
    end type param_model_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(param_model_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%bathy_cell)) then
        dshape(1:1) = shape(this_ptr%p%bathy_cell)
        dloc = loc(this_ptr%p%bathy_cell)
    else
        dloc = 0
    end if
end subroutine f90wrap_param_model__array__bathy_cell

subroutine f90wrap_param_model_initialise(this)
    use m_model, only: param_model
    implicit none
    
    type param_model_ptr_type
        type(param_model), pointer :: p => NULL()
    end type param_model_ptr_type
    type(param_model_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_param_model_initialise

subroutine f90wrap_param_model_finalise(this)
    use m_model, only: param_model
    implicit none
    
    type param_model_ptr_type
        type(param_model), pointer :: p => NULL()
    end type param_model_ptr_type
    type(param_model_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_param_model_finalise

subroutine f90wrap_hydrograph__get__group(this, f90wrap_group)
    use m_model, only: hydrograph
    implicit none
    type hydrograph_ptr_type
        type(hydrograph), pointer :: p => NULL()
    end type hydrograph_ptr_type
    integer, intent(in)   :: this(2)
    type(hydrograph_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_group = this_ptr%p%group
end subroutine f90wrap_hydrograph__get__group

subroutine f90wrap_hydrograph__set__group(this, f90wrap_group)
    use m_model, only: hydrograph
    implicit none
    type hydrograph_ptr_type
        type(hydrograph), pointer :: p => NULL()
    end type hydrograph_ptr_type
    integer, intent(in)   :: this(2)
    type(hydrograph_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%group = f90wrap_group
end subroutine f90wrap_hydrograph__set__group

subroutine f90wrap_hydrograph__array__t(this, nd, dtype, dshape, dloc)
    use m_model, only: hydrograph
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type hydrograph_ptr_type
        type(hydrograph), pointer :: p => NULL()
    end type hydrograph_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(hydrograph_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%t)) then
        dshape(1:1) = shape(this_ptr%p%t)
        dloc = loc(this_ptr%p%t)
    else
        dloc = 0
    end if
end subroutine f90wrap_hydrograph__array__t

subroutine f90wrap_hydrograph__array__q(this, nd, dtype, dshape, dloc)
    use m_model, only: hydrograph
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type hydrograph_ptr_type
        type(hydrograph), pointer :: p => NULL()
    end type hydrograph_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(hydrograph_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%q)) then
        dshape(1:1) = shape(this_ptr%p%q)
        dloc = loc(this_ptr%p%q)
    else
        dloc = 0
    end if
end subroutine f90wrap_hydrograph__array__q

subroutine f90wrap_hydrograph_initialise(this)
    use m_model, only: hydrograph
    implicit none
    
    type hydrograph_ptr_type
        type(hydrograph), pointer :: p => NULL()
    end type hydrograph_ptr_type
    type(hydrograph_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_hydrograph_initialise

subroutine f90wrap_hydrograph_finalise(this)
    use m_model, only: hydrograph
    implicit none
    
    type hydrograph_ptr_type
        type(hydrograph), pointer :: p => NULL()
    end type hydrograph_ptr_type
    type(hydrograph_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_hydrograph_finalise

subroutine f90wrap_hpresc__get__group(this, f90wrap_group)
    use m_model, only: hpresc
    implicit none
    type hpresc_ptr_type
        type(hpresc), pointer :: p => NULL()
    end type hpresc_ptr_type
    integer, intent(in)   :: this(2)
    type(hpresc_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_group = this_ptr%p%group
end subroutine f90wrap_hpresc__get__group

subroutine f90wrap_hpresc__set__group(this, f90wrap_group)
    use m_model, only: hpresc
    implicit none
    type hpresc_ptr_type
        type(hpresc), pointer :: p => NULL()
    end type hpresc_ptr_type
    integer, intent(in)   :: this(2)
    type(hpresc_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%group = f90wrap_group
end subroutine f90wrap_hpresc__set__group

subroutine f90wrap_hpresc__array__t(this, nd, dtype, dshape, dloc)
    use m_model, only: hpresc
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type hpresc_ptr_type
        type(hpresc), pointer :: p => NULL()
    end type hpresc_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(hpresc_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%t)) then
        dshape(1:1) = shape(this_ptr%p%t)
        dloc = loc(this_ptr%p%t)
    else
        dloc = 0
    end if
end subroutine f90wrap_hpresc__array__t

subroutine f90wrap_hpresc__array__h(this, nd, dtype, dshape, dloc)
    use m_model, only: hpresc
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type hpresc_ptr_type
        type(hpresc), pointer :: p => NULL()
    end type hpresc_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(hpresc_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%h)) then
        dshape(1:1) = shape(this_ptr%p%h)
        dloc = loc(this_ptr%p%h)
    else
        dloc = 0
    end if
end subroutine f90wrap_hpresc__array__h

subroutine f90wrap_hpresc_initialise(this)
    use m_model, only: hpresc
    implicit none
    
    type hpresc_ptr_type
        type(hpresc), pointer :: p => NULL()
    end type hpresc_ptr_type
    type(hpresc_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_hpresc_initialise

subroutine f90wrap_hpresc_finalise(this)
    use m_model, only: hpresc
    implicit none
    
    type hpresc_ptr_type
        type(hpresc), pointer :: p => NULL()
    end type hpresc_ptr_type
    type(hpresc_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_hpresc_finalise

subroutine f90wrap_zspresc__get__group(this, f90wrap_group)
    use m_model, only: zspresc
    implicit none
    type zspresc_ptr_type
        type(zspresc), pointer :: p => NULL()
    end type zspresc_ptr_type
    integer, intent(in)   :: this(2)
    type(zspresc_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_group = this_ptr%p%group
end subroutine f90wrap_zspresc__get__group

subroutine f90wrap_zspresc__set__group(this, f90wrap_group)
    use m_model, only: zspresc
    implicit none
    type zspresc_ptr_type
        type(zspresc), pointer :: p => NULL()
    end type zspresc_ptr_type
    integer, intent(in)   :: this(2)
    type(zspresc_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%group = f90wrap_group
end subroutine f90wrap_zspresc__set__group

subroutine f90wrap_zspresc__array__t(this, nd, dtype, dshape, dloc)
    use m_model, only: zspresc
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type zspresc_ptr_type
        type(zspresc), pointer :: p => NULL()
    end type zspresc_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(zspresc_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%t)) then
        dshape(1:1) = shape(this_ptr%p%t)
        dloc = loc(this_ptr%p%t)
    else
        dloc = 0
    end if
end subroutine f90wrap_zspresc__array__t

subroutine f90wrap_zspresc__array__z(this, nd, dtype, dshape, dloc)
    use m_model, only: zspresc
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type zspresc_ptr_type
        type(zspresc), pointer :: p => NULL()
    end type zspresc_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(zspresc_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%z)) then
        dshape(1:1) = shape(this_ptr%p%z)
        dloc = loc(this_ptr%p%z)
    else
        dloc = 0
    end if
end subroutine f90wrap_zspresc__array__z

subroutine f90wrap_zspresc_initialise(this)
    use m_model, only: zspresc
    implicit none
    
    type zspresc_ptr_type
        type(zspresc), pointer :: p => NULL()
    end type zspresc_ptr_type
    type(zspresc_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_zspresc_initialise

subroutine f90wrap_zspresc_finalise(this)
    use m_model, only: zspresc
    implicit none
    
    type zspresc_ptr_type
        type(zspresc), pointer :: p => NULL()
    end type zspresc_ptr_type
    type(zspresc_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_zspresc_finalise

subroutine f90wrap_ratcurve__get__group(this, f90wrap_group)
    use m_model, only: ratcurve
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(in)   :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_group = this_ptr%p%group
end subroutine f90wrap_ratcurve__get__group

subroutine f90wrap_ratcurve__set__group(this, f90wrap_group)
    use m_model, only: ratcurve
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(in)   :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%group = f90wrap_group
end subroutine f90wrap_ratcurve__set__group

subroutine f90wrap_ratcurve__array__h(this, nd, dtype, dshape, dloc)
    use m_model, only: ratcurve
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%h)) then
        dshape(1:1) = shape(this_ptr%p%h)
        dloc = loc(this_ptr%p%h)
    else
        dloc = 0
    end if
end subroutine f90wrap_ratcurve__array__h

subroutine f90wrap_ratcurve__array__q(this, nd, dtype, dshape, dloc)
    use m_model, only: ratcurve
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%q)) then
        dshape(1:1) = shape(this_ptr%p%q)
        dloc = loc(this_ptr%p%q)
    else
        dloc = 0
    end if
end subroutine f90wrap_ratcurve__array__q

subroutine f90wrap_ratcurve__get__z_rat_ref(this, f90wrap_z_rat_ref)
    use m_model, only: ratcurve
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(in)   :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_z_rat_ref
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_z_rat_ref = this_ptr%p%z_rat_ref
end subroutine f90wrap_ratcurve__get__z_rat_ref

subroutine f90wrap_ratcurve__set__z_rat_ref(this, f90wrap_z_rat_ref)
    use m_model, only: ratcurve
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(in)   :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_z_rat_ref
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%z_rat_ref = f90wrap_z_rat_ref
end subroutine f90wrap_ratcurve__set__z_rat_ref

subroutine f90wrap_ratcurve__get__zout(this, f90wrap_zout)
    use m_model, only: ratcurve
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(in)   :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_zout
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_zout = this_ptr%p%zout
end subroutine f90wrap_ratcurve__get__zout

subroutine f90wrap_ratcurve__set__zout(this, f90wrap_zout)
    use m_model, only: ratcurve
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(in)   :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_zout
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%zout = f90wrap_zout
end subroutine f90wrap_ratcurve__set__zout

subroutine f90wrap_ratcurve__get__c1(this, f90wrap_c1)
    use m_model, only: ratcurve
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(in)   :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_c1
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_c1 = this_ptr%p%c1
end subroutine f90wrap_ratcurve__get__c1

subroutine f90wrap_ratcurve__set__c1(this, f90wrap_c1)
    use m_model, only: ratcurve
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(in)   :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_c1
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%c1 = f90wrap_c1
end subroutine f90wrap_ratcurve__set__c1

subroutine f90wrap_ratcurve__get__c2(this, f90wrap_c2)
    use m_model, only: ratcurve
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(in)   :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_c2
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_c2 = this_ptr%p%c2
end subroutine f90wrap_ratcurve__get__c2

subroutine f90wrap_ratcurve__set__c2(this, f90wrap_c2)
    use m_model, only: ratcurve
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(in)   :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_c2
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%c2 = f90wrap_c2
end subroutine f90wrap_ratcurve__set__c2

subroutine f90wrap_ratcurve__array__pow(this, nd, dtype, dshape, dloc)
    use m_model, only: ratcurve
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(ratcurve_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%pow)
    dloc = loc(this_ptr%p%pow)
end subroutine f90wrap_ratcurve__array__pow

subroutine f90wrap_ratcurve_initialise(this)
    use m_model, only: ratcurve
    implicit none
    
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    type(ratcurve_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_ratcurve_initialise

subroutine f90wrap_ratcurve_finalise(this)
    use m_model, only: ratcurve
    implicit none
    
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    type(ratcurve_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_ratcurve_finalise

subroutine f90wrap_rain__get__x_min(this, f90wrap_x_min)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_x_min
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_x_min = this_ptr%p%x_min
end subroutine f90wrap_rain__get__x_min

subroutine f90wrap_rain__set__x_min(this, f90wrap_x_min)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_x_min
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%x_min = f90wrap_x_min
end subroutine f90wrap_rain__set__x_min

subroutine f90wrap_rain__get__x_max(this, f90wrap_x_max)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_x_max
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_x_max = this_ptr%p%x_max
end subroutine f90wrap_rain__get__x_max

subroutine f90wrap_rain__set__x_max(this, f90wrap_x_max)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_x_max
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%x_max = f90wrap_x_max
end subroutine f90wrap_rain__set__x_max

subroutine f90wrap_rain__get__y_min(this, f90wrap_y_min)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_y_min
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_y_min = this_ptr%p%y_min
end subroutine f90wrap_rain__get__y_min

subroutine f90wrap_rain__set__y_min(this, f90wrap_y_min)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_y_min
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%y_min = f90wrap_y_min
end subroutine f90wrap_rain__set__y_min

subroutine f90wrap_rain__get__y_max(this, f90wrap_y_max)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_y_max
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_y_max = this_ptr%p%y_max
end subroutine f90wrap_rain__get__y_max

subroutine f90wrap_rain__set__y_max(this, f90wrap_y_max)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_y_max
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%y_max = f90wrap_y_max
end subroutine f90wrap_rain__set__y_max

subroutine f90wrap_rain__get__tile_index(this, f90wrap_tile_index)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_tile_index
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_tile_index = this_ptr%p%tile_index
end subroutine f90wrap_rain__get__tile_index

subroutine f90wrap_rain__set__tile_index(this, f90wrap_tile_index)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_tile_index
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%tile_index = f90wrap_tile_index
end subroutine f90wrap_rain__set__tile_index

subroutine f90wrap_rain__array__t(this, nd, dtype, dshape, dloc)
    use m_model, only: rain
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(rain_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%t)) then
        dshape(1:1) = shape(this_ptr%p%t)
        dloc = loc(this_ptr%p%t)
    else
        dloc = 0
    end if
end subroutine f90wrap_rain__array__t

subroutine f90wrap_rain__array__q(this, nd, dtype, dshape, dloc)
    use m_model, only: rain
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(rain_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%q)) then
        dshape(1:1) = shape(this_ptr%p%q)
        dloc = loc(this_ptr%p%q)
    else
        dloc = 0
    end if
end subroutine f90wrap_rain__array__q

subroutine f90wrap_rain__get__qin(this, f90wrap_qin)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_qin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_qin = this_ptr%p%qin
end subroutine f90wrap_rain__get__qin

subroutine f90wrap_rain__set__qin(this, f90wrap_qin)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_qin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%qin = f90wrap_qin
end subroutine f90wrap_rain__set__qin

subroutine f90wrap_rain__get__cumul(this, f90wrap_cumul)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_cumul
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_cumul = this_ptr%p%cumul
end subroutine f90wrap_rain__get__cumul

subroutine f90wrap_rain__set__cumul(this, f90wrap_cumul)
    use m_model, only: rain
    implicit none
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in)   :: this(2)
    type(rain_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_cumul
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%cumul = f90wrap_cumul
end subroutine f90wrap_rain__set__cumul

subroutine f90wrap_rain_initialise(this)
    use m_model, only: rain
    implicit none
    
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    type(rain_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_rain_initialise

subroutine f90wrap_rain_finalise(this)
    use m_model, only: rain
    implicit none
    
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    type(rain_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_rain_finalise

subroutine f90wrap_bcs__get__nb(this, f90wrap_nb)
    use m_model, only: bcs
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(in)   :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nb
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nb = this_ptr%p%nb
end subroutine f90wrap_bcs__get__nb

subroutine f90wrap_bcs__set__nb(this, f90wrap_nb)
    use m_model, only: bcs
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(in)   :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nb
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nb = f90wrap_nb
end subroutine f90wrap_bcs__set__nb

subroutine f90wrap_bcs__get__nb_in(this, f90wrap_nb_in)
    use m_model, only: bcs
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(in)   :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nb_in
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nb_in = this_ptr%p%nb_in
end subroutine f90wrap_bcs__get__nb_in

subroutine f90wrap_bcs__set__nb_in(this, f90wrap_nb_in)
    use m_model, only: bcs
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(in)   :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nb_in
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nb_in = f90wrap_nb_in
end subroutine f90wrap_bcs__set__nb_in

subroutine f90wrap_bcs__get__nb_out(this, f90wrap_nb_out)
    use m_model, only: bcs
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(in)   :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nb_out
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nb_out = this_ptr%p%nb_out
end subroutine f90wrap_bcs__get__nb_out

subroutine f90wrap_bcs__set__nb_out(this, f90wrap_nb_out)
    use m_model, only: bcs
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(in)   :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nb_out
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nb_out = f90wrap_nb_out
end subroutine f90wrap_bcs__set__nb_out

subroutine f90wrap_bcs__get__nb_rn(this, f90wrap_nb_rn)
    use m_model, only: bcs
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(in)   :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nb_rn
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nb_rn = this_ptr%p%nb_rn
end subroutine f90wrap_bcs__get__nb_rn

subroutine f90wrap_bcs__set__nb_rn(this, f90wrap_nb_rn)
    use m_model, only: bcs
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(in)   :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nb_rn
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nb_rn = f90wrap_nb_rn
end subroutine f90wrap_bcs__set__nb_rn

subroutine f90wrap_bcs__array__typ(this, nd, dtype, dshape, dloc)
    use m_model, only: bcs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%typ)) then
        dshape(1:3) = (/len(this_ptr%p%typ(1,1)), shape(this_ptr%p%typ)/)
        dloc = loc(this_ptr%p%typ)
    else
        dloc = 0
    end if
end subroutine f90wrap_bcs__array__typ

subroutine f90wrap_bcs__array__grpf(this, nd, dtype, dshape, dloc)
    use m_model, only: bcs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%grpf)) then
        dshape(1:1) = shape(this_ptr%p%grpf)
        dloc = loc(this_ptr%p%grpf)
    else
        dloc = 0
    end if
end subroutine f90wrap_bcs__array__grpf

subroutine f90wrap_bcs__array__inflow(this, nd, dtype, dshape, dloc)
    use m_model, only: bcs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%inflow)) then
        dshape(1:1) = shape(this_ptr%p%inflow)
        dloc = loc(this_ptr%p%inflow)
    else
        dloc = 0
    end if
end subroutine f90wrap_bcs__array__inflow

subroutine f90wrap_bcs__array__outflow(this, nd, dtype, dshape, dloc)
    use m_model, only: bcs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%outflow)) then
        dshape(1:1) = shape(this_ptr%p%outflow)
        dloc = loc(this_ptr%p%outflow)
    else
        dloc = 0
    end if
end subroutine f90wrap_bcs__array__outflow

subroutine f90wrap_bcs__array_getitem__hyd(f90wrap_this, f90wrap_i, hyditem)
    
    use m_model, only: hydrograph, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type hydrograph_ptr_type
        type(hydrograph), pointer :: p => NULL()
    end type hydrograph_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: hyditem(2)
    type(hydrograph_ptr_type) :: hyd_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%hyd)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%hyd)) then
            call f90wrap_abort("array index out of range")
        else
            hyd_ptr%p => this_ptr%p%hyd(f90wrap_i)
            hyditem = transfer(hyd_ptr,hyditem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_bcs__array_getitem__hyd

subroutine f90wrap_bcs__array_setitem__hyd(f90wrap_this, f90wrap_i, hyditem)
    
    use m_model, only: hydrograph, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type hydrograph_ptr_type
        type(hydrograph), pointer :: p => NULL()
    end type hydrograph_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: hyditem(2)
    type(hydrograph_ptr_type) :: hyd_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%hyd)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%hyd)) then
            call f90wrap_abort("array index out of range")
        else
            hyd_ptr = transfer(hyditem,hyd_ptr)
            this_ptr%p%hyd(f90wrap_i) = hyd_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_bcs__array_setitem__hyd

subroutine f90wrap_bcs__array_len__hyd(f90wrap_this, f90wrap_n)
    
    use m_model, only: hydrograph, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type hydrograph_ptr_type
        type(hydrograph), pointer :: p => NULL()
    end type hydrograph_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%hyd)) then
        f90wrap_n = size(this_ptr%p%hyd)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_bcs__array_len__hyd

subroutine f90wrap_bcs__array_getitem__rat(f90wrap_this, f90wrap_i, ratitem)
    
    use m_model, only: ratcurve, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: ratitem(2)
    type(ratcurve_ptr_type) :: rat_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%rat)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%rat)) then
            call f90wrap_abort("array index out of range")
        else
            rat_ptr%p => this_ptr%p%rat(f90wrap_i)
            ratitem = transfer(rat_ptr,ratitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_bcs__array_getitem__rat

subroutine f90wrap_bcs__array_setitem__rat(f90wrap_this, f90wrap_i, ratitem)
    
    use m_model, only: ratcurve, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: ratitem(2)
    type(ratcurve_ptr_type) :: rat_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%rat)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%rat)) then
            call f90wrap_abort("array index out of range")
        else
            rat_ptr = transfer(ratitem,rat_ptr)
            this_ptr%p%rat(f90wrap_i) = rat_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_bcs__array_setitem__rat

subroutine f90wrap_bcs__array_len__rat(f90wrap_this, f90wrap_n)
    
    use m_model, only: ratcurve, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type ratcurve_ptr_type
        type(ratcurve), pointer :: p => NULL()
    end type ratcurve_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%rat)) then
        f90wrap_n = size(this_ptr%p%rat)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_bcs__array_len__rat

subroutine f90wrap_bcs__array_getitem__hpresc(f90wrap_this, f90wrap_i, hprescitem)
    
    use m_model, only: hpresc, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type hpresc_ptr_type
        type(hpresc), pointer :: p => NULL()
    end type hpresc_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: hprescitem(2)
    type(hpresc_ptr_type) :: hpresc_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%hpresc)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%hpresc)) then
            call f90wrap_abort("array index out of range")
        else
            hpresc_ptr%p => this_ptr%p%hpresc(f90wrap_i)
            hprescitem = transfer(hpresc_ptr,hprescitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_bcs__array_getitem__hpresc

subroutine f90wrap_bcs__array_setitem__hpresc(f90wrap_this, f90wrap_i, hprescitem)
    
    use m_model, only: hpresc, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type hpresc_ptr_type
        type(hpresc), pointer :: p => NULL()
    end type hpresc_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: hprescitem(2)
    type(hpresc_ptr_type) :: hpresc_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%hpresc)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%hpresc)) then
            call f90wrap_abort("array index out of range")
        else
            hpresc_ptr = transfer(hprescitem,hpresc_ptr)
            this_ptr%p%hpresc(f90wrap_i) = hpresc_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_bcs__array_setitem__hpresc

subroutine f90wrap_bcs__array_len__hpresc(f90wrap_this, f90wrap_n)
    
    use m_model, only: hpresc, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type hpresc_ptr_type
        type(hpresc), pointer :: p => NULL()
    end type hpresc_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%hpresc)) then
        f90wrap_n = size(this_ptr%p%hpresc)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_bcs__array_len__hpresc

subroutine f90wrap_bcs__array_getitem__zspresc(f90wrap_this, f90wrap_i, zsprescitem)
    
    use m_model, only: zspresc, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type zspresc_ptr_type
        type(zspresc), pointer :: p => NULL()
    end type zspresc_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: zsprescitem(2)
    type(zspresc_ptr_type) :: zspresc_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%zspresc)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%zspresc)) then
            call f90wrap_abort("array index out of range")
        else
            zspresc_ptr%p => this_ptr%p%zspresc(f90wrap_i)
            zsprescitem = transfer(zspresc_ptr,zsprescitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_bcs__array_getitem__zspresc

subroutine f90wrap_bcs__array_setitem__zspresc(f90wrap_this, f90wrap_i, zsprescitem)
    
    use m_model, only: zspresc, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type zspresc_ptr_type
        type(zspresc), pointer :: p => NULL()
    end type zspresc_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: zsprescitem(2)
    type(zspresc_ptr_type) :: zspresc_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%zspresc)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%zspresc)) then
            call f90wrap_abort("array index out of range")
        else
            zspresc_ptr = transfer(zsprescitem,zspresc_ptr)
            this_ptr%p%zspresc(f90wrap_i) = zspresc_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_bcs__array_setitem__zspresc

subroutine f90wrap_bcs__array_len__zspresc(f90wrap_this, f90wrap_n)
    
    use m_model, only: zspresc, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type zspresc_ptr_type
        type(zspresc), pointer :: p => NULL()
    end type zspresc_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%zspresc)) then
        f90wrap_n = size(this_ptr%p%zspresc)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_bcs__array_len__zspresc

subroutine f90wrap_bcs__array_getitem__rain(f90wrap_this, f90wrap_i, rainitem)
    
    use m_model, only: rain, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: rainitem(2)
    type(rain_ptr_type) :: rain_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%rain)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%rain)) then
            call f90wrap_abort("array index out of range")
        else
            rain_ptr%p => this_ptr%p%rain(f90wrap_i)
            rainitem = transfer(rain_ptr,rainitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_bcs__array_getitem__rain

subroutine f90wrap_bcs__array_setitem__rain(f90wrap_this, f90wrap_i, rainitem)
    
    use m_model, only: rain, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: rainitem(2)
    type(rain_ptr_type) :: rain_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%rain)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%rain)) then
            call f90wrap_abort("array index out of range")
        else
            rain_ptr = transfer(rainitem,rain_ptr)
            this_ptr%p%rain(f90wrap_i) = rain_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_bcs__array_setitem__rain

subroutine f90wrap_bcs__array_len__rain(f90wrap_this, f90wrap_n)
    
    use m_model, only: rain, bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type rain_ptr_type
        type(rain), pointer :: p => NULL()
    end type rain_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(bcs_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%rain)) then
        f90wrap_n = size(this_ptr%p%rain)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_bcs__array_len__rain

subroutine f90wrap_bcs__array__sum_mass_flux(this, nd, dtype, dshape, dloc)
    use m_model, only: bcs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(bcs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%sum_mass_flux)) then
        dshape(1:1) = shape(this_ptr%p%sum_mass_flux)
        dloc = loc(this_ptr%p%sum_mass_flux)
    else
        dloc = 0
    end if
end subroutine f90wrap_bcs__array__sum_mass_flux

subroutine f90wrap_bcs_initialise(this)
    use m_model, only: bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type(bcs_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_bcs_initialise

subroutine f90wrap_bcs_finalise(this)
    use m_model, only: bcs
    implicit none
    
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    type(bcs_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_bcs_finalise

subroutine f90wrap_station_obs__get__weight(this, f90wrap_weight)
    use m_model, only: station_obs
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_weight
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_weight = this_ptr%p%weight
end subroutine f90wrap_station_obs__get__weight

subroutine f90wrap_station_obs__set__weight(this, f90wrap_weight)
    use m_model, only: station_obs
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_weight
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%weight = f90wrap_weight
end subroutine f90wrap_station_obs__set__weight

subroutine f90wrap_station_obs__get__length(this, f90wrap_length)
    use m_model, only: station_obs
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_length
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_length = this_ptr%p%length
end subroutine f90wrap_station_obs__get__length

subroutine f90wrap_station_obs__set__length(this, f90wrap_length)
    use m_model, only: station_obs
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_length
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%length = f90wrap_length
end subroutine f90wrap_station_obs__set__length

subroutine f90wrap_station_obs__get__dt_offset(this, f90wrap_dt_offset)
    use m_model, only: station_obs
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_dt_offset
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_dt_offset = this_ptr%p%dt_offset
end subroutine f90wrap_station_obs__get__dt_offset

subroutine f90wrap_station_obs__set__dt_offset(this, f90wrap_dt_offset)
    use m_model, only: station_obs
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_dt_offset
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%dt_offset = f90wrap_dt_offset
end subroutine f90wrap_station_obs__set__dt_offset

subroutine f90wrap_station_obs__get__dt(this, f90wrap_dt)
    use m_model, only: station_obs
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_dt
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_dt = this_ptr%p%dt
end subroutine f90wrap_station_obs__get__dt

subroutine f90wrap_station_obs__set__dt(this, f90wrap_dt)
    use m_model, only: station_obs
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_dt
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%dt = f90wrap_dt
end subroutine f90wrap_station_obs__set__dt

subroutine f90wrap_station_obs__array__dt_obs(this, nd, dtype, dshape, dloc)
    use m_model, only: station_obs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%dt_obs)) then
        dshape(1:1) = shape(this_ptr%p%dt_obs)
        dloc = loc(this_ptr%p%dt_obs)
    else
        dloc = 0
    end if
end subroutine f90wrap_station_obs__array__dt_obs

subroutine f90wrap_station_obs__get__ind_t(this, f90wrap_ind_t)
    use m_model, only: station_obs
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ind_t
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ind_t = this_ptr%p%ind_t
end subroutine f90wrap_station_obs__get__ind_t

subroutine f90wrap_station_obs__set__ind_t(this, f90wrap_ind_t)
    use m_model, only: station_obs
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ind_t
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ind_t = f90wrap_ind_t
end subroutine f90wrap_station_obs__set__ind_t

subroutine f90wrap_station_obs__get__nb_dt(this, f90wrap_nb_dt)
    use m_model, only: station_obs
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nb_dt
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nb_dt = this_ptr%p%nb_dt
end subroutine f90wrap_station_obs__get__nb_dt

subroutine f90wrap_station_obs__set__nb_dt(this, f90wrap_nb_dt)
    use m_model, only: station_obs
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nb_dt
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nb_dt = f90wrap_nb_dt
end subroutine f90wrap_station_obs__set__nb_dt

subroutine f90wrap_station_obs__array__t(this, nd, dtype, dshape, dloc)
    use m_model, only: station_obs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%t)) then
        dshape(1:1) = shape(this_ptr%p%t)
        dloc = loc(this_ptr%p%t)
    else
        dloc = 0
    end if
end subroutine f90wrap_station_obs__array__t

subroutine f90wrap_station_obs__array__h(this, nd, dtype, dshape, dloc)
    use m_model, only: station_obs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%h)) then
        dshape(1:1) = shape(this_ptr%p%h)
        dloc = loc(this_ptr%p%h)
    else
        dloc = 0
    end if
end subroutine f90wrap_station_obs__array__h

subroutine f90wrap_station_obs__array__u(this, nd, dtype, dshape, dloc)
    use m_model, only: station_obs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%u)) then
        dshape(1:1) = shape(this_ptr%p%u)
        dloc = loc(this_ptr%p%u)
    else
        dloc = 0
    end if
end subroutine f90wrap_station_obs__array__u

subroutine f90wrap_station_obs__array__v(this, nd, dtype, dshape, dloc)
    use m_model, only: station_obs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%v)) then
        dshape(1:1) = shape(this_ptr%p%v)
        dloc = loc(this_ptr%p%v)
    else
        dloc = 0
    end if
end subroutine f90wrap_station_obs__array__v

subroutine f90wrap_station_obs__array__q(this, nd, dtype, dshape, dloc)
    use m_model, only: station_obs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%q)) then
        dshape(1:1) = shape(this_ptr%p%q)
        dloc = loc(this_ptr%p%q)
    else
        dloc = 0
    end if
end subroutine f90wrap_station_obs__array__q

subroutine f90wrap_station_obs__array__w(this, nd, dtype, dshape, dloc)
    use m_model, only: station_obs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(station_obs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%w)) then
        dshape(1:1) = shape(this_ptr%p%w)
        dloc = loc(this_ptr%p%w)
    else
        dloc = 0
    end if
end subroutine f90wrap_station_obs__array__w

subroutine f90wrap_station_obs_initialise(this)
    use m_model, only: station_obs
    implicit none
    
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    type(station_obs_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_station_obs_initialise

subroutine f90wrap_station_obs_finalise(this)
    use m_model, only: station_obs
    implicit none
    
    type station_obs_ptr_type
        type(station_obs), pointer :: p => NULL()
    end type station_obs_ptr_type
    type(station_obs_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_station_obs_finalise

subroutine f90wrap_section_obs__get__dt(this, f90wrap_dt)
    use m_model, only: section_obs
    implicit none
    type section_obs_ptr_type
        type(section_obs), pointer :: p => NULL()
    end type section_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(section_obs_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_dt
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_dt = this_ptr%p%dt
end subroutine f90wrap_section_obs__get__dt

subroutine f90wrap_section_obs__set__dt(this, f90wrap_dt)
    use m_model, only: section_obs
    implicit none
    type section_obs_ptr_type
        type(section_obs), pointer :: p => NULL()
    end type section_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(section_obs_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_dt
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%dt = f90wrap_dt
end subroutine f90wrap_section_obs__set__dt

subroutine f90wrap_section_obs__get__dx(this, f90wrap_dx)
    use m_model, only: section_obs
    implicit none
    type section_obs_ptr_type
        type(section_obs), pointer :: p => NULL()
    end type section_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(section_obs_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_dx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_dx = this_ptr%p%dx
end subroutine f90wrap_section_obs__get__dx

subroutine f90wrap_section_obs__set__dx(this, f90wrap_dx)
    use m_model, only: section_obs
    implicit none
    type section_obs_ptr_type
        type(section_obs), pointer :: p => NULL()
    end type section_obs_ptr_type
    integer, intent(in)   :: this(2)
    type(section_obs_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_dx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%dx = f90wrap_dx
end subroutine f90wrap_section_obs__set__dx

subroutine f90wrap_section_obs__array__t(this, nd, dtype, dshape, dloc)
    use m_model, only: section_obs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type section_obs_ptr_type
        type(section_obs), pointer :: p => NULL()
    end type section_obs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(section_obs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%t)) then
        dshape(1:1) = shape(this_ptr%p%t)
        dloc = loc(this_ptr%p%t)
    else
        dloc = 0
    end if
end subroutine f90wrap_section_obs__array__t

subroutine f90wrap_section_obs__array__h(this, nd, dtype, dshape, dloc)
    use m_model, only: section_obs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type section_obs_ptr_type
        type(section_obs), pointer :: p => NULL()
    end type section_obs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(section_obs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%h)) then
        dshape(1:1) = shape(this_ptr%p%h)
        dloc = loc(this_ptr%p%h)
    else
        dloc = 0
    end if
end subroutine f90wrap_section_obs__array__h

subroutine f90wrap_section_obs__array__u(this, nd, dtype, dshape, dloc)
    use m_model, only: section_obs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type section_obs_ptr_type
        type(section_obs), pointer :: p => NULL()
    end type section_obs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(section_obs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%u)) then
        dshape(1:1) = shape(this_ptr%p%u)
        dloc = loc(this_ptr%p%u)
    else
        dloc = 0
    end if
end subroutine f90wrap_section_obs__array__u

subroutine f90wrap_section_obs__array__v(this, nd, dtype, dshape, dloc)
    use m_model, only: section_obs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type section_obs_ptr_type
        type(section_obs), pointer :: p => NULL()
    end type section_obs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(section_obs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%v)) then
        dshape(1:1) = shape(this_ptr%p%v)
        dloc = loc(this_ptr%p%v)
    else
        dloc = 0
    end if
end subroutine f90wrap_section_obs__array__v

subroutine f90wrap_section_obs__array__q(this, nd, dtype, dshape, dloc)
    use m_model, only: section_obs
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type section_obs_ptr_type
        type(section_obs), pointer :: p => NULL()
    end type section_obs_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(section_obs_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%q)) then
        dshape(1:1) = shape(this_ptr%p%q)
        dloc = loc(this_ptr%p%q)
    else
        dloc = 0
    end if
end subroutine f90wrap_section_obs__array__q

subroutine f90wrap_section_obs_initialise(this)
    use m_model, only: section_obs
    implicit none
    
    type section_obs_ptr_type
        type(section_obs), pointer :: p => NULL()
    end type section_obs_ptr_type
    type(section_obs_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_section_obs_initialise

subroutine f90wrap_section_obs_finalise(this)
    use m_model, only: section_obs
    implicit none
    
    type section_obs_ptr_type
        type(section_obs), pointer :: p => NULL()
    end type section_obs_ptr_type
    type(section_obs_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_section_obs_finalise

subroutine f90wrap_soil_data__get__clay(this, f90wrap_clay)
    use m_model, only: soil_data
    implicit none
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    integer, intent(in)   :: this(2)
    type(soil_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_clay
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_clay = this_ptr%p%clay
end subroutine f90wrap_soil_data__get__clay

subroutine f90wrap_soil_data__set__clay(this, f90wrap_clay)
    use m_model, only: soil_data
    implicit none
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    integer, intent(in)   :: this(2)
    type(soil_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_clay
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%clay = f90wrap_clay
end subroutine f90wrap_soil_data__set__clay

subroutine f90wrap_soil_data__get__silt(this, f90wrap_silt)
    use m_model, only: soil_data
    implicit none
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    integer, intent(in)   :: this(2)
    type(soil_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_silt
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_silt = this_ptr%p%silt
end subroutine f90wrap_soil_data__get__silt

subroutine f90wrap_soil_data__set__silt(this, f90wrap_silt)
    use m_model, only: soil_data
    implicit none
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    integer, intent(in)   :: this(2)
    type(soil_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_silt
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%silt = f90wrap_silt
end subroutine f90wrap_soil_data__set__silt

subroutine f90wrap_soil_data__get__sand(this, f90wrap_sand)
    use m_model, only: soil_data
    implicit none
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    integer, intent(in)   :: this(2)
    type(soil_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_sand
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_sand = this_ptr%p%sand
end subroutine f90wrap_soil_data__get__sand

subroutine f90wrap_soil_data__set__sand(this, f90wrap_sand)
    use m_model, only: soil_data
    implicit none
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    integer, intent(in)   :: this(2)
    type(soil_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_sand
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%sand = f90wrap_sand
end subroutine f90wrap_soil_data__set__sand

subroutine f90wrap_soil_data__get__soil_group(this, f90wrap_soil_group)
    use m_model, only: soil_data
    implicit none
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    integer, intent(in)   :: this(2)
    type(soil_data_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_soil_group
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_soil_group = this_ptr%p%soil_group
end subroutine f90wrap_soil_data__get__soil_group

subroutine f90wrap_soil_data__set__soil_group(this, f90wrap_soil_group)
    use m_model, only: soil_data
    implicit none
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    integer, intent(in)   :: this(2)
    type(soil_data_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_soil_group
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%soil_group = f90wrap_soil_group
end subroutine f90wrap_soil_data__set__soil_group

subroutine f90wrap_soil_data_initialise(this)
    use m_model, only: soil_data
    implicit none
    
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    type(soil_data_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_soil_data_initialise

subroutine f90wrap_soil_data_finalise(this)
    use m_model, only: soil_data
    implicit none
    
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    type(soil_data_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_soil_data_finalise

subroutine f90wrap_surface_data__get__imperm(this, f90wrap_imperm)
    use m_model, only: surface_data
    implicit none
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(in)   :: this(2)
    type(surface_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_imperm
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_imperm = this_ptr%p%imperm
end subroutine f90wrap_surface_data__get__imperm

subroutine f90wrap_surface_data__set__imperm(this, f90wrap_imperm)
    use m_model, only: surface_data
    implicit none
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(in)   :: this(2)
    type(surface_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_imperm
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%imperm = f90wrap_imperm
end subroutine f90wrap_surface_data__set__imperm

subroutine f90wrap_surface_data__get__imperm_group(this, f90wrap_imperm_group)
    use m_model, only: surface_data
    implicit none
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(in)   :: this(2)
    type(surface_data_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_imperm_group
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_imperm_group = this_ptr%p%imperm_group
end subroutine f90wrap_surface_data__get__imperm_group

subroutine f90wrap_surface_data__set__imperm_group(this, f90wrap_imperm_group)
    use m_model, only: surface_data
    implicit none
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(in)   :: this(2)
    type(surface_data_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_imperm_group
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%imperm_group = f90wrap_imperm_group
end subroutine f90wrap_surface_data__set__imperm_group

subroutine f90wrap_surface_data__get__Dmax(this, f90wrap_Dmax)
    use m_model, only: surface_data
    implicit none
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(in)   :: this(2)
    type(surface_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_Dmax
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Dmax = this_ptr%p%Dmax
end subroutine f90wrap_surface_data__get__Dmax

subroutine f90wrap_surface_data__set__Dmax(this, f90wrap_Dmax)
    use m_model, only: surface_data
    implicit none
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(in)   :: this(2)
    type(surface_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_Dmax
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Dmax = f90wrap_Dmax
end subroutine f90wrap_surface_data__set__Dmax

subroutine f90wrap_surface_data__get__Dmax_group(this, f90wrap_Dmax_group)
    use m_model, only: surface_data
    implicit none
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(in)   :: this(2)
    type(surface_data_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_Dmax_group
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_Dmax_group = this_ptr%p%Dmax_group
end subroutine f90wrap_surface_data__get__Dmax_group

subroutine f90wrap_surface_data__set__Dmax_group(this, f90wrap_Dmax_group)
    use m_model, only: surface_data
    implicit none
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(in)   :: this(2)
    type(surface_data_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_Dmax_group
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%Dmax_group = f90wrap_Dmax_group
end subroutine f90wrap_surface_data__set__Dmax_group

subroutine f90wrap_surface_data__get__soil_occ_type(this, f90wrap_soil_occ_type)
    use m_model, only: surface_data
    implicit none
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(in)   :: this(2)
    type(surface_data_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_soil_occ_type
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_soil_occ_type = this_ptr%p%soil_occ_type
end subroutine f90wrap_surface_data__get__soil_occ_type

subroutine f90wrap_surface_data__set__soil_occ_type(this, f90wrap_soil_occ_type)
    use m_model, only: surface_data
    implicit none
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(in)   :: this(2)
    type(surface_data_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_soil_occ_type
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%soil_occ_type = f90wrap_soil_occ_type
end subroutine f90wrap_surface_data__set__soil_occ_type

subroutine f90wrap_surface_data_initialise(this)
    use m_model, only: surface_data
    implicit none
    
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    type(surface_data_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_surface_data_initialise

subroutine f90wrap_surface_data_finalise(this)
    use m_model, only: surface_data
    implicit none
    
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    type(surface_data_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_surface_data_finalise

subroutine f90wrap_structure_data__get__s_type(this, f90wrap_s_type)
    use m_model, only: structure_data
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in)   :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_s_type
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_s_type = this_ptr%p%s_type
end subroutine f90wrap_structure_data__get__s_type

subroutine f90wrap_structure_data__set__s_type(this, f90wrap_s_type)
    use m_model, only: structure_data
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in)   :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_s_type
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%s_type = f90wrap_s_type
end subroutine f90wrap_structure_data__set__s_type

subroutine f90wrap_structure_data__array__name(this, nd, dtype, dshape, dloc)
    use m_model, only: structure_data
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 3
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%name)) then
        dshape(1:3) = (/len(this_ptr%p%name(1,1)), shape(this_ptr%p%name)/)
        dloc = loc(this_ptr%p%name)
    else
        dloc = 0
    end if
end subroutine f90wrap_structure_data__array__name

subroutine f90wrap_structure_data__get__C1(this, f90wrap_C1)
    use m_model, only: structure_data
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in)   :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_C1
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_C1 = this_ptr%p%C1
end subroutine f90wrap_structure_data__get__C1

subroutine f90wrap_structure_data__set__C1(this, f90wrap_C1)
    use m_model, only: structure_data
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in)   :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_C1
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%C1 = f90wrap_C1
end subroutine f90wrap_structure_data__set__C1

subroutine f90wrap_structure_data__get__C2(this, f90wrap_C2)
    use m_model, only: structure_data
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in)   :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_C2
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_C2 = this_ptr%p%C2
end subroutine f90wrap_structure_data__get__C2

subroutine f90wrap_structure_data__set__C2(this, f90wrap_C2)
    use m_model, only: structure_data
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in)   :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_C2
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%C2 = f90wrap_C2
end subroutine f90wrap_structure_data__set__C2

subroutine f90wrap_structure_data__get__C3(this, f90wrap_C3)
    use m_model, only: structure_data
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in)   :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_C3
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_C3 = this_ptr%p%C3
end subroutine f90wrap_structure_data__get__C3

subroutine f90wrap_structure_data__set__C3(this, f90wrap_C3)
    use m_model, only: structure_data
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in)   :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_C3
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%C3 = f90wrap_C3
end subroutine f90wrap_structure_data__set__C3

subroutine f90wrap_structure_data__get__true_x(this, f90wrap_true_x)
    use m_model, only: structure_data
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in)   :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_true_x
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_true_x = this_ptr%p%true_x
end subroutine f90wrap_structure_data__get__true_x

subroutine f90wrap_structure_data__set__true_x(this, f90wrap_true_x)
    use m_model, only: structure_data
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in)   :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_true_x
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%true_x = f90wrap_true_x
end subroutine f90wrap_structure_data__set__true_x

subroutine f90wrap_structure_data__get__true_y(this, f90wrap_true_y)
    use m_model, only: structure_data
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in)   :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_true_y
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_true_y = this_ptr%p%true_y
end subroutine f90wrap_structure_data__get__true_y

subroutine f90wrap_structure_data__set__true_y(this, f90wrap_true_y)
    use m_model, only: structure_data
    implicit none
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in)   :: this(2)
    type(structure_data_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_true_y
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%true_y = f90wrap_true_y
end subroutine f90wrap_structure_data__set__true_y

subroutine f90wrap_structure_data_initialise(this)
    use m_model, only: structure_data
    implicit none
    
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    type(structure_data_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_structure_data_initialise

subroutine f90wrap_structure_data_finalise(this)
    use m_model, only: structure_data
    implicit none
    
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    type(structure_data_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_structure_data_finalise

subroutine f90wrap_input_data__array_getitem__soil(f90wrap_this, f90wrap_i, soilitem)
    
    use m_model, only: soil_data, input_data
    implicit none
    
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(input_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: soilitem(2)
    type(soil_data_ptr_type) :: soil_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%soil)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%soil)) then
            call f90wrap_abort("array index out of range")
        else
            soil_ptr%p => this_ptr%p%soil(f90wrap_i)
            soilitem = transfer(soil_ptr,soilitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_input_data__array_getitem__soil

subroutine f90wrap_input_data__array_setitem__soil(f90wrap_this, f90wrap_i, soilitem)
    
    use m_model, only: soil_data, input_data
    implicit none
    
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(input_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: soilitem(2)
    type(soil_data_ptr_type) :: soil_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%soil)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%soil)) then
            call f90wrap_abort("array index out of range")
        else
            soil_ptr = transfer(soilitem,soil_ptr)
            this_ptr%p%soil(f90wrap_i) = soil_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_input_data__array_setitem__soil

subroutine f90wrap_input_data__array_len__soil(f90wrap_this, f90wrap_n)
    
    use m_model, only: soil_data, input_data
    implicit none
    
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type soil_data_ptr_type
        type(soil_data), pointer :: p => NULL()
    end type soil_data_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(input_data_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%soil)) then
        f90wrap_n = size(this_ptr%p%soil)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_input_data__array_len__soil

subroutine f90wrap_input_data__array_getitem__surf(f90wrap_this, f90wrap_i, surfitem)
    
    use m_model, only: surface_data, input_data
    implicit none
    
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(input_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: surfitem(2)
    type(surface_data_ptr_type) :: surf_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%surf)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%surf)) then
            call f90wrap_abort("array index out of range")
        else
            surf_ptr%p => this_ptr%p%surf(f90wrap_i)
            surfitem = transfer(surf_ptr,surfitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_input_data__array_getitem__surf

subroutine f90wrap_input_data__array_setitem__surf(f90wrap_this, f90wrap_i, surfitem)
    
    use m_model, only: surface_data, input_data
    implicit none
    
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(input_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: surfitem(2)
    type(surface_data_ptr_type) :: surf_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%surf)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%surf)) then
            call f90wrap_abort("array index out of range")
        else
            surf_ptr = transfer(surfitem,surf_ptr)
            this_ptr%p%surf(f90wrap_i) = surf_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_input_data__array_setitem__surf

subroutine f90wrap_input_data__array_len__surf(f90wrap_this, f90wrap_n)
    
    use m_model, only: surface_data, input_data
    implicit none
    
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type surface_data_ptr_type
        type(surface_data), pointer :: p => NULL()
    end type surface_data_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(input_data_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%surf)) then
        f90wrap_n = size(this_ptr%p%surf)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_input_data__array_len__surf

subroutine f90wrap_input_data__array_getitem__structures(f90wrap_this, f90wrap_i, structuresitem)
    
    use m_model, only: structure_data, input_data
    implicit none
    
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(input_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: structuresitem(2)
    type(structure_data_ptr_type) :: structures_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%structures)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%structures)) then
            call f90wrap_abort("array index out of range")
        else
            structures_ptr%p => this_ptr%p%structures(f90wrap_i)
            structuresitem = transfer(structures_ptr,structuresitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_input_data__array_getitem__structures

subroutine f90wrap_input_data__array_setitem__structures(f90wrap_this, f90wrap_i, structuresitem)
    
    use m_model, only: structure_data, input_data
    implicit none
    
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(input_data_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: structuresitem(2)
    type(structure_data_ptr_type) :: structures_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%structures)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%structures)) then
            call f90wrap_abort("array index out of range")
        else
            structures_ptr = transfer(structuresitem,structures_ptr)
            this_ptr%p%structures(f90wrap_i) = structures_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_input_data__array_setitem__structures

subroutine f90wrap_input_data__array_len__structures(f90wrap_this, f90wrap_n)
    
    use m_model, only: structure_data, input_data
    implicit none
    
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type structure_data_ptr_type
        type(structure_data), pointer :: p => NULL()
    end type structure_data_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(input_data_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%structures)) then
        f90wrap_n = size(this_ptr%p%structures)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_input_data__array_len__structures

subroutine f90wrap_input_data_initialise(this)
    use m_model, only: input_data
    implicit none
    
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type(input_data_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_input_data_initialise

subroutine f90wrap_input_data_finalise(this)
    use m_model, only: input_data
    implicit none
    
    type input_data_ptr_type
        type(input_data), pointer :: p => NULL()
    end type input_data_ptr_type
    type(input_data_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_input_data_finalise

subroutine f90wrap_input_param__get__mesh_type(this, f90wrap_mesh_type)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_mesh_type
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_mesh_type = this_ptr%p%mesh_type
end subroutine f90wrap_input_param__get__mesh_type

subroutine f90wrap_input_param__set__mesh_type(this, f90wrap_mesh_type)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_mesh_type
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%mesh_type = f90wrap_mesh_type
end subroutine f90wrap_input_param__set__mesh_type

subroutine f90wrap_input_param__get__mesh_name(this, f90wrap_mesh_name)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_mesh_name
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_mesh_name = this_ptr%p%mesh_name
end subroutine f90wrap_input_param__get__mesh_name

subroutine f90wrap_input_param__set__mesh_name(this, f90wrap_mesh_name)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_mesh_name
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%mesh_name = f90wrap_mesh_name
end subroutine f90wrap_input_param__set__mesh_name

subroutine f90wrap_input_param__get__bc_N(this, f90wrap_bc_N)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_bc_N
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_bc_N = this_ptr%p%bc_N
end subroutine f90wrap_input_param__get__bc_N

subroutine f90wrap_input_param__set__bc_N(this, f90wrap_bc_N)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_bc_N
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%bc_N = f90wrap_bc_N
end subroutine f90wrap_input_param__set__bc_N

subroutine f90wrap_input_param__get__bc_S(this, f90wrap_bc_S)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_bc_S
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_bc_S = this_ptr%p%bc_S
end subroutine f90wrap_input_param__get__bc_S

subroutine f90wrap_input_param__set__bc_S(this, f90wrap_bc_S)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_bc_S
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%bc_S = f90wrap_bc_S
end subroutine f90wrap_input_param__set__bc_S

subroutine f90wrap_input_param__get__bc_W(this, f90wrap_bc_W)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_bc_W
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_bc_W = this_ptr%p%bc_W
end subroutine f90wrap_input_param__get__bc_W

subroutine f90wrap_input_param__set__bc_W(this, f90wrap_bc_W)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_bc_W
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%bc_W = f90wrap_bc_W
end subroutine f90wrap_input_param__set__bc_W

subroutine f90wrap_input_param__get__bc_E(this, f90wrap_bc_E)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_bc_E
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_bc_E = this_ptr%p%bc_E
end subroutine f90wrap_input_param__get__bc_E

subroutine f90wrap_input_param__set__bc_E(this, f90wrap_bc_E)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_bc_E
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%bc_E = f90wrap_bc_E
end subroutine f90wrap_input_param__set__bc_E

subroutine f90wrap_input_param__get__bc_rain(this, f90wrap_bc_rain)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_bc_rain
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_bc_rain = this_ptr%p%bc_rain
end subroutine f90wrap_input_param__get__bc_rain

subroutine f90wrap_input_param__set__bc_rain(this, f90wrap_bc_rain)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_bc_rain
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%bc_rain = f90wrap_bc_rain
end subroutine f90wrap_input_param__set__bc_rain

subroutine f90wrap_input_param__get__bc_infil(this, f90wrap_bc_infil)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_bc_infil
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_bc_infil = this_ptr%p%bc_infil
end subroutine f90wrap_input_param__get__bc_infil

subroutine f90wrap_input_param__set__bc_infil(this, f90wrap_bc_infil)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_bc_infil
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%bc_infil = f90wrap_bc_infil
end subroutine f90wrap_input_param__set__bc_infil

subroutine f90wrap_input_param__get__lx(this, f90wrap_lx)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_lx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_lx = this_ptr%p%lx
end subroutine f90wrap_input_param__get__lx

subroutine f90wrap_input_param__set__lx(this, f90wrap_lx)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_lx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lx = f90wrap_lx
end subroutine f90wrap_input_param__set__lx

subroutine f90wrap_input_param__get__ly(this, f90wrap_ly)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_ly
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ly = this_ptr%p%ly
end subroutine f90wrap_input_param__get__ly

subroutine f90wrap_input_param__set__ly(this, f90wrap_ly)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_ly
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ly = f90wrap_ly
end subroutine f90wrap_input_param__set__ly

subroutine f90wrap_input_param__get__nx(this, f90wrap_nx)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nx
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nx = this_ptr%p%nx
end subroutine f90wrap_input_param__get__nx

subroutine f90wrap_input_param__set__nx(this, f90wrap_nx)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nx
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nx = f90wrap_nx
end subroutine f90wrap_input_param__set__nx

subroutine f90wrap_input_param__get__ny(this, f90wrap_ny)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ny
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ny = this_ptr%p%ny
end subroutine f90wrap_input_param__get__ny

subroutine f90wrap_input_param__set__ny(this, f90wrap_ny)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ny
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ny = f90wrap_ny
end subroutine f90wrap_input_param__set__ny

subroutine f90wrap_input_param__get__ts(this, f90wrap_ts)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_ts
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ts = this_ptr%p%ts
end subroutine f90wrap_input_param__get__ts

subroutine f90wrap_input_param__set__ts(this, f90wrap_ts)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_ts
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ts = f90wrap_ts
end subroutine f90wrap_input_param__set__ts

subroutine f90wrap_input_param__get__adapt_dt(this, f90wrap_adapt_dt)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_adapt_dt
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_adapt_dt = this_ptr%p%adapt_dt
end subroutine f90wrap_input_param__get__adapt_dt

subroutine f90wrap_input_param__set__adapt_dt(this, f90wrap_adapt_dt)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_adapt_dt
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%adapt_dt = f90wrap_adapt_dt
end subroutine f90wrap_input_param__set__adapt_dt

subroutine f90wrap_input_param__get__dt(this, f90wrap_dt)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_dt
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_dt = this_ptr%p%dt
end subroutine f90wrap_input_param__get__dt

subroutine f90wrap_input_param__set__dt(this, f90wrap_dt)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_dt
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%dt = f90wrap_dt
end subroutine f90wrap_input_param__set__dt

subroutine f90wrap_input_param__get__cfl(this, f90wrap_cfl)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_cfl
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_cfl = this_ptr%p%cfl
end subroutine f90wrap_input_param__get__cfl

subroutine f90wrap_input_param__set__cfl(this, f90wrap_cfl)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_cfl
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%cfl = f90wrap_cfl
end subroutine f90wrap_input_param__set__cfl

subroutine f90wrap_input_param__get__dtw(this, f90wrap_dtw)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_dtw
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_dtw = this_ptr%p%dtw
end subroutine f90wrap_input_param__get__dtw

subroutine f90wrap_input_param__set__dtw(this, f90wrap_dtw)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_dtw
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%dtw = f90wrap_dtw
end subroutine f90wrap_input_param__set__dtw

subroutine f90wrap_input_param__get__dtp(this, f90wrap_dtp)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_dtp
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_dtp = this_ptr%p%dtp
end subroutine f90wrap_input_param__get__dtp

subroutine f90wrap_input_param__set__dtp(this, f90wrap_dtp)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_dtp
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%dtp = f90wrap_dtp
end subroutine f90wrap_input_param__set__dtp

subroutine f90wrap_input_param__get__dta(this, f90wrap_dta)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_dta
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_dta = this_ptr%p%dta
end subroutine f90wrap_input_param__get__dta

subroutine f90wrap_input_param__set__dta(this, f90wrap_dta)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_dta
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%dta = f90wrap_dta
end subroutine f90wrap_input_param__set__dta

subroutine f90wrap_input_param__get__w_tecplot(this, f90wrap_w_tecplot)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_w_tecplot
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_w_tecplot = this_ptr%p%w_tecplot
end subroutine f90wrap_input_param__get__w_tecplot

subroutine f90wrap_input_param__set__w_tecplot(this, f90wrap_w_tecplot)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_w_tecplot
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%w_tecplot = f90wrap_w_tecplot
end subroutine f90wrap_input_param__set__w_tecplot

subroutine f90wrap_input_param__get__w_vtk(this, f90wrap_w_vtk)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_w_vtk
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_w_vtk = this_ptr%p%w_vtk
end subroutine f90wrap_input_param__get__w_vtk

subroutine f90wrap_input_param__set__w_vtk(this, f90wrap_w_vtk)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_w_vtk
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%w_vtk = f90wrap_w_vtk
end subroutine f90wrap_input_param__set__w_vtk

subroutine f90wrap_input_param__get__w_gnuplot(this, f90wrap_w_gnuplot)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_w_gnuplot
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_w_gnuplot = this_ptr%p%w_gnuplot
end subroutine f90wrap_input_param__get__w_gnuplot

subroutine f90wrap_input_param__set__w_gnuplot(this, f90wrap_w_gnuplot)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_w_gnuplot
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%w_gnuplot = f90wrap_w_gnuplot
end subroutine f90wrap_input_param__set__w_gnuplot

subroutine f90wrap_input_param__get__w_bin(this, f90wrap_w_bin)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_w_bin
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_w_bin = this_ptr%p%w_bin
end subroutine f90wrap_input_param__get__w_bin

subroutine f90wrap_input_param__set__w_bin(this, f90wrap_w_bin)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_w_bin
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%w_bin = f90wrap_w_bin
end subroutine f90wrap_input_param__set__w_bin

subroutine f90wrap_input_param__get__w_exact(this, f90wrap_w_exact)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_w_exact
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_w_exact = this_ptr%p%w_exact
end subroutine f90wrap_input_param__get__w_exact

subroutine f90wrap_input_param__set__w_exact(this, f90wrap_w_exact)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_w_exact
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%w_exact = f90wrap_w_exact
end subroutine f90wrap_input_param__set__w_exact

subroutine f90wrap_input_param__get__w_norm(this, f90wrap_w_norm)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_w_norm
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_w_norm = this_ptr%p%w_norm
end subroutine f90wrap_input_param__get__w_norm

subroutine f90wrap_input_param__set__w_norm(this, f90wrap_w_norm)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_w_norm
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%w_norm = f90wrap_w_norm
end subroutine f90wrap_input_param__set__w_norm

subroutine f90wrap_input_param__get__w_obs(this, f90wrap_w_obs)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_w_obs
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_w_obs = this_ptr%p%w_obs
end subroutine f90wrap_input_param__get__w_obs

subroutine f90wrap_input_param__set__w_obs(this, f90wrap_w_obs)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_w_obs
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%w_obs = f90wrap_w_obs
end subroutine f90wrap_input_param__set__w_obs

subroutine f90wrap_input_param__get__use_obs(this, f90wrap_use_obs)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_use_obs
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_use_obs = this_ptr%p%use_obs
end subroutine f90wrap_input_param__get__use_obs

subroutine f90wrap_input_param__set__use_obs(this, f90wrap_use_obs)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_use_obs
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%use_obs = f90wrap_use_obs
end subroutine f90wrap_input_param__set__use_obs

subroutine f90wrap_input_param__get__spatial_scheme(this, f90wrap_spatial_scheme)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_spatial_scheme
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_spatial_scheme = this_ptr%p%spatial_scheme
end subroutine f90wrap_input_param__get__spatial_scheme

subroutine f90wrap_input_param__set__spatial_scheme(this, f90wrap_spatial_scheme)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_spatial_scheme
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%spatial_scheme = f90wrap_spatial_scheme
end subroutine f90wrap_input_param__set__spatial_scheme

subroutine f90wrap_input_param__get__temp_scheme(this, f90wrap_temp_scheme)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_temp_scheme
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_temp_scheme = this_ptr%p%temp_scheme
end subroutine f90wrap_input_param__get__temp_scheme

subroutine f90wrap_input_param__set__temp_scheme(this, f90wrap_temp_scheme)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_temp_scheme
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%temp_scheme = f90wrap_temp_scheme
end subroutine f90wrap_input_param__set__temp_scheme

subroutine f90wrap_input_param__array__args(this, nd, dtype, dshape, dloc)
    use m_model, only: input_param
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%args)) then
        dshape(1:2) = (/len(this_ptr%p%args(1)), shape(this_ptr%p%args)/)
        dloc = loc(this_ptr%p%args)
    else
        dloc = 0
    end if
end subroutine f90wrap_input_param__array__args

subroutine f90wrap_input_param__get__max_nt_for_direct(this, f90wrap_max_nt_for_direct)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_max_nt_for_direct
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_nt_for_direct = this_ptr%p%max_nt_for_direct
end subroutine f90wrap_input_param__get__max_nt_for_direct

subroutine f90wrap_input_param__set__max_nt_for_direct(this, f90wrap_max_nt_for_direct)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_max_nt_for_direct
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_nt_for_direct = f90wrap_max_nt_for_direct
end subroutine f90wrap_input_param__set__max_nt_for_direct

subroutine f90wrap_input_param__get__max_nt_for_adjoint(this, f90wrap_max_nt_for_adjoint)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_max_nt_for_adjoint
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_nt_for_adjoint = this_ptr%p%max_nt_for_adjoint
end subroutine f90wrap_input_param__get__max_nt_for_adjoint

subroutine f90wrap_input_param__set__max_nt_for_adjoint(this, f90wrap_max_nt_for_adjoint)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_max_nt_for_adjoint
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_nt_for_adjoint = f90wrap_max_nt_for_adjoint
end subroutine f90wrap_input_param__set__max_nt_for_adjoint

subroutine f90wrap_input_param__get__g(this, f90wrap_g)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_g
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_g = this_ptr%p%g
end subroutine f90wrap_input_param__get__g

subroutine f90wrap_input_param__set__g(this, f90wrap_g)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_g
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%g = f90wrap_g
end subroutine f90wrap_input_param__set__g

subroutine f90wrap_input_param__get__heps(this, f90wrap_heps)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_heps
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_heps = this_ptr%p%heps
end subroutine f90wrap_input_param__get__heps

subroutine f90wrap_input_param__set__heps(this, f90wrap_heps)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_heps
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%heps = f90wrap_heps
end subroutine f90wrap_input_param__set__heps

subroutine f90wrap_input_param__get__friction(this, f90wrap_friction)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_friction
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_friction = this_ptr%p%friction
end subroutine f90wrap_input_param__get__friction

subroutine f90wrap_input_param__set__friction(this, f90wrap_friction)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_friction
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%friction = f90wrap_friction
end subroutine f90wrap_input_param__set__friction

subroutine f90wrap_input_param__get__c_manning(this, f90wrap_c_manning)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_c_manning
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_c_manning = this_ptr%p%c_manning
end subroutine f90wrap_input_param__get__c_manning

subroutine f90wrap_input_param__set__c_manning(this, f90wrap_c_manning)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_c_manning
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%c_manning = f90wrap_c_manning
end subroutine f90wrap_input_param__set__c_manning

subroutine f90wrap_input_param__get__c_manning_beta(this, f90wrap_c_manning_beta)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_c_manning_beta
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_c_manning_beta = this_ptr%p%c_manning_beta
end subroutine f90wrap_input_param__get__c_manning_beta

subroutine f90wrap_input_param__set__c_manning_beta(this, f90wrap_c_manning_beta)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_c_manning_beta
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%c_manning_beta = f90wrap_c_manning_beta
end subroutine f90wrap_input_param__set__c_manning_beta

subroutine f90wrap_input_param__get__c_bathy(this, f90wrap_c_bathy)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_c_bathy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_c_bathy = this_ptr%p%c_bathy
end subroutine f90wrap_input_param__get__c_bathy

subroutine f90wrap_input_param__set__c_bathy(this, f90wrap_c_bathy)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_c_bathy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%c_bathy = f90wrap_c_bathy
end subroutine f90wrap_input_param__set__c_bathy

subroutine f90wrap_input_param__get__c_ic(this, f90wrap_c_ic)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_c_ic
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_c_ic = this_ptr%p%c_ic
end subroutine f90wrap_input_param__get__c_ic

subroutine f90wrap_input_param__set__c_ic(this, f90wrap_c_ic)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_c_ic
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%c_ic = f90wrap_c_ic
end subroutine f90wrap_input_param__set__c_ic

subroutine f90wrap_input_param__get__c_hydrograph(this, f90wrap_c_hydrograph)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_c_hydrograph
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_c_hydrograph = this_ptr%p%c_hydrograph
end subroutine f90wrap_input_param__get__c_hydrograph

subroutine f90wrap_input_param__set__c_hydrograph(this, f90wrap_c_hydrograph)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_c_hydrograph
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%c_hydrograph = f90wrap_c_hydrograph
end subroutine f90wrap_input_param__set__c_hydrograph

subroutine f90wrap_input_param__get__c_ratcurve(this, f90wrap_c_ratcurve)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_c_ratcurve
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_c_ratcurve = this_ptr%p%c_ratcurve
end subroutine f90wrap_input_param__get__c_ratcurve

subroutine f90wrap_input_param__set__c_ratcurve(this, f90wrap_c_ratcurve)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_c_ratcurve
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%c_ratcurve = f90wrap_c_ratcurve
end subroutine f90wrap_input_param__set__c_ratcurve

subroutine f90wrap_input_param__get__c_rain(this, f90wrap_c_rain)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_c_rain
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_c_rain = this_ptr%p%c_rain
end subroutine f90wrap_input_param__get__c_rain

subroutine f90wrap_input_param__set__c_rain(this, f90wrap_c_rain)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_c_rain
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%c_rain = f90wrap_c_rain
end subroutine f90wrap_input_param__set__c_rain

subroutine f90wrap_input_param__get__eps_min(this, f90wrap_eps_min)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps_min
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eps_min = this_ptr%p%eps_min
end subroutine f90wrap_input_param__get__eps_min

subroutine f90wrap_input_param__set__eps_min(this, f90wrap_eps_min)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps_min
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps_min = f90wrap_eps_min
end subroutine f90wrap_input_param__set__eps_min

subroutine f90wrap_input_param__get__eps_manning(this, f90wrap_eps_manning)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps_manning
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eps_manning = this_ptr%p%eps_manning
end subroutine f90wrap_input_param__get__eps_manning

subroutine f90wrap_input_param__set__eps_manning(this, f90wrap_eps_manning)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps_manning
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps_manning = f90wrap_eps_manning
end subroutine f90wrap_input_param__set__eps_manning

subroutine f90wrap_input_param__get__eps_bathy(this, f90wrap_eps_bathy)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps_bathy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eps_bathy = this_ptr%p%eps_bathy
end subroutine f90wrap_input_param__get__eps_bathy

subroutine f90wrap_input_param__set__eps_bathy(this, f90wrap_eps_bathy)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps_bathy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps_bathy = f90wrap_eps_bathy
end subroutine f90wrap_input_param__set__eps_bathy

subroutine f90wrap_input_param__get__eps_ic(this, f90wrap_eps_ic)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps_ic
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eps_ic = this_ptr%p%eps_ic
end subroutine f90wrap_input_param__get__eps_ic

subroutine f90wrap_input_param__set__eps_ic(this, f90wrap_eps_ic)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps_ic
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps_ic = f90wrap_eps_ic
end subroutine f90wrap_input_param__set__eps_ic

subroutine f90wrap_input_param__get__eps_hydrograph(this, f90wrap_eps_hydrograph)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps_hydrograph
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eps_hydrograph = this_ptr%p%eps_hydrograph
end subroutine f90wrap_input_param__get__eps_hydrograph

subroutine f90wrap_input_param__set__eps_hydrograph(this, f90wrap_eps_hydrograph)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps_hydrograph
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps_hydrograph = f90wrap_eps_hydrograph
end subroutine f90wrap_input_param__set__eps_hydrograph

subroutine f90wrap_input_param__get__eps_ratcurve(this, f90wrap_eps_ratcurve)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps_ratcurve
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eps_ratcurve = this_ptr%p%eps_ratcurve
end subroutine f90wrap_input_param__get__eps_ratcurve

subroutine f90wrap_input_param__set__eps_ratcurve(this, f90wrap_eps_ratcurve)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps_ratcurve
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps_ratcurve = f90wrap_eps_ratcurve
end subroutine f90wrap_input_param__set__eps_ratcurve

subroutine f90wrap_input_param__get__eps_rain(this, f90wrap_eps_rain)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps_rain
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eps_rain = this_ptr%p%eps_rain
end subroutine f90wrap_input_param__get__eps_rain

subroutine f90wrap_input_param__set__eps_rain(this, f90wrap_eps_rain)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps_rain
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps_rain = f90wrap_eps_rain
end subroutine f90wrap_input_param__set__eps_rain

subroutine f90wrap_input_param__get__eps_Ks(this, f90wrap_eps_Ks)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps_Ks
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eps_Ks = this_ptr%p%eps_Ks
end subroutine f90wrap_input_param__get__eps_Ks

subroutine f90wrap_input_param__set__eps_Ks(this, f90wrap_eps_Ks)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps_Ks
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps_Ks = f90wrap_eps_Ks
end subroutine f90wrap_input_param__set__eps_Ks

subroutine f90wrap_input_param__get__eps_PsiF(this, f90wrap_eps_PsiF)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps_PsiF
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eps_PsiF = this_ptr%p%eps_PsiF
end subroutine f90wrap_input_param__get__eps_PsiF

subroutine f90wrap_input_param__set__eps_PsiF(this, f90wrap_eps_PsiF)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps_PsiF
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps_PsiF = f90wrap_eps_PsiF
end subroutine f90wrap_input_param__set__eps_PsiF

subroutine f90wrap_input_param__get__eps_DeltaTheta(this, f90wrap_eps_DeltaTheta)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps_DeltaTheta
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eps_DeltaTheta = this_ptr%p%eps_DeltaTheta
end subroutine f90wrap_input_param__get__eps_DeltaTheta

subroutine f90wrap_input_param__set__eps_DeltaTheta(this, f90wrap_eps_DeltaTheta)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps_DeltaTheta
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps_DeltaTheta = f90wrap_eps_DeltaTheta
end subroutine f90wrap_input_param__set__eps_DeltaTheta

subroutine f90wrap_input_param__get__eps_lambda(this, f90wrap_eps_lambda)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps_lambda
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eps_lambda = this_ptr%p%eps_lambda
end subroutine f90wrap_input_param__get__eps_lambda

subroutine f90wrap_input_param__set__eps_lambda(this, f90wrap_eps_lambda)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps_lambda
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps_lambda = f90wrap_eps_lambda
end subroutine f90wrap_input_param__set__eps_lambda

subroutine f90wrap_input_param__get__eps_CN(this, f90wrap_eps_CN)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_eps_CN
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_eps_CN = this_ptr%p%eps_CN
end subroutine f90wrap_input_param__get__eps_CN

subroutine f90wrap_input_param__set__eps_CN(this, f90wrap_eps_CN)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_eps_CN
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%eps_CN = f90wrap_eps_CN
end subroutine f90wrap_input_param__set__eps_CN

subroutine f90wrap_input_param__get__regul_manning(this, f90wrap_regul_manning)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_regul_manning
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_regul_manning = this_ptr%p%regul_manning
end subroutine f90wrap_input_param__get__regul_manning

subroutine f90wrap_input_param__set__regul_manning(this, f90wrap_regul_manning)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_regul_manning
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%regul_manning = f90wrap_regul_manning
end subroutine f90wrap_input_param__set__regul_manning

subroutine f90wrap_input_param__get__regul_bathy(this, f90wrap_regul_bathy)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_regul_bathy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_regul_bathy = this_ptr%p%regul_bathy
end subroutine f90wrap_input_param__get__regul_bathy

subroutine f90wrap_input_param__set__regul_bathy(this, f90wrap_regul_bathy)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_regul_bathy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%regul_bathy = f90wrap_regul_bathy
end subroutine f90wrap_input_param__set__regul_bathy

subroutine f90wrap_input_param__get__regul_ic(this, f90wrap_regul_ic)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_regul_ic
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_regul_ic = this_ptr%p%regul_ic
end subroutine f90wrap_input_param__get__regul_ic

subroutine f90wrap_input_param__set__regul_ic(this, f90wrap_regul_ic)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_regul_ic
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%regul_ic = f90wrap_regul_ic
end subroutine f90wrap_input_param__set__regul_ic

subroutine f90wrap_input_param__get__regul_hydrograph(this, f90wrap_regul_hydrograph)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_regul_hydrograph
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_regul_hydrograph = this_ptr%p%regul_hydrograph
end subroutine f90wrap_input_param__get__regul_hydrograph

subroutine f90wrap_input_param__set__regul_hydrograph(this, f90wrap_regul_hydrograph)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_regul_hydrograph
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%regul_hydrograph = f90wrap_regul_hydrograph
end subroutine f90wrap_input_param__set__regul_hydrograph

subroutine f90wrap_input_param__get__regul_ratcurve(this, f90wrap_regul_ratcurve)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_regul_ratcurve
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_regul_ratcurve = this_ptr%p%regul_ratcurve
end subroutine f90wrap_input_param__get__regul_ratcurve

subroutine f90wrap_input_param__set__regul_ratcurve(this, f90wrap_regul_ratcurve)
    use m_model, only: input_param
    implicit none
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    integer, intent(in)   :: this(2)
    type(input_param_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_regul_ratcurve
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%regul_ratcurve = f90wrap_regul_ratcurve
end subroutine f90wrap_input_param__set__regul_ratcurve

subroutine f90wrap_input_param_initialise(this)
    use m_model, only: input_param
    implicit none
    
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    type(input_param_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_input_param_initialise

subroutine f90wrap_input_param_finalise(this)
    use m_model, only: input_param
    implicit none
    
    type input_param_ptr_type
        type(input_param), pointer :: p => NULL()
    end type input_param_ptr_type
    type(input_param_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_input_param_finalise

subroutine f90wrap_default_values
    use m_model, only: default_values
    implicit none
    
    call default_values()
end subroutine f90wrap_default_values

subroutine f90wrap_alloc_dof(dof, mesh)
    use m_mesh, only: msh
    use m_model, only: unk, alloc_dof
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(unk_ptr_type) :: dof_ptr
    integer, intent(out), dimension(2) :: dof
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    mesh_ptr = transfer(mesh, mesh_ptr)
    allocate(dof_ptr%p)
    call alloc_dof(dof=dof_ptr%p, mesh=mesh_ptr%p)
    dof = transfer(dof_ptr, dof)
end subroutine f90wrap_alloc_dof

subroutine f90wrap_dealloc_dof(dof)
    use m_model, only: unk, dealloc_dof
    implicit none
    
    type unk_ptr_type
        type(unk), pointer :: p => NULL()
    end type unk_ptr_type
    type(unk_ptr_type) :: dof_ptr
    integer, intent(in), dimension(2) :: dof
    dof_ptr = transfer(dof, dof_ptr)
    call dealloc_dof(dof=dof_ptr%p)
end subroutine f90wrap_dealloc_dof

subroutine f90wrap_dealloc_model
    use m_model, only: dealloc_model
    implicit none
    
    call dealloc_model()
end subroutine f90wrap_dealloc_model

subroutine f90wrap_m_model__get__sw_nb(f90wrap_sw_nb)
    use m_model, only: m_model_sw_nb => sw_nb
    implicit none
    integer(4), intent(out) :: f90wrap_sw_nb
    
    f90wrap_sw_nb = m_model_sw_nb
end subroutine f90wrap_m_model__get__sw_nb

subroutine f90wrap_m_model__set__sw_nb(f90wrap_sw_nb)
    use m_model, only: m_model_sw_nb => sw_nb
    implicit none
    integer(4), intent(in) :: f90wrap_sw_nb
    
    m_model_sw_nb = f90wrap_sw_nb
end subroutine f90wrap_m_model__set__sw_nb

subroutine f90wrap_m_model__array__bathy_node(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_time_screen
    use m_model, only: m_model_bathy_node => bathy_node
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_model_bathy_node)) then
        dshape(1:1) = shape(m_model_bathy_node)
        dloc = loc(m_model_bathy_node)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_model__array__bathy_node

subroutine f90wrap_m_model__array__bathy_cell(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_time_screen
    use m_model, only: m_model_bathy_cell => bathy_cell
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_model_bathy_cell)) then
        dshape(1:1) = shape(m_model_bathy_cell)
        dloc = loc(m_model_bathy_cell)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_model__array__bathy_cell

subroutine f90wrap_m_model__array__manning(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_time_screen
    use m_model, only: m_model_manning => manning
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_model_manning)) then
        dshape(1:1) = shape(m_model_manning)
        dloc = loc(m_model_manning)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_model__array__manning

subroutine f90wrap_m_model__array__manning_beta(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_time_screen
    use m_model, only: m_model_manning_beta => manning_beta
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(m_model_manning_beta)) then
        dshape(1:1) = shape(m_model_manning_beta)
        dloc = loc(m_model_manning_beta)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_model__array__manning_beta

subroutine f90wrap_m_model__get__nland(f90wrap_nland)
    use m_model, only: m_model_nland => nland
    implicit none
    integer(4), intent(out) :: f90wrap_nland
    
    f90wrap_nland = m_model_nland
end subroutine f90wrap_m_model__get__nland

subroutine f90wrap_m_model__set__nland(f90wrap_nland)
    use m_model, only: m_model_nland => nland
    implicit none
    integer(4), intent(in) :: f90wrap_nland
    
    m_model_nland = f90wrap_nland
end subroutine f90wrap_m_model__set__nland

subroutine f90wrap_m_model__array__land(dummy_this, nd, dtype, dshape, dloc)
    use m_common
    use m_linear_algebra
    use m_mesh
    use m_time_screen
    use m_model, only: m_model_land => land
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    if (allocated(m_model_land)) then
        dshape(1:1) = shape(m_model_land)
        dloc = loc(m_model_land)
    else
        dloc = 0
    end if
end subroutine f90wrap_m_model__array__land

subroutine f90wrap_m_model__get__mass_cut(f90wrap_mass_cut)
    use m_model, only: m_model_mass_cut => mass_cut
    implicit none
    real(8), intent(out) :: f90wrap_mass_cut
    
    f90wrap_mass_cut = m_model_mass_cut
end subroutine f90wrap_m_model__get__mass_cut

subroutine f90wrap_m_model__set__mass_cut(f90wrap_mass_cut)
    use m_model, only: m_model_mass_cut => mass_cut
    implicit none
    real(8), intent(in) :: f90wrap_mass_cut
    
    m_model_mass_cut = f90wrap_mass_cut
end subroutine f90wrap_m_model__set__mass_cut

subroutine f90wrap_m_model__get__manning_data_glob(f90wrap_manning_data_glob)
    use m_model, only: m_model_manning_data_glob => manning_data_glob
    implicit none
    integer(4), intent(out) :: f90wrap_manning_data_glob
    
    f90wrap_manning_data_glob = m_model_manning_data_glob
end subroutine f90wrap_m_model__get__manning_data_glob

subroutine f90wrap_m_model__set__manning_data_glob(f90wrap_manning_data_glob)
    use m_model, only: m_model_manning_data_glob => manning_data_glob
    implicit none
    integer(4), intent(in) :: f90wrap_manning_data_glob
    
    m_model_manning_data_glob = f90wrap_manning_data_glob
end subroutine f90wrap_m_model__set__manning_data_glob

subroutine f90wrap_m_model__get__feedback_inflow(f90wrap_feedback_inflow)
    use m_model, only: m_model_feedback_inflow => feedback_inflow
    implicit none
    integer(4), intent(out) :: f90wrap_feedback_inflow
    
    f90wrap_feedback_inflow = m_model_feedback_inflow
end subroutine f90wrap_m_model__get__feedback_inflow

subroutine f90wrap_m_model__set__feedback_inflow(f90wrap_feedback_inflow)
    use m_model, only: m_model_feedback_inflow => feedback_inflow
    implicit none
    integer(4), intent(in) :: f90wrap_feedback_inflow
    
    m_model_feedback_inflow = f90wrap_feedback_inflow
end subroutine f90wrap_m_model__set__feedback_inflow

subroutine f90wrap_m_model__get__coef_feedback(f90wrap_coef_feedback)
    use m_model, only: m_model_coef_feedback => coef_feedback
    implicit none
    real(8), intent(out) :: f90wrap_coef_feedback
    
    f90wrap_coef_feedback = m_model_coef_feedback
end subroutine f90wrap_m_model__get__coef_feedback

subroutine f90wrap_m_model__set__coef_feedback(f90wrap_coef_feedback)
    use m_model, only: m_model_coef_feedback => coef_feedback
    implicit none
    real(8), intent(in) :: f90wrap_coef_feedback
    
    m_model_coef_feedback = f90wrap_coef_feedback
end subroutine f90wrap_m_model__set__coef_feedback

subroutine f90wrap_m_model__get__bc(f90wrap_bc)
    use m_model, only: m_model_bc => bc, bcs
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(out) :: f90wrap_bc(2)
    type(bcs_ptr_type) :: bc_ptr
    
    bc_ptr%p => m_model_bc
    f90wrap_bc = transfer(bc_ptr,f90wrap_bc)
end subroutine f90wrap_m_model__get__bc

subroutine f90wrap_m_model__set__bc(f90wrap_bc)
    use m_model, only: m_model_bc => bc, bcs
    implicit none
    type bcs_ptr_type
        type(bcs), pointer :: p => NULL()
    end type bcs_ptr_type
    integer, intent(in) :: f90wrap_bc(2)
    type(bcs_ptr_type) :: bc_ptr
    
    bc_ptr = transfer(f90wrap_bc,bc_ptr)
    m_model_bc = bc_ptr%p
end subroutine f90wrap_m_model__set__bc

subroutine f90wrap_m_model__get__g(f90wrap_g)
    use m_model, only: m_model_g => g
    implicit none
    real(8), intent(out) :: f90wrap_g
    
    f90wrap_g = m_model_g
end subroutine f90wrap_m_model__get__g

subroutine f90wrap_m_model__set__g(f90wrap_g)
    use m_model, only: m_model_g => g
    implicit none
    real(8), intent(in) :: f90wrap_g
    
    m_model_g = f90wrap_g
end subroutine f90wrap_m_model__set__g

subroutine f90wrap_m_model__get__heps(f90wrap_heps)
    use m_model, only: m_model_heps => heps
    implicit none
    real(8), intent(out) :: f90wrap_heps
    
    f90wrap_heps = m_model_heps
end subroutine f90wrap_m_model__get__heps

subroutine f90wrap_m_model__set__heps(f90wrap_heps)
    use m_model, only: m_model_heps => heps
    implicit none
    real(8), intent(in) :: f90wrap_heps
    
    m_model_heps = f90wrap_heps
end subroutine f90wrap_m_model__set__heps

subroutine f90wrap_m_model__get__friction(f90wrap_friction)
    use m_model, only: m_model_friction => friction
    implicit none
    integer(4), intent(out) :: f90wrap_friction
    
    f90wrap_friction = m_model_friction
end subroutine f90wrap_m_model__get__friction

subroutine f90wrap_m_model__set__friction(f90wrap_friction)
    use m_model, only: m_model_friction => friction
    implicit none
    integer(4), intent(in) :: f90wrap_friction
    
    m_model_friction = f90wrap_friction
end subroutine f90wrap_m_model__set__friction

subroutine f90wrap_m_model__get__c_manning(f90wrap_c_manning)
    use m_model, only: m_model_c_manning => c_manning
    implicit none
    integer(4), intent(out) :: f90wrap_c_manning
    
    f90wrap_c_manning = m_model_c_manning
end subroutine f90wrap_m_model__get__c_manning

subroutine f90wrap_m_model__set__c_manning(f90wrap_c_manning)
    use m_model, only: m_model_c_manning => c_manning
    implicit none
    integer(4), intent(in) :: f90wrap_c_manning
    
    m_model_c_manning = f90wrap_c_manning
end subroutine f90wrap_m_model__set__c_manning

subroutine f90wrap_m_model__get__c_manning_beta(f90wrap_c_manning_beta)
    use m_model, only: m_model_c_manning_beta => c_manning_beta
    implicit none
    integer(4), intent(out) :: f90wrap_c_manning_beta
    
    f90wrap_c_manning_beta = m_model_c_manning_beta
end subroutine f90wrap_m_model__get__c_manning_beta

subroutine f90wrap_m_model__set__c_manning_beta(f90wrap_c_manning_beta)
    use m_model, only: m_model_c_manning_beta => c_manning_beta
    implicit none
    integer(4), intent(in) :: f90wrap_c_manning_beta
    
    m_model_c_manning_beta = f90wrap_c_manning_beta
end subroutine f90wrap_m_model__set__c_manning_beta

subroutine f90wrap_m_model__get__c_bathy(f90wrap_c_bathy)
    use m_model, only: m_model_c_bathy => c_bathy
    implicit none
    integer(4), intent(out) :: f90wrap_c_bathy
    
    f90wrap_c_bathy = m_model_c_bathy
end subroutine f90wrap_m_model__get__c_bathy

subroutine f90wrap_m_model__set__c_bathy(f90wrap_c_bathy)
    use m_model, only: m_model_c_bathy => c_bathy
    implicit none
    integer(4), intent(in) :: f90wrap_c_bathy
    
    m_model_c_bathy = f90wrap_c_bathy
end subroutine f90wrap_m_model__set__c_bathy

subroutine f90wrap_m_model__get__c_ic(f90wrap_c_ic)
    use m_model, only: m_model_c_ic => c_ic
    implicit none
    integer(4), intent(out) :: f90wrap_c_ic
    
    f90wrap_c_ic = m_model_c_ic
end subroutine f90wrap_m_model__get__c_ic

subroutine f90wrap_m_model__set__c_ic(f90wrap_c_ic)
    use m_model, only: m_model_c_ic => c_ic
    implicit none
    integer(4), intent(in) :: f90wrap_c_ic
    
    m_model_c_ic = f90wrap_c_ic
end subroutine f90wrap_m_model__set__c_ic

subroutine f90wrap_m_model__get__c_hydrograph(f90wrap_c_hydrograph)
    use m_model, only: m_model_c_hydrograph => c_hydrograph
    implicit none
    integer(4), intent(out) :: f90wrap_c_hydrograph
    
    f90wrap_c_hydrograph = m_model_c_hydrograph
end subroutine f90wrap_m_model__get__c_hydrograph

subroutine f90wrap_m_model__set__c_hydrograph(f90wrap_c_hydrograph)
    use m_model, only: m_model_c_hydrograph => c_hydrograph
    implicit none
    integer(4), intent(in) :: f90wrap_c_hydrograph
    
    m_model_c_hydrograph = f90wrap_c_hydrograph
end subroutine f90wrap_m_model__set__c_hydrograph

subroutine f90wrap_m_model__get__c_ratcurve(f90wrap_c_ratcurve)
    use m_model, only: m_model_c_ratcurve => c_ratcurve
    implicit none
    integer(4), intent(out) :: f90wrap_c_ratcurve
    
    f90wrap_c_ratcurve = m_model_c_ratcurve
end subroutine f90wrap_m_model__get__c_ratcurve

subroutine f90wrap_m_model__set__c_ratcurve(f90wrap_c_ratcurve)
    use m_model, only: m_model_c_ratcurve => c_ratcurve
    implicit none
    integer(4), intent(in) :: f90wrap_c_ratcurve
    
    m_model_c_ratcurve = f90wrap_c_ratcurve
end subroutine f90wrap_m_model__set__c_ratcurve

subroutine f90wrap_m_model__get__c_rain(f90wrap_c_rain)
    use m_model, only: m_model_c_rain => c_rain
    implicit none
    integer(4), intent(out) :: f90wrap_c_rain
    
    f90wrap_c_rain = m_model_c_rain
end subroutine f90wrap_m_model__get__c_rain

subroutine f90wrap_m_model__set__c_rain(f90wrap_c_rain)
    use m_model, only: m_model_c_rain => c_rain
    implicit none
    integer(4), intent(in) :: f90wrap_c_rain
    
    m_model_c_rain = f90wrap_c_rain
end subroutine f90wrap_m_model__set__c_rain

subroutine f90wrap_m_model__get__c_Ks(f90wrap_c_Ks)
    use m_model, only: m_model_c_Ks => c_Ks
    implicit none
    integer(4), intent(out) :: f90wrap_c_Ks
    
    f90wrap_c_Ks = m_model_c_Ks
end subroutine f90wrap_m_model__get__c_Ks

subroutine f90wrap_m_model__set__c_Ks(f90wrap_c_Ks)
    use m_model, only: m_model_c_Ks => c_Ks
    implicit none
    integer(4), intent(in) :: f90wrap_c_Ks
    
    m_model_c_Ks = f90wrap_c_Ks
end subroutine f90wrap_m_model__set__c_Ks

subroutine f90wrap_m_model__get__c_PsiF(f90wrap_c_PsiF)
    use m_model, only: m_model_c_PsiF => c_PsiF
    implicit none
    integer(4), intent(out) :: f90wrap_c_PsiF
    
    f90wrap_c_PsiF = m_model_c_PsiF
end subroutine f90wrap_m_model__get__c_PsiF

subroutine f90wrap_m_model__set__c_PsiF(f90wrap_c_PsiF)
    use m_model, only: m_model_c_PsiF => c_PsiF
    implicit none
    integer(4), intent(in) :: f90wrap_c_PsiF
    
    m_model_c_PsiF = f90wrap_c_PsiF
end subroutine f90wrap_m_model__set__c_PsiF

subroutine f90wrap_m_model__get__c_DeltaTheta(f90wrap_c_DeltaTheta)
    use m_model, only: m_model_c_DeltaTheta => c_DeltaTheta
    implicit none
    integer(4), intent(out) :: f90wrap_c_DeltaTheta
    
    f90wrap_c_DeltaTheta = m_model_c_DeltaTheta
end subroutine f90wrap_m_model__get__c_DeltaTheta

subroutine f90wrap_m_model__set__c_DeltaTheta(f90wrap_c_DeltaTheta)
    use m_model, only: m_model_c_DeltaTheta => c_DeltaTheta
    implicit none
    integer(4), intent(in) :: f90wrap_c_DeltaTheta
    
    m_model_c_DeltaTheta = f90wrap_c_DeltaTheta
end subroutine f90wrap_m_model__set__c_DeltaTheta

subroutine f90wrap_m_model__get__c_lambda(f90wrap_c_lambda)
    use m_model, only: m_model_c_lambda => c_lambda
    implicit none
    integer(4), intent(out) :: f90wrap_c_lambda
    
    f90wrap_c_lambda = m_model_c_lambda
end subroutine f90wrap_m_model__get__c_lambda

subroutine f90wrap_m_model__set__c_lambda(f90wrap_c_lambda)
    use m_model, only: m_model_c_lambda => c_lambda
    implicit none
    integer(4), intent(in) :: f90wrap_c_lambda
    
    m_model_c_lambda = f90wrap_c_lambda
end subroutine f90wrap_m_model__set__c_lambda

subroutine f90wrap_m_model__get__c_CN(f90wrap_c_CN)
    use m_model, only: m_model_c_CN => c_CN
    implicit none
    integer(4), intent(out) :: f90wrap_c_CN
    
    f90wrap_c_CN = m_model_c_CN
end subroutine f90wrap_m_model__get__c_CN

subroutine f90wrap_m_model__set__c_CN(f90wrap_c_CN)
    use m_model, only: m_model_c_CN => c_CN
    implicit none
    integer(4), intent(in) :: f90wrap_c_CN
    
    m_model_c_CN = f90wrap_c_CN
end subroutine f90wrap_m_model__set__c_CN

subroutine f90wrap_m_model__get__eps_manning(f90wrap_eps_manning)
    use m_model, only: m_model_eps_manning => eps_manning
    implicit none
    real(8), intent(out) :: f90wrap_eps_manning
    
    f90wrap_eps_manning = m_model_eps_manning
end subroutine f90wrap_m_model__get__eps_manning

subroutine f90wrap_m_model__set__eps_manning(f90wrap_eps_manning)
    use m_model, only: m_model_eps_manning => eps_manning
    implicit none
    real(8), intent(in) :: f90wrap_eps_manning
    
    m_model_eps_manning = f90wrap_eps_manning
end subroutine f90wrap_m_model__set__eps_manning

subroutine f90wrap_m_model__get__eps_bathy(f90wrap_eps_bathy)
    use m_model, only: m_model_eps_bathy => eps_bathy
    implicit none
    real(8), intent(out) :: f90wrap_eps_bathy
    
    f90wrap_eps_bathy = m_model_eps_bathy
end subroutine f90wrap_m_model__get__eps_bathy

subroutine f90wrap_m_model__set__eps_bathy(f90wrap_eps_bathy)
    use m_model, only: m_model_eps_bathy => eps_bathy
    implicit none
    real(8), intent(in) :: f90wrap_eps_bathy
    
    m_model_eps_bathy = f90wrap_eps_bathy
end subroutine f90wrap_m_model__set__eps_bathy

subroutine f90wrap_m_model__get__eps_ic(f90wrap_eps_ic)
    use m_model, only: m_model_eps_ic => eps_ic
    implicit none
    real(8), intent(out) :: f90wrap_eps_ic
    
    f90wrap_eps_ic = m_model_eps_ic
end subroutine f90wrap_m_model__get__eps_ic

subroutine f90wrap_m_model__set__eps_ic(f90wrap_eps_ic)
    use m_model, only: m_model_eps_ic => eps_ic
    implicit none
    real(8), intent(in) :: f90wrap_eps_ic
    
    m_model_eps_ic = f90wrap_eps_ic
end subroutine f90wrap_m_model__set__eps_ic

subroutine f90wrap_m_model__get__eps_hydrograph(f90wrap_eps_hydrograph)
    use m_model, only: m_model_eps_hydrograph => eps_hydrograph
    implicit none
    real(8), intent(out) :: f90wrap_eps_hydrograph
    
    f90wrap_eps_hydrograph = m_model_eps_hydrograph
end subroutine f90wrap_m_model__get__eps_hydrograph

subroutine f90wrap_m_model__set__eps_hydrograph(f90wrap_eps_hydrograph)
    use m_model, only: m_model_eps_hydrograph => eps_hydrograph
    implicit none
    real(8), intent(in) :: f90wrap_eps_hydrograph
    
    m_model_eps_hydrograph = f90wrap_eps_hydrograph
end subroutine f90wrap_m_model__set__eps_hydrograph

subroutine f90wrap_m_model__get__eps_ratcurve(f90wrap_eps_ratcurve)
    use m_model, only: m_model_eps_ratcurve => eps_ratcurve
    implicit none
    real(8), intent(out) :: f90wrap_eps_ratcurve
    
    f90wrap_eps_ratcurve = m_model_eps_ratcurve
end subroutine f90wrap_m_model__get__eps_ratcurve

subroutine f90wrap_m_model__set__eps_ratcurve(f90wrap_eps_ratcurve)
    use m_model, only: m_model_eps_ratcurve => eps_ratcurve
    implicit none
    real(8), intent(in) :: f90wrap_eps_ratcurve
    
    m_model_eps_ratcurve = f90wrap_eps_ratcurve
end subroutine f90wrap_m_model__set__eps_ratcurve

subroutine f90wrap_m_model__get__eps_rain(f90wrap_eps_rain)
    use m_model, only: m_model_eps_rain => eps_rain
    implicit none
    real(8), intent(out) :: f90wrap_eps_rain
    
    f90wrap_eps_rain = m_model_eps_rain
end subroutine f90wrap_m_model__get__eps_rain

subroutine f90wrap_m_model__set__eps_rain(f90wrap_eps_rain)
    use m_model, only: m_model_eps_rain => eps_rain
    implicit none
    real(8), intent(in) :: f90wrap_eps_rain
    
    m_model_eps_rain = f90wrap_eps_rain
end subroutine f90wrap_m_model__set__eps_rain

subroutine f90wrap_m_model__get__eps_Ks(f90wrap_eps_Ks)
    use m_model, only: m_model_eps_Ks => eps_Ks
    implicit none
    real(8), intent(out) :: f90wrap_eps_Ks
    
    f90wrap_eps_Ks = m_model_eps_Ks
end subroutine f90wrap_m_model__get__eps_Ks

subroutine f90wrap_m_model__set__eps_Ks(f90wrap_eps_Ks)
    use m_model, only: m_model_eps_Ks => eps_Ks
    implicit none
    real(8), intent(in) :: f90wrap_eps_Ks
    
    m_model_eps_Ks = f90wrap_eps_Ks
end subroutine f90wrap_m_model__set__eps_Ks

subroutine f90wrap_m_model__get__eps_PsiF(f90wrap_eps_PsiF)
    use m_model, only: m_model_eps_PsiF => eps_PsiF
    implicit none
    real(8), intent(out) :: f90wrap_eps_PsiF
    
    f90wrap_eps_PsiF = m_model_eps_PsiF
end subroutine f90wrap_m_model__get__eps_PsiF

subroutine f90wrap_m_model__set__eps_PsiF(f90wrap_eps_PsiF)
    use m_model, only: m_model_eps_PsiF => eps_PsiF
    implicit none
    real(8), intent(in) :: f90wrap_eps_PsiF
    
    m_model_eps_PsiF = f90wrap_eps_PsiF
end subroutine f90wrap_m_model__set__eps_PsiF

subroutine f90wrap_m_model__get__eps_DeltaTheta(f90wrap_eps_DeltaTheta)
    use m_model, only: m_model_eps_DeltaTheta => eps_DeltaTheta
    implicit none
    real(8), intent(out) :: f90wrap_eps_DeltaTheta
    
    f90wrap_eps_DeltaTheta = m_model_eps_DeltaTheta
end subroutine f90wrap_m_model__get__eps_DeltaTheta

subroutine f90wrap_m_model__set__eps_DeltaTheta(f90wrap_eps_DeltaTheta)
    use m_model, only: m_model_eps_DeltaTheta => eps_DeltaTheta
    implicit none
    real(8), intent(in) :: f90wrap_eps_DeltaTheta
    
    m_model_eps_DeltaTheta = f90wrap_eps_DeltaTheta
end subroutine f90wrap_m_model__set__eps_DeltaTheta

subroutine f90wrap_m_model__get__eps_lambda(f90wrap_eps_lambda)
    use m_model, only: m_model_eps_lambda => eps_lambda
    implicit none
    real(8), intent(out) :: f90wrap_eps_lambda
    
    f90wrap_eps_lambda = m_model_eps_lambda
end subroutine f90wrap_m_model__get__eps_lambda

subroutine f90wrap_m_model__set__eps_lambda(f90wrap_eps_lambda)
    use m_model, only: m_model_eps_lambda => eps_lambda
    implicit none
    real(8), intent(in) :: f90wrap_eps_lambda
    
    m_model_eps_lambda = f90wrap_eps_lambda
end subroutine f90wrap_m_model__set__eps_lambda

subroutine f90wrap_m_model__get__eps_CN(f90wrap_eps_CN)
    use m_model, only: m_model_eps_CN => eps_CN
    implicit none
    real(8), intent(out) :: f90wrap_eps_CN
    
    f90wrap_eps_CN = m_model_eps_CN
end subroutine f90wrap_m_model__get__eps_CN

subroutine f90wrap_m_model__set__eps_CN(f90wrap_eps_CN)
    use m_model, only: m_model_eps_CN => eps_CN
    implicit none
    real(8), intent(in) :: f90wrap_eps_CN
    
    m_model_eps_CN = f90wrap_eps_CN
end subroutine f90wrap_m_model__set__eps_CN

subroutine f90wrap_m_model__get__regul_manning(f90wrap_regul_manning)
    use m_model, only: m_model_regul_manning => regul_manning
    implicit none
    real(8), intent(out) :: f90wrap_regul_manning
    
    f90wrap_regul_manning = m_model_regul_manning
end subroutine f90wrap_m_model__get__regul_manning

subroutine f90wrap_m_model__set__regul_manning(f90wrap_regul_manning)
    use m_model, only: m_model_regul_manning => regul_manning
    implicit none
    real(8), intent(in) :: f90wrap_regul_manning
    
    m_model_regul_manning = f90wrap_regul_manning
end subroutine f90wrap_m_model__set__regul_manning

subroutine f90wrap_m_model__get__regul_bathy(f90wrap_regul_bathy)
    use m_model, only: m_model_regul_bathy => regul_bathy
    implicit none
    real(8), intent(out) :: f90wrap_regul_bathy
    
    f90wrap_regul_bathy = m_model_regul_bathy
end subroutine f90wrap_m_model__get__regul_bathy

subroutine f90wrap_m_model__set__regul_bathy(f90wrap_regul_bathy)
    use m_model, only: m_model_regul_bathy => regul_bathy
    implicit none
    real(8), intent(in) :: f90wrap_regul_bathy
    
    m_model_regul_bathy = f90wrap_regul_bathy
end subroutine f90wrap_m_model__set__regul_bathy

subroutine f90wrap_m_model__get__regul_ic(f90wrap_regul_ic)
    use m_model, only: m_model_regul_ic => regul_ic
    implicit none
    real(8), intent(out) :: f90wrap_regul_ic
    
    f90wrap_regul_ic = m_model_regul_ic
end subroutine f90wrap_m_model__get__regul_ic

subroutine f90wrap_m_model__set__regul_ic(f90wrap_regul_ic)
    use m_model, only: m_model_regul_ic => regul_ic
    implicit none
    real(8), intent(in) :: f90wrap_regul_ic
    
    m_model_regul_ic = f90wrap_regul_ic
end subroutine f90wrap_m_model__set__regul_ic

subroutine f90wrap_m_model__get__regul_hydrograph(f90wrap_regul_hydrograph)
    use m_model, only: m_model_regul_hydrograph => regul_hydrograph
    implicit none
    real(8), intent(out) :: f90wrap_regul_hydrograph
    
    f90wrap_regul_hydrograph = m_model_regul_hydrograph
end subroutine f90wrap_m_model__get__regul_hydrograph

subroutine f90wrap_m_model__set__regul_hydrograph(f90wrap_regul_hydrograph)
    use m_model, only: m_model_regul_hydrograph => regul_hydrograph
    implicit none
    real(8), intent(in) :: f90wrap_regul_hydrograph
    
    m_model_regul_hydrograph = f90wrap_regul_hydrograph
end subroutine f90wrap_m_model__set__regul_hydrograph

subroutine f90wrap_m_model__get__regul_ratcurve(f90wrap_regul_ratcurve)
    use m_model, only: m_model_regul_ratcurve => regul_ratcurve
    implicit none
    real(8), intent(out) :: f90wrap_regul_ratcurve
    
    f90wrap_regul_ratcurve = m_model_regul_ratcurve
end subroutine f90wrap_m_model__get__regul_ratcurve

subroutine f90wrap_m_model__set__regul_ratcurve(f90wrap_regul_ratcurve)
    use m_model, only: m_model_regul_ratcurve => regul_ratcurve
    implicit none
    real(8), intent(in) :: f90wrap_regul_ratcurve
    
    m_model_regul_ratcurve = f90wrap_regul_ratcurve
end subroutine f90wrap_m_model__set__regul_ratcurve

subroutine f90wrap_m_model__get__fix_time_step_serie(f90wrap_fix_time_step_serie)
    use m_model, only: m_model_fix_time_step_serie => fix_time_step_serie
    implicit none
    integer(4), intent(out) :: f90wrap_fix_time_step_serie
    
    f90wrap_fix_time_step_serie = m_model_fix_time_step_serie
end subroutine f90wrap_m_model__get__fix_time_step_serie

subroutine f90wrap_m_model__set__fix_time_step_serie(f90wrap_fix_time_step_serie)
    use m_model, only: m_model_fix_time_step_serie => fix_time_step_serie
    implicit none
    integer(4), intent(in) :: f90wrap_fix_time_step_serie
    
    m_model_fix_time_step_serie = f90wrap_fix_time_step_serie
end subroutine f90wrap_m_model__set__fix_time_step_serie

! End of module m_model defined in file m_sw_mono.f90

