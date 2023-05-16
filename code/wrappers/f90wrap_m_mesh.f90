! Module m_mesh defined in file m_mesh.f90

subroutine f90wrap_nodetype__array__cell(this, nd, dtype, dshape, dloc)
    use m_mesh, only: nodetype
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(nodetype_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%cell)) then
        dshape(1:1) = shape(this_ptr%p%cell)
        dloc = loc(this_ptr%p%cell)
    else
        dloc = 0
    end if
end subroutine f90wrap_nodetype__array__cell

subroutine f90wrap_nodetype__array__edge(this, nd, dtype, dshape, dloc)
    use m_mesh, only: nodetype
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(nodetype_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%edge)) then
        dshape(1:1) = shape(this_ptr%p%edge)
        dloc = loc(this_ptr%p%edge)
    else
        dloc = 0
    end if
end subroutine f90wrap_nodetype__array__edge

subroutine f90wrap_nodetype__get__coord(this, f90wrap_coord)
    use m_mesh, only: nodetype
    use m_linear_algebra, only: vec2d
    implicit none
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(nodetype_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_coord(2)
    type(vec2d_ptr_type) :: coord_ptr
    
    this_ptr = transfer(this, this_ptr)
    coord_ptr%p => this_ptr%p%coord
    f90wrap_coord = transfer(coord_ptr,f90wrap_coord)
end subroutine f90wrap_nodetype__get__coord

subroutine f90wrap_nodetype__set__coord(this, f90wrap_coord)
    use m_mesh, only: nodetype
    use m_linear_algebra, only: vec2d
    implicit none
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(nodetype_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_coord(2)
    type(vec2d_ptr_type) :: coord_ptr
    
    this_ptr = transfer(this, this_ptr)
    coord_ptr = transfer(f90wrap_coord,coord_ptr)
    this_ptr%p%coord = coord_ptr%p
end subroutine f90wrap_nodetype__set__coord

subroutine f90wrap_nodetype__get__boundary(this, f90wrap_boundary)
    use m_mesh, only: nodetype
    implicit none
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    integer, intent(in)   :: this(2)
    type(nodetype_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_boundary
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_boundary = this_ptr%p%boundary
end subroutine f90wrap_nodetype__get__boundary

subroutine f90wrap_nodetype__set__boundary(this, f90wrap_boundary)
    use m_mesh, only: nodetype
    implicit none
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    integer, intent(in)   :: this(2)
    type(nodetype_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_boundary
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%boundary = f90wrap_boundary
end subroutine f90wrap_nodetype__set__boundary

subroutine f90wrap_nodetype__get__lim(this, f90wrap_lim)
    use m_mesh, only: nodetype
    implicit none
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    integer, intent(in)   :: this(2)
    type(nodetype_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_lim
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_lim = this_ptr%p%lim
end subroutine f90wrap_nodetype__get__lim

subroutine f90wrap_nodetype__set__lim(this, f90wrap_lim)
    use m_mesh, only: nodetype
    implicit none
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    integer, intent(in)   :: this(2)
    type(nodetype_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_lim
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lim = f90wrap_lim
end subroutine f90wrap_nodetype__set__lim

subroutine f90wrap_nodetype_initialise(node, dim1, dim2)
    use m_mesh, only: nodetype_initialise, nodetype
    implicit none
    
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    type(nodetype_ptr_type) :: node_ptr
    integer, intent(out), dimension(2) :: node
    integer, intent(in) :: dim1
    integer, intent(in) :: dim2
    allocate(node_ptr%p)
    call nodetype_initialise(node=node_ptr%p, dim1=dim1, dim2=dim2)
    node = transfer(node_ptr, node)
end subroutine f90wrap_nodetype_initialise

subroutine f90wrap_nodetype_finalise(this)
    use m_mesh, only: nodetype
    implicit none
    
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    type(nodetype_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_nodetype_finalise

subroutine f90wrap_nodetypelim__get__ind(this, f90wrap_ind)
    use m_mesh, only: nodetypelim
    implicit none
    type nodetypelim_ptr_type
        type(nodetypelim), pointer :: p => NULL()
    end type nodetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(nodetypelim_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ind
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ind = this_ptr%p%ind
end subroutine f90wrap_nodetypelim__get__ind

subroutine f90wrap_nodetypelim__set__ind(this, f90wrap_ind)
    use m_mesh, only: nodetypelim
    implicit none
    type nodetypelim_ptr_type
        type(nodetypelim), pointer :: p => NULL()
    end type nodetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(nodetypelim_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ind
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ind = f90wrap_ind
end subroutine f90wrap_nodetypelim__set__ind

subroutine f90wrap_nodetypelim__get__typlim(this, f90wrap_typlim)
    use m_mesh, only: nodetypelim
    implicit none
    type nodetypelim_ptr_type
        type(nodetypelim), pointer :: p => NULL()
    end type nodetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(nodetypelim_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_typlim
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_typlim = this_ptr%p%typlim
end subroutine f90wrap_nodetypelim__get__typlim

subroutine f90wrap_nodetypelim__set__typlim(this, f90wrap_typlim)
    use m_mesh, only: nodetypelim
    implicit none
    type nodetypelim_ptr_type
        type(nodetypelim), pointer :: p => NULL()
    end type nodetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(nodetypelim_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_typlim
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%typlim = f90wrap_typlim
end subroutine f90wrap_nodetypelim__set__typlim

subroutine f90wrap_nodetypelim__get__group(this, f90wrap_group)
    use m_mesh, only: nodetypelim
    implicit none
    type nodetypelim_ptr_type
        type(nodetypelim), pointer :: p => NULL()
    end type nodetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(nodetypelim_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_group = this_ptr%p%group
end subroutine f90wrap_nodetypelim__get__group

subroutine f90wrap_nodetypelim__set__group(this, f90wrap_group)
    use m_mesh, only: nodetypelim
    implicit none
    type nodetypelim_ptr_type
        type(nodetypelim), pointer :: p => NULL()
    end type nodetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(nodetypelim_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%group = f90wrap_group
end subroutine f90wrap_nodetypelim__set__group

subroutine f90wrap_nodetypelim_initialise(nodelim, ind, group)
    use m_mesh, only: nodetypelim, nodetypelim_initialise
    implicit none
    
    type nodetypelim_ptr_type
        type(nodetypelim), pointer :: p => NULL()
    end type nodetypelim_ptr_type
    type(nodetypelim_ptr_type) :: nodelim_ptr
    integer, intent(out), dimension(2) :: nodelim
    integer, intent(in) :: ind
    integer, intent(in) :: group
    allocate(nodelim_ptr%p)
    call nodetypelim_initialise(nodelim=nodelim_ptr%p, ind=ind, group=group)
    nodelim = transfer(nodelim_ptr, nodelim)
end subroutine f90wrap_nodetypelim_initialise

subroutine f90wrap_nodetypelim_finalise(this)
    use m_mesh, only: nodetypelim
    implicit none
    
    type nodetypelim_ptr_type
        type(nodetypelim), pointer :: p => NULL()
    end type nodetypelim_ptr_type
    type(nodetypelim_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_nodetypelim_finalise

subroutine f90wrap_celltype__array__node(this, nd, dtype, dshape, dloc)
    use m_mesh, only: celltype
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(celltype_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%node)
    dloc = loc(this_ptr%p%node)
end subroutine f90wrap_celltype__array__node

subroutine f90wrap_celltype__array__cell(this, nd, dtype, dshape, dloc)
    use m_mesh, only: celltype
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(celltype_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%cell)
    dloc = loc(this_ptr%p%cell)
end subroutine f90wrap_celltype__array__cell

subroutine f90wrap_celltype__array__edge(this, nd, dtype, dshape, dloc)
    use m_mesh, only: celltype
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(celltype_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%edge)
    dloc = loc(this_ptr%p%edge)
end subroutine f90wrap_celltype__array__edge

subroutine f90wrap_celltype__get__nbed(this, f90wrap_nbed)
    use m_mesh, only: celltype
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nbed
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nbed = this_ptr%p%nbed
end subroutine f90wrap_celltype__get__nbed

subroutine f90wrap_celltype__set__nbed(this, f90wrap_nbed)
    use m_mesh, only: celltype
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nbed
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nbed = f90wrap_nbed
end subroutine f90wrap_celltype__set__nbed

subroutine f90wrap_celltype__get__boundary(this, f90wrap_boundary)
    use m_mesh, only: celltype
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_boundary
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_boundary = this_ptr%p%boundary
end subroutine f90wrap_celltype__get__boundary

subroutine f90wrap_celltype__set__boundary(this, f90wrap_boundary)
    use m_mesh, only: celltype
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_boundary
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%boundary = f90wrap_boundary
end subroutine f90wrap_celltype__set__boundary

subroutine f90wrap_celltype__get__surf(this, f90wrap_surf)
    use m_mesh, only: celltype
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_surf
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_surf = this_ptr%p%surf
end subroutine f90wrap_celltype__get__surf

subroutine f90wrap_celltype__set__surf(this, f90wrap_surf)
    use m_mesh, only: celltype
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_surf
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%surf = f90wrap_surf
end subroutine f90wrap_celltype__set__surf

subroutine f90wrap_celltype__get__invsurf(this, f90wrap_invsurf)
    use m_mesh, only: celltype
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_invsurf
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_invsurf = this_ptr%p%invsurf
end subroutine f90wrap_celltype__get__invsurf

subroutine f90wrap_celltype__set__invsurf(this, f90wrap_invsurf)
    use m_mesh, only: celltype
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_invsurf
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%invsurf = f90wrap_invsurf
end subroutine f90wrap_celltype__set__invsurf

subroutine f90wrap_celltype__get__peri(this, f90wrap_peri)
    use m_mesh, only: celltype
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_peri
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_peri = this_ptr%p%peri
end subroutine f90wrap_celltype__get__peri

subroutine f90wrap_celltype__set__peri(this, f90wrap_peri)
    use m_mesh, only: celltype
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_peri
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%peri = f90wrap_peri
end subroutine f90wrap_celltype__set__peri

subroutine f90wrap_celltype__get__grav(this, f90wrap_grav)
    use m_mesh, only: celltype
    use m_linear_algebra, only: vec2d
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_grav(2)
    type(vec2d_ptr_type) :: grav_ptr
    
    this_ptr = transfer(this, this_ptr)
    grav_ptr%p => this_ptr%p%grav
    f90wrap_grav = transfer(grav_ptr,f90wrap_grav)
end subroutine f90wrap_celltype__get__grav

subroutine f90wrap_celltype__set__grav(this, f90wrap_grav)
    use m_mesh, only: celltype
    use m_linear_algebra, only: vec2d
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_grav(2)
    type(vec2d_ptr_type) :: grav_ptr
    
    this_ptr = transfer(this, this_ptr)
    grav_ptr = transfer(f90wrap_grav,grav_ptr)
    this_ptr%p%grav = grav_ptr%p
end subroutine f90wrap_celltype__set__grav

subroutine f90wrap_celltype__get__rain(this, f90wrap_rain)
    use m_mesh, only: celltype
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_rain
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_rain = this_ptr%p%rain
end subroutine f90wrap_celltype__get__rain

subroutine f90wrap_celltype__set__rain(this, f90wrap_rain)
    use m_mesh, only: celltype
    implicit none
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in)   :: this(2)
    type(celltype_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_rain
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%rain = f90wrap_rain
end subroutine f90wrap_celltype__set__rain

subroutine f90wrap_celltype_initialise(this)
    use m_mesh, only: celltype
    implicit none
    
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    type(celltype_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_celltype_initialise

subroutine f90wrap_celltype_finalise(this)
    use m_mesh, only: celltype
    implicit none
    
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    type(celltype_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_celltype_finalise

subroutine f90wrap_celltypelim__get__ind(this, f90wrap_ind)
    use m_mesh, only: celltypelim
    implicit none
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(celltypelim_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ind
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ind = this_ptr%p%ind
end subroutine f90wrap_celltypelim__get__ind

subroutine f90wrap_celltypelim__set__ind(this, f90wrap_ind)
    use m_mesh, only: celltypelim
    implicit none
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(celltypelim_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ind
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ind = f90wrap_ind
end subroutine f90wrap_celltypelim__set__ind

subroutine f90wrap_celltypelim__get__typlim(this, f90wrap_typlim)
    use m_mesh, only: celltypelim
    implicit none
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(celltypelim_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_typlim
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_typlim = this_ptr%p%typlim
end subroutine f90wrap_celltypelim__get__typlim

subroutine f90wrap_celltypelim__set__typlim(this, f90wrap_typlim)
    use m_mesh, only: celltypelim
    implicit none
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(celltypelim_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_typlim
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%typlim = f90wrap_typlim
end subroutine f90wrap_celltypelim__set__typlim

subroutine f90wrap_celltypelim__get__group(this, f90wrap_group)
    use m_mesh, only: celltypelim
    implicit none
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(celltypelim_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_group = this_ptr%p%group
end subroutine f90wrap_celltypelim__get__group

subroutine f90wrap_celltypelim__set__group(this, f90wrap_group)
    use m_mesh, only: celltypelim
    implicit none
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(celltypelim_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%group = f90wrap_group
end subroutine f90wrap_celltypelim__set__group

subroutine f90wrap_celltypelim__get__cell(this, f90wrap_cell)
    use m_mesh, only: celltypelim
    implicit none
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(celltypelim_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_cell
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_cell = this_ptr%p%cell
end subroutine f90wrap_celltypelim__get__cell

subroutine f90wrap_celltypelim__set__cell(this, f90wrap_cell)
    use m_mesh, only: celltypelim
    implicit none
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(celltypelim_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_cell
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%cell = f90wrap_cell
end subroutine f90wrap_celltypelim__set__cell

subroutine f90wrap_celltypelim__get__grav(this, f90wrap_grav)
    use m_mesh, only: celltypelim
    use m_linear_algebra, only: vec2d
    implicit none
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(celltypelim_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_grav(2)
    type(vec2d_ptr_type) :: grav_ptr
    
    this_ptr = transfer(this, this_ptr)
    grav_ptr%p => this_ptr%p%grav
    f90wrap_grav = transfer(grav_ptr,f90wrap_grav)
end subroutine f90wrap_celltypelim__get__grav

subroutine f90wrap_celltypelim__set__grav(this, f90wrap_grav)
    use m_mesh, only: celltypelim
    use m_linear_algebra, only: vec2d
    implicit none
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(celltypelim_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_grav(2)
    type(vec2d_ptr_type) :: grav_ptr
    
    this_ptr = transfer(this, this_ptr)
    grav_ptr = transfer(f90wrap_grav,grav_ptr)
    this_ptr%p%grav = grav_ptr%p
end subroutine f90wrap_celltypelim__set__grav

subroutine f90wrap_celltypelim_initialise(this)
    use m_mesh, only: celltypelim
    implicit none
    
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    type(celltypelim_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_celltypelim_initialise

subroutine f90wrap_celltypelim_finalise(this)
    use m_mesh, only: celltypelim
    implicit none
    
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    type(celltypelim_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_celltypelim_finalise

subroutine f90wrap_edgetype__array__node(this, nd, dtype, dshape, dloc)
    use m_mesh, only: edgetype
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%node)
    dloc = loc(this_ptr%p%node)
end subroutine f90wrap_edgetype__array__node

subroutine f90wrap_edgetype__array__cell(this, nd, dtype, dshape, dloc)
    use m_mesh, only: edgetype
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%cell)
    dloc = loc(this_ptr%p%cell)
end subroutine f90wrap_edgetype__array__cell

subroutine f90wrap_edgetype__get__cell1D2D(this, f90wrap_cell1D2D)
    use m_mesh, only: edgetype
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_cell1D2D
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_cell1D2D = this_ptr%p%cell1D2D
end subroutine f90wrap_edgetype__get__cell1D2D

subroutine f90wrap_edgetype__set__cell1D2D(this, f90wrap_cell1D2D)
    use m_mesh, only: edgetype
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_cell1D2D
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%cell1D2D = f90wrap_cell1D2D
end subroutine f90wrap_edgetype__set__cell1D2D

subroutine f90wrap_edgetype__get__boundary(this, f90wrap_boundary)
    use m_mesh, only: edgetype
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_boundary
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_boundary = this_ptr%p%boundary
end subroutine f90wrap_edgetype__get__boundary

subroutine f90wrap_edgetype__set__boundary(this, f90wrap_boundary)
    use m_mesh, only: edgetype
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_boundary
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%boundary = f90wrap_boundary
end subroutine f90wrap_edgetype__set__boundary

subroutine f90wrap_edgetype__get__subdomain(this, f90wrap_subdomain)
    use m_mesh, only: edgetype
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_subdomain
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_subdomain = this_ptr%p%subdomain
end subroutine f90wrap_edgetype__get__subdomain

subroutine f90wrap_edgetype__set__subdomain(this, f90wrap_subdomain)
    use m_mesh, only: edgetype
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_subdomain
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%subdomain = f90wrap_subdomain
end subroutine f90wrap_edgetype__set__subdomain

subroutine f90wrap_edgetype__get__lim(this, f90wrap_lim)
    use m_mesh, only: edgetype
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_lim
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_lim = this_ptr%p%lim
end subroutine f90wrap_edgetype__get__lim

subroutine f90wrap_edgetype__set__lim(this, f90wrap_lim)
    use m_mesh, only: edgetype
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_lim
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lim = f90wrap_lim
end subroutine f90wrap_edgetype__set__lim

subroutine f90wrap_edgetype__get__length(this, f90wrap_length)
    use m_mesh, only: edgetype
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_length
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_length = this_ptr%p%length
end subroutine f90wrap_edgetype__get__length

subroutine f90wrap_edgetype__set__length(this, f90wrap_length)
    use m_mesh, only: edgetype
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_length
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%length = f90wrap_length
end subroutine f90wrap_edgetype__set__length

subroutine f90wrap_edgetype__get__center(this, f90wrap_center)
    use m_mesh, only: edgetype
    use m_linear_algebra, only: vec2d
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_center(2)
    type(vec2d_ptr_type) :: center_ptr
    
    this_ptr = transfer(this, this_ptr)
    center_ptr%p => this_ptr%p%center
    f90wrap_center = transfer(center_ptr,f90wrap_center)
end subroutine f90wrap_edgetype__get__center

subroutine f90wrap_edgetype__set__center(this, f90wrap_center)
    use m_mesh, only: edgetype
    use m_linear_algebra, only: vec2d
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_center(2)
    type(vec2d_ptr_type) :: center_ptr
    
    this_ptr = transfer(this, this_ptr)
    center_ptr = transfer(f90wrap_center,center_ptr)
    this_ptr%p%center = center_ptr%p
end subroutine f90wrap_edgetype__set__center

subroutine f90wrap_edgetype__get__normal(this, f90wrap_normal)
    use m_mesh, only: edgetype
    use m_linear_algebra, only: vec2d
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_normal(2)
    type(vec2d_ptr_type) :: normal_ptr
    
    this_ptr = transfer(this, this_ptr)
    normal_ptr%p => this_ptr%p%normal
    f90wrap_normal = transfer(normal_ptr,f90wrap_normal)
end subroutine f90wrap_edgetype__get__normal

subroutine f90wrap_edgetype__set__normal(this, f90wrap_normal)
    use m_mesh, only: edgetype
    use m_linear_algebra, only: vec2d
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_normal(2)
    type(vec2d_ptr_type) :: normal_ptr
    
    this_ptr = transfer(this, this_ptr)
    normal_ptr = transfer(f90wrap_normal,normal_ptr)
    this_ptr%p%normal = normal_ptr%p
end subroutine f90wrap_edgetype__set__normal

subroutine f90wrap_edgetype__get__tangent(this, f90wrap_tangent)
    use m_mesh, only: edgetype
    use m_linear_algebra, only: vec2d
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_tangent(2)
    type(vec2d_ptr_type) :: tangent_ptr
    
    this_ptr = transfer(this, this_ptr)
    tangent_ptr%p => this_ptr%p%tangent
    f90wrap_tangent = transfer(tangent_ptr,f90wrap_tangent)
end subroutine f90wrap_edgetype__get__tangent

subroutine f90wrap_edgetype__set__tangent(this, f90wrap_tangent)
    use m_mesh, only: edgetype
    use m_linear_algebra, only: vec2d
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_tangent(2)
    type(vec2d_ptr_type) :: tangent_ptr
    
    this_ptr = transfer(this, this_ptr)
    tangent_ptr = transfer(f90wrap_tangent,tangent_ptr)
    this_ptr%p%tangent = tangent_ptr%p
end subroutine f90wrap_edgetype__set__tangent

subroutine f90wrap_edgetype__get__vcell(this, f90wrap_vcell)
    use m_mesh, only: edgetype
    use m_linear_algebra, only: vec2d
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_vcell(2)
    type(vec2d_ptr_type) :: vcell_ptr
    
    this_ptr = transfer(this, this_ptr)
    vcell_ptr%p => this_ptr%p%vcell
    f90wrap_vcell = transfer(vcell_ptr,f90wrap_vcell)
end subroutine f90wrap_edgetype__get__vcell

subroutine f90wrap_edgetype__set__vcell(this, f90wrap_vcell)
    use m_mesh, only: edgetype
    use m_linear_algebra, only: vec2d
    implicit none
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_vcell(2)
    type(vec2d_ptr_type) :: vcell_ptr
    
    this_ptr = transfer(this, this_ptr)
    vcell_ptr = transfer(f90wrap_vcell,vcell_ptr)
    this_ptr%p%vcell = vcell_ptr%p
end subroutine f90wrap_edgetype__set__vcell

subroutine f90wrap_edgetype__array_getitem__v_edge_cell(f90wrap_this, f90wrap_i, v_edge_cellitem)
    
    use m_mesh, only: edgetype
    use m_linear_algebra, only: vec2d
    implicit none
    
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: v_edge_cellitem(2)
    type(vec2d_ptr_type) :: v_edge_cell_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%v_edge_cell)) then
        call f90wrap_abort("array index out of range")
    else
        v_edge_cell_ptr%p => this_ptr%p%v_edge_cell(f90wrap_i)
        v_edge_cellitem = transfer(v_edge_cell_ptr,v_edge_cellitem)
    endif
end subroutine f90wrap_edgetype__array_getitem__v_edge_cell

subroutine f90wrap_edgetype__array_setitem__v_edge_cell(f90wrap_this, f90wrap_i, v_edge_cellitem)
    
    use m_mesh, only: edgetype
    use m_linear_algebra, only: vec2d
    implicit none
    
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(edgetype_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: v_edge_cellitem(2)
    type(vec2d_ptr_type) :: v_edge_cell_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%v_edge_cell)) then
        call f90wrap_abort("array index out of range")
    else
        v_edge_cell_ptr = transfer(v_edge_cellitem,v_edge_cell_ptr)
        this_ptr%p%v_edge_cell(f90wrap_i) = v_edge_cell_ptr%p
    endif
end subroutine f90wrap_edgetype__array_setitem__v_edge_cell

subroutine f90wrap_edgetype__array_len__v_edge_cell(f90wrap_this, f90wrap_n)
    
    use m_mesh, only: edgetype
    use m_linear_algebra, only: vec2d
    implicit none
    
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(edgetype_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    f90wrap_n = size(this_ptr%p%v_edge_cell)
end subroutine f90wrap_edgetype__array_len__v_edge_cell

subroutine f90wrap_edgetype_initialise(this)
    use m_mesh, only: edgetype
    implicit none
    
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type(edgetype_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_edgetype_initialise

subroutine f90wrap_edgetype_finalise(this)
    use m_mesh, only: edgetype
    implicit none
    
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    type(edgetype_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_edgetype_finalise

subroutine f90wrap_edgetypelim__get__ind(this, f90wrap_ind)
    use m_mesh, only: edgetypelim
    implicit none
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetypelim_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ind
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ind = this_ptr%p%ind
end subroutine f90wrap_edgetypelim__get__ind

subroutine f90wrap_edgetypelim__set__ind(this, f90wrap_ind)
    use m_mesh, only: edgetypelim
    implicit none
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetypelim_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ind
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ind = f90wrap_ind
end subroutine f90wrap_edgetypelim__set__ind

subroutine f90wrap_edgetypelim__get__typlim(this, f90wrap_typlim)
    use m_mesh, only: edgetypelim
    implicit none
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetypelim_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_typlim
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_typlim = this_ptr%p%typlim
end subroutine f90wrap_edgetypelim__get__typlim

subroutine f90wrap_edgetypelim__set__typlim(this, f90wrap_typlim)
    use m_mesh, only: edgetypelim
    implicit none
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetypelim_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_typlim
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%typlim = f90wrap_typlim
end subroutine f90wrap_edgetypelim__set__typlim

subroutine f90wrap_edgetypelim__get__group(this, f90wrap_group)
    use m_mesh, only: edgetypelim
    implicit none
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetypelim_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_group = this_ptr%p%group
end subroutine f90wrap_edgetypelim__get__group

subroutine f90wrap_edgetypelim__set__group(this, f90wrap_group)
    use m_mesh, only: edgetypelim
    implicit none
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetypelim_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_group
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%group = f90wrap_group
end subroutine f90wrap_edgetypelim__set__group

subroutine f90wrap_edgetypelim__get__perio(this, f90wrap_perio)
    use m_mesh, only: edgetypelim
    implicit none
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetypelim_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_perio
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_perio = this_ptr%p%perio
end subroutine f90wrap_edgetypelim__get__perio

subroutine f90wrap_edgetypelim__set__perio(this, f90wrap_perio)
    use m_mesh, only: edgetypelim
    implicit none
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    integer, intent(in)   :: this(2)
    type(edgetypelim_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_perio
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%perio = f90wrap_perio
end subroutine f90wrap_edgetypelim__set__perio

subroutine f90wrap_edgetypelim_initialise(this)
    use m_mesh, only: edgetypelim
    implicit none
    
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    type(edgetypelim_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_edgetypelim_initialise

subroutine f90wrap_edgetypelim_finalise(this)
    use m_mesh, only: edgetypelim
    implicit none
    
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    type(edgetypelim_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_edgetypelim_finalise

subroutine f90wrap_msh__get__nn(this, f90wrap_nn)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nn
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nn = this_ptr%p%nn
end subroutine f90wrap_msh__get__nn

subroutine f90wrap_msh__set__nn(this, f90wrap_nn)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nn
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nn = f90wrap_nn
end subroutine f90wrap_msh__set__nn

subroutine f90wrap_msh__get__nnb(this, f90wrap_nnb)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nnb
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nnb = this_ptr%p%nnb
end subroutine f90wrap_msh__get__nnb

subroutine f90wrap_msh__set__nnb(this, f90wrap_nnb)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nnb
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nnb = f90wrap_nnb
end subroutine f90wrap_msh__set__nnb

subroutine f90wrap_msh__get__nc(this, f90wrap_nc)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_nc
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nc = this_ptr%p%nc
end subroutine f90wrap_msh__get__nc

subroutine f90wrap_msh__set__nc(this, f90wrap_nc)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_nc
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nc = f90wrap_nc
end subroutine f90wrap_msh__set__nc

subroutine f90wrap_msh__get__ncb(this, f90wrap_ncb)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ncb
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ncb = this_ptr%p%ncb
end subroutine f90wrap_msh__get__ncb

subroutine f90wrap_msh__set__ncb(this, f90wrap_ncb)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ncb
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ncb = f90wrap_ncb
end subroutine f90wrap_msh__set__ncb

subroutine f90wrap_msh__get__ne(this, f90wrap_ne)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_ne
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_ne = this_ptr%p%ne
end subroutine f90wrap_msh__get__ne

subroutine f90wrap_msh__set__ne(this, f90wrap_ne)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_ne
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%ne = f90wrap_ne
end subroutine f90wrap_msh__set__ne

subroutine f90wrap_msh__get__neb(this, f90wrap_neb)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_neb
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_neb = this_ptr%p%neb
end subroutine f90wrap_msh__get__neb

subroutine f90wrap_msh__set__neb(this, f90wrap_neb)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_neb
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%neb = f90wrap_neb
end subroutine f90wrap_msh__set__neb

subroutine f90wrap_msh__array_getitem__node(f90wrap_this, f90wrap_i, nodeitem)
    
    use m_mesh, only: msh, nodetype
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: nodeitem(2)
    type(nodetype_ptr_type) :: node_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%node)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%node)) then
            call f90wrap_abort("array index out of range")
        else
            node_ptr%p => this_ptr%p%node(f90wrap_i)
            nodeitem = transfer(node_ptr,nodeitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_msh__array_getitem__node

subroutine f90wrap_msh__array_setitem__node(f90wrap_this, f90wrap_i, nodeitem)
    
    use m_mesh, only: msh, nodetype
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: nodeitem(2)
    type(nodetype_ptr_type) :: node_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%node)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%node)) then
            call f90wrap_abort("array index out of range")
        else
            node_ptr = transfer(nodeitem,node_ptr)
            this_ptr%p%node(f90wrap_i) = node_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_msh__array_setitem__node

subroutine f90wrap_msh__array_len__node(f90wrap_this, f90wrap_n)
    
    use m_mesh, only: msh, nodetype
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type nodetype_ptr_type
        type(nodetype), pointer :: p => NULL()
    end type nodetype_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%node)) then
        f90wrap_n = size(this_ptr%p%node)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_msh__array_len__node

subroutine f90wrap_msh__array_getitem__nodeb(f90wrap_this, f90wrap_i, nodebitem)
    
    use m_mesh, only: msh, nodetypelim
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type nodetypelim_ptr_type
        type(nodetypelim), pointer :: p => NULL()
    end type nodetypelim_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: nodebitem(2)
    type(nodetypelim_ptr_type) :: nodeb_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%nodeb)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%nodeb)) then
            call f90wrap_abort("array index out of range")
        else
            nodeb_ptr%p => this_ptr%p%nodeb(f90wrap_i)
            nodebitem = transfer(nodeb_ptr,nodebitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_msh__array_getitem__nodeb

subroutine f90wrap_msh__array_setitem__nodeb(f90wrap_this, f90wrap_i, nodebitem)
    
    use m_mesh, only: msh, nodetypelim
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type nodetypelim_ptr_type
        type(nodetypelim), pointer :: p => NULL()
    end type nodetypelim_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: nodebitem(2)
    type(nodetypelim_ptr_type) :: nodeb_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%nodeb)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%nodeb)) then
            call f90wrap_abort("array index out of range")
        else
            nodeb_ptr = transfer(nodebitem,nodeb_ptr)
            this_ptr%p%nodeb(f90wrap_i) = nodeb_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_msh__array_setitem__nodeb

subroutine f90wrap_msh__array_len__nodeb(f90wrap_this, f90wrap_n)
    
    use m_mesh, only: msh, nodetypelim
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type nodetypelim_ptr_type
        type(nodetypelim), pointer :: p => NULL()
    end type nodetypelim_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%nodeb)) then
        f90wrap_n = size(this_ptr%p%nodeb)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_msh__array_len__nodeb

subroutine f90wrap_msh__array_getitem__cell(f90wrap_this, f90wrap_i, cellitem)
    
    use m_mesh, only: msh, celltype
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: cellitem(2)
    type(celltype_ptr_type) :: cell_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%cell)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%cell)) then
            call f90wrap_abort("array index out of range")
        else
            cell_ptr%p => this_ptr%p%cell(f90wrap_i)
            cellitem = transfer(cell_ptr,cellitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_msh__array_getitem__cell

subroutine f90wrap_msh__array_setitem__cell(f90wrap_this, f90wrap_i, cellitem)
    
    use m_mesh, only: msh, celltype
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: cellitem(2)
    type(celltype_ptr_type) :: cell_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%cell)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%cell)) then
            call f90wrap_abort("array index out of range")
        else
            cell_ptr = transfer(cellitem,cell_ptr)
            this_ptr%p%cell(f90wrap_i) = cell_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_msh__array_setitem__cell

subroutine f90wrap_msh__array_len__cell(f90wrap_this, f90wrap_n)
    
    use m_mesh, only: msh, celltype
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type celltype_ptr_type
        type(celltype), pointer :: p => NULL()
    end type celltype_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%cell)) then
        f90wrap_n = size(this_ptr%p%cell)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_msh__array_len__cell

subroutine f90wrap_msh__array_getitem__cellb(f90wrap_this, f90wrap_i, cellbitem)
    
    use m_mesh, only: msh, celltypelim
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: cellbitem(2)
    type(celltypelim_ptr_type) :: cellb_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%cellb)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%cellb)) then
            call f90wrap_abort("array index out of range")
        else
            cellb_ptr%p => this_ptr%p%cellb(f90wrap_i)
            cellbitem = transfer(cellb_ptr,cellbitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_msh__array_getitem__cellb

subroutine f90wrap_msh__array_setitem__cellb(f90wrap_this, f90wrap_i, cellbitem)
    
    use m_mesh, only: msh, celltypelim
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: cellbitem(2)
    type(celltypelim_ptr_type) :: cellb_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%cellb)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%cellb)) then
            call f90wrap_abort("array index out of range")
        else
            cellb_ptr = transfer(cellbitem,cellb_ptr)
            this_ptr%p%cellb(f90wrap_i) = cellb_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_msh__array_setitem__cellb

subroutine f90wrap_msh__array_len__cellb(f90wrap_this, f90wrap_n)
    
    use m_mesh, only: msh, celltypelim
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type celltypelim_ptr_type
        type(celltypelim), pointer :: p => NULL()
    end type celltypelim_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%cellb)) then
        f90wrap_n = size(this_ptr%p%cellb)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_msh__array_len__cellb

subroutine f90wrap_msh__array_getitem__edge(f90wrap_this, f90wrap_i, edgeitem)
    
    use m_mesh, only: msh, edgetype
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: edgeitem(2)
    type(edgetype_ptr_type) :: edge_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%edge)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%edge)) then
            call f90wrap_abort("array index out of range")
        else
            edge_ptr%p => this_ptr%p%edge(f90wrap_i)
            edgeitem = transfer(edge_ptr,edgeitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_msh__array_getitem__edge

subroutine f90wrap_msh__array_setitem__edge(f90wrap_this, f90wrap_i, edgeitem)
    
    use m_mesh, only: msh, edgetype
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: edgeitem(2)
    type(edgetype_ptr_type) :: edge_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%edge)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%edge)) then
            call f90wrap_abort("array index out of range")
        else
            edge_ptr = transfer(edgeitem,edge_ptr)
            this_ptr%p%edge(f90wrap_i) = edge_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_msh__array_setitem__edge

subroutine f90wrap_msh__array_len__edge(f90wrap_this, f90wrap_n)
    
    use m_mesh, only: msh, edgetype
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type edgetype_ptr_type
        type(edgetype), pointer :: p => NULL()
    end type edgetype_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%edge)) then
        f90wrap_n = size(this_ptr%p%edge)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_msh__array_len__edge

subroutine f90wrap_msh__array_getitem__edgeb(f90wrap_this, f90wrap_i, edgebitem)
    
    use m_mesh, only: msh, edgetypelim
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: edgebitem(2)
    type(edgetypelim_ptr_type) :: edgeb_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%edgeb)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%edgeb)) then
            call f90wrap_abort("array index out of range")
        else
            edgeb_ptr%p => this_ptr%p%edgeb(f90wrap_i)
            edgebitem = transfer(edgeb_ptr,edgebitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_msh__array_getitem__edgeb

subroutine f90wrap_msh__array_setitem__edgeb(f90wrap_this, f90wrap_i, edgebitem)
    
    use m_mesh, only: msh, edgetypelim
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: edgebitem(2)
    type(edgetypelim_ptr_type) :: edgeb_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%edgeb)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%edgeb)) then
            call f90wrap_abort("array index out of range")
        else
            edgeb_ptr = transfer(edgebitem,edgeb_ptr)
            this_ptr%p%edgeb(f90wrap_i) = edgeb_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_msh__array_setitem__edgeb

subroutine f90wrap_msh__array_len__edgeb(f90wrap_this, f90wrap_n)
    
    use m_mesh, only: msh, edgetypelim
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type edgetypelim_ptr_type
        type(edgetypelim), pointer :: p => NULL()
    end type edgetypelim_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(msh_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%edgeb)) then
        f90wrap_n = size(this_ptr%p%edgeb)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_msh__array_len__edgeb

subroutine f90wrap_msh__get__file_name(this, f90wrap_file_name)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_file_name
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_file_name = this_ptr%p%file_name
end subroutine f90wrap_msh__get__file_name

subroutine f90wrap_msh__set__file_name(this, f90wrap_file_name)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_file_name
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%file_name = f90wrap_file_name
end subroutine f90wrap_msh__set__file_name

subroutine f90wrap_msh__get__scal(this, f90wrap_scal)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_scal
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_scal = this_ptr%p%scal
end subroutine f90wrap_msh__get__scal

subroutine f90wrap_msh__set__scal(this, f90wrap_scal)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_scal
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%scal = f90wrap_scal
end subroutine f90wrap_msh__set__scal

subroutine f90wrap_msh__get__surf(this, f90wrap_surf)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    real(8), intent(out) :: f90wrap_surf
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_surf = this_ptr%p%surf
end subroutine f90wrap_msh__get__surf

subroutine f90wrap_msh__set__surf(this, f90wrap_surf)
    use m_mesh, only: msh
    implicit none
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    integer, intent(in)   :: this(2)
    type(msh_ptr_type) :: this_ptr
    real(8), intent(in) :: f90wrap_surf
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%surf = f90wrap_surf
end subroutine f90wrap_msh__set__surf

subroutine f90wrap_msh_initialise(mesh)
    use m_mesh, only: msh, msh_initialise
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(out), dimension(2) :: mesh
    allocate(mesh_ptr%p)
    call msh_initialise(mesh=mesh_ptr%p)
    mesh = transfer(mesh_ptr, mesh)
end subroutine f90wrap_msh_initialise

subroutine f90wrap_msh_finalise(mesh)
    use m_mesh, only: msh, msh_finalise
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    mesh_ptr = transfer(mesh, mesh_ptr)
    call msh_finalise(mesh=mesh_ptr%p)
    deallocate(mesh_ptr%p)
end subroutine f90wrap_msh_finalise

subroutine f90wrap_point_in_mesh__get__coord(this, f90wrap_coord)
    use m_mesh, only: point_in_mesh
    use m_linear_algebra, only: vec2d
    implicit none
    type point_in_mesh_ptr_type
        type(point_in_mesh), pointer :: p => NULL()
    end type point_in_mesh_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(point_in_mesh_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_coord(2)
    type(vec2d_ptr_type) :: coord_ptr
    
    this_ptr = transfer(this, this_ptr)
    coord_ptr%p => this_ptr%p%coord
    f90wrap_coord = transfer(coord_ptr,f90wrap_coord)
end subroutine f90wrap_point_in_mesh__get__coord

subroutine f90wrap_point_in_mesh__set__coord(this, f90wrap_coord)
    use m_mesh, only: point_in_mesh
    use m_linear_algebra, only: vec2d
    implicit none
    type point_in_mesh_ptr_type
        type(point_in_mesh), pointer :: p => NULL()
    end type point_in_mesh_ptr_type
    type vec2d_ptr_type
        type(vec2d), pointer :: p => NULL()
    end type vec2d_ptr_type
    integer, intent(in)   :: this(2)
    type(point_in_mesh_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_coord(2)
    type(vec2d_ptr_type) :: coord_ptr
    
    this_ptr = transfer(this, this_ptr)
    coord_ptr = transfer(f90wrap_coord,coord_ptr)
    this_ptr%p%coord = coord_ptr%p
end subroutine f90wrap_point_in_mesh__set__coord

subroutine f90wrap_point_in_mesh__get__cell(this, f90wrap_cell)
    use m_mesh, only: point_in_mesh
    implicit none
    type point_in_mesh_ptr_type
        type(point_in_mesh), pointer :: p => NULL()
    end type point_in_mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(point_in_mesh_ptr_type) :: this_ptr
    integer(4), intent(out) :: f90wrap_cell
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_cell = this_ptr%p%cell
end subroutine f90wrap_point_in_mesh__get__cell

subroutine f90wrap_point_in_mesh__set__cell(this, f90wrap_cell)
    use m_mesh, only: point_in_mesh
    implicit none
    type point_in_mesh_ptr_type
        type(point_in_mesh), pointer :: p => NULL()
    end type point_in_mesh_ptr_type
    integer, intent(in)   :: this(2)
    type(point_in_mesh_ptr_type) :: this_ptr
    integer(4), intent(in) :: f90wrap_cell
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%cell = f90wrap_cell
end subroutine f90wrap_point_in_mesh__set__cell

subroutine f90wrap_point_in_mesh_initialise(this)
    use m_mesh, only: point_in_mesh
    implicit none
    
    type point_in_mesh_ptr_type
        type(point_in_mesh), pointer :: p => NULL()
    end type point_in_mesh_ptr_type
    type(point_in_mesh_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_point_in_mesh_initialise

subroutine f90wrap_point_in_mesh_finalise(this)
    use m_mesh, only: point_in_mesh
    implicit none
    
    type point_in_mesh_ptr_type
        type(point_in_mesh), pointer :: p => NULL()
    end type point_in_mesh_ptr_type
    type(point_in_mesh_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_point_in_mesh_finalise

subroutine f90wrap_display_mesh_cell(mesh)
    use m_mesh, only: msh, display_mesh_cell
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    mesh_ptr = transfer(mesh, mesh_ptr)
    call display_mesh_cell(mesh=mesh_ptr%p)
end subroutine f90wrap_display_mesh_cell

subroutine f90wrap_calc_cells_connectivity(mesh)
    use m_mesh, only: msh, calc_cells_connectivity
    implicit none
    
    type msh_ptr_type
        type(msh), pointer :: p => NULL()
    end type msh_ptr_type
    type(msh_ptr_type) :: mesh_ptr
    integer, intent(in), dimension(2) :: mesh
    mesh_ptr = transfer(mesh, mesh_ptr)
    call calc_cells_connectivity(mesh=mesh_ptr%p)
end subroutine f90wrap_calc_cells_connectivity

subroutine f90wrap_m_mesh__get__maxed(f90wrap_maxed)
    use m_mesh, only: m_mesh_maxed => maxed
    implicit none
    integer, intent(out) :: f90wrap_maxed
    
    f90wrap_maxed = m_mesh_maxed
end subroutine f90wrap_m_mesh__get__maxed

! End of module m_mesh defined in file m_mesh.f90

