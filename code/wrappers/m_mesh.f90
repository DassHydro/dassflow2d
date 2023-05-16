MODULE m_mesh
   USE m_common
   USE m_linear_algebra
   implicit none
   !===================================================================================================================!
   ! Maximum number of Nodes/Edges connecting the Cell
   !===================================================================================================================!
   integer, parameter :: maxed = 4 ! Maximum number of Nodes/Edges connecting the Cell
   !===================================================================================================================!
   ! Node Structure
   !===================================================================================================================!
    !> Node Structure.
   TYPE NodeType
      integer(ip), dimension(:), allocatable :: cell !< Connected Cells Indexes
      integer(ip), dimension(:), allocatable :: edge !< Connected Edges Indexes
      type(vec2d) :: coord !< Node coordinates
      logical :: boundary !< true boolean if Node is a at Mesh Boundary
      integer(ip) :: lim !< Boundary Node Index
   END TYPE NodeType
   !===================================================================================================================!
   ! Boundary Node Structure
   !===================================================================================================================!
    !> boundary Node Structure.
   TYPE NodeTypeLim
      integer(ip) :: ind !< Global Node Index
      character(len=lchar) :: typlim !< Type of the Boundary condition
      integer(ip) :: group !< Group Number in case of Multi Boundary Conditions
   END TYPE NodeTypeLim
   !===================================================================================================================!
   ! Cell Structure
   !===================================================================================================================!
    !> Cell Structure
   TYPE CellType
      integer(ip), dimension(maxed) :: node !< Suited list of Nodes connecting the Cell
      integer(ip), dimension(maxed) :: cell !< Neighboring Cells Index (referenced to upper Nodes list)
      integer(ip), dimension(maxed) :: edge !< Edges Indexes
      integer(ip) :: nbed !< Number of Neighboring Cells/Edges
      logical :: boundary !< True boolean if Cell is a at Mesh Boundary
      real(rp) :: surf !< Cell surface
      real(rp) :: invsurf !< Inverse of Cell surface
      real(rp) :: peri !< Cell perimeter
      type(vec2d) :: grav !< Cell gravity center
      integer(ip) :: rain !< Rain type
   END TYPE CellType
   !===================================================================================================================!
   ! Boundary Cell Structure
   !===================================================================================================================!
    !> Boundary Cell Structure
   TYPE CellTypeLim
      integer(ip) :: ind !< Global Cell Index
      character(len=lchar) :: typlim !< Type of the Boundary Condition
      integer(ip) :: group !< Group Number in case of Multi Boundary Conditions
      integer(ip) :: cell !< Neighboring Cell
      type(vec2d) :: grav !< Cell gravity center
   END TYPE CellTypeLim
   !===================================================================================================================!
   ! Edge Structure
   !===================================================================================================================!
    !> Edge Structure
   TYPE EdgeType
      integer(ip) :: node(2) !< Node Indexes connecting the Edge
      integer(ip) :: cell(2) !< Cell Indexes linked to the Edge, 3rd is index is for 1D2D connections
      integer(ip) :: cell1D2D
      logical :: boundary !< True boolean if Edge is a at mesh Boundary
      logical :: subdomain !< True boolean if Edge is a at sub-domain Boundary
      integer(ip) :: lim !< Boundary Edge Index
      real(rp) :: length !< Edge length
      type(vec2d) :: center !< Edge center
      type(vec2d) :: normal !< Edge normal (oriented from 1 to 2 linked cells)
      type(vec2d) :: tangent !< Edge tangent (oriented from 1 to 2 linked nodes)
      type(vec2d) :: vcell !< Cells(1&2) gravity center vector (oriented from 1 to 2)
      type(vec2d) :: v_edge_cell(2) !< Vector going from the edge center to the cell(1&2) gravity center
   END TYPE EdgeType
   !===================================================================================================================!
   ! Boundary Edge Structure
   !===================================================================================================================!
    ! Boundary Edge Structure
   TYPE EdgeTypeLim
      integer(ip) :: ind !< Global Edge Index
      character(len=lchar) :: typlim !< Type of the Boundary condition
      integer(ip) :: group !< Group Number in case of Multi Boundary Conditions
      integer(ip) :: perio !< Symetric Edge Index in case of Periodic Boundary Condition
   END TYPE EdgeTypeLim
   !===================================================================================================================!
   ! Global Mesh Structure
   !===================================================================================================================!
    !> Global Mesh Structure
   TYPE msh
      integer(ip) :: nn !< Number of mesh Nodes
      integer(ip) :: nnb !< Number of mesh Nodes at Boundary
      integer(ip) :: nc !< Number of mesh Cells
      integer(ip) :: ncb !< Number of mesh Cells at Boundary
      integer(ip) :: ne !< Number of mesh Edges
      integer(ip) :: neb !< Number of mesh Edges at Boundary
      type(NodeType) , dimension(:), allocatable :: node !< Node Structure
      type(NodeTypeLim), dimension(:), allocatable :: nodeb !< Boundary Node Structure
      type(CellType) , dimension(:), allocatable :: cell !< Cell Structure
      type(CellTypeLim), dimension(:), allocatable :: cellb !< Boundary Cell Structure
      type(EdgeType) , dimension(:), allocatable :: edge !< Edge Structure
      type(EdgeTypeLim), dimension(:), allocatable :: edgeb !< Boundary Edge Structure
      character(len=lchar) :: file_name !< Name of the mesh file
      real(rp) :: scal !< Mesh scaling factor
      real(rp) :: surf !< Mesh total surface
   END TYPE msh
   !===================================================================================================================!
   ! Point in Mesh Structure
   !===================================================================================================================!
    !< Point in Mesh Structure
    !! \details correspondance between id cell and coordinates ?????
   TYPE point_in_mesh
      type(vec2d) :: coord !< coordinates (of center cell ???)
      integer(ip) :: cell !< id of the cell ()
   END TYPE point_in_mesh
CONTAINS
 SUBROUTINE dealloc_mesh( mesh )
      implicit none
      type( msh ), intent(out) :: mesh
      if ( allocated( mesh%cell ) ) deallocate( mesh%cell )
      if ( allocated( mesh%cellb ) ) deallocate( mesh%cellb )
      if ( allocated( mesh%edge ) ) deallocate( mesh%edge )
      if ( allocated( mesh%edgeb ) ) deallocate( mesh%edgeb )
      if ( allocated( mesh%node ) ) deallocate( mesh%node )
      if ( allocated( mesh%nodeb ) ) deallocate( mesh%nodeb )
   END SUBROUTINE dealloc_mesh
    !> \brief TO VISUALIZE FOR IMPORTING DATA
   subroutine display_mesh_cell(mesh)
      implicit none
      type(msh), intent(in) :: mesh
      open(10, file = 'XYcells.txt')
      write(10,*) "index ", "x ", "y ", "invsurf ", "surf ", "boundary"
      do i=1,mesh%nc
         write(10,*) i, mesh%cell(i)%grav%x, mesh%cell(i)%grav%y, mesh%cell(i)%invsurf, mesh%cell(i)%surf, mesh%cell(i)%boundary
      end do
      close(10)
    end subroutine
   !===================================================================================================================!
   ! Nodes variables (boundary included in global indexes)
   !===================================================================================================================!
    !> \brief Allocation of one Cell variable
   SUBROUTINE alloc_cell( var , mesh )
      implicit none
      type( msh ), intent(in) :: mesh
      real(rp), dimension(:), allocatable, intent(out) :: var
      allocate( var( mesh%nc ) )
      var(:) = 0._rp
   END SUBROUTINE alloc_cell
    !> \brief Allocation of one Edge variable
   SUBROUTINE alloc_edge( var , mesh )
      implicit none
      type( msh ), intent(in) :: mesh
      real(rp), dimension(:), allocatable, intent(out) :: var
      allocate( var( mesh%ne ) )
      var(:) = 0._rp
   END SUBROUTINE alloc_edge
    !> \brief Allocation of one Node variable
   SUBROUTINE alloc_node( var , mesh )
      implicit none
      type( msh ), intent(in) :: mesh
      real(rp), dimension(:), allocatable, intent(out) :: var
      allocate( var( mesh%nn ) )
      var(:) = 0._rp
   END SUBROUTINE alloc_node
    !> \brief Increase or Decrease Allocated Mesh (CELL) Array Memory
   SUBROUTINE reallocate_cell( var , new )
      implicit none
      type( CellType ), dimension(:), allocatable, intent(inout) :: var
      integer(ip), intent(in) :: new
      integer(ip) :: old
      type( CellType ), dimension(:), allocatable :: temp
      intrinsic move_alloc
      old = size(var)
      if ( new == old ) then
         return
      else if ( new < old ) then
         allocate( temp( new ) )
         temp( 1 : new ) = var( 1 : new )
         call move_alloc( temp , var )
      else
         allocate( temp( new ) )
         temp( 1 : old ) = var( : )
         call move_alloc( temp , var )
      end if
   END SUBROUTINE reallocate_cell
    !> \brief Increase or Decrease Allocated Mesh (EDGE) Array Memory
   SUBROUTINE reallocate_edge( var , new )
      implicit none
      type( EdgeType ), dimension(:), allocatable, intent(inout) :: var
      integer(ip), intent(in) :: new
      integer(ip) :: old
      type( EdgeType ), dimension(:), allocatable :: temp
      intrinsic move_alloc
      old = size(var)
      if ( new == old ) then
         return
      else if ( new < old ) then
         allocate( temp( new ) )
         temp( 1 : new ) = var( 1 : new )
         call move_alloc( temp , var )
      else
         allocate( temp( new ) )
         temp( 1 : old ) = var( : )
         call move_alloc( temp , var )
      end if
   END SUBROUTINE reallocate_edge
    !> \brief Increase or Decrease Allocated Mesh (NODE) Array Memory
   SUBROUTINE reallocate_node( var , new )
      implicit none
      type( NodeType ), dimension(:), allocatable, intent(inout) :: var
      integer(ip), intent(in) :: new
      integer(ip) :: old
      type( NodeType ), dimension(:), allocatable :: temp
      intrinsic move_alloc
      old = size(var)
      if ( new == old ) then
         return
      else if ( new < old ) then
         allocate( temp( new ) )
         temp( 1 : new ) = var( 1 : new )
         call move_alloc( temp , var )
      else
         allocate( temp( new ) )
         temp( 1 : old ) = var( : )
         call move_alloc( temp , var )
      end if
   END SUBROUTINE reallocate_node
    !> \brief Increase or Decrease Allocated Mesh (CELL BOUNDARY) Array Memory
   SUBROUTINE reallocate_cellb( var , new )
      implicit none
      type( CellTypeLim ), dimension(:), allocatable, intent(inout) :: var
      integer(ip), intent(in) :: new
      integer(ip) :: old
      type( CellTypeLim ), dimension(:), allocatable :: temp
      intrinsic move_alloc
      old = size(var)
      if ( new == old ) then
         return
      else if ( new < old ) then
         allocate( temp( new ) )
         temp( 1 : new ) = var( 1 : new )
         call move_alloc( temp , var )
      else
         allocate( temp( new ) )
         temp( 1 : old ) = var( : )
         call move_alloc( temp , var )
      end if
   END SUBROUTINE reallocate_cellb
    !> \brief Increase or Decrease Allocated Mesh (EDGE BOUNDARY) Array Memory
   SUBROUTINE reallocate_edgeb( var , new )
      implicit none
      type( EdgeTypeLim ), dimension(:), allocatable, intent(inout) :: var
      integer(ip), intent(in) :: new
      integer(ip) :: old
      type( EdgeTypeLim ), dimension(:), allocatable :: temp
      intrinsic move_alloc
      old = size(var)
      if ( new == old ) then
         return
      else if ( new < old ) then
         allocate( temp( new ) )
         temp( 1 : new ) = var( 1 : new )
         call move_alloc( temp , var )
      else
         allocate( temp( new ) )
         temp( 1 : old ) = var( : )
         call move_alloc( temp , var )
      end if
   END SUBROUTINE reallocate_edgeb
    !> \brief Increase or Decrease Allocated Mesh (NODE BOUNDARY) Array Memory
   SUBROUTINE reallocate_nodeb( var , new )
      implicit none
      type( NodeTypeLim ), dimension(:), allocatable, intent(inout) :: var
      integer(ip), intent(in) :: new
      integer(ip) :: old
      type( NodeTypeLim ), dimension(:), allocatable :: temp
      intrinsic move_alloc
      old = size(var)
      if ( new == old ) then
         return
      else if ( new < old ) then
         allocate( temp( new ) )
         temp( 1 : new ) = var( 1 : new )
         call move_alloc( temp , var )
      else
         allocate( temp( new ) )
         temp( 1 : old ) = var( : )
         call move_alloc( temp , var )
      end if
   END SUBROUTINE reallocate_nodeb
    !> \brief Calculation of cells connectivity
   SUBROUTINE calc_cells_connectivity( mesh )
      implicit none
      type(msh), intent(inout) :: mesh
      integer(ip) :: side1 , side2 , tri1 , tri2 , col( 4 , maxed * mesh%nc ) , node( maxed ) , n1 , n2
      mesh%cell(:)%nbed = maxed
      j = 0
write(*,*)"INTO calc_cells_connectivity"
  ! following loop compute the table col
  ! 1. first and second lines of col contain the global indices of nodes that are in the same edge
  ! (small indice always at first to make sure that the same edge with be adjacent after sorting)
  ! 2. third line contains the local indice of the edge in the cell
  ! 3. forth line contains the global indices of the cells
      do i = 1,mesh%nc ! compute the neighbor cells
         mesh%cell(i)%cell(:) = 0
         node(:) = mesh%cell(i)%node(:)
         do k = 1,maxed
            n1 = k
            n2 = mod(k,maxed) + 1
            if ( node(n1) == node(n2) ) then ! triangular case
               mesh%cell(i)%cell(k) = -1
               mesh%cell(i)%nbed = mesh%cell(i)%nbed - 1
            else if ( node(n1) < node(n2) ) then
               j = j + 1 ; col(1:4,j) = (/ node(n1), node(n2), k , i /)
            else
               j = j + 1 ; col(1:4,j) = (/ node(n2), node(n1), k , i /)
            end if
         end do
      end do
write(*,*)"STEP 1 DONE"
      ! after sorting, column(s) that contain the same edge will be adjacent in col
      ! the same pair of indices (an edge) appears 2 times or 1 time in col
      ! 2 times => the edge connecting previous pair of indices is inside the domain
      ! 1 time => the edge is on the boundary
      call i4col_sort_a( 4 , j , col )
      mesh%ne = 0
      mesh%neb = 0
      i = 1
      do
         if ( j <= i ) then
            if ( i == j ) then
    mesh%neb = mesh%neb + 1
    mesh%ne = mesh%ne + 1 ! SHENYUAN
   end if
            exit
         end if
         mesh%ne = mesh%ne + 1 ! always add one more edges
         if ( col(1,i) /= col(1,i+1) .or. col(2,i) /= col(2,i+1) ) then !if the current edge doesnot repeat in col
            i = i + 1
            mesh%neb = mesh%neb + 1 ! there is one more boundary edge.
            cycle
         end if
         side1 = col(3,i )
         side2 = col(3,i+1)
         tri1 = col(4,i )
         tri2 = col(4,i+1)
         mesh%cell( tri1 )%cell( side1 ) = tri2
         mesh%cell( tri2 )%cell( side2 ) = tri1
         i = i + 2 ! jump to the next edge
      end do
write(*,*)"OUT calc_cells_connectivity"
   END SUBROUTINE calc_cells_connectivity
    !> FOR PYTHON WRAPPING : initialise NodeType
   subroutine NodeType_initialise(node,dim1,dim2)
      implicit none
      type(NodeType), intent(out) :: node
      integer, intent(in) :: dim1
      integer, intent(in) :: dim2
      allocate(node%cell(dim1))
      allocate(node%edge(dim2))
      node%boundary = .true.
      node%lim = 0
      node%coord%x = 0
      node%coord%y = 0
   end subroutine NodeType_initialise
        !> FOR PYTHON WRAPPING : initialise NodeType boundary
   subroutine NodeTypeLim_initialise(nodelim,ind,group)
      implicit none
      integer, intent(in) :: ind
      integer, intent(in) :: group
      type(NodeTypeLim), intent(out) :: nodelim
      nodelim%ind = ind
      nodelim%group = group
      nodelim%typlim = 'NULL'
   end subroutine NodeTypeLim_initialise
        !> FOR PYTHON WRAPPING : initialise mesh
   subroutine msh_initialise(mesh)
      implicit none
      type(msh), intent(inout) :: mesh
      call Mesh_Input(mesh)
   end subroutine msh_initialise
         !> FOR PYTHON WRAPPING : finalise mesh (deallocate mesh infos)
 SUBROUTINE msh_finalise(mesh)
      implicit none
      type(msh), intent(inout) :: mesh
      if ( allocated( mesh%cell ) ) deallocate( mesh%cell )
      if ( allocated( mesh%cellb ) ) deallocate( mesh%cellb )
      if ( allocated( mesh%edge ) ) deallocate( mesh%edge )
      if ( allocated( mesh%edgeb ) ) deallocate( mesh%edgeb )
      if ( allocated( mesh%node ) ) deallocate( mesh%node )
      if ( allocated( mesh%nodeb ) ) deallocate( mesh%nodeb )
   END SUBROUTINE msh_finalise
END MODULE m_mesh
