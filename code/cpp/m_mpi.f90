MODULE m_mpi
   USE m_common
   USE m_linear_algebra
   USE m_mesh
   implicit none
   integer(ip) :: np !> number of processes
   integer(ip) :: proc !> process number from 0 to np-1
   integer(ip) :: code !> ???
   integer(ip) :: nneighb !> number of neighboring proc
    integer(ip), dimension(:), allocatable :: swap_index , inv_swap_index !> ??? !> ?????
   integer(ip), dimension(: ), allocatable :: part !> ???
   integer(ip), dimension(: ), allocatable :: part_size !> ???
   integer(ip), dimension(:,:), allocatable :: part_neighbs !> ???
   real(rp) :: val_tmp_r !> ???
   integer(ip) :: val_tmp_i !> ???
CONTAINS
   SUBROUTINE Init_MPI
      implicit none
         proc = 0 ! compatible execution value for sequential Dassflow version
         np = 1 ! compatible execution value for sequential Dassflow version
   END SUBROUTINE Init_MPI
   SUBROUTINE End_MPI
      implicit none
      if ( proc == 0 ) write(6,'(A)')
      call exit()
   END SUBROUTINE End_MPI
   SUBROUTINE Mesh_Partition_Scotch( mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(inout) :: mesh
         allocate( swap_index( mesh%nc ) ) ; swap_index(:) = (/ ( i , i=1,mesh%nc ) /)
         allocate( part( mesh%nc ) ) ; part (:) = 0! (/ ( i , i=1,mesh%nc ) /)
         allocate( inv_swap_index( mesh%nc ) ) ; inv_swap_index(:) = (/ ( i , i=1,mesh%nc ) /)
         mesh%edge(:)%subdomain = mesh%edge(:)%boundary
   END SUBROUTINE Mesh_Partition_Scotch
   SUBROUTINE swap_mesh( swap , inv_swap , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(inout) :: mesh
      integer(ip), allocatable, dimension(:), intent(inout) :: swap
      integer(ip), allocatable, dimension(:), intent(inout) :: inv_swap
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      type(msh) :: mesh_tmp
      !integer(ip), dimension( mesh%nn ) :: swap_index_node , inv_swap_index_node
      integer(ip), dimension( mesh%ne ) :: swap_index_edge , inv_swap_index_edge
      !================================================================================================================!
      ! Creation of a temporary mesh to swap indexes
      !================================================================================================================!
      allocate( mesh_tmp%node ( mesh%nn ) )
      allocate( mesh_tmp%edge ( mesh%ne ) )
      allocate( mesh_tmp%cell ( mesh%nc ) )
      allocate( mesh_tmp%edgeb( mesh%neb ) )
      allocate( mesh_tmp%cellb( mesh%ncb ) )
      mesh_tmp = mesh
      !================================================================================================================!
      ! Swapping Cells
      !================================================================================================================!
      j = 0
      do k = 1 , part_size(proc) + part_neighbs(proc,np)
         mesh%cell( k )%nbed = mesh_tmp%cell( swap(k) )%nbed
         mesh%cell( k )%node(:) = mesh_tmp%cell( swap(k) )%node(:)
         mesh%cell( k )%edge(:) = mesh_tmp%cell( swap(k) )%edge(:)
         do ke = 1,maxed
            i = mesh_tmp%cell( swap(k) )%cell( ke )
            if ( i <= 0 ) then
               mesh%cell( k )%cell( ke ) = i
            else if ( i <= mesh%nc ) then
               mesh%cell( k )%cell( ke ) = inv_swap( i )
            else
               j = j + 1
               mesh%cell( k )%cell( ke ) = part_size( proc ) + part_neighbs( proc , np ) + j
            end if
         end do
      end do
      call reallocate_cell( mesh%cell , part_size( proc ) + part_neighbs( proc , np ) ) !reduce the size of local (proc wise) cells list, nc = part_size(proc) ncb = part_neighbs()
      !================================================================================================================!
      ! Swapping Edges
      !================================================================================================================!
      ie = 0
      swap_index_edge(:) = 0
      do k = 1,part_size( proc ) + part_neighbs( proc , np )
         do ke = 1,maxed
            if ( mesh%cell(k)%cell(ke) > k ) then
               ie = ie + 1
               je = mesh%cell(k)%edge(ke)
                   swap_index_edge( ie ) = je
               inv_swap_index_edge( je ) = ie
               mesh%edge( ie ) = mesh_tmp%edge( je )
               mesh%edge( ie )%cell(1) = k
               mesh%edge( ie )%cell(2) = mesh%cell(k)%cell(ke)
            end if
         end do
      end do
      do k = 1,part_size(proc)
         do ke = 1,maxed
            je = mesh%cell(k)%edge(ke)
            if ( je > 0 ) mesh%cell(k)%edge(ke) = inv_swap_index_edge( je )
         end do
      end do
      j = 0
      ib = 0
      jb = 0
      mesh%edge(:)%subdomain = .false.
      do i = 1,ie
         if ( mesh%edge(i)%boundary ) then
            j = j + 1
            mesh%edge ( i )%lim = j
            mesh%edgeb( j )%ind = i
            mesh%cellb( j )%cell = mesh%edge(i)%cell(1)
         end if
         if ( minval( mesh%edge(i)%cell(:) ) <= part_size( proc ) ) then
            ib = ib + 1
            if ( mesh%edge(i)%boundary ) then
               jb = jb + 1
            else if ( maxval( mesh%edge(i)%cell(:) ) > part_size( proc ) ) then
               mesh%edge(i)%subdomain = .true.
            end if
         end if
      end do
      call reallocate_edge( mesh%edge , ib )
      !================================================================================================================!
      ! New local Cell number
      !================================================================================================================!
      mesh%nc = part_size ( proc )
      mesh%ncb = part_neighbs( proc , np ) + jb
      !================================================================================================================!
      ! New local Edge number
      !================================================================================================================!
      mesh%ne = ib
      mesh%neb = jb
      !================================================================================================================!
      ! Swapping Mesh Boundaries
      !================================================================================================================!
      call reallocate_i( swap , mesh%nc + mesh%ncb )
      do i = 1,mesh%neb
         ie = swap_index_edge( mesh%edgeb(i)%ind )
         mesh%edgeb(i)%typlim = mesh_tmp%edgeb( mesh_tmp%edge( ie )%lim )%typlim
         mesh%edgeb(i)%group = mesh_tmp%edgeb( mesh_tmp%edge( ie )%lim )%group
         swap_index( mesh%nc + part_neighbs( proc , np ) + i ) = mesh_tmp%edge( ie )%cell(2)
      end do
      !================================================================================================================!
      ! Deallocating temporary mesh to swap indexes
      !================================================================================================================!
      call dealloc_mesh( mesh_tmp )
   END SUBROUTINE swap_mesh
   SUBROUTINE fill_swap_lists( mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(inout) :: mesh
      !================================================================================================================!
      ! Constructing base inv_swap_index
      !================================================================================================================!
      do i = 1,mesh%nc
         if ( swap_index(i) /= 0 ) then
            inv_swap_index( swap_index(i) ) = i
         else
            exit
         end if
      end do
      !================================================================================================================!
      ! Filling swap_index using gaps in inv_swap_index
      !================================================================================================================!
      do j = 1,mesh%nc
         if ( inv_swap_index(j) == 0 ) then
            swap_index(i) = j
            inv_swap_index(j) = i
            i = i + 1
         end if
      end do
   END SUBROUTINE fill_swap_lists
   SUBROUTINE fill_swap_index( mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(inout) :: mesh
      !================================================================================================================!
      ! Filling swap_index vector
      !================================================================================================================!
      do j = 1,mesh%nc
         do i = 1,mesh%nc
            if ( swap_index(i) == j ) exit
            if ( swap_index(i) == 0 ) then
               swap_index(i) = j
               exit
            end if
         end do
      end do
   END SUBROUTINE fill_swap_index
   SUBROUTINE fill_inv_swap_index( mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(inout) :: mesh
      !================================================================================================================!
      ! Filling inv_swap_index vector
      !================================================================================================================!
      do j = 1,mesh%nc
         do i = 1,mesh%nc
            if ( swap_index(i) == j ) then
               inv_swap_index(j) = i
               exit
            end if
         end do
      end do
   END SUBROUTINE fill_inv_swap_index
   SUBROUTINE com_var_i( var , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(in) :: mesh
      integer(ip), dimension( mesh%nc + mesh%ncb ), intent(inout) :: var
      !================================================================================================================!
      !
      !================================================================================================================!
   END SUBROUTINE com_var_i
   SUBROUTINE com_var_r( var , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(in) :: mesh
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(inout) :: var
      !================================================================================================================!
      !
      !================================================================================================================!
   END SUBROUTINE com_var_r
   SUBROUTINE mpi_send_recv_scal_i( to_send , to_recv , proc_send , proc_recv )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      integer(ip), intent(in ) :: to_send
      integer(ip), intent(out) :: to_recv
      integer(ip), intent(in) :: proc_send , proc_recv
      !================================================================================================================!
      !
      !================================================================================================================!
         to_recv = to_send
   END SUBROUTINE mpi_send_recv_scal_i
   SUBROUTINE mpi_send_recv_scal_r( to_send , to_recv , proc_send , proc_recv )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      real(rp), intent(in ) :: to_send
      real(rp), intent(out) :: to_recv
      integer(ip), intent(in) :: proc_send , proc_recv
      !================================================================================================================!
      !
      !================================================================================================================!
         to_recv = to_send
   END SUBROUTINE mpi_send_recv_scal_r
   SUBROUTINE mpi_send_recv_array_i( to_send , to_recv , proc_send , proc_recv )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      integer(ip), dimension(:), intent(in ) :: to_send
      integer(ip), dimension(:), intent(out) :: to_recv
      integer(ip), intent(in) :: proc_send , proc_recv
      !================================================================================================================!
      !
      !================================================================================================================!
         to_recv = to_send
   END SUBROUTINE mpi_send_recv_array_i
   SUBROUTINE mpi_send_recv_array_r( to_send , to_recv , proc_send , proc_recv )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      real(rp), dimension(:), intent(in ) :: to_send
      real(rp), dimension(:), intent(out) :: to_recv
      integer(ip), intent(in) :: proc_send , proc_recv
      !================================================================================================================!
      !
      !================================================================================================================!
         to_recv = to_send
   END SUBROUTINE mpi_send_recv_array_r
   SUBROUTINE mpi_sum_r( val )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      real(rp), intent(inout) :: val
      !================================================================================================================!
      !
      !================================================================================================================!
      call mpi_wait_all
   END SUBROUTINE mpi_sum_r
   SUBROUTINE mpi_sum_i( val )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      integer(ip), intent(inout) :: val
      !================================================================================================================!
      !
      !================================================================================================================!
   END SUBROUTINE mpi_sum_i
   SUBROUTINE mpi_max_r( val )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      real(rp), intent(inout) :: val
      !================================================================================================================!
      !
      !================================================================================================================!
   END SUBROUTINE mpi_max_r
   SUBROUTINE mpi_max_i( val )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      integer(ip), intent(inout) :: val
      !================================================================================================================!
      !
      !================================================================================================================!
   END SUBROUTINE mpi_max_i
   SUBROUTINE mpi_min_r( val )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      real(rp), intent(inout) :: val
      !================================================================================================================!
      !
      !================================================================================================================!
   END SUBROUTINE mpi_min_r
   SUBROUTINE mpi_min_i( val )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      integer(ip), intent(inout) :: val
      !================================================================================================================!
      !
      !================================================================================================================!
   END SUBROUTINE mpi_min_i
   SUBROUTINE mpi_bcast_r( val , pr )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      real(rp), intent(inout) :: val
      integer(ip), intent(in) :: pr
      !================================================================================================================!
      !
      !================================================================================================================!
   END SUBROUTINE mpi_bcast_r
   SUBROUTINE mpi_bcast_i( val , pr )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      integer(ip), intent(inout) :: val
      integer(ip), intent(in) :: pr
      !================================================================================================================!
      !
      !================================================================================================================!
   END SUBROUTINE mpi_bcast_i
   SUBROUTINE mpi_allgather_r( val , temp )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      real(rp), intent(inout) :: val
      real(rp), dimension(np), intent(inout) :: temp
      !================================================================================================================!
      !
      !================================================================================================================!
         temp(1) = val
   END SUBROUTINE mpi_allgather_r
   SUBROUTINE mpi_allgather_i( val , temp )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      integer(ip), intent(inout) :: val
      integer(ip), dimension(np), intent(inout) :: temp
      !================================================================================================================!
      !
      !================================================================================================================!
         temp(1) = val
   END SUBROUTINE mpi_allgather_i
   SUBROUTINE mpi_wait_all
      implicit none
   END SUBROUTINE mpi_wait_all
   SUBROUTINE Stopping_Program_Sub( comment )
      implicit none
      character(len=*), intent(in) :: comment
      if ( proc == 0 ) then
         write(6,*)
         write(6,'(A80)') '================================================================================'
         write(6,'(A,A)') '*  ' , comment
         write(6,'(A80)') '================================================================================'
      end if
      if ( proc == 0 ) then
         write(6,*)
         write(6,'(A)') '================================================================================'
         write(6,'(A)') '*  STOPPING DASSFLOW RUN'
         write(6,'(A)') '================================================================================'
         write(6,*)
      end if
      ! call f90wrap_abort(comment) ! lilian choix douteux sur conseil de fran\U000000e7ois
   END SUBROUTINE Stopping_Program_Sub
END MODULE m_mpi
