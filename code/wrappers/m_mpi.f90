MODULE m_mpi
   USE m_common
   USE m_linear_algebra
   USE m_mesh
   use mpi
   implicit none
   !include 'mpif.h'
   !
      include 'scotchf.h'
   integer(ip) :: np !> number of processes
   integer(ip) :: proc !> process number from 0 to np-1
   integer(ip) :: code !> ???
   integer(ip) :: nneighb !> number of neighboring proc
    integer(ip), dimension(:), allocatable :: swap_index , inv_swap_index !> ??? !> ?????
      integer(ip), dimension(:), allocatable :: type_com !> ????
      integer(ip), parameter :: realtype = mpi_double_precision !> double precision for mpi ?
      integer(ip), parameter :: inttype = mpi_integer !> integer precision for mpi ?
      integer(ip), dimension( mpi_status_size ) :: stat !> ???
      real(rp), dimension( Scotch_graphdim ) :: Scotchgraph !> ???
      real(rp), dimension( Scotch_stratdim ) :: Scotchstrat !> ???
   integer(ip), dimension(: ), allocatable :: part !> ???
   integer(ip), dimension(: ), allocatable :: part_size !> ???
   integer(ip), dimension(:,:), allocatable :: part_neighbs !> ???
   real(rp) :: val_tmp_r !> ???
   integer(ip) :: val_tmp_i !> ???
CONTAINS
   SUBROUTINE Init_MPI
      implicit none
         call MPI_COMM_SIZE( MPI_COMM_WORLD , np , code ) ! np -> number of processes
         call MPI_COMM_RANK( MPI_COMM_WORLD , proc , code ) ! proc -> process number from 0 to np-1
write(*,*) "mpi in fortran:np, proc", np, proc
   END SUBROUTINE Init_MPI
   SUBROUTINE End_MPI
      implicit none
      if ( proc == 0 ) write(6,'(A)')
        call MPI_BARRIER( MPI_COMM_WORLD , code )
        call MPI_FINALIZE( code )
      call exit()
   END SUBROUTINE End_MPI
   SUBROUTINE Mesh_Partition_Scotch( mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(inout) :: mesh
         !=============================================================================================================!
         ! Local Variables
         !=============================================================================================================!
         integer(ip) :: file_u , ierr , ireq
         character(lchar) :: filename_xyz , filename_grf , filename_map
         character(256) :: command
         integer(ip), dimension( 2*np ) :: req
         integer(ip), dimension( mpi_status_size , 2*np ) :: status_mpi
         integer(ip), dimension(:), allocatable :: long_bloc , depla_bloc
         !=============================================================================================================!
         ! Begin Subroutine
         !=============================================================================================================!
         if ( proc == 0 ) then
            write(6,'(A)')
            write(6,'(A)') '================================================================================'
            write(6,'(A)') '*  Mesh Partitionning'
            write(6,'(A)') '================================================================================'
         end if
         !=============================================================================================================!
         ! Creating Files Names for each process
         !=============================================================================================================!
         call system('mkdir -p msh')
         write(filename_xyz,'("msh/cells_",I4.4,".xyz")') proc + 1
         write(filename_map,'("msh/cells_",I4.4,".map")') proc + 1
         write(filename_grf,'("msh/cells_",I4.4,".grf")') proc + 1
         !=============================================================================================================!
         ! Generating cells xyz file
         !=============================================================================================================!
         open(10+proc,file=filename_xyz,status='replace')
         write(10+proc,'(I1)') 2
         write(10+proc,'(I8)') mesh%nc
         do i = 1,mesh%nc
            if ( mesh%cell(i)%node(4) == mesh%cell(i)%node(1) ) then
               mesh%cell(i)%grav%x = sum( mesh%node( mesh%cell(i)%node(1:3) )%coord%x ) * d1p3
               mesh%cell(i)%grav%y = sum( mesh%node( mesh%cell(i)%node(1:3) )%coord%y ) * d1p3
            else
               mesh%cell(i)%grav%x = sum( mesh%node( mesh%cell(i)%node(1:4) )%coord%x ) * d1p4
               mesh%cell(i)%grav%y = sum( mesh%node( mesh%cell(i)%node(1:4) )%coord%y ) * d1p4
            end if
            write(10+proc,'(I8,2ES15.7)') i-1 , mesh%cell(i)%grav%x , mesh%cell(i)%grav%y
         end do
         close(10+proc)
         !=============================================================================================================!
         ! Generating cells grf file
         !=============================================================================================================!
         k = 0
         do i = 1,mesh%nc
            k = k + count( mesh%cell(i)%cell(:) >= 1 .and. mesh%cell(i)%cell(:) <= mesh%nc )
         end do
         open(10+proc,file=filename_grf,status='replace')
         write(10+proc,'(1I10    )') 0
         write(10+proc,'(2I10    )') mesh%nc , k
         write(10+proc,'(1I10,A10)') 1 , '000'
         do i = 1,mesh%nc
            write(10+proc,'(20I10)') count( mesh%cell(i)%cell(:) >= 1 .and. &
                                            mesh%cell(i)%cell(:) <= mesh%nc ) , &
                                     pack ( mesh%cell(i)%cell(:) , mesh%cell(i)%cell(:) >= 1 .and. &
                                                                   mesh%cell(i)%cell(:) <= mesh%nc )
         end do
         close(10+proc)
         !=============================================================================================================!
         ! Graph and Paritionning Strategy Initialization
         !=============================================================================================================!
         call ScotchfGraphInit( Scotchgraph(1) , ierr )
         if ( ierr /= 0 ) call Stopping_Program_Sub( 'Scotch Error : cannot initialize strat' )
         call ScotchfStratInit( Scotchstrat(1) , ierr )
         if ( ierr /= 0 ) call Stopping_Program_Sub( 'Scotch Error : cannot initialize graph' )
         !=============================================================================================================!
         ! Graph File Loading
         !=============================================================================================================!
         open(10+proc,file=filename_grf)
            file_u = FNUM (10+proc)
         call ScotchfGraphLoad( Scotchgraph(1) , file_u , 1 , 0 , ierr )
         if ( ierr /= 0 ) call Stopping_Program_Sub( 'Scotch error : cannot load graph' )
         close(10+proc)
         !=============================================================================================================!
         ! Scotch Graph File Checking
         !=============================================================================================================!
         call ScotchfGraphCheck( Scotchgraph(1) , ierr )
         if ( ierr /= 0 ) call Stopping_Program_Sub( 'Scotch Error : invalid graph file checking' )
         !=============================================================================================================!
         ! Scotch Partitionning SubRoutine Calling
         !=============================================================================================================!
         allocate( part( mesh%nc ) )
         call ScotchfGraphPart( Scotchgraph(1) , np , Scotchstrat(1) , part(1) , ierr )
         if ( ierr /= 0 ) call Stopping_Program_Sub( 'Scotch Error : cannot partitioning graph' )
         open(10+proc,file=filename_map)
         write(10+proc,*) mesh%nc
         do i = 1,mesh%nc
            write(10+proc,'(2I10)') i-1 , part(i)
         end do
         close(10+proc)
         !=============================================================================================================!
         ! Generating eps and pdf Graphic Mesh Partitionning Files
         !=============================================================================================================!
         if ( proc == 0 ) then
            write(command,'("../../libs/scotch_5.1.12_esmumps/bin/gout -mn -Op{''c,d,e,s''} ",a," ",a," > ",a)' ) &
            trim(filename_grf),trim(filename_xyz),"msh/mesh_elt.eps"
            call system(command)
            write(command,'("../../libs/scotch_5.1.12_esmumps/bin/gout     -Op{''c,d,e,s''} ",a," ",a," ",a," > ",a)') &
            trim(filename_grf),trim(filename_xyz),trim(filename_map),"msh/mesh_elt_part.eps"
            call system(command)
         end if
         !=============================================================================================================!
         ! Partitionning Operations
         !=============================================================================================================!
         allocate( swap_index( 1 : mesh%nc ) ) ; swap_index(:) = 0
         allocate( part_size ( 0 : np - 1 ) ) ; part_size (: ) = 0
         allocate( part_neighbs( 0 : np - 1 , 0 : np ) ) ; part_neighbs(:,:) = 0 ! column np represents the overall number of neighers
         do k = 1,mesh%nc
            part_size( part(k) ) = part_size( part(k) ) + 1
            if ( part(k) == proc ) swap_index( part_size( proc ) ) = k
         end do
         i = 0
         do j = 0,np-1
            if ( j == proc ) cycle !don't do the next steps if proc == j
            do k = 1,mesh%nc
               if ( part( k ) /= j ) cycle ! if the current cell doesnot belong to proc j
               do ke = 1,maxed
                  ie = mesh%cell(k)%cell(ke) ! neighbor cell index to cell k
                  if ( ie >= 1 .and. ie <= mesh%nc ) then
                     if ( part( ie ) == proc ) then
                        part_neighbs( proc , j ) = part_neighbs( proc , j ) + 1 ! add proc j as one more neighbor to proc
                        i = i + 1
                        swap_index( part_size( proc ) + i ) = k !? add cell k to be the "ghost cell" of local mesh
                        exit
                     end if
                  end if
               end do
            end do
         end do
         !=============================================================================================================!
         ! All processes have neighbors account of all other ones
         !=============================================================================================================!
         do i = 0,np-1
            do j = 0,np-1
               call MPI_BCAST( part_neighbs( i , j ) , 1 , inttype , i , MPI_COMM_WORLD , code )
            end do
         end do
         part_neighbs( : , np ) = sum( part_neighbs( : , 0 : np-1 ) , dim = 2 )
         nneighb = count( part_neighbs( proc , 0 : np-1 ) > 0 )
         !=============================================================================================================!
         ! Sending at the end of swap_index indexes that proc will send to other ones
         !=============================================================================================================!
         ireq = 0
         do k = 0,np-1
            if ( k == proc ) cycle
            ireq = ireq + 1
            call MPI_ISEND( swap_index( part_size( proc ) + 1 + sum( part_neighbs(proc,0:k-1) ) ) , &
                            part_neighbs( proc , k ) , inttype , k , k , MPI_COMM_WORLD , req(ireq) , code )
         end do
         do k = 0,np-1
            if ( k == proc ) cycle
            ireq = ireq + 1
            call MPI_IRECV( swap_index( part_size( proc ) + 1 + part_neighbs(proc,np) + sum( part_neighbs(0:k-1,proc) ) ) , &
                            part_neighbs( k , proc ) , inttype , k , proc , MPI_COMM_WORLD , req(ireq) , code )
         end do
         call MPI_WAITALL( ireq , req , status_mpi , code )
         !=============================================================================================================!
         ! Communication Data Array Construction
         !=============================================================================================================!
         if (allocated (type_com)) deallocate (type_com)
         allocate( long_bloc ( 1 : sum( part_neighbs( 0 : np-1 , proc ) ) ) )
         allocate( depla_bloc ( 1 : sum( part_neighbs( 0 : np-1 , proc ) ) ) )
         allocate( type_com ( 0 : np-1 ) )
         long_bloc (:) = 1
         depla_bloc(:) = 0
         do j = 0,np-1
            if ( part_neighbs( j , proc ) == 0 ) cycle
            do i = 1 + sum( part_neighbs( 0 : j-1 , proc ) ) , sum( part_neighbs( 0 : j , proc ) )
               do k = 1,part_size( proc )
                  if ( swap_index(k) == swap_index( part_size( proc ) + part_neighbs(proc,np) + i ) ) then
                     depla_bloc(i) = k - 1
                  end if
               end do
            end do
            call MPI_TYPE_INDEXED( part_neighbs( j , proc ) , long_bloc ( 1 + sum( part_neighbs( 0 : j-1 , proc ) ) : &
                                                                              sum( part_neighbs( 0 : j , proc ) ) ) , &
                                                              depla_bloc( 1 + sum( part_neighbs( 0 : j-1 , proc ) ) : &
                                                                              sum( part_neighbs( 0 : j , proc ) ) ) &
                                                            , realtype , type_com(j) ,code)
            call MPI_TYPE_COMMIT( type_com(j) , code )
         end do
         deallocate( depla_bloc )
         deallocate( long_bloc )
         swap_index( part_size( proc ) + 1 + part_neighbs(proc,np) : mesh%nc ) = 0
      !=============================================================================================================!
    !Get swap index without border cells, to send to python
         !=============================================================================================================!
         ! Array Index Reorganisation
         !=============================================================================================================!
         allocate( inv_swap_index( mesh%nc ) ) ; inv_swap_index(:) = 0
         call fill_swap_lists( mesh )
         call swap_mesh( swap_index , inv_swap_index , mesh ) !some pb in this function ?
         !=============================================================================================================!
         ! Freeing Graph and Paritionning Strategy
         !=============================================================================================================!
         call ScotchfgraphExit ( Scotchgraph(1) )
         call ScotchfstratExit ( Scotchstrat(1) )
         call MPI_BARRIER( MPI_COMM_WORLD , code )
         if ( proc == 0 ) then
            write(6,'(A)') '================================================================================'
            write(6,'(A)') '*  Mesh Partitionning Finished'
            write(6,'(A)') '================================================================================'
         end if
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
         integer(ip), dimension( 2*np ) :: req
         integer(ip), dimension( mpi_status_size , 2*np ) :: status_mpi
         integer(ip) :: ireq
         ireq = 0
         do k = 0,np-1
            if ( part_neighbs( k , proc ) == 0 ) cycle
            ireq = ireq + 1
            call MPI_ISEND( var( 1 ) , 1 , type_com(k) , k , k , MPI_COMM_WORLD , req(ireq) , code )
         end do
         do k = 0,np-1
            if ( part_neighbs( proc , k ) == 0 ) cycle
            ireq = ireq + 1
            call MPI_IRECV( var( part_size( proc ) + 1 + sum( part_neighbs( proc , 0:k-1 ) ) ) , &
                            part_neighbs( proc , k ) , inttype , k , proc , MPI_COMM_WORLD , req(ireq) , code )
         end do
         call MPI_WAITALL( ireq , req , status_mpi , code )
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
         integer(ip), dimension( 2*np ) :: req
         integer(ip), dimension( mpi_status_size , 2*np ) :: status_mpi
         integer(ip) :: ireq
         ireq = 0
         do k = 0,np-1
            if ( part_neighbs( k , proc ) == 0 ) cycle
            ireq = ireq + 1
            call MPI_ISEND( var( 1 ) , 1 , type_com(k) , k , k , MPI_COMM_WORLD , req(ireq) , code )
         end do
         do k = 0,np-1
            if ( part_neighbs( proc , k ) == 0 ) cycle
            ireq = ireq + 1
            call MPI_IRECV( var( part_size( proc ) + 1 + sum( part_neighbs( proc , 0:k-1 ) ) ) , &
                            part_neighbs( proc , k ) , realtype , k , proc , MPI_COMM_WORLD , req(ireq) , code )
         end do
         call MPI_WAITALL( ireq , req , status_mpi , code )
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
         integer(ip), dimension( mpi_status_size ) :: status_mpi
         if ( proc == proc_send ) then
            call MPI_SEND( to_send , 1 , inttype , proc_recv , 1 , MPI_COMM_WORLD , code )
         else if ( proc == proc_recv ) then
            call MPI_RECV( to_recv , 1 , inttype , proc_send , 1 , MPI_COMM_WORLD , status_mpi , code )
         end if
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
         integer(ip), dimension( mpi_status_size ) :: status_mpi
         if ( proc == proc_send ) then
            call MPI_SEND( to_send , 1 , realtype , proc_recv , 1 , MPI_COMM_WORLD , code )
         else if ( proc == proc_recv ) then
            call MPI_RECV( to_recv , 1 , realtype , proc_send , 1 , MPI_COMM_WORLD , status_mpi , code )
         end if
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
         integer(ip), dimension( mpi_status_size ) :: status_mpi
         if ( proc == proc_send ) then
            call MPI_SEND( to_send(1) , size( to_send ) , inttype , proc_recv , 1 , MPI_COMM_WORLD , code )
         else if ( proc == proc_recv ) then
            call MPI_RECV( to_recv(1) , size( to_recv ) , inttype , proc_send , 1 , MPI_COMM_WORLD , status_mpi , code )
         end if
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
         integer(ip), dimension( mpi_status_size ) :: status_mpi
         if ( proc == proc_send ) then
            call MPI_SEND( to_send(1) , size( to_send ) , realtype , proc_recv , 1 , MPI_COMM_WORLD , code )
         else if ( proc == proc_recv ) then
            call MPI_RECV( to_recv(1) , size( to_recv ) , realtype , proc_send , 1 , MPI_COMM_WORLD , status_mpi , code )
         end if
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
         val_tmp_r = val
         call MPI_ALLREDUCE( val_tmp_r , val , 1 , realtype , MPI_SUM , MPI_COMM_WORLD , code )
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
         val_tmp_i = val
         call MPI_ALLREDUCE( val_tmp_i , val , 1 , inttype , MPI_SUM , MPI_COMM_WORLD , code )
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
         val_tmp_r = val
         call MPI_ALLREDUCE( val_tmp_r , val , 1 , realtype , MPI_MAX , MPI_COMM_WORLD , code )
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
         val_tmp_i = val
         call MPI_ALLREDUCE( val_tmp_i , val , 1 , inttype , MPI_MAX , MPI_COMM_WORLD , code )
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
         val_tmp_r = val
         call MPI_ALLREDUCE( val_tmp_r , val , 1 , realtype , MPI_MIN , MPI_COMM_WORLD , code )
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
         val_tmp_i = val
         call MPI_ALLREDUCE( val_tmp_i , val , 1 , inttype , MPI_MIN , MPI_COMM_WORLD , code )
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
         call MPI_BCAST( val , 1 , realtype , pr , MPI_COMM_WORLD , code )
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
         call MPI_BCAST( val , 1 , inttype , pr , MPI_COMM_WORLD , code )
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
         call MPI_ALLGATHER( val , 1 , realtype , temp , 1 , realtype , MPI_COMM_WORLD , code )
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
         call MPI_ALLGATHER( val , 1 , inttype , temp , 1 , inttype , MPI_COMM_WORLD , code )
   END SUBROUTINE mpi_allgather_i
   SUBROUTINE mpi_wait_all
      implicit none
        call MPI_BARRIER( MPI_COMM_WORLD , code )
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
         call MPI_BARRIER( MPI_COMM_WORLD , code )
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
