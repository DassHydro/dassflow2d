SUBROUTINE Read_Input(filename)
      USE m_model
   implicit none
   character(len=*), intent(in) :: filename
   !===================================================================================================================!
   !
   !===================================================================================================================!
   call Default_values
   !===================================================================================================================!
   !
   !===================================================================================================================!
   if ( proc == 0 ) then
      write(buffer,'(A)') "sed -e '/^!/d;/^$/d;s/!.*//g;s/[ \t]//g' "//trim(filename)//" > input.post"
      call system( buffer )
   end if
   open(10,file='input.post',form='formatted',status='old')
   read(10,nml=list_input)
   close(10)
END SUBROUTINE
Subroutine print_all() ! afficher toutes les valeurs lues
 Use m_common
 Use m_model
 implicit none
 write(*, nml=list_input)
 print *, "boundary conditions : "
 print *, bc%nb, bc%nb_in, bc%nb_out
End Subroutine
SUBROUTINE Mesh_Input(mesh)
   USE m_common
   USE m_mesh
   USE m_time_screen
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ) , intent(inout) :: mesh
   !===================================================================================================================!
   ! Calling appropriate Subroutine
   !===================================================================================================================!
   select case (mesh_type)
      case default
    ! call Stopping_Program( 'unknow_mesh' )
      case( 'basic' )
         call Create_Cartesian_Mesh( mesh )
         case( 'dassflow' )
            call Read_Dass_Mesh( mesh )
         case( 'gmsh' )
            call Read_Gmsh_Mesh( mesh )
   end select
write(*,*) " call Read_Dass_Mesh( mesh ) done "
   !===================================================================================================================!
   ! Mesh Partition
   !===================================================================================================================!
write(*,*) " call Mesh_Partition_Scotch( mesh ) "
   ! call Mesh_Partition_Scotch( mesh )
   !===================================================================================================================!
   ! Compute geometrical mesh properties Schemes needed
   !===================================================================================================================!
write(*,*) " call Mesh_Geometric_Properties( mesh ) "
   call Mesh_Geometric_Properties( mesh )
write(*,*) " done mesh input "
do i=1,mesh%neb
  write(*,*) mesh%edgeb(i)%typlim
end do
END SUBROUTINE Mesh_Input
SUBROUTINE Create_Cartesian_Mesh( mesh )
   USE m_common
   USE m_mesh
   USE m_time_screen
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(inout) :: mesh
   !===================================================================================================================!
   ! Number of Cells, Nodes and Edges
   !===================================================================================================================!
   mesh%nn = nx * ny
   mesh%ne = nx * (ny-1) + ny * (nx-1)
   mesh%nc = (nx-1) * (ny-1)
   mesh%nnb = 2 * ( (nx-2) + (ny-2) ) + 4
   mesh%neb = 2 * ( (nx-1) + (ny-1) )
   mesh%ncb = mesh%neb
   allocate( mesh%node( mesh%nn ) )
   allocate( mesh%edge( mesh%ne ) )
   allocate( mesh%cell( mesh%nc ) )
   allocate( mesh%nodeb( mesh%nnb ) )
   allocate( mesh%edgeb( mesh%neb ) )
   allocate( mesh%cellb( mesh%ncb ) )
   !===================================================================================================================!
   ! Nodes variables (boundary included in global indexes)
   !===================================================================================================================!
   dx = lx / float( nx - 1 )
   dy = ly / float( ny - 1 )
   k = 0
   kb = 0
   do j = 1,ny
      do i = 1,nx
         k = k + 1
         mesh%node(k)%coord%x = float( i - 1 ) * dx
         mesh%node(k)%coord%y = float( j - 1 ) * dy
         if ( i == 1 .or. i == nx .or. j == 1 .or. j == ny ) then
            kb = kb + 1
            mesh%node(k)%boundary = .true.
            mesh%node(k)%lim = kb
            mesh%nodeb(kb)%ind = k
         else
            mesh%node(k)%boundary = .false.
            mesh%node(k)%lim = 0
         end if
      end do
   end do
   !===================================================================================================================!
   ! Edges variables (boundary included in global indexes)
   !===================================================================================================================!
   k = 0
   kb = 0
   do j = 1,ny-1
      do i = 1,nx
         k = k + 1
         mesh%edge(k)%node(1) = i + (j-1) * nx
         mesh%edge(k)%node(2) = i + j * nx
         if ( i == 1 ) then
            kb = kb + 1
            mesh%edge(k)%boundary = .true.
            mesh%edge(k)%lim = kb
            mesh%edgeb(kb)%ind = k
            mesh%edgeb(kb)%typlim = bc_W
            mesh%edgeb(kb)%group = 1
            mesh%edge(k)%cell(1) = i + (j-1) * (nx-1)
            mesh%edge(k)%cell(2) = mesh%nc + kb
            if ( bc_W == 'periodic' ) then
               mesh%cellb(kb)%cell = mesh%edge(k)%cell(1) + nx - 2
               mesh%edgeb(kb)%perio = k + nx - 1
            else
               mesh%cellb(kb)%cell = mesh%edge(k)%cell(1)
               mesh%edgeb(kb)%perio = -1
            end if
         else if ( i == nx ) then
            kb = kb + 1
            mesh%edge(k)%boundary = .true.
            mesh%edge(k)%lim = kb
            mesh%edgeb(kb)%ind = k
            mesh%edgeb(kb)%typlim = bc_E
            mesh%edgeb(kb)%group = 2
            mesh%edge(k)%cell(1) = i-1 + (j-1) * (nx-1)
            mesh%edge(k)%cell(2) = mesh%nc + kb
            if ( bc_E == 'periodic' ) then
               mesh%cellb(kb)%cell = mesh%edge(k)%cell(1) - nx + 2
               mesh%edgeb(kb)%perio = k - nx + 1
            else
               mesh%cellb(kb)%cell = mesh%edge(k)%cell(1)
               mesh%edgeb(kb)%perio = -1
            end if
         else
            mesh%edge(k)%boundary = .false.
            mesh%edge(k)%lim = 0
            mesh%edge(k)%cell(1) = i-1 + (j-1) * (nx-1)
            mesh%edge(k)%cell(2) = i + (j-1) * (nx-1)
         end if
      end do
   end do
   do i = 1,nx-1
      do j = 1,ny
         k = k + 1
         mesh%edge(k)%node(1) = i + (j-1) * nx
         mesh%edge(k)%node(2) = i+1 + (j-1) * nx
         if ( j == 1 ) then
            kb = kb + 1
            mesh%edge(k)%boundary = .true.
            mesh%edge(k)%lim = kb
            mesh%edgeb(kb)%ind = k
            mesh%edgeb(kb)%typlim = bc_S
            mesh%edgeb(kb)%group = 3
            mesh%edge(k)%cell(1) = i + j - 1
            mesh%edge(k)%cell(2) = mesh%nc + kb
            if ( bc_S == 'periodic' ) then
               mesh%cellb(kb)%cell = mesh%edge(k)%cell(1) + ( nx - 1 ) * ( ny - 2 )
               mesh%edgeb(kb)%perio = k + ny - 1
            else
               mesh%cellb(kb)%cell = mesh%edge(k)%cell(1)
               mesh%edgeb(kb)%perio = -1
            end if
         else if ( j == ny ) then
            kb = kb + 1
            mesh%edge(k)%boundary = .true.
            mesh%edge(k)%lim = kb
            mesh%edgeb(kb)%ind = k
            mesh%edgeb(kb)%typlim = bc_N
            mesh%edgeb(kb)%group = 4
            mesh%edge(k)%cell(1) = i + (j-2) * (nx-1)
            mesh%edge(k)%cell(2) = mesh%nc + kb
            if ( bc_N == 'periodic' ) then
               mesh%cellb(kb)%cell = mesh%edge(k)%cell(1) - ( nx - 1 ) * ( ny - 2 )
               mesh%edgeb(kb)%perio = k - ny + 1
            else
               mesh%cellb(kb)%cell = mesh%edge(k)%cell(1)
               mesh%edgeb(kb)%perio = -1
            end if
         else
            mesh%edge(k)%boundary = .false.
            mesh%edge(k)%lim = 0
            mesh%edge(k)%cell(1) = i + (j-2) * (nx-1)
            mesh%edge(k)%cell(2) = i + (j-1) * (nx-1)
         end if
      end do
   end do
   !===================================================================================================================!
   ! Cells variables (boundary excluded from global indexes, see below)
   !===================================================================================================================!
   k = 0
   do j = 1,ny-1
      do i = 1,nx-1
         k = k + 1
         mesh%cell(k)%node(1) = k + j - 1
         mesh%cell(k)%node(2) = k + j
         mesh%cell(k)%node(3) = k + nx + j
         mesh%cell(k)%node(4) = k + nx + j - 1
         mesh%cell(k)%cell(1) = k - (nx-1)
         mesh%cell(k)%cell(2) = k + 1
         mesh%cell(k)%cell(3) = k + (nx-1)
         mesh%cell(k)%cell(4) = k - 1
         mesh%cell(k)%edge(1) = j + nx * ( ny - 1 ) + ny * ( i - 1 )
         mesh%cell(k)%edge(2) = k + j
         mesh%cell(k)%edge(3) = mesh%cell(k)%edge(1) + 1
         mesh%cell(k)%edge(4) = mesh%cell(k)%edge(2) - 1
         if ( i == 1 .or. i == nx-1 .or. j == 1 .or. j == ny-1 ) then
            mesh%cell(k)%boundary = .true.
            if ( j == 1 ) mesh%cell(k)%cell(1) = mesh%edge( mesh%cell(k)%edge(1) )%cell(2)
            if ( i == nx-1 ) mesh%cell(k)%cell(2) = mesh%edge( mesh%cell(k)%edge(2) )%cell(2)
            if ( j == ny-1 ) mesh%cell(k)%cell(3) = mesh%edge( mesh%cell(k)%edge(3) )%cell(2)
            if ( i == 1 ) mesh%cell(k)%cell(4) = mesh%edge( mesh%cell(k)%edge(4) )%cell(2)
         else
            mesh%cell(k)%boundary = .false.
         end if
         mesh%cell(k)%nbed = 4
      end do
   end do
END SUBROUTINE Create_Cartesian_Mesh
SUBROUTINE Read_Dass_Mesh( mesh )
   USE m_common
   USE m_mesh
   USE m_model
   USE m_time_screen
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(out) :: mesh
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   integer(ip) :: neighb(4) , n1 , n2 , nc_bc , nb_grp_in , nb_grp_out , num_grp , bc_number
   real(rp) :: ghost_cell_bathy
   character(6) :: bc_type
   real(rp) :: tmp ! temporary variable to store useless data that are read (here, it is cell's nland value, that is read and applied in initialize subroutine)
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   call Print_Screen( 'start_mesh' )
   mesh%file_name = mesh_name
   write(*,*) " mesh%file_name", mesh%file_name
   open(10,file=mesh%file_name,status='old')
   !===================================================================================================================!
   ! Reading header
   !===================================================================================================================!
   read(10,*)
   read(10,*) mesh%nn , mesh%nc , mesh%scal
   allocate( mesh%cell( mesh%nc ) )
   allocate( mesh%node( mesh%nn ) )
   allocate( bathy_cell( mesh%nc ) )
   allocate( bathy_node( mesh%nn ) )
   allocate( land( mesh%nc ) )
   inquire( file = 'land_uses.txt' , exist = file_exist(1) )
   if ( file_exist(1) ) then
      open(20,file='land_uses.txt',status='old',form='formatted')
      read(20,*) ! comment line
      read(20,*) ! comment line
      read(20,*) ! comment line
      read(20,*) nland
      read(20,*) ! comment line
      read(20,*) ! comment line
      read(20,*) ! comment line
      allocate( manning( nland ) )
      allocate( manning_beta( nland ) )
      do i = 1,nland
         read(20,*) k , manning(k), manning_beta(k)
      end do
      close(20)
   else ! then file not exist
    write(*,*) "mesh type is basic and no mesh is provided, we stop here"
    STOP
   end if ! end if file exist
      manning_data_glob = 1 ! ??
   !===================================================================================================================!
   ! Reading Nodes
   !===================================================================================================================!
   read(10,*)
   do i = 1,mesh%nn
      read(10,*) k , mesh%node(k)%coord%x , mesh%node(k)%coord%y , bathy_node(k)
   end do
   !===================================================================================================================!
   ! Reading Nodes connectivity
   !===================================================================================================================!
   read(10,*)
   do i = 1,mesh%nc
      read(10,*) k , &
                   mesh%cell(k)%node(1) , &
                   mesh%cell(k)%node(2) , &
                   mesh%cell(k)%node(3) , &
                   mesh%cell(k)%node(4) , &
                   land(k) , & ! land(k) is defined in src/sw_mono/initialization.f90/my_friction_2_fortran(my_friction)
                   bathy_cell(k)
      if ( mesh%cell(k)%node(4) == 0 .or. &
           mesh%cell(k)%node(4) == mesh%cell(k)%node(3) ) mesh%cell(k)%node(4) = mesh%cell(k)%node(1)
   end do
   !===================================================================================================================!
   ! Calculating Cells connectivity + Edges count
   !===================================================================================================================!
   write(buffer,'(A,A,A)') 'msh/' , trim(mesh%file_name) , '_cells_connect.bin'
   inquire( file = trim(buffer) , exist = file_exist(1) )
   if ( .not. file_exist(1) ) then
      call calc_cells_connectivity( mesh )
      if ( proc == 0 ) then
         call system('mkdir -p msh')
         open(20,file=trim(buffer),form='unformatted',status='replace')
         do i = 1,mesh%nc
            write(20) mesh%cell(i)%nbed , mesh%cell(i)%cell(:)
         end do
         write(20) mesh%ne , mesh%neb
      end if
   else
      open(20,file=trim(buffer),form='unformatted',status='old')
      do i = 1,mesh%nc
         read(20) mesh%cell(i)%nbed , mesh%cell(i)%cell(:)
      end do
      read(20) mesh%ne , mesh%neb
      close(20)
   end if
   !===================================================================================================================!
   ! Edges variables
   !===================================================================================================================!
   allocate( mesh%edge ( mesh%ne ) )
   allocate( mesh%edgeb( mesh%neb ) )
   i = 0
   j = 0
   do k = 1,mesh%nc
      neighb(:) = mesh%cell(k)%cell(:)
      do ke = 1,4
         if ( neighb(ke) == 0 ) then
            i = i + 1
            j = j + 1
            mesh%edge(i)%boundary = .true.
            mesh%edge(i)%lim = j
            mesh%edgeb(j)%ind = i
            mesh%edge(i)%cell(1) = k
            mesh%edge(i)%cell(2) = mesh%nc + j
            mesh%edge(i)%node(1) = mesh%cell(k)%node( ke )
            mesh%edge(i)%node(2) = mesh%cell(k)%node( mod(ke,4) + 1 )
            mesh%cell(k)%edge(ke) = i
            mesh%cell(k)%cell(ke) = mesh%edge(i)%cell(2)
         else if ( neighb(ke) > k .and. neighb(ke) <= mesh%nc ) then
            i = i + 1
            mesh%edge(i)%boundary = .false.
            mesh%edge(i)%cell(1) = k
            mesh%edge(i)%cell(2) = neighb(ke)
            mesh%edge(i)%node(1) = mesh%cell(k)%node( ke )
            mesh%edge(i)%node(2) = mesh%cell(k)%node( mod(ke,4) + 1 )
            mesh%cell(k)%edge(ke) = i
         else if ( neighb(ke) == -1 ) then
            mesh%cell(k)%edge(ke) = -1
         else
            do ie = 1,4
               if ( mesh%cell( neighb(ke) )%cell(ie) == k ) exit
            end do
            mesh%cell(k)%edge(ke) = mesh%cell( neighb(ke) )%edge(ie)
         end if
      end do
   end do
   !===================================================================================================================!
   ! Creating Neighboring Boundary Cells
   !===================================================================================================================!
   mesh%ncb = mesh%neb
   allocate( mesh%cellb( mesh%ncb ) )
   do i = 1,mesh%ncb
      mesh%cellb(i)%cell = mesh%edge( mesh%edgeb(i)%ind )%cell(1)
   end do
   !===================================================================================================================!
   ! Initialization of Ghost bathymetry and Boundary Condition Type
   !===================================================================================================================!
   mesh%edgeb(:)%typlim = 'wall'
   call reallocate_r( bathy_cell , mesh%nc + mesh%ncb )
   do i = 1,mesh%ncb
      bathy_cell( mesh%nc + i ) = bathy_cell( mesh%cellb(i)%cell )
   end do
   !===================================================================================================================!
   ! Reading Boundary Conditions Type File
   !===================================================================================================================!
   call read_bc_file
   !===================================================================================================================!
   ! INLET
   !===================================================================================================================!
   read(10,*) buffer
   read(10,*) bc_type , nc_bc , nb_grp_in
   if ( nb_grp_in == 0 ) then
      do i = 1,nc_bc
         read(10,*) k , j , bc_number , ghost_cell_bathy
         mesh%edgeb( mesh%edge(ie)%lim )%typlim = bc%typ( 1 , 1 )
         mesh%edgeb( mesh%edge(ie)%lim )%group = 1
         bathy_cell( mesh%nc + mesh%edge(ie)%lim ) = ghost_cell_bathy
      end do
   else
      do i = 1,nc_bc
         read(10,*) k , j , bc_number , ghost_cell_bathy, num_grp
         ie = mesh%cell(k)%edge(j)
         mesh%edgeb( mesh%edge(ie)%lim )%typlim = bc%typ( num_grp , 1 )
         mesh%edgeb( mesh%edge(ie)%lim )%group = num_grp
         bathy_cell( mesh%nc + mesh%edge(ie)%lim ) = ghost_cell_bathy
      end do
   end if
   !===================================================================================================================!
   ! OUTLET
   !===================================================================================================================!
   read(10,*) bc_type , nc_bc , nb_grp_out
   if ( nb_grp_out == 0 ) then
      do i = 1,nc_bc
         read(10,*) k , j , bc_number , ghost_cell_bathy
         ie = mesh%cell(k)%edge(j)
         mesh%edgeb( mesh%edge(ie)%lim )%typlim = bc%typ( 2 , 1 )
         mesh%edgeb( mesh%edge(ie)%lim )%group = 2
         bathy_cell( mesh%nc + mesh%edge(ie)%lim ) = ghost_cell_bathy
      end do
   else
      do i = 1,nc_bc
         read(10,*) k , j , bc_number , ghost_cell_bathy , num_grp
         ie = mesh%cell(k)%edge(j)
         mesh%edgeb( mesh%edge(ie)%lim )%typlim = bc%typ( nb_grp_in + num_grp , 1 )
         mesh%edgeb( mesh%edge(ie)%lim )%group = nb_grp_in + num_grp
         bathy_cell( mesh%nc + mesh%edge(ie)%lim ) = ghost_cell_bathy
      end do
   end if
   !===================================================================================================================!
   ! Finishing Subroutine
   !===================================================================================================================!
   close(10)
   close(20)
   call Print_Screen( 'end_mesh' )
END SUBROUTINE Read_Dass_Mesh
SUBROUTINE Read_Gmsh_Mesh( mesh )
   USE m_common
   USE m_mesh
   USE m_model
   USE m_time_screen
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(out) :: mesh
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   integer(ip) :: typ , nb_tags , tags(4) , nb_elt
   integer(ip) :: neighb(4) , n1 , n2 , n3 , n4 , nc_bc , nb_group , bc_number
   real(rp) :: ghost_cell_bathy
   character(6) :: bc_type
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   call Print_Screen( 'start_mesh' )
   mesh%file_name = mesh_name
   open(10,file=mesh%file_name,status='old')
   !===================================================================================================================!
   ! Reading Number of Boundary Conditions
   !===================================================================================================================!
   buffer = ''
   do while ( trim(buffer) /= '$PhysicalNames' ) ; read(10,*) buffer ; end do
   read(10,*) nc_bc
   !===================================================================================================================!
   ! Reading Boundary Conditions Type File
   !===================================================================================================================!
   call read_bc_file
   !===================================================================================================================!
   ! Reading Nodes + eventually bathy_node (z coordinate)
   !===================================================================================================================!
   do while ( trim(buffer) /= '$Nodes' ) ; read(10,*) buffer ; end do
   read(10,*) mesh%nn
   allocate( mesh%node ( mesh%nn ) )
   allocate( bathy_node( mesh%nn ) )
   do i = 1,mesh%nn
      read(10,*) k , mesh%node(k)%coord%x , mesh%node(k)%coord%y , bathy_node(k)
   end do
   !===================================================================================================================!
   ! Reading Gmsh Elements
   !===================================================================================================================!
   do while ( trim(buffer) /= '$Elements' ) ; read(10,*) buffer ; end do
   read(10,*) nb_elt
   allocate( mesh%cell ( nb_elt ) )
   allocate( mesh%cellb( nb_elt ) )
   allocate( mesh%edge ( maxed * nb_elt ) )
   allocate( mesh%edgeb( nb_elt ) )
   mesh%nc = 0
   mesh%ne = 0
   mesh%neb = 0
   do i = 1,nb_elt
      mesh%cell(i)%node(:) = -1
      mesh%edge(i)%node(:) = -1
      mesh%edge(i)%lim = -1
      read(10,'(a)') buffer ; read(buffer,*) k , typ , nb_tags
      if ( typ == 1 ) then
         mesh%ne = mesh%ne + 1
         mesh%neb = mesh%neb + 1
         read(buffer,*) k , typ , nb_tags , tags(1:nb_tags) , mesh%edge(mesh%ne)%node(1:2)
         mesh%edge( mesh%ne )%boundary = .true.
         mesh%edge( mesh%ne )%lim = mesh%neb
         mesh%edgeb( mesh%neb )%ind = mesh%ne
         mesh%edgeb( mesh%neb )%group = tags(1)
         mesh%edgeb( mesh%neb )%typlim = bc%typ(tags(1),1)
      else if ( typ == 2 ) then
         mesh%nc = mesh%nc + 1
         mesh%cell(mesh%nc)%nbed = 3
         read(buffer,*) k , typ , nb_tags , tags(1:nb_tags) , mesh%cell(mesh%nc)%node(1:3)
      else if ( typ == 3 ) then
         mesh%nc = mesh%nc + 1
         mesh%cell(mesh%nc)%nbed = 4
         read(buffer,*) k , typ , nb_tags , tags(1:nb_tags) , mesh%cell(mesh%nc)%node(1:4)
      else
  ! call Stopping_Program_Sub( 'Unknow Element Type in Gmsh mesh format' )
      end if
   end do
   !===================================================================================================================!
   ! Cells/Edges/Nodes Connectivities
   !===================================================================================================================!
   k = 0
   do i = 1,mesh%nc
      mesh%cell(i)%cell(:) = -1
      mesh%cell(i)%edge(:) = -1
      do ie = 1,mesh%cell(i)%nbed
         n1 = mesh%cell(i)%node( ie )
         n2 = mesh%cell(i)%node( mod( ie , mesh%cell(i)%nbed ) + 1 )
         loop_inner: do j = 1,mesh%nc
            do je = 1,mesh%cell(j)%nbed
               if ( n2 == mesh%cell(j)%node(je) ) then
                  if ( n1 == mesh%cell(j)%node( mod( je , mesh%cell(j)%nbed ) + 1 ) ) then
                     mesh%cell(i)%cell(ie) = j
                     if ( j > i ) then
                        mesh%ne = mesh%ne + 1
                        mesh%edge( mesh%ne )%boundary = .false.
                        mesh%edge( mesh%ne )%cell(1) = i
                        mesh%edge( mesh%ne )%cell(2) = j
                        mesh%edge( mesh%ne )%node(1) = n1
                        mesh%edge( mesh%ne )%node(2) = n2
                        mesh%cell(i)%edge(ie) = mesh%ne
                     else
                        mesh%cell(i)%edge(ie) = mesh%cell(j)%edge(je)
                     end if
                     exit loop_inner
                  end if
               end if
            end do
         end do loop_inner
         if ( mesh%cell(i)%cell(ie) == -1 ) then
            k = k + 1
            logic_test = .false.
            do je = 1,mesh%neb
               if ( ( n1 == mesh%edge(je)%node(1) .and. n2 == mesh%edge(je)%node(2) ) .or.&
                    ( n1 == mesh%edge(je)%node(2) .and. n2 == mesh%edge(je)%node(1) ) ) then
                  logic_test = .true.
                  exit
               end if
            end do
            if ( logic_test ) then
               mesh%edge(je)%cell(1) = i
               mesh%edge(je)%cell(2) = mesh%nc + k
               mesh%cell(i)%cell(ie) = mesh%nc + k
               mesh%cell(i)%edge(ie) = je
            else
               mesh%ne = mesh%ne + 1
               mesh%neb = mesh%neb + 1
               mesh%edge( mesh%ne )%cell(1) = i
               mesh%edge( mesh%ne )%cell(2) = mesh%nc + k
               mesh%edge( mesh%ne )%node(1) = n1
               mesh%edge( mesh%ne )%node(2) = n2
               mesh%cell(i)%cell(ie) = mesh%nc + k
               mesh%cell(i)%edge(ie) = mesh%ne
               mesh%edge( mesh%ne )%boundary = .true.
               mesh%edge( mesh%ne )%lim = mesh%neb
               mesh%edgeb( mesh%neb )%ind = mesh%ne
               mesh%edgeb( mesh%neb )%group = nc_bc + 1
               mesh%edgeb( mesh%neb )%typlim = 'wall'
            end if
         end if
      end do
      if ( mesh%cell(i)%nbed == 3 ) mesh%cell(i)%node(4) = mesh%cell(i)%node(1) !FIX
   end do
   mesh%ncb = mesh%neb
   call reallocate_cell ( mesh%cell , mesh%nc )
   call reallocate_cellb( mesh%cellb , mesh%ncb )
   call reallocate_edge ( mesh%edge , mesh%ne )
   call reallocate_edgeb( mesh%edgeb , mesh%neb )
   close(10)
   call Print_Screen( 'end_mesh' )
END SUBROUTINE Read_Gmsh_Mesh
SUBROUTINE read_bc_file
   USE m_common
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   integer(ip) :: n_hydrograph , n_ratcurve, n_internal, n_zspresc, n_hpresc! n_internal_discharg, n_internal_ratcurve
   character(20) :: temp
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   open(20,file='bc.txt',form='formatted',status='old')
   read(20,*)
   read(20,*)
   read(20,*)
   read(20,*) bc%nb
   read(20,*)
   read(20,*)
   read(20,*)
   allocate( bc%typ ( bc%nb , 3 ) )! 3 to account for new internal frontier, 3rd item in line does not matter for non-internal borders
   allocate( bc%grpf( bc%nb ) )
   bc%typ (:,:) = 'NAN'
   bc%grpf(: ) = 0
   n_hydrograph = 0
   n_ratcurve = 0
   n_zspresc = 0
   n_hpresc = 0
   n_internal = 0
   do i = 1,bc%nb
      read(20,'(a)') buffer
      read(buffer,*) k , temp
      if ( temp(1:8) == 'internal' ) then
        read(buffer,*) k , bc%typ(k,1:3) !item 3 of the line is read here only, it contains the index of the corresponding internal 1D2D interface
        n_internal = n_internal + 1
        bc%grpf(k) = n_internal
     endif
      if ( temp(1:6) /= 'transm' ) read(buffer,*) k , bc%typ(k,1:2)
      if ( temp(1:7) /= 'neumann' ) read(buffer,*) k , bc%typ(k,1:2)
      if ( temp(1:8) == 'discharg' ) then
         n_hydrograph = n_hydrograph + 1
         bc%grpf(k) = n_hydrograph
      end if
      if ( temp(1:8) == 'ratcurve' ) then
         n_ratcurve = n_ratcurve + 1
         bc%grpf(k) = n_ratcurve
      end if
      if ( temp(1:7) == 'zspresc' ) then
         n_zspresc = n_zspresc + 1
         bc%grpf(k) = n_zspresc
      end if
      if ( temp(1:6) == 'hpresc' ) then
      n_hpresc = n_hpresc +1
         bc%grpf(k) = n_hpresc
      end if
   end do
 close(20)
END SUBROUTINE
