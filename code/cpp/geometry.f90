SUBROUTINE Mesh_Geometric_Properties( mesh )
   USE m_common
   USE m_mesh
   USE m_model
   implicit none
   type( msh ), intent(inout) :: mesh
   type(vec2d) :: node_coord(4)
   integer(ip) :: nbed, connected_num_bc
   integer(ip), dimension(size(mesh%node)) :: connected_to_node
WRITE(*,*) "INTO Mesh_Geometric_Properties"
   mesh%surf = 0._rp
   connected_to_node(:) = 0
   do i = 1,size(mesh%cell)
      nbed = mesh%cell(i)%nbed
      node_coord(1:nbed) = mesh%node( mesh%cell(i)%node(1:nbed) )%coord
      if ( nbed == 3 ) then
         mesh%cell(i)%surf = 0.5_rp * abs( ( node_coord(3)%x - node_coord(2)%x ) * node_coord(1)%y + &
                                             ( node_coord(1)%x - node_coord(3)%x ) * node_coord(2)%y + &
                                             ( node_coord(2)%x - node_coord(1)%x ) * node_coord(3)%y )
      else if ( nbed == 4 ) then
         mesh%cell(i)%surf = 0.5_rp * abs( ( node_coord(4)%x - node_coord(2)%x ) * node_coord(1)%y + &
                                             ( node_coord(1)%x - node_coord(3)%x ) * node_coord(2)%y + &
                                             ( node_coord(2)%x - node_coord(4)%x ) * node_coord(3)%y + &
                                             ( node_coord(3)%x - node_coord(1)%x ) * node_coord(4)%y )
     else
  ! call Stopping_Program_Sub( 'nbed problem Mesh_Geometric_Properties' )
      end if
      mesh%surf = mesh%surf + mesh%cell(i)%surf
      mesh%cell(i)%invsurf = 1._rp / mesh%cell(i)%surf
      mesh%cell(i)%grav%x = sum( node_coord(1:nbed)%x ) / real(nbed,rp)
      mesh%cell(i)%grav%y = sum( node_coord(1:nbed)%y ) / real(nbed,rp)
      connected_to_node( mesh%cell(i)%node(1:nbed) ) = &
      connected_to_node( mesh%cell(i)%node(1:nbed) ) + 1
   end do
   do ie = 1,mesh%neb
      connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(:) ) = &
      connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(:) ) + 1
   end do
   do i = 1,size(mesh%node)
      allocate( mesh%node(i)%cell( connected_to_node(i) ) )
   end do
   connected_to_node(:) = 0
   do i = 1,size(mesh%cell)
      nbed = mesh%cell(i)%nbed
      connected_to_node( mesh%cell(i)%node(1:nbed) ) = &
      connected_to_node( mesh%cell(i)%node(1:nbed) ) + 1
      do j = 1,nbed
         mesh%node( mesh%cell(i)%node(j) )%cell( connected_to_node( mesh%cell(i)%node(j) ) ) = i
      end do
  enddo
do i = 1, mesh%ne
      if ( mesh%edge(i)%boundary ) then
        if ( mesh%edgeb(mesh%edge(i)%lim)%typlim == 'internal_2D' ) then !then change connectivity to connected 1D-like cell
                read(bc%typ (mesh%edgeb(mesh%edge(i)%lim)%group, 3 ),'(i3)') connected_num_bc !Get connectivity to 1D-like cell from bc.txt
                do ib = 1,mesh%neb
                    if ( mesh%edgeb(ib)%group == connected_num_bc ) then
                        mesh%edge(i)%cell1D2D = mesh%edge( mesh%edgeb( ib )%ind )%cell(1) !Get id of the single 1D-like cell with interface in the connected bc number
                    endif
                enddo
        endif
      endif
   end do
   do ie = 1,mesh%neb
      connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(:) ) = &
      connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(:) ) + 1
      if ( mesh%edge( mesh%edgeb(ie)%ind )%cell(2) > mesh%nc ) then
         mesh%node( mesh%edge( mesh%edgeb(ie)%ind )%node(1) )%cell( connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(1) ) ) = mesh%edge( mesh%edgeb(ie)%ind )%cell(2)
         mesh%node( mesh%edge( mesh%edgeb(ie)%ind )%node(2) )%cell( connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(2) ) ) = mesh%edge( mesh%edgeb(ie)%ind )%cell(2)
      else
         mesh%node( mesh%edge( mesh%edgeb(ie)%ind )%node(1) )%cell( connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(1) ) ) = mesh%edge( mesh%edgeb(ie)%ind )%lim + mesh%nc
         mesh%node( mesh%edge( mesh%edgeb(ie)%ind )%node(2) )%cell( connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(2) ) ) = mesh%edge( mesh%edgeb(ie)%ind )%lim + mesh%nc
      end if
   end do
   mesh%cell(:)%peri = 0._rp
   connected_to_node(:) = 0
   do ie = 1,size(mesh%edge)
      node_coord(1:2) = mesh%node( mesh%edge(ie)%node(1:2) )%coord
      mesh%edge(ie)%center%x = demi * sum( node_coord(1:2)%x )
      mesh%edge(ie)%center%y = demi * sum( node_coord(1:2)%y )
      mesh%edge(ie)%length = sqrt( ( node_coord(2)%x - node_coord(1)%x )**2 + &
                                     ( node_coord(2)%y - node_coord(1)%y )**2 )
      mesh%edge(ie)%normal%x = ( node_coord(1)%y - node_coord(2)%y ) / mesh%edge(ie)%length
      mesh%edge(ie)%normal%y = ( node_coord(2)%x - node_coord(1)%x ) / mesh%edge(ie)%length
      mesh%edge(ie)%tangent%x = ( node_coord(2)%x - node_coord(1)%x ) / mesh%edge(ie)%length
      mesh%edge(ie)%tangent%y = ( node_coord(2)%y - node_coord(1)%y ) / mesh%edge(ie)%length
      if ( .not. mesh%edge(ie)%boundary ) then
         mesh%cell( mesh%edge(ie)%cell(1:2) )%peri = &
         mesh%cell( mesh%edge(ie)%cell(1:2) )%peri + mesh%edge(ie)%length
         mesh%edge(ie)%vcell%x = mesh%cell( mesh%edge(ie)%cell(2) )%grav%x - &
                                 mesh%cell( mesh%edge(ie)%cell(1) )%grav%x
         mesh%edge(ie)%vcell%y = mesh%cell( mesh%edge(ie)%cell(2) )%grav%y - &
                                 mesh%cell( mesh%edge(ie)%cell(1) )%grav%y
         if ( ( mesh%edge(ie)%vcell .dotprod. mesh%edge(ie)%normal ) < 0._rp ) then
            mesh%edge(ie)%normal%x = - mesh%edge(ie)%normal%x
            mesh%edge(ie)%normal%y = - mesh%edge(ie)%normal%y
         end if
         mesh%edge(ie)%v_edge_cell(1)%x = mesh%edge(ie)%center%x - &
                                          mesh%cell( mesh%edge(ie)%cell(1) )%grav%x
         mesh%edge(ie)%v_edge_cell(1)%y = mesh%edge(ie)%center%y - &
                                          mesh%cell( mesh%edge(ie)%cell(1) )%grav%y
         mesh%edge(ie)%v_edge_cell(2)%x = mesh%edge(ie)%center%x - &
                                          mesh%cell( mesh%edge(ie)%cell(2) )%grav%x
         mesh%edge(ie)%v_edge_cell(2)%y = mesh%edge(ie)%center%y - &
                                          mesh%cell( mesh%edge(ie)%cell(2) )%grav%y
      else
         mesh%cell( mesh%edge(ie)%cell(1) )%boundary = .true.
         mesh%cell( mesh%edge(ie)%cell(1) )%peri = &
         mesh%cell( mesh%edge(ie)%cell(1) )%peri + mesh%edge(ie)%length
         mesh%edge(ie)%vcell%x = 2._rp * ( mesh%edge(ie)%center%x - mesh%cell( mesh%edge(ie)%cell(1) )%grav%x )
         mesh%edge(ie)%vcell%y = 2._rp * ( mesh%edge(ie)%center%y - mesh%cell( mesh%edge(ie)%cell(1) )%grav%y )
         if ( ( mesh%edge(ie)%vcell .dotprod. mesh%edge(ie)%normal ) < 0._rp ) then
            mesh%edge(ie)%normal%x = - mesh%edge(ie)%normal%x
            mesh%edge(ie)%normal%y = - mesh%edge(ie)%normal%y
         end if
         mesh%edge(ie)%v_edge_cell(1)%x = mesh%edge(ie)%center%x - &
                                          mesh%cell( mesh%edge(ie)%cell(1) )%grav%x
         mesh%edge(ie)%v_edge_cell(1)%y = mesh%edge(ie)%center%y - &
                                          mesh%cell( mesh%edge(ie)%cell(1) )%grav%y
         ib = mesh%edge(ie)%lim
         mesh%cellb( ib )%grav%x = mesh%cell( mesh%edge(ie)%cell(1) )%grav%x + mesh%edge(ie)%vcell%x
         mesh%cellb( ib )%grav%y = mesh%cell( mesh%edge(ie)%cell(1) )%grav%y + mesh%edge(ie)%vcell%y
      end if
      connected_to_node( mesh%edge(ie)%node(1:2) ) = &
      connected_to_node( mesh%edge(ie)%node(1:2) ) + 1
   end do
   do ib = 1,size(mesh%edgeb)
      if ( mesh%edgeb(ib)%typlim == 'periodic' ) then
         mesh%edge( mesh%edgeb(ib)%ind )%v_edge_cell(2)%x = mesh%edge( mesh%edgeb(ib)%perio )%center%x - &
                                                            mesh%cell( mesh%cellb(ib)%cell )%grav%x
         mesh%edge( mesh%edgeb(ib)%ind )%v_edge_cell(2)%y = mesh%edge( mesh%edgeb(ib)%perio )%center%y - &
                                                            mesh%cell( mesh%cellb(ib)%cell )%grav%y
      end if
   end do
   do i = 1,size(mesh%node)
      allocate( mesh%node(i)%edge( connected_to_node(i) ) )
   end do
   connected_to_node(:) = 0
   do ie = 1,size(mesh%edge)
      connected_to_node( mesh%edge(ie)%node(1:2) ) = &
      connected_to_node( mesh%edge(ie)%node(1:2) ) + 1
      mesh%node( mesh%edge(ie)%node(1) )%edge( connected_to_node( mesh%edge(ie)%node(1) ) ) = ie
      mesh%node( mesh%edge(ie)%node(2) )%edge( connected_to_node( mesh%edge(ie)%node(2) ) ) = ie
   end do
WRITE(*,*) "OUT Mesh_Geometric_Properties"
END SUBROUTINE Mesh_Geometric_Properties
real(kind=rp) FUNCTION muscl_interp( var , grad_var , point , i_cell , mesh )
   USE m_mesh
   implicit none
   TYPE(msh), intent(in) :: mesh
   real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
   TYPE(vec2d), dimension( mesh%nc ), intent(in) :: grad_var
   TYPE(vec2d), intent(in) :: point
   integer(ip), intent(in) :: i_cell
   muscl_interp = var(i_cell) + ( grad_var(i_cell)%x * ( point%x - mesh%cell(i_cell)%grav%x ) + &
                                    grad_var(i_cell)%y * ( point%y - mesh%cell(i_cell)%grav%y ) )
END FUNCTION muscl_interp
real(kind=rp) FUNCTION Shepard_interp( var , point , i_cell , mesh )
   USE m_mesh
   implicit none
   TYPE( msh ), intent(in) :: mesh
   real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
   TYPE( vec2d ), intent(in) :: point
   integer(ip), intent(in) :: i_cell
   real(rp), parameter :: shepard_coef = 1.5_rp ! Power Law Basis Function ( > 0.5 and <= 1.5 )
   real(rp) :: basis_func_val , sum_d , sum_n
   integer(ip) :: neib , i_celle , i_celln
   basis_func_val = one / ( zerom + ( point%x - mesh%cell(i_cell)%grav%x )**2 + &
                                    ( point%y - mesh%cell(i_cell)%grav%y )**2 )**shepard_coef
   sum_n = var(i_cell) * basis_func_val
   sum_d = basis_func_val
   do i_celle = 1,maxed
      i_celln = mesh%cell(i_cell)%cell(i_celle)
      basis_func_val = one / ( zerom + ( point%x - mesh%cell(i_celln)%grav%x )**2 + &
                                       ( point%y - mesh%cell(i_celln)%grav%y )**2 )**shepard_coef
      sum_n = sum_n + var(i_celln) * basis_func_val
      sum_d = sum_d + basis_func_val
   end do
   Shepard_interp = sum_n / sum_d
END FUNCTION Shepard_interp
logical FUNCTION point_in_cell( point , cell_ind , mesh )
   USE m_mesh
   implicit none
   TYPE( vec2d ), intent(in) :: point
   integer(ip), intent(in) :: cell_ind
   TYPE( msh ), intent(in) :: mesh
   integer(ip) :: cnt , n1 , n2
   logical :: ray_inter
   type( vec2d ) :: point_offset
   cnt = 0
   point_offset = point * ( 1._rp - 10._rp * zerom )
   do ie = 1,maxed
      n1 = mesh%cell( cell_ind )%node( ie )
      n2 = mesh%cell( cell_ind )%node( mod( ie , maxed ) + 1 )
      if ( point_offset%x == mesh%node(n1)%coord%x .and. point_offset%y == mesh%node(n1)%coord%y ) then
         point_in_cell = .true. ; return
      else if ( ray_inter( point_offset , mesh%node(n1)%coord , mesh%node(n2)%coord ) ) then
         cnt = cnt + 1
      end if
   end do
   if ( mod( cnt , 2 ) == 0 ) then
      point_in_cell = .false.
   else
      point_in_cell = .true.
   end if
END FUNCTION point_in_cell
logical FUNCTION ray_inter( p , a0 , b0 )
   USE m_mesh
   implicit none
   type(vec2d), intent(in) :: p , a0 , b0
   type(vec2d) :: a , b
   real(rp) :: m_red , m_blue
   ray_inter = .false.
   if ( ( p%y /= a0%y ) .and. ( p%y /= b0%y ) ) then
      if ( a0%y > b0%y ) then
         b = a0
         a = b0
      else
         a = a0
         b = b0
      end if
      if ( ( p%y > b%y ) .or. ( p%y < a%y ) ) return
      if ( p%x > max( a%x , b%x ) ) return
      if ( p%x < min( a%x , b%x ) ) then
         ray_inter = .true.
      else
         if ( abs( a%x - b%x ) > tiny( a%x ) ) then
            m_red = ( b%y - a%y ) / ( b%x - a%x )
         else
            m_red = huge( m_red )
         end if
         if ( abs( a%x - p%x ) > tiny( a%x ) ) then
            m_blue = ( p%y - a%y ) / ( p%x - a%x )
         else
            m_blue = huge( m_blue )
         end if
         if ( m_blue >= m_red ) then
            ray_inter = .true.
         else
            ray_inter = .false.
         end if
      end if
   else if ( ( p%x /= a0%x ) .and. ( p%x /= b0%x ) ) then
      if ( a0%x > b0%x ) then
         b = a0
         a = b0
      else
         a = a0
         b = b0
      end if
      if ( ( p%x > b%x ) .or. ( p%x < a%x ) ) return
      if ( p%y > max( a%y , b%y ) ) return
      if ( p%y < min( a%y , b%y ) ) then
         ray_inter = .true.
      else
         if ( abs( a%y - b%y ) > tiny( a%y ) ) then
            m_red = ( b%x - a%x ) / ( b%y - a%y )
         else
            m_red = huge( m_red )
         end if
         if ( abs( a%y - p%y ) > tiny( a%y ) ) then
            m_blue = ( p%x - a%x ) / ( p%y - a%y )
         else
            m_blue = huge( m_blue )
         end if
         if ( m_blue >= m_red ) then
            ray_inter = .true.
         else
            ray_inter = .false.
         end if
      end if
   end if
END FUNCTION ray_inter
FUNCTION search_cell_inc_point( mesh , point ) RESULT( cell )
   USE m_mesh
   implicit none
   type( msh ), intent(in) :: mesh
   type( vec2d ), intent(in) :: point
   integer(ip) :: cell , icell
   logical :: point_in_cell
   cell = -1
   do icell = 1,mesh%nc
      if ( point_in_cell( point , icell , mesh ) ) then
         cell = icell
         exit
      end if
   end do
END FUNCTION search_cell_inc_point
