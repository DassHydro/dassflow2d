FUNCTION vf4_normal_discrete_gradient_scal( var , mesh ) RESULT( grad_var )
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(in) :: mesh
   real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
   real(rp), dimension( mesh%ne ) :: grad_var
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   integer(ip) :: iL , iR
   real(rp) :: dist
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   do ie = 1,mesh%ne
      iL = mesh%edge(ie)%cell(1)
      iR = mesh%edge(ie)%cell(2)
      dist = .norm. mesh%edge(ie)%vcell
      grad_var(ie) = ( var( iR ) - var( iL ) ) / dist
   end do
END FUNCTION vf4_normal_discrete_gradient_scal
FUNCTION vf4_normal_discrete_gradient_vec2d( var , mesh ) RESULT( grad_var )
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(in) :: mesh
   type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
   type( vec2d ), dimension( mesh%ne ) :: grad_var
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   integer(ip) :: iL , iR
   real(rp) :: dist
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   do ie = 1,mesh%ne
      iL = mesh%edge(ie)%cell(1)
      iR = mesh%edge(ie)%cell(2)
      dist = .norm. mesh%edge(ie)%vcell
      grad_var(ie) = ( var( iR ) - var( iL ) ) / dist
   end do
END FUNCTION vf4_normal_discrete_gradient_vec2d
FUNCTION diamond_discrete_gradient_scal( var , mesh ) RESULT( grad_var )
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(in) :: mesh
   real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
   type( vec2d ), dimension( mesh%ne ) :: grad_var
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   real(rp), dimension( mesh%nn ) :: var_node
   real(rp) :: dist(2) , correction , grad(2)
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   var_node = Least_Square_Cell_to_Node( var , mesh )
   do ie = 1,mesh%ne
      dist(1) = .norm. mesh%edge(ie)%vcell
      correction = ( mesh%edge(ie)%vcell .dotprod. mesh%edge(ie)%normal ) / dist(1)
      if ( correction + 10._rp * zerom >= one ) then
         correction = 1._rp
      else
         correction = acos( correction )
      end if
      dist(1) = dist(1) * correction
      dist(2) = mesh%edge(ie)%length * correction
      grad(1) = ( var ( mesh%edge(ie)%cell(2) ) - var ( mesh%edge(ie)%cell(1) ) ) / dist(1)
      grad(2) = ( var_node( mesh%edge(ie)%node(2) ) - var_node( mesh%edge(ie)%node(1) ) ) / dist(2)
      grad_var(ie) = grad(1) * mesh%edge(ie)%normal + &
                       grad(2) * mesh%edge(ie)%tangent
   end do
END FUNCTION diamond_discrete_gradient_scal
FUNCTION diamond_discrete_gradient_vec2d( var , mesh ) RESULT( grad_var )
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(in) :: mesh
   type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
   type( tens2d ), dimension( mesh%ne ) :: grad_var
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   type(vec2d), dimension( mesh%nn ) :: var_node
   real(rp) :: dist(2) , correction
   type(vec2d) :: grad(2)
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   var_node%x = Least_Square_Cell_to_Node( var%x , mesh )
   var_node%y = Least_Square_Cell_to_Node( var%y , mesh )
   do ie = 1,mesh%ne
      dist(1) = .norm. mesh%edge(ie)%vcell
      correction = ( mesh%edge(ie)%vcell .dotprod. mesh%edge(ie)%normal ) / dist(1)
      if ( correction + 10._rp * zerom >= one ) then
         correction = 1._rp
      else
         correction = acos( correction )
      end if
      dist(1) = dist(1) * correction
      dist(2) = mesh%edge(ie)%length * correction
      grad(1) = ( var ( mesh%edge(ie)%cell(2) ) - var ( mesh%edge(ie)%cell(1) ) ) / dist(1)
      grad(2) = ( var_node( mesh%edge(ie)%node(2) ) - var_node( mesh%edge(ie)%node(1) ) ) / dist(2)
      grad_var(ie) = ( grad(1) .tensprod. mesh%edge(ie)%normal ) + &
                       ( grad(2) .tensprod. mesh%edge(ie)%tangent )
   end do
END FUNCTION diamond_discrete_gradient_vec2d
FUNCTION Least_Square_Cell_to_Node( cell , mesh ) RESULT( node )
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(in) :: mesh
   real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: cell
   real(rp), dimension( mesh%nn ) :: node
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   do k = 1,mesh%nn
      node(k) = sum( cell_to_node(k)%weights(:) * cell( mesh%node(k)%cell(:) ) )
   end do
END FUNCTION Least_Square_Cell_to_Node
SUBROUTINE Init_Cell_to_Node( mesh )
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(in) :: mesh
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   integer(ip) :: nb_cell
   real(rp) :: Rx , Ry , Ixx , Iyy , Ixy , D , Lx , Ly , xK , yK
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   allocate( cell_to_node( size( mesh%node ) ) )
   do k = 1,size(mesh%node)
      nb_cell = size( mesh%node(k)%cell(:) )
      allocate( cell_to_node(k)%weights( nb_cell ) )
      if ( mesh%node(k)%boundary .and. nb_cell == 3 ) then
         cell_to_node(k)%weights(1) = 1._rp
         cell_to_node(k)%weights(2) = 0._rp
         cell_to_node(k)%weights(3) = 0._rp
      else
         Rx = 0._rp
         Ry = 0._rp
         Ixx = 0._rp
         Iyy = 0._rp
         Ixy = 0._rp
         do i = 1,nb_cell
            iK = mesh%node(k)%cell(i)
            if ( iK <= mesh%nc ) then
               xK = mesh%node(k)%coord%x - mesh%cell(iK)%grav%x
               yK = mesh%node(k)%coord%y - mesh%cell(iK)%grav%y
            else
               xK = mesh%node(k)%coord%x - mesh%cellb(iK-mesh%nc)%grav%x
               yK = mesh%node(k)%coord%y - mesh%cellb(iK-mesh%nc)%grav%y
            end if
            Rx = Rx + xK
            Ry = Ry + yK
            Ixx = Ixx + xK**2
            Iyy = Iyy + yK**2
            Ixy = Ixy + xK * yK
         end do
         if ( abs( Rx ) < zerom ) Rx = 0._rp
         if ( abs( Ry ) < zerom ) Ry = 0._rp
         if ( abs( Ixy ) < zerom ) Ixy = 0._rp
         D = Ixx * Iyy - Ixy**2
         Lx = ( Ixy * Ry - Iyy * Rx ) / D
         Ly = ( Ixy * Rx - Ixx * Ry ) / D
         do i = 1,nb_cell
            iK = mesh%node(k)%cell(i)
            if ( iK <= mesh%nc ) then
               xK = mesh%node(k)%coord%x - mesh%cell(iK)%grav%x
               yK = mesh%node(k)%coord%y - mesh%cell(iK)%grav%y
            else
               xK = mesh%node(k)%coord%x - mesh%cellb(iK-mesh%nc)%grav%x
               yK = mesh%node(k)%coord%y - mesh%cellb(iK-mesh%nc)%grav%y
            end if
            cell_to_node(k)%weights(i) = ( 1._rp + Lx * xK + Ly * yK ) / ( real(nb_cell,rp) + Lx * Rx + Ly * Ry )
         end do
      end if
   end do! ; STOP
END SUBROUTINE Init_Cell_to_Node
