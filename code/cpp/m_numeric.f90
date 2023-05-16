MODULE m_numeric
   USE m_common
   USE m_linear_algebra
   USE m_mesh
   USE m_time_screen !NOADJ
   implicit none
   !===================================================================================================================!
   ! Weights to perform MUSCL Scheme arising from Least Square Problem
   !===================================================================================================================!
   type( weights ), dimension(:), allocatable :: muscl
   !===================================================================================================================!
   ! Weights to reconstruct a node value from connected cell values
   ! ( Resulting from the resolution of a Least Square Problem )
   !===================================================================================================================!
   type( weights ), dimension(:), allocatable :: cell_to_node
CONTAINS
SUBROUTINE dealloc_m_numeric
      ! m_numeric.f90
      if ( allocated( muscl ) ) deallocate(muscl)
      if ( allocated( cell_to_node ) ) deallocate(cell_to_node)
END SUBROUTINE dealloc_m_numeric
   !> Linear interpolation
   !!
   !! \details linear interpolation.
   !! \param[in] tx Array
   !! \param[in] ty Array
   !! \param[in] x Value of x
   !! \return Value of y
   real(rp) FUNCTION linear_interp( tx , ty , x )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      real(rp), dimension(:), intent(in) :: tx , ty
      real(rp) , intent(in) :: x
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      real(rp) :: alpha
      integer(ip) :: ind , imin, imax
      !================================================================================================================!
      ! Begin Function
      !================================================================================================================!
      imin = lbound( tx(:) , dim=1 )
      imax = ubound( tx(:) , dim=1 )
      if ( x + zerom < tx( imin ) ) then
         alpha = ( x - tx( imin ) ) / &
                 ( tx( imin + 1 ) - tx( imin ) )
         linear_interp = ( one - alpha ) * ty( imin ) + &
                                 alpha * ty( imin + 1 )
      else if ( x + zerom < tx( imax ) ) then
         ind = imin + 1
         do while ( tx( ind ) < x + zerom ) ; ind = ind + 1 ; end do
         alpha = ( x - tx( ind - 1 ) ) / &
                 ( tx( ind ) - tx( ind - 1 ) )
         linear_interp = ( one - alpha ) * ty( ind - 1 ) + &
                                 alpha * ty( ind )
      else
         alpha = ( x - tx( imax - 1 ) ) / &
                 ( tx( imax ) - tx( imax - 1 ) )
         linear_interp = ( one - alpha ) * ty( imax - 1 ) + &
                                 alpha * ty( imax )
      end if
   END FUNCTION linear_interp
   !> \brief Calculate the Cell Gradient of a variable
   !! \details Calculate the Cell Gradient of a variable using Green Formula and basic interpolation at edge
   SUBROUTINE FV_Cell_Grad( grad_var , var , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in) :: mesh
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
      type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(out) :: grad_var
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      integer(ip) :: iL , iR ! Left and Right cells indexes to edge
      !================================================================================================================!
      ! Begin Subroutine
      !================================================================================================================!
      grad_var(:)%x = 0._rp
      grad_var(:)%y = 0._rp
      do ie = 1,mesh%ne
         iL = mesh%edge(ie)%cell(1)
         iR = mesh%edge(ie)%cell(2)
         grad_var( iL )%x = grad_var( iL )%x + mesh%edge(ie)%length * mesh%edge(ie)%normal%x * &
                                                   demi * ( var( iL ) + var( iR ) )
         grad_var( iL )%y = grad_var( iL )%y + mesh%edge(ie)%length * mesh%edge(ie)%normal%y * &
                                                   demi * ( var( iL ) + var( iR ) )
         if ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain ) then
            grad_var( iR )%x = grad_var( iR )%x - mesh%edge(ie)%length * mesh%edge(ie)%normal%x * &
                                                      demi * ( var( iL ) + var( iR ) )
            grad_var( iR )%y = grad_var( iR )%y - mesh%edge(ie)%length * mesh%edge(ie)%normal%y * &
                                                      demi * ( var( iL ) + var( iR ) )
         end if
      end do
      do i = 1,mesh%nc
         grad_var(i)%x = grad_var(i)%x * mesh%cell(i)%invsurf
         grad_var(i)%y = grad_var(i)%y * mesh%cell(i)%invsurf
      end do
   END SUBROUTINE FV_Cell_Grad
                                                                                                                 !<NOADJ
   include 'discrete_operators.inc'
   !> \brief Initialization of Scheme
   !! \details Arrays Initialization of Geometrical Mesh Properties to Construct Schemes
   SUBROUTINE Init_Schemes( mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in) :: mesh
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      real(rp) :: sum_for_muscl( 3 ) , muscl_det
      !================================================================================================================!
      ! Identity Matrix
      !================================================================================================================!
      Id2d%xx = 1._rp
      Id2d%yy = 1._rp
      Id2d%xy = 0._rp
      Id2d%yx = 0._rp
      Id3d%xx = 1._rp
      Id3d%yy = 1._rp
      Id3d%zz = 1._rp
      Id3d%xy = 0._rp
      Id3d%xz = 0._rp
      Id3d%yx = 0._rp
      Id3d%yz = 0._rp
      Id3d%zx = 0._rp
      Id3d%zy = 0._rp
      !================================================================================================================!
      ! Filling Weights to perform MUSCL Scheme arising from Least Square Problem
      !================================================================================================================!
      ! if ( spatial_scheme(1:5) == 'muscl' ) then ! ALWAYS DO THIS
         allocate( muscl( mesh%nc ) )
         do i = 1,mesh%nc
            k = mesh%cell(i)%nbed
            allocate( muscl( i )%weights( 2 * k ) )
            sum_for_muscl(:) = 0._rp
            do j = 1,k
               ie = mesh%cell(i)%edge(j)
               sum_for_muscl(1) = sum_for_muscl(1) + mesh%edge(ie)%vcell%x * mesh%edge(ie)%vcell%x
               sum_for_muscl(2) = sum_for_muscl(2) + mesh%edge(ie)%vcell%y * mesh%edge(ie)%vcell%y
               sum_for_muscl(3) = sum_for_muscl(3) + mesh%edge(ie)%vcell%x * mesh%edge(ie)%vcell%y
            end do
            muscl_det = sum_for_muscl(1) * sum_for_muscl(2) - sum_for_muscl(3)**2
          ! if ( muscl_det < zerom ) call Stopping_Program_Sub( 'Init_Schemes -> Matrice Determinant Problem' )
            muscl(i)%weights(:) = 0._rp
            do j = 1,k
               ie = mesh%cell(i)%edge(j)
               if ( mesh%edge(ie)%cell(1) == i ) then
                  muscl(i)%weights(j ) = sum_for_muscl(2) * mesh%edge(ie)%vcell%x - &
                                            sum_for_muscl(3) * mesh%edge(ie)%vcell%y
                  muscl(i)%weights(k+j) = sum_for_muscl(1) * mesh%edge(ie)%vcell%y - &
                                            sum_for_muscl(3) * mesh%edge(ie)%vcell%x
               else
                  muscl(i)%weights(j ) = sum_for_muscl(3) * mesh%edge(ie)%vcell%y - &
                                            sum_for_muscl(2) * mesh%edge(ie)%vcell%x
                  muscl(i)%weights(k+j) = sum_for_muscl(3) * mesh%edge(ie)%vcell%x - &
                                            sum_for_muscl(1) * mesh%edge(ie)%vcell%y
               end if
            end do
            muscl(i)%weights(:) = muscl(i)%weights(:) / muscl_det
         end do
    ! end if ! end if ( spatial_scheme(1:5) == 'muscl' ) then ! ALWAYS DO THIS
   END SUBROUTINE
   !> \brief Least Square Cell Slope
   !! \details Least Square Cell Slope given mesh and var from which we calculate the slope
   SUBROUTINE Least_Square_Slope( var_slope , var , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in) :: mesh
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
      type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(out) :: var_slope
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      real(rp) :: var_cell
      real(rp), dimension(maxed) :: var_diff
      !================================================================================================================!
      ! Begin Subroutine
      !================================================================================================================!
      do i = 1,mesh%nc
         k = mesh%cell(i)%nbed
         var_cell = var(i)
         var_diff(1:k) = var( (/ mesh%cell(i)%cell(1:k) /) ) - var_cell
         var_slope(i)%x = sum( muscl(i)%weights(1 : k) * var_diff(1:k) )
         var_slope(i)%y = sum( muscl(i)%weights(1+k:2*k) * var_diff(1:k) )
      end do
   END SUBROUTINE Least_Square_Slope
   !> \brief Least Square Cell Slope with Filter
   !! \details Least Square Cell Slope with Filter given mesh and var from which we calculate the slope + filter and cutoff
   SUBROUTINE Least_Square_Slope_Filter( var_slope , var , filter , cutoff , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in) :: mesh
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var , filter
      type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(out) :: var_slope
      real(rp), intent(in) :: cutoff
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      real(rp) :: var_cell
      real(rp), dimension(maxed) :: var_diff
      !================================================================================================================!
      ! Begin Subroutine
      !================================================================================================================!
      do i = 1,mesh%nc
         if ( filter(i) > cutoff ) then
            k = mesh%cell(i)%nbed
            var_cell = var(i)
            var_diff(1:k) = var( (/ mesh%cell(i)%cell(1:k) /) ) - var_cell
            var_slope(i)%x = sum( muscl(i)%weights(1 : k) * var_diff(1:k) )
            var_slope(i)%y = sum( muscl(i)%weights(1+k:2*k) * var_diff(1:k) )
         else
            var_slope(i)%x = 0._rp
            var_slope(i)%y = 0._rp
         end if
      end do
   END SUBROUTINE Least_Square_Slope_Filter
   !> \brief Least Square Cell Slope with Filter
   !! \details Least Square Cell Slope with DOUBLE Filter given mesh and var from which we calculate the slope + filter and cutoff
   SUBROUTINE Least_Square_Slope_Double_Filter( var_slope , var , filter , cutoff , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in) :: mesh
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var , filter
      type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(out) :: var_slope
      real(rp), intent(in) :: cutoff
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      real(rp) :: var_cell
      real(rp), dimension(maxed) :: var_diff
      !================================================================================================================!
      ! Begin Subroutine
      !================================================================================================================!
      do i = 1,mesh%nc
         if ( filter(i) > cutoff ) then
            k = mesh%cell(i)%nbed
            var_cell = var(i)
            var_diff(1:k) = var( (/ mesh%cell(i)%cell(1:k) /) ) - var_cell
            var_slope(i)%x = sum( muscl(i)%weights(1 : k) * var_diff(1:k) , &
                                    filter( (/ mesh%cell(i)%cell(1:k) /) ) > cutoff )
            var_slope(i)%y = sum( muscl(i)%weights(1+k:2*k) * var_diff(1:k) , &
                                    filter( (/ mesh%cell(i)%cell(1:k) /) ) > cutoff )
         else
            var_slope(i)%x = 0._rp
            var_slope(i)%y = 0._rp
         end if
      end do
   END SUBROUTINE Least_Square_Slope_Double_Filter
   !> \brief Least Square Cell Slope with Filter
   !! \details Least Square Cell Slope with SUPER Filter given mesh and var from which we calculate the slope + filter and cutoff
   SUBROUTINE Least_Square_Slope_Super_Filter( var_slope , var , filter , cutoff , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in) :: mesh
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var , filter
      type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(out) :: var_slope
      real(rp), intent(in) :: cutoff
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      real(rp) :: var_cell
      real(rp), dimension(maxed) :: var_diff
      !================================================================================================================!
      ! Begin Subroutine
      !================================================================================================================!
      do i = 1,mesh%nc
         k = mesh%cell(i)%nbed
         if ( all( filter( (/ i , mesh%cell(i)%cell(1:k) /) ) > cutoff ) ) then
            var_cell = var(i)
            var_diff(1:k) = var( (/ mesh%cell(i)%cell(1:k) /) ) - var_cell
            var_slope(i)%x = sum( muscl(i)%weights(1 : k) * var_diff(1:k) )
            var_slope(i)%y = sum( muscl(i)%weights(1+k:2*k) * var_diff(1:k) )
         else
            var_slope(i)%x = 0._rp
            var_slope(i)%y = 0._rp
         end if
      end do
   END SUBROUTINE Least_Square_Slope_Super_Filter
   !> \brief Dimensional Cell Slope for cartesian mesh
   SUBROUTINE Dimensional_Slope( var_slope , var , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in) :: mesh
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
      type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(out) :: var_slope
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      real(rp), dimension( mesh%nc ) :: gam
      integer(ip) :: i_p , i_m , ie_p , ie_m
      real(rp) :: var_p , var_m , var_c , r
      !================================================================================================================!
      ! Begin Subroutine
      !================================================================================================================!
  ! if ( mesh_type /= 'basic' ) call Stopping_Program_Sub( 'Can t use dimensional splitting MUSCL scheme' )
      do i = 1,mesh%nc
         i_m = mesh%cell(i)%cell(4)
         i_p = mesh%cell(i)%cell(2)
         ie_m = mesh%cell(i)%edge(4)
         ie_p = mesh%cell(i)%edge(2)
         var_m = var( i_m )
         var_p = var( i_p )
         var_c = var( i )
         if ( var_c == var_p .or. &
              var_c == var_m ) then
            var_slope(i)%x = 0._rp
         else
            r = ( var_c - var_m ) / ( var_p - var_c )
            var_slope(i)%x = ( var_p - var_c ) * Van_Leer_Limiter( r ) / dx
         end if
         i_m = mesh%cell(i)%cell(1)
         i_p = mesh%cell(i)%cell(3)
         ie_m = mesh%cell(i)%edge(1)
         ie_p = mesh%cell(i)%edge(3)
         var_m = var( i_m )
         var_p = var( i_p )
         var_c = var( i )
         if ( var_c == var_p .or. &
              var_c == var_m ) then
            var_slope(i)%y = 0._rp
         else
            r = ( var_c - var_m ) / ( var_p - var_c )
            var_slope(i)%y = ( var_p - var_c ) * Van_Leer_Limiter( r ) / dy
         end if
      end do
   END SUBROUTINE Dimensional_Slope
   !> \brief Barth Slope Limiter
   !! \details ...
   SUBROUTINE Slope_Limiter_Barth( var_slope , var , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in) :: mesh
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
      type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(inout) :: var_slope
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      real(rp), dimension(maxed) :: t1 , t2 , limiter
      real(rp) :: cell_limiter
      type( vec2d ) :: vec2d_loc
      !================================================================================================================!
      ! Begin Subroutine
      !================================================================================================================!
      do i = 1,mesh%nc
         vec2d_loc%x = var_slope(i)%x
         vec2d_loc%y = var_slope(i)%y
         if ( .norm. vec2d_loc < zerom ) cycle
         limiter(:) = one
         do j = 1,mesh%cell(i)%nbed
            ie = mesh%cell(i)%edge(j)
            t1(j) = var( mesh%cell(i)%cell(j) ) - var(i)
            if ( mesh%edge(ie)%cell(1) == i ) then
               t2(j) = vec2d_loc .dotprod. mesh%edge( mesh%cell(i)%edge(j) )%v_edge_cell(1)
            else
               t2(j) = vec2d_loc .dotprod. mesh%edge( mesh%cell(i)%edge(j) )%v_edge_cell(2)
            end if
            if ( t1(j) * t2(j) < zerom ) then
               limiter(j) = zero
            else
               limiter(j) = t1(j) / t2(j)
            end if
         end do
         cell_limiter = min( one , minval( limiter(:) ) )
         var_slope(i)%x = cell_limiter * vec2d_loc%x
         var_slope(i)%y = cell_limiter * vec2d_loc%y
      end do
   END SUBROUTINE Slope_Limiter_Barth
   !> \brief MP_property
   !! \details ...
   real(rp) FUNCTION MP_property( var , left , right )
      USE m_common
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      real(rp), intent(in) :: var , left , right
      !================================================================================================================!
      ! Begin Function
      !================================================================================================================!
      MP_property = min( max( left , right ) , max( min( left , right ) , var ) )
   END FUNCTION MP_property
   !> \brief TVD Property (left)
   !! \details
   real(rp) FUNCTION TVD_left_property( var , left , right )
      USE m_common
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      real(rp), intent(in) :: var , left , right
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      real(rp) :: middle
      !================================================================================================================!
      ! Begin Function
      !================================================================================================================!
      middle = 0.5_rp * ( left + right )
      TVD_left_property = min( max( left , middle ) , max( min( left , middle ) , var ) )
   END FUNCTION TVD_left_property
   !> \brief TVD Property (right)
   !! \details
   real(rp) FUNCTION TVD_right_property( var , left , right )
      USE m_common
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      real(rp), intent(in) :: var , left , right
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      real(rp) :: middle
      !================================================================================================================!
      ! Begin Function
      !================================================================================================================!
      middle = 0.5_rp * ( left + right )
      TVD_right_property = min( max( middle , right ) , max( min( middle , right ) , var ) )
   END FUNCTION TVD_right_property
   !> \brief Least Square Cell Slope with Filter (old)
   !! \details ...
   SUBROUTINE Old_Least_Square_Slope_Filter( var_slope , var , filter , cutoff , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in) :: mesh
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var , filter
      type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(out) :: var_slope
      real(rp), intent(in) :: cutoff
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      real(rp), dimension(2,3) :: coef
      real(rp) :: det
      !================================================================================================================!
      ! Begin Subroutine
      !================================================================================================================!
      do i = 1,mesh%nc
         coef(:,:) = 0._rp
         if ( filter(i) > cutoff ) then
            do j = 1,mesh%cell(i)%nbed
               ie = mesh%cell(i)%edge(j)
               if ( mesh%edge(ie)%cell(1) == i ) then
                  coef(1,1) = coef(1,1) + mesh%edge(ie)%vcell%x * mesh%edge(ie)%vcell%x
                  coef(1,2) = coef(1,2) + mesh%edge(ie)%vcell%x * mesh%edge(ie)%vcell%y
                  coef(1,3) = coef(1,3) + mesh%edge(ie)%vcell%x * ( var( mesh%cell(i)%cell(j) ) - var(i) )
                  coef(2,1) = coef(2,1) + mesh%edge(ie)%vcell%y * mesh%edge(ie)%vcell%x
                  coef(2,2) = coef(2,2) + mesh%edge(ie)%vcell%y * mesh%edge(ie)%vcell%y
                  coef(2,3) = coef(2,3) + mesh%edge(ie)%vcell%y * ( var( mesh%cell(i)%cell(j) ) - var(i) )
               else
                  coef(1,1) = coef(1,1) + mesh%edge(ie)%vcell%x * mesh%edge(ie)%vcell%x
                  coef(1,2) = coef(1,2) + mesh%edge(ie)%vcell%x * mesh%edge(ie)%vcell%y
                  coef(1,3) = coef(1,3) - mesh%edge(ie)%vcell%x * ( var( mesh%cell(i)%cell(j) ) - var(i) )
                  coef(2,1) = coef(2,1) + mesh%edge(ie)%vcell%y * mesh%edge(ie)%vcell%x
                  coef(2,2) = coef(2,2) + mesh%edge(ie)%vcell%y * mesh%edge(ie)%vcell%y
                  coef(2,3) = coef(2,3) - mesh%edge(ie)%vcell%y * ( var( mesh%cell(i)%cell(j) ) - var(i) )
               end if
            end do
            det = coef(1,1) * coef(2,2) - coef(1,2) * coef(2,1)
            if ( abs( det ) < zerom ) then
               var_slope(i)%x = 0._rp
               var_slope(i)%y = 0._rp
            else
               var_slope(i)%x = ( coef(1,3) * coef(2,2) - coef(1,2) * coef(2,3) ) / det
               var_slope(i)%y = ( coef(1,1) * coef(2,3) - coef(2,1) * coef(1,3) ) / det
            end if
         else
            var_slope(i)%x = 0._rp
            var_slope(i)%y = 0._rp
         end if
      end do
   END SUBROUTINE Old_Least_Square_Slope_Filter
   !> \brief Calculate the Cell Gradient of a variable using Green Formula and basic interpolation at edge
   !! \details i don't know the difference with first one
   SUBROUTINE FV_Cell_Grad2( grad_var , var , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in) :: mesh
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
      type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(out) :: grad_var
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      integer(ip) :: iL , iR ! Left and Right cells indexes to edge
      !================================================================================================================!
      ! Begin Subroutine
      !================================================================================================================!
      grad_var(:)%x = 0._rp
      grad_var(:)%y = 0._rp
      do ie = 1,mesh%ne
         iL = mesh%edge(ie)%cell(1)
         iR = mesh%edge(ie)%cell(2)
         grad_var( iL )%x = grad_var( iL )%x + mesh%edge(ie)%length * mesh%edge(ie)%normal%x * &
                                                   d1p4 * ( var( iL ) + var( iR ) ) ** 2
         grad_var( iL )%y = grad_var( iL )%y + mesh%edge(ie)%length * mesh%edge(ie)%normal%y * &
                                                   d1p4 * ( var( iL ) + var( iR ) ) ** 2
         if ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain ) then
            grad_var( iR )%x = grad_var( iR )%x - mesh%edge(ie)%length * mesh%edge(ie)%normal%x * &
                                                      d1p4 * ( var( iL ) + var( iR ) ) ** 2
            grad_var( iR )%y = grad_var( iR )%y - mesh%edge(ie)%length * mesh%edge(ie)%normal%y * &
                                                      d1p4 * ( var( iL ) + var( iR ) ) ** 2
         end if
      end do
      do i = 1,mesh%nc
         grad_var(i)%x = grad_var(i)%x * mesh%cell(i)%invsurf
         grad_var(i)%y = grad_var(i)%y * mesh%cell(i)%invsurf
      end do
   END SUBROUTINE FV_Cell_Grad2
   !> \brief Van_Leer_Limiter
   !! \details ...
   real(rp) FUNCTION Van_Leer_Limiter( r )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      real(rp), intent(in) :: r
      !================================================================================================================!
      ! Begin Function
      !================================================================================================================!
      Van_Leer_Limiter = max( 0._rp , min( 1._rp , 2._rp * r ) , min( r , 2._rp ) )
   END FUNCTION Van_Leer_Limiter
   !> \brief Solving ODE using RK4 method
   !! \details STEP
   SUBROUTINE rk4_step( f , t , dtr , y )
      USE m_common
      implicit none
      real(rp), intent(in ) :: t , dtr
      real(rp), intent(inout) :: y
      real(rp) :: f ; external f
      integer(ip) :: i_rk
      real(rp) :: k_rk(4)
      real(rp), dimension(4) :: c = (/ 0.0_8 , &
                                         0.5_8 , &
                                         0.5_8 , &
                                         1.0_8 /)
      real(rp), dimension(4) :: b = (/ 1._8 / 6._8 , &
                                         1._8 / 3._8 , &
                                         1._8 / 3._8 , &
                                         1._8 / 6._8 /)
      real(rp), dimension(4,4) :: a = reshape( (/ 0.0_8 , 0.5_8 , 0.0_8 , 0.0_8 , &
                                                    0.0_8 , 0.0_8 , 0.5_8 , 0.0_8 , &
                                                    0.0_8 , 0.0_8 , 0.0_8 , 1.0_8 , &
                                                    0.0_8 , 0.0_8 , 0.0_8 , 0.0_8 /) , (/ 4,4 /) )
      k_rk(:) = 0._rp
      do i_rk = 1,4
         k_rk(i_rk) = f( t + dtr * c(i_rk) , y + dtr * sum( a(i_rk,:) * k_rk(:) ) )
      end do
      y = y + dtr * sum( b(:) * k_rk(:) )
      return
   END SUBROUTINE rk4_step
   !> \brief Solving ODE using RK4 method
   !! \details POINT
   SUBROUTINE ode_rk4_point(t_start , t_end , y0 , y1, dtr , f)
      USE m_common
      implicit none
      real(rp), intent(in) :: t_start , t_end , y0
      real(rp), intent(out) :: y1
      real(rp), intent(inout) :: dtr
      real(rp) :: f ; external f
      real(rp) :: t
      t = t_start
      y1 = y0
      dtr = sign(1._rp , t_end - t_start) * abs(dtr)
      do while(t /= t_end)
         call rk4_step(f , t , dtr , y1 )
         t = t + dtr
         if ( sign( 1._rp , t_end - t_start ) * t > t_end ) t = t_end
      end do
   END SUBROUTINE ode_rk4_point
   !> \brief Solving ODE using RK4 method
   !! \details TOTAL
   SUBROUTINE ode_rk4_total( t_start , t_end , y0 , y , dtr , f )
      USE m_common
      implicit none
      real(rp), intent(in ) :: t_start , t_end , y0
      real(rp), intent(inout) :: dtr
      real(rp), dimension(:,:), allocatable, intent(out) :: y
      real(rp) :: f ; external f
      integer(ip) :: ir , n
      n = int( ( t_end - t_start ) / abs( dtr ) )
      allocate( y( 2 , 0 : n ) )
      dtr = sign( 1._rp , t_end - t_start ) * abs( dtr )
      ir = 0
      y(1,ir) = t_start
      y(2,ir) = y0
      do while( y(1,ir) /= t_end .and. ir <= n )
         ir = ir + 1
         y(2,ir) = y(2,ir-1)
         call rk4_step( f , y(1,ir-1) , dtr , y(2,ir) )
         y(1,ir) = y(1,ir-1) + dtr
      end do
   END SUBROUTINE ode_rk4_total
   !> \brief Numerical integration using Simpson rule
   !! \details USED for macdonnald's benchmark with perturbated topography
   real(rp) FUNCTION Simpson( f , xmin , xmax , nbp )
      USE m_common
      implicit none
      real(rp) :: f ; external f
      real(rp), intent(in) :: xmin , xmax
      integer(ip), intent(in) :: nbp
      real(rp) :: step , sum_s , x_stage
      integer(ip) :: i_stage , nb_stage
      if ( modulo( nbp , 2 ) == 0 ) then
         nb_stage = nbp
      else
         nb_stage = nbp + 1
      end if
      step = ( xmax - xmin ) / real( nb_stage , 8 )
      sum_s = f( xmin )
      do i_stage = 1,nb_stage,2
         x_stage = xmin + real( i_stage , 8 ) * step
         sum_s = sum_s + 4._rp * f( x_stage )
      end do
      do i_stage = 2,nb_stage-1,2
         x_stage = xmin + real( i_stage , 8 ) * step
         sum_s = sum_s + 2._rp * f( x_stage )
      end do
      sum_s = sum_s + f( xmax )
      Simpson = step * sum_s / 3._rp
   END FUNCTION Simpson
   !> \brief Numerical integration over a cell using Gauss rule
   !! \details USED for ???
   real(rp) FUNCTION Gauss_rule_cell( var , it , mesh )
      USE m_common
         implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in) :: mesh
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
      integer(ip), intent(in) :: it
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      integer(ip), dimension(maxed) :: i_cell_nodes
      !================================================================================================================!
      ! Begin Function
      !================================================================================================================!
      Gauss_rule_cell = 0._rp
   END FUNCTION Gauss_rule_cell !>NOADJ
END MODULE m_numeric
