MODULE m_obs
   USE m_common
   USE m_model
   implicit none
   TYPE innovation_obs
      integer(ip) :: nb_dt , nb_dx , ind_t
      real(rp), dimension(:), allocatable :: diff
   END TYPE
   type( innovation_obs ), dimension(:), allocatable :: innovation
   integer(ip) :: nb_obs , nb_grp , iobs
   real(rp) :: w_mean
CONTAINS
   SUBROUTINE calc_cost_function( cost , mesh )
      USE m_numeric
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in) :: mesh
      real(rp), intent(inout) :: cost
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      real(rp) :: cost_part(3) , filtered(4)
      integer(ip) :: idiff
      type( vec2d ), dimension( mesh%nc + mesh%ncb ) :: grad_var
      !================================================================================================================!
      ! Loop on observations and quadratic norm of innovation vector / Regularization terms
      !================================================================================================================!
         if ( use_obs == 1 ) then
            !==========================================================================================================!
            ! Initialisation to zero of each part
            !==========================================================================================================!
            cost_part(:) = 0._rp
            !==========================================================================================================!
            ! Loop on observations in Time/Space
            !==========================================================================================================!
            do iobs = 1,size( station )
               do idiff = 1,size( innovation( iobs )%diff )
                  cost_part(1) = cost_part(1) + station( iobs )%weight * innovation ( iobs )%diff( idiff )**2
               end do
            end do
            do i = 1,mesh%nc
               cost_part(2) = cost_part(2) + ( ( grad_var(i)%x )**2 + ( grad_var(i)%y )**2 )
            end do
            cost_part(2) = cost_part(2) * regul_bathy
             do k = 1,bc%nb_in
                filtered(1) = bc%hyd(k)%q(1)
                do i = 2,size( bc%hyd(k)%q(:) )-3
                   filtered(2) = filtered(1) + 0.2_rp * ( bc%hyd(k)%q(i ) - filtered(1) )
                   filtered(3) = filtered(2) + 0.2_rp * ( bc%hyd(k)%q(i+1) - filtered(2) )
                   filtered(4) = filtered(3) + 0.2_rp * ( bc%hyd(k)%q(i+2) - filtered(3) )
                   cost_part(3) = cost_part(3) + ( bc%hyd(k)%q(i) - filtered(4) )**2
                   filtered(1) = filtered(4)
                end do
             end do
             cost_part(3) = cost_part(3) * regul_hydrograph
           cost = sum( cost_part )
         end if
      !================================================================================================================!
      ! Fake operation for Tapenade Automatic Differentiation (Last operation ...)
      !================================================================================================================!
      cost = sqrt( cost**2 )
   END SUBROUTINE calc_cost_function
   SUBROUTINE update_cost_function( dof , cost )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( unk ), intent(in) :: dof
      real(rp), intent(inout) :: cost
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      integer(ip) :: cell , pt
      real(rp) :: h_mean
      !================================================================================================================!
      ! Begin
      !================================================================================================================!
        if (.not. allocated( station ) ) return
         do iobs = 1,size( station )
            if ( .not. test_dt_just_after( station( iobs )%dt ) ) cycle
            h_mean = 0._rp
            do pt = 1,size( station( iobs )%pt )
               cell = station( iobs )%pt( pt )%cell
               if ( cell < 0 ) cycle
               h_mean = h_mean + dof%h( cell )
            end do
            h_mean = h_mean / real( size( station( iobs )%pt ) , 8 )
            cost = cost + station( iobs )%weight * h_mean**2
         end do
   END SUBROUTINE update_cost_function
END MODULE m_obs
