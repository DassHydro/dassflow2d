MODULE m_obs
   USE m_common
   USE m_mpi
      USE m_model
   implicit none
   TYPE innovation_obs
      integer(ip) :: nb_dt , nb_dx , ind_t
      real(rp), dimension(:), allocatable :: diff
   END TYPE
   type( innovation_obs ), dimension(:), allocatable :: innovation, innovW, innovQ, innovUV
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
      integer(ip) :: idiff, temp, sizeQ
      type( vec2d ), dimension( mesh%nc + mesh%ncb ) :: grad_var
      real(rp) , dimension( mesh%nc + mesh%ncb ) :: bathy_temp
      !================================================================================================================!
      ! Loop on observations and quadratic norm of innovation vector / Regularization terms
      !================================================================================================================!
            if ( use_obs == 1 ) then
             if (( use_Zobs == 0 ) .and. ( use_Qobs == 0 ).and. &
                 ( use_UVobs == 0 ) .and. ( use_Qobs_gr4 == 0 )) then
                call Stopping_Program_Sub( 'You need to pick at least one observation type to be used with use_obs = 1. Please set any of use_Zobs, use_Qobs, use_Qobs_gr4 or use_UVobs to 1 and provide the appropriate observation files.')
            endif
            !==========================================================================================================!
            ! Initialisation to zero of each part
            !==========================================================================================================!
            cost_part(:) = 0._rp
            !==========================================================================================================!
            ! Loop on observations in Time/Space
            !==========================================================================================================!
            if ( use_Zobs == 1 ) then
                do iobs = 1,size( station )
                  do idiff = 1,size( innovation( iobs )%diff )
                    cost_part(1) = cost_part(1) + station( iobs )%weight * innovation ( iobs )%diff( idiff )**2
                  enddo
                enddo
            endif
write(*,*) "cost from use_Zobs", cost_part(1)!NOADJ
            if ( use_UVobs == 1 ) then
                do iobs = 1,size( station )
                  do idiff = 1,size( innovUV( iobs )%diff )
                    cost_part(1) = cost_part(1) + station( iobs )%weight * innovUV ( iobs )%diff( idiff )**2
                  enddo
                enddo
            endif
            !!!!!!!!!!!!!
            ! InnovW
            !!!!!!!!!!!!!
            !!!!!!!!!!!!!
            ! InnovW
            !!!!!!!!!!!!!
            !do iobs = 1,size( station )
            ! do idiff = 1,size( innovW( iobs )%diff )
            ! cost_part(1) = cost_part(1) + station( iobs )%weight * innovW ( iobs )%diff( idiff )**2
            ! end do
            !end do
            if (( use_Qobs == 1 ) .or. ( use_Qobs_gr4 == 1 )) then
                if ( use_NSE == 0 ) then
                 ! InnovQ with RMSE
                  do iobs = 1,size( stationQ )
                    do idiff = 1,size( innovQ( iobs )%diff )
                    !write(*,*) proc, idiff, stationQ( iobs )%weight, innovQ ( iobs )%diff( idiff )**2!NOADJ
                        cost_part(1) = cost_part(1) + stationQ( iobs )%weight * innovQ ( iobs )%diff( idiff )**2
                    end do
                  end do
                else
                  ! InnovQ with NSE
                  sizeQ = size( stationQ )
                  do iobs = 1,size( stationQ )
                    cost_part(3) = 0._rp
                    cost_part(2) = 0._rp
                      do idiff = 1,size( innovQ( iobs )%diff )
                        cost_part(2) = cost_part(2) + (innovQ ( iobs )%diff( idiff ))**2
                        cost_part(3) = cost_part(3) + (innovQ ( iobs + sizeQ )%diff( idiff ))**2
                      end do
                    cost_part(1) = cost_part(1) + cost_part(2) / cost_part(3)
                  enddo
                endif
            endif
            !~ ! Regularize with bathymetry gradient (either full gradient, or only the free distributed bathymetry
            if (regul_bathy_grad == 1) then
                call FV_Cell_Grad( grad_var , bathy_cell , mesh )
            elseif ((regul_bathy_grad .eq. 2) .and. (use_xsshp .eq. 1)) then
                bathy_temp(:) = 0._rp
                if (xsshp_along_y == 1) then
                    do ie = 1, mesh%nc
                    if (mesh%cell(ie)%grav%x <= XSshape(1)%xcenter) then
                        bathy_temp( ie ) = bathy_cell(ie) - ( &
                        XSshape(1)%topz + min( 0._rp, - abs(XSshape(1)%hmax) *&
                        (1 - ( abs((mesh%cell(ie)%grav%x - XSshape(1)%xleft) -&
                        (XSshape(1)%xcenter-XSshape(1)%xleft))/&
                        ((XSshape(1)%xcenter-XSshape(1)%xleft)) )**XSshape(1)%s)) &
                        + slope_y(1) * mesh%cell(ie)%grav%y &
                        + slope_x(1) * mesh%cell(ie)%grav%x )
                    else
                        bathy_temp( ie ) = bathy_cell(ie) - ( &
                        XSshape(1)%topz + &
                        min( 0._rp, - abs(XSshape(1)%hmax) *&
                        (1 - ( abs((mesh%cell(ie)%grav%x - XSshape(1)%xleft) -&
                        (XSshape(1)%xcenter-XSshape(1)%xleft))/&
                        ((XSshape(1)%xright-XSshape(1)%xcenter)) )**XSshape(1)%s)) &
                        + slope_y(1) * mesh%cell(ie)%grav%y &
                        + slope_x(1) * mesh%cell(ie)%grav%x )
                    endif
                    enddo
                elseif (xsshp_along_x == 1) then
                    do ie = 1, mesh%nc
                    if (mesh%cell(ie)%grav%y <= XSshape(1)%xcenter) then
                        bathy_temp( ie ) = bathy_cell(ie) - ( &
                        XSshape(1)%topz + min( 0._rp, - abs(XSshape(1)%hmax) *&
                        (1 - ( abs((mesh%cell(ie)%grav%y - XSshape(1)%xleft) -&
                        (XSshape(1)%xcenter-XSshape(1)%xleft))/&
                        ((XSshape(1)%xcenter-XSshape(1)%xleft)) )**XSshape(1)%s)) &
                        + slope_y(1) * mesh%cell(ie)%grav%y &
                        + slope_x(1) * mesh%cell(ie)%grav%x )
                    else
                        bathy_temp( ie ) = bathy_cell(ie) - ( &
                        XSshape(1)%topz + &
                        min( 0._rp, - abs(XSshape(1)%hmax) *&
                        (1 - ( abs((mesh%cell(ie)%grav%y - XSshape(1)%xleft) -&
                        (XSshape(1)%xcenter-XSshape(1)%xleft))/&
                        ((XSshape(1)%xright-XSshape(1)%xcenter)) )**XSshape(1)%s)) &
                        + slope_y(1) * mesh%cell(ie)%grav%y &
                        + slope_x(1) * mesh%cell(ie)%grav%x )
                    endif
                    enddo
                else
                    write(*,*) "ERROR: You must define a direction (x or y-axis) for the cross-section parameterization, with xsshp_along_x or xsshp_along_y." !NOADJ
                endif
                call FV_Cell_Grad( grad_var , bathy_temp , mesh )
            elseif (regul_bathy_grad .eq. 2) then
                write(*,*) "ERROR: regul_bathy_grad is equal to 2 but use_xsshp is not equal to 1. You need to use cross-section parameterization to use this bathymetry gradient regularization." !NOADJ
            endif
            if (regul_bathy_grad .ne. 0) then
                do i = 1,mesh%nc
                    cost_part(2) = cost_part(2) + ( ( grad_var(i)%x )**2 + ( grad_var(i)%y )**2 )
                end do
                call mpi_sum_r( cost_part(2) )
            endif
        !~ ! Regularize with free bathymetry misfit to XS shape and slope
            if ((regul_bathy_shape .eq. 1) .and. (use_xsshp .eq. 1)) then
              if (xsshp_along_x == 1) then
                do ie= 1,mesh%nc
                  if (mesh%cell(ie)%grav%x <= XSshape(1)%xcenter) then
                    cost_part(2) = cost_part(2) + ( bathy_cell(ie) - ( &
                    XSshape(1)%topz + min( 0._rp, - abs(XSshape(1)%hmax) *&
                    (1 - ( abs((mesh%cell(ie)%grav%x - XSshape(1)%xleft) -&
                    (XSshape(1)%xcenter-XSshape(1)%xleft))/&
                    ((XSshape(1)%xcenter-XSshape(1)%xleft)) )**XSshape(1)%s)) &
                    + slope_y(1) * mesh%cell(ie)%grav%y &
                    + slope_x(1) * mesh%cell(ie)%grav%x ) ) ** 2
                  else
                    cost_part(2) = cost_part(2) + ( bathy_cell(ie) - ( &
                    XSshape(1)%topz + &
                    min( 0._rp, - abs(XSshape(1)%hmax) *&
                    (1 - ( abs((mesh%cell(ie)%grav%x - XSshape(1)%xleft) -&
                    (XSshape(1)%xcenter-XSshape(1)%xleft))/&
                    ((XSshape(1)%xright-XSshape(1)%xcenter)) )**XSshape(1)%s)) &
                    + slope_y(1) * mesh%cell(ie)%grav%y &
                    + slope_x(1) * mesh%cell(ie)%grav%x ) ) **2
                  endif
                enddo
                call mpi_sum_r( cost_part(2) )
              elseif (xsshp_along_y == 1) then
                do ie= 1,mesh%nc
                  if (mesh%cell(ie)%grav%y <= XSshape(1)%xcenter) then
                    cost_part(2) = cost_part(2) + ( bathy_cell(ie) - ( &
                    XSshape(1)%topz + min( 0._rp, - abs(XSshape(1)%hmax) *&
                    (1 - ( abs((mesh%cell(ie)%grav%y - XSshape(1)%xleft) -&
                    (XSshape(1)%xcenter-XSshape(1)%xleft))/&
                    ((XSshape(1)%xcenter-XSshape(1)%xleft)) )**XSshape(1)%s)) &
                    + slope_y(1) * mesh%cell(ie)%grav%y &
                    + slope_x(1) * mesh%cell(ie)%grav%x ) ) ** 2
                  else
                    cost_part(2) = cost_part(2) + ( bathy_cell(ie) - ( &
                    XSshape(1)%topz + &
                    min( 0._rp, - abs(XSshape(1)%hmax) *&
                    (1 - ( abs((mesh%cell(ie)%grav%y - XSshape(1)%xleft) -&
                    (XSshape(1)%xcenter-XSshape(1)%xleft))/&
                    ((XSshape(1)%xright-XSshape(1)%xcenter)) )**XSshape(1)%s)) &
                    + slope_y(1) * mesh%cell(ie)%grav%y &
                    + slope_x(1) * mesh%cell(ie)%grav%x ) ) **2
                  endif
                end do
              endif
            endif
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
endif
      !================================================================================================================!
      ! Fake operation for Tapenade Automatic Differentiation (Last operation ...)
      !================================================================================================================!
      cost = sqrt( cost**2 )
      write(*,*) "cost", cost
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
            call mpi_sum_r( h_mean )
            h_mean = h_mean / real( size( station( iobs )%pt ) , 8 )
            cost = cost + station( iobs )%weight * h_mean**2
         end do
   END SUBROUTINE update_cost_function
END MODULE m_obs
