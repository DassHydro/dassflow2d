SUBROUTINE euler_time_step_first_b1( dof , mesh )
   USE m_common
   USE m_mesh
   USE m_time_screen !NOADJ
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(inout) :: mesh
   type( unk ), intent(inout) :: dof
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   integer(ip) :: iL , iR ! Left and Right cells indexes to edge
   real(rp) :: hL(2) , uL(2) , vL(2) , zL ! Left State in edge cell normal coordinates
   real(rp) :: hR(2) , uR(2) , vR(2) , zR ! Right State in edge cell normal coordinates
   real(rp), dimension( sw_nb ) :: nflux ! Finite Volume normal edge flux
   real(rp), dimension( sw_nb ) :: lflux ! Finite Volume edge flux in (x,y) coordinates
   real(rp), dimension( sw_nb , mesh%nc ) :: tflux ! Finite Volume total flux for each cell
   real(rp) :: h , u , v ! Temporal primitive variables
   !Infiltration variables
   real(rp) :: S ! potential maximal retention
   real(rp) :: Fn1 ! Temporal Fn+1
   real(rp) :: aFn1 , bFn1
   real(rp) :: h_infil ! local variable of infil calculated depth 
   real(rp) :: vel ! Velocity norm
   real(rp) :: sfl ! Manning
   real(rp) :: madd ! mass rain >TGADJ
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
 !write(*,*) "/\/\/\ WITHIN euler_time_step_first_b1"
   tflux(:,:) = 0._rp
   do ie = 1,mesh%ne
 !write(*,*) "--- edge=",ie
      !================================================================================================================!
      ! Calculate Left and Right States
      !================================================================================================================!
      iL = mesh%edge(ie)%cell(1)
      iR = mesh%edge(ie)%cell(2) !Left cell id for a normal cell
    if ( mesh%edge(ie)%boundary ) then !Check if bounfary first so typlim exists
        if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_1D' ) cycle
        if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' ) then !then change connectivity to connected 1D-like cell
            iR = mesh%edge(ie)%cell1D2D !Get id of the single 1D-like cell with interface in the connected bc number => this should be done once!
        endif
    endif
      hL(1) = dof%h( iL )
      hR(1) = dof%h( iR )
       if ( hL(1) > heps .or. hR(1) > heps ) then
         zL = bathy_cell( iL )! + global_bathy_shift(1)
         zR = bathy_cell( iR )! + global_bathy_shift(1)
         uL(1) = dof%u( iL )
         vL(1) = dof%v( iL )
         uL(2) = mesh%edge(ie)%normal%x * uL(1) + mesh%edge(ie)%normal%y * vL(1)
         vL(2) = mesh%edge(ie)%normal%x * vL(1) - mesh%edge(ie)%normal%y * uL(1)
         if ( mesh%edge(ie)%boundary) then
             if (.not. ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' )) then !do not call boundary calculations for internal BCs
                call calc_boundary_state( mesh , hL(1) , zL , uL(2) , vL(2) , &
                                                 hR(1) , zR , uR(2) , vR(2) )
             else
                uR(1) = dof%u( iR )
                vR(1) = dof%v( iR )
                uR(2) = mesh%edge(ie)%normal%x * uR(1) + mesh%edge(ie)%normal%y * vR(1)
                vR(2) = mesh%edge(ie)%normal%x * vR(1) - mesh%edge(ie)%normal%y * uR(1)
            endif
         else
            uR(1) = dof%u( iR )
            vR(1) = dof%v( iR )
            uR(2) = mesh%edge(ie)%normal%x * uR(1) + mesh%edge(ie)%normal%y * vR(1)
            vR(2) = mesh%edge(ie)%normal%x * vR(1) - mesh%edge(ie)%normal%y * uR(1)
         end if
         !=============================================================================================================!
         ! New reconstructed well balanced water depth
         !=============================================================================================================!
         hL(2) = max( 0._rp , hL(1) + zL - max( zL , zR ) )
         hR(2) = max( 0._rp , hR(1) + zR - max( zL , zR ) )
         !=============================================================================================================!
         ! Calling the balanced HLLC Solver dedicated to Shallow-Water Equations
         !=============================================================================================================!
         call sw_hllc( hL(2) , uL(2) , vL(2) , &
                       hR(2) , uR(2) , vR(2) , nflux )
         !=============================================================================================================!
         ! Boundary post treatment :
         ! - Feedback control of bathy_cell in ghost cells to properly control the Qin imposed
         ! - Calculation of nflux sum for each inflow
         !=============================================================================================================!
         if ( mesh%edge(ie)%boundary ) call boundary_post( nflux(1) , iR , mesh )
         !=============================================================================================================!
         ! Flux rotation and summation (as antisymmetric part to save time computation)
         !=============================================================================================================!
         lflux(1) = nflux(1)
         lflux(2) = mesh%edge(ie)%normal%x * nflux(2) - mesh%edge(ie)%normal%y * nflux(3)
         lflux(3) = mesh%edge(ie)%normal%y * nflux(2) + mesh%edge(ie)%normal%x * nflux(3)
         lflux(1:3) = lflux(1:3) * mesh%edge(ie)%length
         tflux( 1 , iL ) = tflux( 1 , iL ) + lflux(1)
         tflux( 2 , iL ) = tflux( 2 , iL ) + lflux(2)
         tflux( 3 , iL ) = tflux( 3 , iL ) + lflux(3)
         tflux( 2 , iL ) = tflux( 2 , iL ) + mesh%edge(ie)%normal%x * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                 ( hL(1)**2 - hL(2)**2 ) )
         tflux( 3 , iL ) = tflux( 3 , iL ) + mesh%edge(ie)%normal%y * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                 ( hL(1)**2 - hL(2)**2 ) )
         if ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain ) then
            tflux( 1 , iR ) = tflux( 1 , iR ) - lflux(1)
            tflux( 2 , iR ) = tflux( 2 , iR ) - lflux(2)
            tflux( 3 , iR ) = tflux( 3 , iR ) - lflux(3)
            tflux( 2 , iR ) = tflux( 2 , iR ) - mesh%edge(ie)%normal%x * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                    ( hR(1)**2 - hR(2)**2 ) )
            tflux( 3 , iR ) = tflux( 3 , iR ) - mesh%edge(ie)%normal%y * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                    ( hR(1)**2 - hR(2)**2 ) )
         end if
         if ( mesh%edge(ie)%boundary ) then
            if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' ) then
                tflux( 1 , iR ) = tflux( 1 , iR ) - lflux(1)
                tflux( 2 , iR ) = tflux( 2 , iR ) - lflux(2)
                tflux( 3 , iR ) = tflux( 3 , iR ) - lflux(3)
                tflux( 2 , iR ) = tflux( 2 , iR ) - mesh%edge(ie)%normal%x * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                        ( hR(1)**2 - hR(2)**2 ) )
                tflux( 3 , iR ) = tflux( 3 , iR ) - mesh%edge(ie)%normal%y * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                        ( hR(1)**2 - hR(2)**2 ) )
            endif
        endif
      end if
   end do
   !===================================================================================================================!
   ! Cumulative rain Calculation
   !===================================================================================================================!
   do k=1,bc%nb_rn
      bc%rain(k)%cumul = bc%rain(k)%cumul + dt*bc%rain(k)%qin
   end do
   !===================================================================================================================!
   ! Euler Time Step
   !===================================================================================================================!
   do i = 1,mesh%nc
      h = dof%h(i)
      u = dof%u(i)
      v = dof%v(i)
      dof%h(i) = max( 0._rp , h - dt * tflux(1,i) * mesh%cell(i)%invsurf )
    if (bc_rain == 1) then
      k = mesh%cell(i)%rain !Get rain group for current cell
      if (k > 0) then
         !SCS coefficient
         if (( bc_infil == 2) .and. (infil%land(i) /= 0)) then
            S = 25.4_rp * ( 1000._rp / abs(infil%SCS( infil%land(i) )%CN) - 10._rp ) / 1000._rp ! infil%land(i) is infiltration group for current cell
            if ( bc%rain(k)%cumul > abs(infil%SCS( infil%land(i) )%lambda) * S ) then
               Fn1 = abs(infil%SCS( infil%land(i) )%lambda) * S + &
                     S*( bc%rain(k)%cumul - abs(infil%SCS( infil%land(i) )%lambda) * S ) &
                     / ( bc%rain(k)%cumul + (1 - abs(infil%SCS( infil%land(i) )%lambda) ) * S )
            else
               Fn1 = dof%infil(i) + dt*bc%rain(k)%qin
            endif
            dof%h(i) = dof%h(i) + dt*bc%rain(k)%qin - Fn1 + dof%infil(i) !SCS-modified rain
            dof%infil(i) = Fn1
         else !Unmodified rain
            dof%h(i) = dof%h(i) + dt*bc%rain(k)%qin
         endif
       endif
    endif
     if ((bc_infil == 1) ) then
       if (infil%land(i) /= 0) then
         aFn1 = dof%infil(i) + dt * infil%GA( infil%land(i) )%Ks * ( 1._rp - infil%GA( infil%land(i) )%DeltaTheta ) !CHECK FORMULATION !CHECK FORMULATION
         bFn1 = infil%GA( infil%land(i) )%Ks*dt*infil%GA( infil%land(i) )%DeltaTheta*( dof%infil(i) + dof%h(i) &
               + infil%GA( infil%land(i) )%PsiF ) ! infil%land(i) is infiltration group for current cell
         Fn1 = ( aFn1 + sqrt( aFn1**2._rp + 4._rp*bFn1 ) ) / 2._rp !CHECK FORMULATION
         !-------original code-----
         !dof%h(i) = dof%h(i) + dof%infil(i) - Fn1
   h_infil = dof%h(i) + dof%infil(i) - Fn1
   !-------original code-------
         !if (dof%h(i) < 0._rp ) then
   ! dof%h(i) = 0._rp
         ! Fn1 = dof%h(i) + dof%infil(i)
         !endif
   if (h_infil < 0._rp ) then
             Fn1 = dof%h(i) + dof%infil(i)
    h_infil = 0.000001_rp
         endif
   dof%h(i) = h_infil !replace the local variable h_infil
         dof%infil(i) = Fn1
       endif
      endif
      !================================================================================================================!
      ! Positivity cut-off
      !================================================================================================================!
      if ( dof%h(i) <= heps ) then
         dof%u(i) = 0._rp
         dof%v(i) = 0._rp
      else
         dof%u(i) = ( h * u - dt * ( tflux(2,i) * mesh%cell(i)%invsurf ) ) / dof%h(i)
         dof%v(i) = ( h * v - dt * ( tflux(3,i) * mesh%cell(i)%invsurf ) ) / dof%h(i)
         !=============================================================================================================!
         ! Semi-Implicit Treatment of Friction Source Term (Manning/Strickler Formula)
         !=============================================================================================================!
         if ( friction == 1 ) then
            vel = sqrt( dof%u(i)**2 + dof%v(i)**2 )
            sfl = dof%h(i)**d2p3 + sqrt( dof%h(i)**d4p3 + 4._rp * dt * g * &
                   (manning(land(i))*dof%h(i)**manning_beta(land(i)))**2 * vel )
            sfl = 2._rp * dof%h(i)**d2p3 / sfl
         else if ( friction == 2 ) then
            sfl = one - dt * manning( land(i) )
         else
            sfl = 1._rp
         end if
         dof%u(i) = dof%u(i) * sfl
         dof%v(i) = dof%v(i) * sfl
      end if
   end do
END SUBROUTINE euler_time_step_first_b1
