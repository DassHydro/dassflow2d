SUBROUTINE muscl_aud_flux_n( dof , dofn , mesh )
   USE m_common
   USE m_mesh
   USE m_time_screen
   USE m_model
   USE m_numeric
   implicit none
   type( msh ), intent(in ) :: mesh
   type( unk ), intent(in ) :: dofn
   type( unk ), intent(inout) :: dof
   integer(ip) :: iL , iR, connected_num_bc ! Left and Right cells indexes to edge + connected interface index for 1D2D interfaces
   real(rp) :: hL(3) , uL(3) , vL(3) , zL(2) ! Left States in edge cell normal coordinates
   real(rp) :: hR(3) , uR(3) , vR(3) , zR(2) ! Right States in edge cell normal coordinates
   real(rp), dimension( sw_nb ) :: nflux ! Finite Volume normal edge flux
   real(rp), dimension( sw_nb ) :: lflux ! Finite Volume edge flux in (x,y) coordinates
   real(rp), dimension( sw_nb , mesh%nc ) :: tflux ! Finite Volume total flux for each cell
   tflux(:,:) = 0._rp
   call fill_bc( dof , mesh )
   call Time_Init_Part(2)
   call Least_Square_Slope_Filter ( dof%grad_h , dof%h , dof%h , heps , mesh )
   call Least_Square_Slope_Filter ( dof%grad_u , dof%u , dof%h , heps , mesh )
   call Least_Square_Slope_Filter ( dof%grad_v , dof%v , dof%h , heps , mesh )
   call Least_Square_Slope_Double_Filter( dof%grad_z , bathy_cell + dof%h , dof%h , heps , mesh )
   call Time_End_Part(2)
   call Time_Init_Part(3)
   if ( spatial_scheme == 'muscl_b1_b' ) call Slope_Limiter_Barth( dof%grad_h , dof%h , mesh )
   call Time_End_Part(3)
   call Time_Init_Part(4)
   do ie = 1,mesh%ne
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
         uL(1) = dof%u( iL )
         uR(1) = dof%u( iR )
         vL(1) = dof%v( iL )
         vR(1) = dof%v( iR )
         zL(1) = bathy_cell( iL )
         zR(1) = bathy_cell( iR )
         hL(2) = hL(1) + ( dof%grad_h( iL )%x * mesh%edge(ie)%v_edge_cell(1)%x + &
                                     dof%grad_h( iL )%y * mesh%edge(ie)%v_edge_cell(1)%y )
         hL(2) = MP_property( hL(2) , hL(1) , hR(1) )
         zL(2) = hL(1) + zL(1) + ( dof%grad_z( iL )%x * mesh%edge(ie)%v_edge_cell(1)%x + &
                                     dof%grad_z( iL )%y * mesh%edge(ie)%v_edge_cell(1)%y ) &
                                 - hL(2)
         zL(2) = MP_property( zL(2) , zL(1) , zR(1) )
         uL(2) = uL(1) + ( dof%grad_u( iL )%x * mesh%edge(ie)%v_edge_cell(1)%x + &
                             dof%grad_u( iL )%y * mesh%edge(ie)%v_edge_cell(1)%y )
         vL(2) = vL(1) + ( dof%grad_v( iL )%x * mesh%edge(ie)%v_edge_cell(1)%x + &
                             dof%grad_v( iL )%y * mesh%edge(ie)%v_edge_cell(1)%y )
         uL(2) = MP_property( uL(2) , uL(1) , uR(1) )
         vL(2) = MP_property( vL(2) , vL(1) , vR(1) )
         uL(3) = mesh%edge(ie)%normal%x * uL(2) + mesh%edge(ie)%normal%y * vL(2)
         vL(3) = mesh%edge(ie)%normal%x * vL(2) - mesh%edge(ie)%normal%y * uL(2)
         if ( mesh%edge(ie)%boundary ) then
             if (.not. ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' )) then !do not call boundary calculations for internal BCs
                    zR(2) = zL(2)
                    call calc_boundary_state( mesh , hL(2) , zL(2) , uL(3) , vL(3) , &
                                                    hR(2) , zR(2) , uR(3) , vR(3) )
             else
                hR(2) = hR(1) + ( dof%grad_h( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                                            dof%grad_h( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )
                hR(2) = MP_property( hR(2) , hL(1) , hR(1) )
                zR(2) = hR(1) + zR(1) + ( dof%grad_z( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                                            dof%grad_z( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y ) &
                                        - hR(2)
                zR(2) = MP_property( zR(2) , zL(1) , zR(1) )
                uR(2) = uR(1) + ( dof%grad_u( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                                    dof%grad_u( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )
                vR(2) = vR(1) + ( dof%grad_v( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                                    dof%grad_v( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )
                uR(2) = MP_property( uR(2) , uR(1) , uL(1) )
                vR(2) = MP_property( vR(2) , vR(1) , vL(1) )
                uR(3) = mesh%edge(ie)%normal%x * uR(2) + mesh%edge(ie)%normal%y * vR(2)
                vR(3) = mesh%edge(ie)%normal%x * vR(2) - mesh%edge(ie)%normal%y * uR(2)
             endif
         else
                hR(2) = hR(1) + ( dof%grad_h( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                                            dof%grad_h( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )
                hR(2) = MP_property( hR(2) , hL(1) , hR(1) )
                zR(2) = hR(1) + zR(1) + ( dof%grad_z( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                                            dof%grad_z( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y ) &
                                        - hR(2)
                zR(2) = MP_property( zR(2) , zL(1) , zR(1) )
                uR(2) = uR(1) + ( dof%grad_u( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                                    dof%grad_u( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )
                vR(2) = vR(1) + ( dof%grad_v( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                                    dof%grad_v( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )
                uR(2) = MP_property( uR(2) , uR(1) , uL(1) )
                vR(2) = MP_property( vR(2) , vR(1) , vL(1) )
                uR(3) = mesh%edge(ie)%normal%x * uR(2) + mesh%edge(ie)%normal%y * vR(2)
                vR(3) = mesh%edge(ie)%normal%x * vR(2) - mesh%edge(ie)%normal%y * uR(2)
         endif
         !=============================================================================================================!
         ! New reconstructed well balanced water depth
         !=============================================================================================================!
         hL(3) = max( 0._rp , hL(2) + zL(2) - max( zL(2) , zR(2) ) )
         hR(3) = max( 0._rp , hR(2) + zR(2) - max( zL(2) , zR(2) ) )
         !=============================================================================================================!
         ! Calling the balanced HLLC Solver dedicated to Shallow-Water Equations
         !=============================================================================================================!
         call sw_hllc( hL(3) , uL(3) , vL(3) , hR(3) , uR(3) , vR(3) , nflux )
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
                                                 ( hL(2)**2 - hL(3)**2 ) + ( hL(1) + hL(2) ) * ( zL(2) - zL(1) ) )
         tflux( 3 , iL ) = tflux( 3 , iL ) + mesh%edge(ie)%normal%y * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                 ( hL(2)**2 - hL(3)**2 ) + ( hL(1) + hL(2) ) * ( zL(2) - zL(1) ) )
         if ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain ) then
            tflux( 1 , iR ) = tflux( 1 , iR ) - lflux(1)
            tflux( 2 , iR ) = tflux( 2 , iR ) - lflux(2)
            tflux( 3 , iR ) = tflux( 3 , iR ) - lflux(3)
            tflux( 2 , iR ) = tflux( 2 , iR ) - mesh%edge(ie)%normal%x * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                    ( hR(2)**2 - hR(3)**2 ) + ( hR(1) + hR(2) ) * ( zR(2) - zR(1) ) )
            tflux( 3 , iR ) = tflux( 3 , iR ) - mesh%edge(ie)%normal%y * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                    ( hR(2)**2 - hR(3)**2 ) + ( hR(1) + hR(2) ) * ( zR(2) - zR(1) ) )
         end if
         if ( mesh%edge(ie)%boundary ) then
            if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' ) then !same as for (.not. mesh%edge(ie)%boundary) above for the particular case of internal BCs
                    tflux( 1 , iR ) = tflux( 1 , iR ) - lflux(1)
                    tflux( 2 , iR ) = tflux( 2 , iR ) - lflux(2)
                    tflux( 3 , iR ) = tflux( 3 , iR ) - lflux(3)
                    tflux( 2 , iR ) = tflux( 2 , iR ) - mesh%edge(ie)%normal%x * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                            ( hR(2)**2 - hR(3)**2 ) + ( hR(1) + hR(2) ) * ( zR(2) - zR(1) ) )
                    tflux( 3 , iR ) = tflux( 3 , iR ) - mesh%edge(ie)%normal%y * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                            ( hR(2)**2 - hR(3)**2 ) + ( hR(1) + hR(2) ) * ( zR(2) - zR(1) ) )
            endif
        endif
      end if
   end do
   call Time_End_Part(4)
   call Time_Init_Part(5)
   !===================================================================================================================!
   ! Euler Time Step
   !===================================================================================================================!
   do i = 1,mesh%nc
      dof%h(i) = dofn%h(i) - dt * tflux(1,i) * mesh%cell(i)%invsurf
      if ( dof%h(i) < 0._rp ) call spread_mass_added( dof , mesh )
      !================================================================================================================!
      ! Positivity cut-off
      !================================================================================================================!
      if ( dof%h(i) <= heps ) then
         dof%u(i) = 0._rp
         dof%v(i) = 0._rp
      else
         dof%u(i) = ( dofn%h(i) * dofn%u(i) - dt * tflux(2,i) * mesh%cell(i)%invsurf ) / dof%h(i)
         dof%v(i) = ( dofn%h(i) * dofn%v(i) - dt * tflux(3,i) * mesh%cell(i)%invsurf ) / dof%h(i)
      end if
   end do
   call Time_End_Part(5)
   !===================================================================================================================!
   ! filling ghost cells (+mpi if not disabled)
   !===================================================================================================================!
END SUBROUTINE muscl_aud_flux_n
