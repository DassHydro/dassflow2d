!======================================================================================================================!
!
!                    DassFlow Version 3.0
!
!======================================================================================================================!
!
!  Copyright University of Toulouse-INSA, Univ. of Strasbourg, INRAE & CNRS (France)
!
!  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
!  DassFlow is a computational software whose purpose is to simulate geophysical free surface flows,
!  designed for variational sensitivities and data assimilation (4D-var). Inverse capabilities are
!  based on the adjoint code generation by a source-to-source algorithmic differentiation (Tapenade software used).
!
!  DassFlow software includes few mostly independent "modules" with common architectures and structures.
!  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
!
!  Many people have contributed to the DassFlow development from the initial version to the latest ones.
!  Current contributions:
!               L. Pujol (PhD Unistra)
!               L. Villenave (PhD student)
!               P.-A. Garambois (INRAE Aix-en-Provence)
!               J. Monnier (INSA & Mathematics Institute of Toulouse IMT).
!               K. Larnier (CS group - IMT-INSA).
!  Former scientific or programming contributions of:
!               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
!               J.-P. Vila (INSA & Mathematics Institute of Toulouse IMT).
!               R. Madec   (Mathematics Institute of Toulouse IMT).
!  plus less recent other developers (M. Honnorat and J. Marin).
!
!  Contact : see the DassFlow webpage
!
!  This software is governed by the CeCILL license under French law and abiding by the rules of distribution
!  of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license
!  as circulated by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
!
!  As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the
!  license, users are provided only with a limited warranty and the software's author, the holder of the economic
!  rights, and the successive licensors have only limited liability.
!
!  In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or
!  developing or reproducing the software by the user in light of its specific status of free software, that may
!  mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
!  experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
!  software's suitability as regards their requirements in conditions enabling the security of their systems and/or
!  data to be ensured and, more generally, to use and operate it in the same conditions as regards security.
!
!  The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you
!  accept its terms.
!
!======================================================================================================================!
!> \file euler_time_step_first_b1.f90
!! \brief This file includes euler_time_step_first_b1 routine.
!! \details The file includes only euler_time_step_first_b1 routine (see doc euler_time_step_first_b1 routine).

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Perform Euler Time Step dedicated to Shallow-Water Equations
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief  Perform Euler Time Step dedicated to Shallow-Water Equations
!! \return dof updated after this new timestep
SUBROUTINE euler_time_step_first_b1( dof , mesh )

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_time_screen                                                                                              !NOADJ
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(inout)  ::  mesh

   type( unk ), intent(inout)  ::  dof

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  iL , iR   ! Left and Right cells indexes to edge

   real(rp)  ::  hL(2) , uL(2) , vL(2) , zL              ! Left  State in edge cell normal coordinates
   real(rp)  ::  hR(2) , uR(2) , vR(2) , zR              ! Right State in edge cell normal coordinates

   real(rp), dimension( sw_nb )  ::  nflux               ! Finite Volume normal edge flux
   real(rp), dimension( sw_nb )  ::  lflux               ! Finite Volume edge flux in (x,y) coordinates

   real(rp), dimension( sw_nb , mesh%nc )  ::  tflux     ! Finite Volume total flux for each cell

   real(rp)  ::  h , u , v                               ! Temporal primitive variables

   !Infiltration variables
   real(rp)  :: S                                        ! potential maximal retention
   real(rp)  :: Fn1                                      ! Temporal Fn+1
   real(rp)  :: aFn1 , bFn1
   real(rp)  :: h_infil                                  !local variable of infil calculated depth
   real(rp)  ::  vel                                     ! Velocity norm
   real(rp)  ::  sfl                                     ! Manning
   real(rp)  :: madd                                     ! mass rain >TGADJ
   
   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!


   tflux(:,:)  =  0._rp
   
   do ie = 1,mesh%ne

      !================================================================================================================!
      !  Calculate Left and Right States
      !================================================================================================================!

      iL  =  mesh%edge(ie)%cell(1)
      iR  =  mesh%edge(ie)%cell(2) !Left cell id for a normal cell

!     Get right cell id for internal BC cells
    if ( mesh%edge(ie)%boundary ) then !Check if bounfary first so typlim exists

        if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_1D' ) cycle

        if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' ) then !then change connectivity to connected 1D-like cell
            iR  =  mesh%edge(ie)%cell1D2D !Get id of the single 1D-like cell with interface in the connected bc number => this should be done once!
        endif

    endif

      hL(1)  =  dof%h( iL )
      hR(1)  =  dof%h( iR )

       if ( hL(1) > heps .or. hR(1) > heps ) then

         !zL  =  bathy_cell( iL )! + global_bathy_shift(1)
         !zR  =  bathy_cell( iR )! + global_bathy_shift(1)

         uL(1)  =  dof%u( iL )
         vL(1)  =  dof%v( iL )

         uL(2)  =  mesh%edge(ie)%normal%x * uL(1) + mesh%edge(ie)%normal%y * vL(1)
         vL(2)  =  mesh%edge(ie)%normal%x * vL(1) - mesh%edge(ie)%normal%y * uL(1)

         if ( mesh%edge(ie)%boundary) then

         zL  =  bathy_cell( iL )


          !================= TEMP FOR ANDROMEDE
            if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'zspresc') then
                zR  =  bathy_cell( iL ) !&
                        !- slope_y(1) * mesh%cell( mesh%edge(iL)%cell(1) )%surf / mesh%edge(iL)%length &
                       !- slope_x(1) * mesh%cell( mesh%edge(iL)%cell(1) )%surf / mesh%edge(iL)%length
            else if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'discharg1') then
                zR  =  bathy_cell( iL ) &
                 + slope_y(1) * mesh%cell( mesh%edge(iL)%cell(1) )%surf / mesh%edge(iL)%length &
                 + slope_x(1) * mesh%cell( mesh%edge(iL)%cell(1) )%surf / mesh%edge(iL)%length
            else if  ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'wall') then
                zR  =  bathy_cell( iL )
            endif
          !================= END TEMP
!          write(*,*) mesh%edgeb(mesh%edge(ie)%lim)%typlim, iL, iR, zR, zL !NOADJ

             if (.not. ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' )) then !do not call boundary calculations for internal BCs

                call calc_boundary_state( mesh , hL(1) , zL , uL(2) , vL(2) , &
                                                 hR(1) , zR , uR(2) , vR(2) )

             else

                uR(1)  =  dof%u( iR )
                vR(1)  =  dof%v( iR )

                uR(2)  =  mesh%edge(ie)%normal%x * uR(1) + mesh%edge(ie)%normal%y * vR(1)
                vR(2)  =  mesh%edge(ie)%normal%x * vR(1) - mesh%edge(ie)%normal%y * uR(1)

            endif

         else

            zL  =  bathy_cell( iL )
            zR  =  bathy_cell( iR )

            uR(1)  =  dof%u( iR )
            vR(1)  =  dof%v( iR )

            uR(2)  =  mesh%edge(ie)%normal%x * uR(1) + mesh%edge(ie)%normal%y * vR(1)
            vR(2)  =  mesh%edge(ie)%normal%x * vR(1) - mesh%edge(ie)%normal%y * uR(1)

         end if

         !=============================================================================================================!
         !   New reconstructed well balanced water depth
         !=============================================================================================================!

         hL(2)  =  max( 0._rp , hL(1) + zL - max( zL , zR ) )
         hR(2)  =  max( 0._rp , hR(1) + zR - max( zL , zR ) )
    
         !=============================================================================================================!
         !  Calling the balanced HLLC Solver dedicated to Shallow-Water Equations
         !=============================================================================================================!

         call sw_hllc( hL(2) , uL(2) , vL(2) , &
                       hR(2) , uR(2) , vR(2) , nflux )

         !=============================================================================================================!
         !  Boundary post treatment :
         !    - Feedback control of bathy_cell in ghost cells to properly control the Qin imposed
         !    - Calculation of nflux sum for each inflow
         !=============================================================================================================!

         if ( mesh%edge(ie)%boundary ) call boundary_post( nflux(1) , iR , mesh )

         !=============================================================================================================!
         !  Flux rotation and summation (as antisymmetric part to save time computation)
         !=============================================================================================================!

         lflux(1)  =                           nflux(1)
         lflux(2)  =  mesh%edge(ie)%normal%x * nflux(2)  -  mesh%edge(ie)%normal%y * nflux(3)
         lflux(3)  =  mesh%edge(ie)%normal%y * nflux(2)  +  mesh%edge(ie)%normal%x * nflux(3)

         lflux(1:3)  =  lflux(1:3)  *  mesh%edge(ie)%length

         tflux( 1 , iL )  =  tflux( 1 , iL )  +  lflux(1)
         tflux( 2 , iL )  =  tflux( 2 , iL )  +  lflux(2)
         tflux( 3 , iL )  =  tflux( 3 , iL )  +  lflux(3)

         tflux( 2 , iL )  =  tflux( 2 , iL )  +  mesh%edge(ie)%normal%x * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                 ( hL(1)**2 - hL(2)**2 ) )

         tflux( 3 , iL )  =  tflux( 3 , iL )  +  mesh%edge(ie)%normal%y * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                 ( hL(1)**2 - hL(2)**2 ) )

         if ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain ) then

            tflux( 1 , iR )  =  tflux( 1 , iR )  -  lflux(1)
            tflux( 2 , iR )  =  tflux( 2 , iR )  -  lflux(2)
            tflux( 3 , iR )  =  tflux( 3 , iR )  -  lflux(3)

            tflux( 2 , iR )  =  tflux( 2 , iR )  -  mesh%edge(ie)%normal%x * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                    ( hR(1)**2 - hR(2)**2 ) )

            tflux( 3 , iR )  =  tflux( 3 , iR )  -  mesh%edge(ie)%normal%y * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                    ( hR(1)**2 - hR(2)**2 ) )

         end if

         if ( mesh%edge(ie)%boundary ) then
            if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' ) then

                tflux( 1 , iR )  =  tflux( 1 , iR )  -  lflux(1)
                tflux( 2 , iR )  =  tflux( 2 , iR )  -  lflux(2)
                tflux( 3 , iR )  =  tflux( 3 , iR )  -  lflux(3)

                tflux( 2 , iR )  =  tflux( 2 , iR )  -  mesh%edge(ie)%normal%x * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                        ( hR(1)**2 - hR(2)**2 ) )

                tflux( 3 , iR )  =  tflux( 3 , iR )  -  mesh%edge(ie)%normal%y * mesh%edge(ie)%length * 0.5_rp * g * ( &
                                                        ( hR(1)**2 - hR(2)**2 ) )

            endif
        endif

      end if

   end do

   !===================================================================================================================!
   !  Cumulative rain Calculation
   !===================================================================================================================!

   do k=1,bc%nb_rn

      bc%rain(k)%cumul = bc%rain(k)%cumul + dt*bc%rain(k)%qin

   end do

   !===================================================================================================================!
   !  Euler Time Step
   !===================================================================================================================!

   do i = 1,mesh%nc

      h  =  dof%h(i)
      u  =  dof%u(i)
      v  =  dof%v(i)

      dof%h(i)  =  max( 0._rp , h  -  dt * tflux(1,i) * mesh%cell(i)%invsurf )

    ! Add rain source term
    if (bc_rain == 1) then

      k = bc%rain_land(i)!mesh%cell(i)%rain !Get rain group for current cell

      if (k > 0) then !If the cell does have a rain value attributed

         if ( bc_infil == 2 ) then ! If SCS-type infiltration is selected and this cell does have an infiltration value attributed

            S = 25.4_rp * ( 1000._rp / abs(infil%SCS( infil%land( i ) )%CN) - 10._rp ) / 1000._rp

            if ( bc%rain(k)%cumul > abs(infil%SCS( infil%land(i) )%lambdacn) * S ) then
               Fn1 = S * abs( infil%SCS( infil%land( i ) )%lambdacn ) + &
                     S * ( bc%rain( k )%cumul -      abs(infil%SCS( infil%land( i ) )%lambdacn)   * S ) / &
                         ( bc%rain( k )%cumul + (1 - abs(infil%SCS( infil%land( i ) )%lambdacn) ) * S )

            else

               Fn1 = dof%infil(i) + dt*bc%rain( k )%qin

            endif

            dof%h( i     ) = dof%h( i ) + dt * bc%rain( k )%qin - Fn1 + dof%infil( i ) ! Output SCS-modified rain
            dof%infil( i ) = Fn1

         else !Unmodified rain

            dof%h( i     ) = dof%h( i ) + dt * bc%rain( k )%qin

         endif

       endif

    endif

     if ( bc_infil == 1 ) then !If Green-Ampt infiltration is selected

        if (infil%land(i) .ne. 0) then !If the current cell does have an infiltration value attributed

         aFn1 = dof%infil(i) + dt * infil%GA( infil%land( i ) )%Ks * ( 1._rp - infil%GA( infil%land( i ) )%DeltaTheta )

         bFn1 = infil%GA( infil%land( i ) )%Ks * dt * infil%GA( infil%land(i) )%DeltaTheta * &
                ( dof%infil(i) + dof%h(i) + infil%GA( infil%land(i) )%PsiF )

         Fn1 = ( aFn1 + sqrt( aFn1**2._rp + 4._rp * bFn1 ) ) / 2._rp

		 h_infil = dof%h(i) + dof%infil(i) - Fn1

		 if (h_infil  <  0._rp ) then

             Fn1 = dof%h( i ) + dof%infil( i )
			 h_infil = 0._rp

         endif

		 dof%h( i ) = h_infil  !Replace the local variable h_infil
         dof%infil( i ) = Fn1

        endif

      endif


      !================================================================================================================!
      !   Positivity cut-off
      !================================================================================================================!

      if ( dof%h(i) <= heps ) then

         dof%u(i)  =  0._rp
         dof%v(i)  =  0._rp

      else

         dof%u(i)  =  (  h * u  -  dt * ( tflux(2,i) * mesh%cell(i)%invsurf )  )  /  dof%h(i)
         dof%v(i)  =  (  h * v  -  dt * ( tflux(3,i) * mesh%cell(i)%invsurf )  )  /  dof%h(i)

         !=============================================================================================================!
         !   Semi-Implicit Treatment of Friction Source Term (Manning/Strickler Formula)
         !=============================================================================================================!

         if ( friction == 1 ) then

            vel  =  sqrt( dof%u( i )**2 + dof%v( i )**2 )

            sfl  =  dof%h( i )**d2p3 + sqrt( dof%h(  i)**d4p3 + 4._rp * dt * g * &
                   ( manning( land( i ) ) * dof%h( i )**manning_beta( land( i ) ))**2 * vel )

            sfl  =  2._rp * dof%h( i )**d2p3 / sfl

         else if ( friction == 2 ) then

            sfl  =  one - dt * manning( land( i ) )

         else

            sfl  =  1._rp

         end if

         dof%u( i )  =  dof%u( i ) * sfl
         dof%v( i )  =  dof%v( i ) * sfl

      end if

   end do

   !===================================================================================================================!
   !  Calling MPI and filling ghost cells
   !===================================================================================================================!

   call com_dof( dof , mesh )

   call com_var_r( bathy_cell , mesh )                   ! Required MPI Communication due to inverse variable dependency

END SUBROUTINE euler_time_step_first_b1

SUBROUTINE euler_time_step_first_b1_HB( dof , mesh , local_slopes ) !

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_time_screen                                                                                              !NOADJ
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(inout)  ::  mesh

   type( unk ), intent(inout)  ::  dof

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!


   double precision, dimension(mesh%nc,2) , intent(in) :: local_slopes
   
   integer(ip)  ::  iL , iR

   real(rp)  ::  hL(2) , uL(2) , vL(2) , zL

   real(rp)  ::  hR(2) , uR(2) , vR(2) , zR

   real(rp), dimension( sw_nb )  ::  nflux               ! Finite Volume normal edge flux
   real(rp), dimension( sw_nb )  ::  lflux               ! Finite Volume edge flux in (x,y) coordinates

   real(rp), dimension( sw_nb , mesh%nc )  ::  tflux     ! Finite Volume total flux for each cell

   real(rp)  ::  h , u , v                             ! Temporal primitive variables

   integer(ip)  ::  index_edge, index_neighbour, index_current, &
   valid_point_count, aux, i_cell  ! Left and Right cells indexes to edge

   real(rp)  :: SxL(2), SyL(2), coeff_x_L(2), coeff_y_L(2), &
               corrective_term_x_L(2), corrective_term_y_L(2)! Left  State in edge cell normal coordinates

   real(rp)  :: SxR(2), SyR(2), coeff_x_R(2), coeff_y_R(2), &
               corrective_term_x_R(2), corrective_term_y_R(2)

                           ! Right State in edge cell normal coordinates
   !real(rp)  ::  x_neighbour(4), y_neighbour(4), z_neighbour(4), h_neighbour(4)
   real(rp) :: xedge, yedge, distance_x, distance_y
   real :: modified_properties(5), arr1(5), epsilon_duplicates
   double precision, allocatable, dimension(:,:) :: properties_for_gradient
   
   !double precision, dimension(mesh%nc,2) :: WS_gradient
   double precision :: dHdx_L, dHdy_L, dbdx_L, dbdy_L
   double precision :: c_coef, d_coef
   real(rp) :: n_powerlaw

   n_powerlaw = 1/m_powerlaw_index
   C_m =  1._rp / ((2._rp * m_powerlaw_index + 3._rp) * (m_powerlaw_index + 2._rp)**2) 			     ! expression involving Power law index

   !real(rp)  ::  local_edge_slope, local_WS_gradient, dist_grav                            ! slope of the terrain from left state to right state
   !real(rp), dimension( mesh%nc )  ::  slopes
   !real(rp), dimension( mesh%nc )  ::  iL_array
   !real(rp), dimension( mesh%nc )  ::  ie_array

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   tflux(:,:)  =  0._rp

   do ie = 1,mesh%ne    ! loop for each edge of the mesh

      !================================================================================================================!
      !  Calculate Left (always the cell in analysis) and Right States
      !================================================================================================================!

      iL  =  mesh%edge(ie)%cell(1) ! id of the current cell for the edge
      iR  =  mesh%edge(ie)%cell(2) ! id of the neighbor cell for the edge

      hL(1)  =  dof%h( iL )
      hR(1)  =  dof%h( iR )

       if ( hL(1) > heps .or. hR(1) > heps ) then

         zL  =  bathy_cell( iL ) !+ global_bathy_shift(1)
         zR  =  bathy_cell( iR ) !+ global_bathy_shift(1)

         uL(1)  =  dof%u( iL )  ! velocity in x-direction of current cell
         vL(1)  =  dof%v( iL )  ! velocity in y-direction of current cell

         SxL(1) = local_slopes( iL, 1)    ! x-slope of the current cell
         SyL(1) = 0.000001    ! y-slope of the current cell

         !! velocity magnitude and slope in the normal direction of edge for current cell

         uL(2)  =  mesh%edge(ie)%normal%x * uL(1) + mesh%edge(ie)%normal%y * vL(1)
         vL(2)  =  mesh%edge(ie)%normal%x * vL(1) - mesh%edge(ie)%normal%y * uL(1)


         if ( mesh%edge(ie)%boundary) then

               !! for dambreak, hL should be read from the file. so hL(1) = value read
               !! but this loop refers too all boundaries. i need to modify just the inlet
               !! but before i would need an interpolation i think

                call calc_boundary_state( mesh , hL(1) , zL , uL(2) , vL(2) , &
                                                 hR(1) , zR , uR(2) , vR(2) )

         else  ! for internal cells

            uR(1)  =  dof%u( iR ) ! velocity in x-direction of neighbor cell
            vR(1)  =  dof%v( iR ) ! velocity in y-direction of neighbor cell


            SxR(1)  =  local_slopes( iR, 1 ) ! x-slope of neighbor cell
            SyR(1)  =  0.000001 ! y-slope of neighbor cell

            !! velocity magnitude in the normal direction of edge for neighbor cell

            uR(2)  =  mesh%edge(ie)%normal%x * uR(1) + mesh%edge(ie)%normal%y * vR(1)
            vR(2)  =  mesh%edge(ie)%normal%x * vR(1) - mesh%edge(ie)%normal%y * uR(1)


         end if

         !=============================================================================================================!
         !   New reconstructed well balanced water depth
         !=============================================================================================================!

         hL(2)  =  max( 0._rp , hL(1) + zL - max( zL , zR ) )
         hR(2)  =  max( 0._rp , hR(1) + zR - max( zL , zR ) )

         corrective_term_x_L(1) = C_m * ((rho*g/K_index)**2 * sin(ABS(SxL(1))) * &
                           sin(ABS(SxL(1))))** m_powerlaw_index * &
                           hL(2)**(2 * m_powerlaw_index + 3)

         corrective_term_y_L(1) = C_m * ((rho*g/K_index)**2 * sin(ABS(SyL(1))) * &
                           sin(ABS(SyL(1))))** m_powerlaw_index * &
                           hL(2)**(2 * m_powerlaw_index + 3)

         corrective_term_x_L(2)  =  mesh%edge(ie)%normal%x * corrective_term_x_L(1) + &
                     mesh%edge(ie)%normal%y * corrective_term_y_L(1)

         corrective_term_y_L(2)  =  mesh%edge(ie)%normal%x * corrective_term_y_L(1) - &
                     mesh%edge(ie)%normal%y * corrective_term_x_L(1)

         !===================================================================!

         corrective_term_x_R(1) = C_m * ((rho*g/K_index)**2 * sin(ABS(SxR(1))) * &
                           sin(ABS(SxR(1))))** m_powerlaw_index * &
                           hR(2)**(2 * m_powerlaw_index + 3)

         corrective_term_y_R(1) = C_m * ((rho*g/K_index)**2 * sin(ABS(SyR(1))) * &
                           sin(ABS(SyR(1))))** m_powerlaw_index * &
                           hR(2)**(2 * m_powerlaw_index + 3)

         corrective_term_x_R(2)  =  mesh%edge(ie)%normal%x * corrective_term_x_R(1) + &
                     mesh%edge(ie)%normal%y * corrective_term_y_R(1)

         corrective_term_y_R(2)  =  mesh%edge(ie)%normal%x * corrective_term_y_R(1) - &
                     mesh%edge(ie)%normal%y * corrective_term_x_R(1)

         !===================================================================!


         coeff_x_L(1) = ((2*n_powerlaw+1)/(3*n_powerlaw+2)) * &
         ((2*rho*g*(hL(2))*sin(ABS(SxL(1)))*(n_powerlaw+1)**2 + &
         tau_c*n_powerlaw*(4*n_powerlaw + 3)) / &
         (rho*g*(hL(2))*sin(ABS(SxL(1)))*(n_powerlaw+1)**2 + &
         2*n_powerlaw*tau_c*(n_powerlaw+1) + &
         (n_powerlaw**2*tau_c**2)/(rho*g*(hL(2))*sin(ABS(SxL(1))))))

         coeff_y_L(1) = ((2*n_powerlaw+1)/(3*n_powerlaw+2)) * &
         ((2*rho*g*(hL(2))*sin(ABS(SyL(1)))*(n_powerlaw+1)**2 + &
         tau_c*n_powerlaw*(4*n_powerlaw + 3)) / &
         (rho*g*(hL(2))*sin(ABS(SyL(1)))*(n_powerlaw+1)**2 + &
         2*n_powerlaw*tau_c*(n_powerlaw+1) + &
         (n_powerlaw**2*tau_c**2)/(rho*g*(hL(2))*sin(ABS(SyL(1))))))

         coeff_x_L(2)  =  mesh%edge(ie)%normal%x * coeff_x_L(1) + &
                     mesh%edge(ie)%normal%y * coeff_y_L(1)

         coeff_y_L(2)  =  mesh%edge(ie)%normal%x * coeff_y_L(1) - &
                     mesh%edge(ie)%normal%y * coeff_x_L(1)

         !===================================================================!

         coeff_x_R(1) = ((2*n_powerlaw+1)/(3*n_powerlaw+2)) * &
            ((2*rho*g*(hR(2))*sin(ABS(SxR(1)))*(n_powerlaw+1)**2 + &
            tau_c*n_powerlaw*(4*n_powerlaw + 3)) / &
            (rho*g*(hR(2))*sin(ABS(SxR(1)))*(n_powerlaw+1)**2 + &
            2*n_powerlaw*tau_c*(n_powerlaw+1) + &
            (n_powerlaw**2*tau_c**2)/(rho*g*(hR(2))*sin(ABS(SxR(1))))))

         coeff_y_R(1) = ((2*n_powerlaw+1)/(3*n_powerlaw+2)) * &
            ((2*rho*g*(hR(2))*sin(ABS(SyR(1)))*(n_powerlaw+1)**2 + &
            tau_c*n_powerlaw*(4*n_powerlaw + 3)) / &
            (rho*g*(hR(2))*sin(ABS(SyR(1)))*(n_powerlaw+1)**2 + &
            2*n_powerlaw*tau_c*(n_powerlaw+1) + &
            (n_powerlaw**2*tau_c**2)/(rho*g*(hR(2))*sin(ABS(SyR(1))))))

         coeff_x_R(2)  =  mesh%edge(ie)%normal%x * coeff_x_R(1) + &
                     mesh%edge(ie)%normal%y * coeff_y_R(1)

         coeff_y_R(2)  =  mesh%edge(ie)%normal%x * coeff_y_R(1) - &
                     mesh%edge(ie)%normal%y * coeff_x_R(1)
                     
         !=============================================================================================================!
         !  Calling the balanced HLLC Solver dedicated to Shallow-Water Equations
         !  (it receives the left and right states of 1D Riemann problem (h, u and v) and returns the nflux).
         !  Please note that the problem is already normalized for each edge (index 2 of h, u, v)
         !=============================================================================================================!

         ! i need to identify the WS gradients related to both cells (current and neighbor)
         ! i have the indices iL and iR, so i can assign them to the matrix WS_gradient, which has the size (mesh%nc, 2)
         ! so, for dHdx of current cell: WS_gradient(iL,1)
         !     for dHdx of neighbor cell: WS_gradient(iR,1)

         dHdx_L = 0!WS_gradient(iL,1)
         dHdy_L = 0!WS_gradient(iL,2)

         dbdx_L = local_slopes(iL, 1) ! already in rad
         dbdy_L = 0.000001!local_slopes(iL, 2) ! already in rad

         call sw_hllc_HB( hL(2) , uL(2) , vL(2) , &
                       hR(2) , uR(2) , vR(2) , &
                       dHdx_L, dHdy_L, dbdx_L, dbdy_L, &
                       coeff_x_L(2) , coeff_y_L(2) , &
                       coeff_x_R(2) , coeff_y_R(2) , &
                       corrective_term_x_L(2), corrective_term_y_L(2), &
                       corrective_term_x_R(2), corrective_term_x_R(2), nflux )

         !=============================================================================================================!
         !  Boundary post treatment :
         !    - Feedback control of bathy_cell in ghost cells to properly control the Qin imposed
         !    - Calculation of nflux sum for each inflow
         !=============================================================================================================!

         if ( mesh%edge(ie)%boundary ) call boundary_post( nflux(1) , iR , mesh )

         !=============================================================================================================!
         !  Flux rotation and summation (as antisymmetric part to save time computation)
         !=============================================================================================================!

         ! these fluxes correspond to the continuity, x-momentum and y-momentum
         ! nflux is the flux in the normal (local) coordinate of the cell
         ! since lflux(1) is scalar, it doesn't need to rotate/scale
         ! however, nflux(2) and nflux(3) are vectors in normal coordinates and need to be rotated and scaled
         ! the rotation/scaling takes them back to the global coordinates of the system (not normal), represented by lflux(2) and lflux(3).
         ! so, lflux(2) and lflux(3) are in the global coordinates x and y.


         lflux(1)  =                           nflux(1)
         lflux(2)  =  mesh%edge(ie)%normal%x * nflux(2)  -  mesh%edge(ie)%normal%y * nflux(3)
         lflux(3)  =  mesh%edge(ie)%normal%y * nflux(2)  +  mesh%edge(ie)%normal%x * nflux(3)

         lflux(1:3)  =  lflux(1:3)  *  mesh%edge(ie)%length

         tflux( 1 , iL )  =  tflux( 1 , iL )  +  lflux(1)
         tflux( 2 , iL )  =  tflux( 2 , iL )  +  lflux(2)
         tflux( 3 , iL )  =  tflux( 3 , iL )  +  lflux(3)


         tflux( 2 , iL )  =  tflux( 2 , iL )  +  mesh%edge(ie)%normal%x * mesh%edge(ie)%length * &
                      0.5_rp * g * cos(dbdx_L)* ( hL(1)**2 - hL(2)**2 )  &
                   + mesh%edge(ie)%normal%x * mesh%edge(ie)%length**2  * 0.5_rp * g * (( hL(1) + hL(2) )* sin(dbdx_L) )

         tflux( 3 , iL )  =  tflux( 3 , iL )  +  mesh%edge(ie)%normal%y * mesh%edge(ie)%length * &
                  0.5_rp * g * cos(dbdy_L) * ( hL(1)**2 - hL(2)**2 )  &
                   +mesh%edge(ie)%normal%y * mesh%edge(ie)%length**2  * 0.5_rp * g * (( hL(1) + hL(2) )* sin(dbdy_L) )


         if ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain ) then

            tflux( 1 , iR )  =  tflux( 1 , iR )  -  lflux(1)
            tflux( 2 , iR )  =  tflux( 2 , iR )  -  lflux(2)
            tflux( 3 , iR )  =  tflux( 3 , iR )  -  lflux(3)

         tflux( 2 , iR )  =  tflux( 2 , iR )  -  mesh%edge(ie)%normal%x * mesh%edge(ie)%length * &
                  0.5_rp * g * cos(dbdx_L)* ( hR(1)**2 - hR(2)**2 )  &
                   - mesh%edge(ie)%normal%x * mesh%edge(ie)%length**2  * 0.5_rp * g * (( hR(1) + hR(2) )* sin(dbdx_L) )

         tflux( 3 , iR )  =  tflux( 3 , iR )  -  mesh%edge(ie)%normal%y * mesh%edge(ie)%length * &
                  0.5_rp * g * cos(dbdy_L)* ( hR(1)**2 - hR(2)**2 ) &
                   - mesh%edge(ie)%normal%y * mesh%edge(ie)%length**2 * 0.5_rp * g * (( hR(1) + hR(2) )* sin(dbdy_L) )
         end if
      end if
   end do

   !===================================================================================================================!
   !  Euler Time Step
   !===================================================================================================================!

   do i = 1,mesh%nc

      h  =  dof%h(i)
      u  =  dof%u(i)
      v  =  dof%v(i)

      dof%h(i)  =  max( 0._rp , h  -  dt * tflux(1,i) * mesh%cell(i)%invsurf )


      !================================================================================================================!
      !   Positivity cut-off
      !================================================================================================================!

      if ( dof%h(i) <= heps ) then

         dof%u(i)  =  0._rp
         dof%v(i)  =  0._rp

      else

         !! equation A8
         dof%u(i)  =  (  h * u  -  dt * ( tflux(2,i) * mesh%cell(i)%invsurf )  )  /  dof%h(i)
         dof%v(i)  =  (  h * v  -  dt * ( tflux(3,i) * mesh%cell(i)%invsurf )  )  /  dof%h(i)
         !print *, 'u:', dof%u(i)
         !print *, 'h:', dof%h(i)
         !! now, we need to use the calculated h, u and v into the friction algorithm (splitted).

         call friction_euler_HB( dof , mesh, local_slopes ) ! local_slopes , WS_gradient

      end if

   end do

   !===================================================================================================================!
   !  Calling MPI and filling ghost cells
   !===================================================================================================================!

   call com_dof( dof , mesh  )

   call com_var_r( bathy_cell , mesh )                   ! Required MPI Communication due to inverse variable dependency

END SUBROUTINE euler_time_step_first_b1_HB
