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

SUBROUTINE euler_time_step_first_b1( dof , mesh, poro_unit, it)
   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_time_screen !NOADJ
   USE m_model

   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(inout) :: mesh
   type( unk ), intent(inout) :: dof
   INTEGER, INTENT(IN) :: poro_unit
   INTEGER, INTENT(IN) :: it
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   INTEGER, PARAMETER :: write_frequency = 1000
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
   real(rp) :: h_infil !local variable of infil calculated depth
   real(rp) :: vel ! Velocity norm
   real(rp) :: sfl ! Manning
   real(rp) :: madd ! mass rain >TGADJ

   !Porosity variables
   real(rp)  ::  phiL , phiR                             ! Left/Right porosity in cell
   real(rp)  ::  s2L  , s2R                              ! Left/Right term source

   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   tflux(:,:) = 0._rp
   

#ifdef USE_PORO
   call update_all_porosities(dof, mesh)

IF (MOD(it, write_frequency) == 0 .OR. it == 1) THEN
   WRITE(poro_unit, *) it, SPorosity%phi(:)
END IF
#endif


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

         !zL = bathy_cell( iL )! + global_bathy_shift(1) !global_bathy_shift is not maintained
         !zR = bathy_cell( iR )! + global_bathy_shift(1)

         uL(1) = dof%u( iL )
         vL(1) = dof%v( iL )

         uL(2) = mesh%edge(ie)%normal%x * uL(1) + mesh%edge(ie)%normal%y * vL(1)
         vL(2) = mesh%edge(ie)%normal%x * vL(1) - mesh%edge(ie)%normal%y * uL(1)

         if ( mesh%edge(ie)%boundary) then

          zL = bathy_cell( iL )

            !================= Temporary modifications for some Andromede cases
            if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'zspresc') then
               zR = bathy_cell( iL ) !&
                     !- slope_y(1) * mesh%cell( mesh%edge(iL)%cell(1) )%surf / mesh%edge(iL)%length &
                     !- slope_x(1) * mesh%cell( mesh%edge(iL)%cell(1) )%surf / mesh%edge(iL)%length
            else if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'discharg1') then
               zR = bathy_cell( iL ) &
               + slope_y(1) * mesh%cell( mesh%edge(iL)%cell(1) )%surf / mesh%edge(iL)%length &
               + slope_x(1) * mesh%cell( mesh%edge(iL)%cell(1) )%surf / mesh%edge(iL)%length
            else if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'wall') then
               zR = bathy_cell( iL )
            endif
            !================= Temporary modifications for some Andromede cases

             if (.not. ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' )) then !do not call boundary calculations for internal BCs
                call calc_boundary_state( mesh , hL(1) , zL , uL(2) , vL(2) , &
                                                 hR(1) , zR , uR(2) , vR(2) )
             else

                uR(1) = dof%u( iR )
                vR(1) = dof%v( iR )

                uR(2) = mesh%edge(ie)%normal%x * uR(1) + mesh%edge(ie)%normal%y * vR(1)
                vR(2) = mesh%edge(ie)%normal%x * vR(1) - mesh%edge(ie)%normal%y * uR(1)

            endif
#ifdef USE_PORO
            phiL  =  SPorosity%Phi( SPorosity%land(iL) )
            phiR  =  phiL
#endif
         else

            zL = bathy_cell( iL )
            zR = bathy_cell( iR )

            uR(1) = dof%u( iR )
            vR(1) = dof%v( iR )

            uR(2) = mesh%edge(ie)%normal%x * uR(1) + mesh%edge(ie)%normal%y * vR(1)
            vR(2) = mesh%edge(ie)%normal%x * vR(1) - mesh%edge(ie)%normal%y * uR(1)

#ifdef USE_PORO
            phiL  =  SPorosity%Phi( SPorosity%land(iL) )
            phiR  =  SPorosity%Phi( SPorosity%land(iR) )
#endif

         end if

         !=============================================================================================================!
         ! New reconstructed well balanced water depth
         !=============================================================================================================!

         hL(2) = max( 0._rp , hL(1) + zL - max( zL , zR ) )
         hR(2) = max( 0._rp , hR(1) + zR - max( zL , zR ) )
         
         !=============================================================================================================!
         ! Calling the balanced HLLC Solver dedicated to Shallow-Water Equations
         !=============================================================================================================!
#ifdef USE_PORO
        if ( phiL > zerom .or. phiR > zerom ) then

            call sw_hllc_SP( hL(1) , uL(2) , vL(2) , zL , phiL , s2L , &
                             hR(1) , uR(2) , vR(2) , zR , phiR , s2R , nflux )

         else

            call sw_hllc_Impervious_SP( hL(2) , uL(2) , vL(2) , zL , phiL , s2L , &
                                        hR(2) , uR(2) , vR(2) , zR , phiR , s2R , nflux )

         end if

#else
         call sw_hllc( hL(2) , uL(2) , vL(2) , &
                       hR(2) , uR(2) , vR(2) , nflux )
#endif
         !=============================================================================================================!
         ! Boundary post treatment :
         ! - Feedback control of bathy_cell in ghost cells to properly control the Qin imposed
         ! - Calculation of nflux sum for each inflow
         !=============================================================================================================!
         
         if ( mesh%edge(ie)%boundary ) then
            call boundary_post( nflux(1) , iR , mesh )
         endif
         !=============================================================================================================!
         ! Flux rotation and summation (as antisymmetric part to save time computation)
         !=============================================================================================================!

         lflux(1) = nflux(1)
         lflux(2) = mesh%edge(ie)%normal%x * nflux(2) - mesh%edge(ie)%normal%y * nflux(3)
         lflux(3) = mesh%edge(ie)%normal%y * nflux(2) + mesh%edge(ie)%normal%x * nflux(3)

         lflux(1:3) = lflux(1:3) * mesh%edge(ie)%length

#ifdef USE_PORO
         tflux( 1 , iL )  =  tflux( 1 , iL )  +  lflux(1)
         tflux( 2 , iL )  =  tflux( 2 , iL )  +  lflux(2)  -  mesh%edge(ie)%normal%x * s2L * mesh%edge(ie)%length
         tflux( 3 , iL )  =  tflux( 3 , iL )  +  lflux(3)  -  mesh%edge(ie)%normal%y * s2L * mesh%edge(ie)%length

         if ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain ) then

            tflux( 1 , iR )  =  tflux( 1 , iR )  -  lflux(1)
            tflux( 2 , iR )  =  tflux( 2 , iR )  -  lflux(2)  -  mesh%edge(ie)%normal%x * s2R * mesh%edge(ie)%length
            tflux( 3 , iR )  =  tflux( 3 , iR )  -  lflux(3)  -  mesh%edge(ie)%normal%y * s2R * mesh%edge(ie)%length

         end if

         if ( mesh%edge(ie)%boundary ) then
            if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' ) then

               tflux( 1 , iR )  =  tflux( 1 , iR )  -  lflux(1)
               tflux( 2 , iR )  =  tflux( 2 , iR )  -  lflux(2)  -  mesh%edge(ie)%normal%x * s2R * mesh%edge(ie)%length
               tflux( 3 , iR )  =  tflux( 3 , iR )  -  lflux(3)  -  mesh%edge(ie)%normal%y * s2R * mesh%edge(ie)%length
            endif
         endif
#else
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
         endif

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

#endif
      endif
   enddo

   !===================================================================================================================!
   ! Cumulative rain Calculation
   !===================================================================================================================!
   
#ifdef USE_INFIL
   do k=1,bc%nb_rn
      bc%rain(k)%cumul = bc%rain(k)%cumul + dt*bc%rain(k)%qin
   end do
#endif

   !===================================================================================================================!
   ! Euler Time Step
   !===================================================================================================================!
   
   do i = 1,mesh%nc

      h = dof%h(i)
      u = dof%u(i)
      v = dof%v(i)

#ifdef USE_PORO
      dof%h(i)  =  max( 0._rp , h  -  dt / SPorosity%Phi( SPorosity%land(i) ) * tflux(1,i) * mesh%cell(i)%invsurf )
#else
      dof%h(i) = max( 0._rp , h - dt * tflux(1,i) * mesh%cell(i)%invsurf )
#endif

#ifdef USE_INFIL
    !===================================================================================================================!
    ! Adding rain source term for GA or SCS method
    !===================================================================================================================!

    if (bc_rain == 1) then ! If rain is accounted for

      k = bc%rain_land(i) !Get rain group for current cell

      if (k > 0) then !If the cell does have a rain value attributed

         if ( bc_infil == 2 ) then ! If SCS-type infiltration is selected and this cell does have an infiltration value attributed
         !Then substract SCS retention from raw rain data and inject the rest in the hydraulic model

            S = 25.4_rp * ( 1000._rp / abs(infil%SCS( infil%land( i ) )%CN) - 10._rp ) / 1000._rp
            if ( bc%rain(k)%cumul > abs(infil%SCS( infil%land(i) )%lambdacn) * S ) then
               Fn1 = S * abs( infil%SCS( infil%land( i ) )%lambdacn ) + &
                     S * ( bc%rain( k )%cumul - abs(infil%SCS( infil%land( i ) )%lambdacn) * S ) / &
                         ( bc%rain( k )%cumul + (1 - abs(infil%SCS( infil%land( i ) )%lambdacn) ) * S )
            else
               Fn1 = dof%infil(i) + dt*bc%rain( k )%qin
            endif
            dof%h( i ) = dof%h( i ) + dt * bc%rain( k )%qin - Fn1 + dof%infil( i ) ! Output SCS-modified rain
            dof%infil( i ) = Fn1

         else ! Raw rain data injected

            dof%h( i ) = dof%h( i ) + dt * bc%rain( k )%qin

         endif
       endif
    endif

     if ( bc_infil == 1 ) then !If Green-Ampt infiltration is selected

        if (infil%land(i) .ne. 0) then !If the current cell does have an infiltration value attributed
        ! Then use compute infiltrated depth at each cell using GA solving method on quadratic equation from Ni et al. (2020) : 10.1002/hyp.13722

          if (dof%infil(i) < infil%h_infil_max( infil%land(i)) ) then ! If current infiltrated depth is lower than the max soil depth, then apply computed infiltration (else do not infiltrate)

            aFn1 = dof%infil(i) + dt * infil%GA( infil%land( i ) )%Ks * ( 1._rp - infil%GA( infil%land( i ) )%DeltaTheta )
            bFn1 = infil%GA( infil%land( i ) )%Ks * dt * infil%GA( infil%land(i) )%DeltaTheta * &
                    ( dof%infil(i) + dof%h(i) + infil%GA( infil%land(i) )%PsiF )

            Fn1 = ( aFn1 + sqrt( aFn1**2._rp + 4._rp * bFn1 ) ) / 2._rp ! Infiltrated depth over current time step (F^(n+1))


            h_infil = dof%h(i)       + dof%infil(i) - Fn1
!           h^(n+1) = (h^n + r * Dt) + F^n          - F^(n+1)


            if (h_infil < 0._rp ) then ! If more than the available water depth should be infiltrated, infiltrate only the total depth (h^n + r * Dt)
                      Fn1 = dof%h( i ) + dof%infil( i )
                      h_infil = 0._rp
            endif

            dof%h( i ) = h_infil ! Pass h^(n+1) to dof%h (h^n)
            dof%infil( i ) = Fn1 ! Pass F^(n+1) to dof%infil (F^n)

        endif
      endif
    endif
#endif

      !================================================================================================================!
      ! Positivity cut-off
      !================================================================================================================!
      if ( dof%h(i) <= heps ) then
         dof%u(i) = 0._rp
         dof%v(i) = 0._rp
      else

#ifdef USE_PORO
         dof%u(i)  =  (  h * u  -  dt / SPorosity%Phi( SPorosity%land(i) ) * ( tflux(2,i) * mesh%cell(i)%invsurf )  )  /  dof%h(i)
         dof%v(i)  =  (  h * v  -  dt / SPorosity%Phi( SPorosity%land(i) ) * ( tflux(3,i) * mesh%cell(i)%invsurf )  )  /  dof%h(i)
#else
         dof%u(i) = ( h * u - dt * ( tflux(2,i) * mesh%cell(i)%invsurf ) ) / dof%h(i)
         dof%v(i) = ( h * v - dt * ( tflux(3,i) * mesh%cell(i)%invsurf ) ) / dof%h(i)
#endif

         !=============================================================================================================!
         ! Semi-Implicit Treatment of Friction Source Term (Manning/Strickler Formula)
         !=============================================================================================================!

         if ( friction == 1 ) then
            vel = sqrt( dof%u( i )**2 + dof%v( i )**2 )
            sfl = dof%h( i )**d2p3 + sqrt( dof%h( i)**d4p3 + 4._rp * dt * g * &
                   ( manning( land( i ) ) * dof%h( i )**manning_beta( land( i ) ))**2 * vel )
            sfl = 2._rp * dof%h( i )**d2p3 / sfl
         else if ( friction == 2 ) then
            sfl = one - dt * manning( land( i ) )
         else
            sfl = 1._rp
         end if

         dof%u( i ) = dof%u( i ) * sfl
         dof%v( i ) = dof%v( i ) * sfl
      end if
   end do

   !===================================================================================================================!
   ! Calling MPI and filling ghost cells
   !===================================================================================================================!

   call com_dof( dof , mesh )
   call com_var_r( bathy_cell , mesh ) ! Required MPI Communication due to inverse variable dependency



CONTAINS


    !==================================================================================================================!
    ! SUBROUTINE 1 : updates SPorosity (calculate the porosity of all cells)
    !==================================================================================================================!
    SUBROUTINE update_all_porosities(dof, mesh)
        IMPLICIT NONE
        TYPE(unk), INTENT(IN)  :: dof
        TYPE(msh), INTENT(IN)  :: mesh
        
        INTEGER  :: icell
        REAL(rp) :: H_k, wetted_area, phi_K_new, macro_area
        REAL(rp) :: W
        
        DO icell = 1, mesh%nc
            W = SPorosity%width(icell)
            H_k = dof%h(icell) + bathy_cell(icell)
            wetted_area = calculate_wetted_area(icell, H_k)
            macro_area = W * dof%h(icell)
            
            IF (macro_area > 1.0E-9_rp) THEN
                phi_K_new = wetted_area / macro_area
            ELSE
                phi_K_new = 1.0_rp 
            END IF
            SPorosity%phi(icell) = phi_K_new
        END DO
    END SUBROUTINE update_all_porosities

    
    !===============================================================================================================!
    ! FUNCTION 2 : calculates yN of a cell (half the width occupied by water)
    !===============================================================================================================!
    FUNCTION calculate_yn(H_k, a, beta, c) RESULT(yN)
        IMPLICIT NONE
        REAL(rp), INTENT(IN) :: H_k
        REAL(rp), INTENT(IN) :: a         !parabola parameter
        REAL(rp), INTENT(IN) :: beta      !parabola parameter
        REAL(rp), INTENT(IN) :: c         !parabola parameter
        REAL(rp) :: yN
        
        IF ((H_k - c) < 0.0_rp .OR. a <= 0.0_rp) THEN
            yN = 0.0_rp
            RETURN
        END IF
 
        yN = ((H_k - c) / a)**(1.0_rp / beta)
    END FUNCTION calculate_yn

    !===============================================================================================================!
    ! FUNCTION 3 : calculates the wetted area of a cell assuming that the bathymetry is a parabola 
    ! ay^beta + bathy_cell, a and beta are parameters fixed by the user
    !===============================================================================================================!
    FUNCTION calculate_wetted_area(icell, H_k) RESULT(area)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: icell
        REAL(rp), INTENT(IN) :: H_k
        REAL(rp) :: area
        REAL(rp) :: a, c, beta, yN, Hbanks       !parabola parameters + half the width occupied by water

        a     = SPorosity%a(icell)
        beta  = SPorosity%beta(icell)
        c     = bathy_cell(icell)
        Hbanks = SPorosity%hbanks(icell) 
        
        IF ((H_k - c) < 0.0_rp .OR. a <= 0.0_rp) THEN
            area = 0.0_rp
            RETURN
        ELSE
            yN = calculate_yn(H_k, a, beta, c)
            IF (H_k < Hbanks) THEN 
                area = (H_k - c) * yN - (a / (beta + 1.0_rp)) * yN**(beta + 1.0_rp)
                area = 2.0_rp * area
                area = MAX(0.0_rp, area)

            ELSE
                area = (Hbanks - c) * yN - (a / (beta + 1.0_rp)) * yN**(beta + 1.0_rp)
                area = 2.0_rp * area
                area = MAX(0.0_rp, area)
                area = MAX(0.0_rp, area + 2*yN*(H_k-Hbanks))
            END IF
        END IF 
    END FUNCTION calculate_wetted_area


END SUBROUTINE euler_time_step_first_b1



!> \brief  Perform Euler Time Step dedicated to Shallow-Water Equations
!! \return dof updated after this new timestep
!
! SUBROUTINE euler_time_step_first_b1_porosity( dof , mesh )
!
!    USE m_common
!    USE m_mesh
!    USE m_mpi
!    USE m_time_screen                                                                                              !NOADJ
!    USE m_model
!
!    implicit none
!
!    !===================================================================================================================!
!    !  Interface Variables
!    !===================================================================================================================!
!
!    type( msh ), intent(inout)  ::  mesh
!
!    type( unk ), intent(inout)  ::  dof
!
!    !===================================================================================================================!
!    !  Local Variables
!    !===================================================================================================================!
!
!    integer(ip)  ::  iL , iR   ! Left and Right cells indexes to edge
!
!    real(rp)  ::  hL(2) , uL(2) , vL(2) , zL              ! Left  State in edge cell normal coordinates
!    real(rp)  ::  hR(2) , uR(2) , vR(2) , zR              ! Right State in edge cell normal coordinates
!
!    real(rp)  ::  phiL , phiR                             ! Left/Right porosity in cell
!    real(rp)  ::  s2L  , s2R                              ! Left/Right term source
!
!    real(rp), dimension( sw_nb )  ::  nflux               ! Finite Volume normal edge flux
!    real(rp), dimension( sw_nb )  ::  lflux               ! Finite Volume edge flux in (x,y) coordinates
!
!    real(rp), dimension( sw_nb , mesh%nc )  ::  tflux     ! Finite Volume total flux for each cell
!
!    real(rp)  ::  h , u , v                               ! Temporal primitive variables
!
!    !Infiltration variables
!    real(rp)  :: S                                        ! potential maximal retention
!    real(rp)  :: Fn1                                      ! Temporal Fn+1
!    real(rp)  :: aFn1 , bFn1
!    real(rp)  :: h_infil                                  !local variable of infil calculated depth
!    real(rp)  ::  vel                                     ! Velocity norm
!    real(rp)  ::  sfl                                     ! Manning
!    real(rp)  :: madd                                     ! mass rain >TGADJ
!    !===================================================================================================================!
!    !  Begin Subroutine
!    !===================================================================================================================!
!
!
!    tflux(:,:)  =  0._rp
!
!    do ie = 1,mesh%ne
!
!       !================================================================================================================!
!       !  Calculate Left and Right States
!       !================================================================================================================!
!
!       iL  =  mesh%edge(ie)%cell(1)
!       iR  =  mesh%edge(ie)%cell(2) !Left cell id for a normal cell
!
!       ! Get right cell id for internal BC cells
!       if ( mesh%edge(ie)%boundary ) then !Check if bounfary first so typlim exists
!
!          if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_1D' ) cycle
!
!          if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' ) then !then change connectivity to connected 1D-like cell
!             iR  =  mesh%edge(ie)%cell1D2D !Get id of the single 1D-like cell with interface in the connected bc number => this should be done once!
!          endif
!
!       endif
!
!       hL(1)  =  dof%h( iL )
!       hR(1)  =  dof%h( iR )
!
!       if ( hL(1) > heps .or. hR(1) > heps ) then
!
!          !zL  =  bathy_cell( iL )! + global_bathy_shift(1)
!          !zR  =  bathy_cell( iR )! + global_bathy_shift(1)
!
!          uL(1)  =  dof%u( iL )
!          vL(1)  =  dof%v( iL )
!
!          uL(2)  =  mesh%edge(ie)%normal%x * uL(1) + mesh%edge(ie)%normal%y * vL(1)
!          vL(2)  =  mesh%edge(ie)%normal%x * vL(1) - mesh%edge(ie)%normal%y * uL(1)
!
!          if ( mesh%edge(ie)%boundary) then
!
!             zL  =  bathy_cell( iL )
!
!
!             !================= TEMP FOR ANDROMEDE
!             if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'zspresc') then
!
!                zR  =  bathy_cell( iL ) !&
!                      !- slope_y(1) * mesh%cell( mesh%edge(iL)%cell(1) )%surf / mesh%edge(iL)%length &
!                      !- slope_x(1) * mesh%cell( mesh%edge(iL)%cell(1) )%surf / mesh%edge(iL)%length
!
!             else if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'discharg1') then
!
!                zR  =  bathy_cell( iL ) &
!                + slope_y(1) * mesh%cell( mesh%edge(iL)%cell(1) )%surf / mesh%edge(iL)%length &
!                + slope_x(1) * mesh%cell( mesh%edge(iL)%cell(1) )%surf / mesh%edge(iL)%length
!
!             else if  ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'wall') then
!
!                zR  =  bathy_cell( iL )
!
!             endif
!             !================= END TEMP
!             !          write(*,*) mesh%edgeb(mesh%edge(ie)%lim)%typlim, iL, iR, zR, zL !NOADJ
!
!             if (.not. ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' )) then !do not call boundary calculations for internal BCs
!
!                call calc_boundary_state( mesh , hL(1) , zL , uL(2) , vL(2) , &
!                                                 hR(1) , zR , uR(2) , vR(2) )
!
!             else
!
!                uR(1)  =  dof%u( iR )
!                vR(1)  =  dof%v( iR )
!
!                uR(2)  =  mesh%edge(ie)%normal%x * uR(1) + mesh%edge(ie)%normal%y * vR(1)
!                vR(2)  =  mesh%edge(ie)%normal%x * vR(1) - mesh%edge(ie)%normal%y * uR(1)
!
!             endif
!
!             phiL  =  SPorosity%Phi( SPorosity%land(iL) )
!             phiR  =  phiL
!
!          else
!
!             zL  =  bathy_cell( iL )
!             zR  =  bathy_cell( iR )
!
!             uR(1)  =  dof%u( iR )
!             vR(1)  =  dof%v( iR )
!
!             uR(2)  =  mesh%edge(ie)%normal%x * uR(1) + mesh%edge(ie)%normal%y * vR(1)
!             vR(2)  =  mesh%edge(ie)%normal%x * vR(1) - mesh%edge(ie)%normal%y * uR(1)
!
!             phiL  =  SPorosity%Phi( SPorosity%land(iL) )
!             phiR  =  SPorosity%Phi( SPorosity%land(iR) )
!
!          end if
!
!          !=============================================================================================================!
!          !   New reconstructed well balanced water depth
!          !=============================================================================================================!
!
! !          hL(2)  =   hL(1) ! max( 0._rp , hL(1) + zL - max( zL , zR ) )
! !          hR(2)  =   hR(1) ! max( 0._rp , hR(1) + zR - max( zL , zR ) )
!
!          !=============================================================================================================!
!          !  Calling the balanced HLLC Solver dedicated to Shallow-Water Equations
!          !=============================================================================================================!
!
!          if ( phiL > zerom .or. phiR > zerom ) then
!
!             call sw_hllc_SP( hL(1) , uL(2) , vL(2) , zL , phiL , s2L , &
!                              hR(1) , uR(2) , vR(2) , zR , phiR , s2R , nflux )
!
!          else
!
!             call sw_hllc_Impervious_SP( hL(2) , uL(2) , vL(2) , zL , phiL , s2L , &
!                                         hR(2) , uR(2) , vR(2) , zR , phiR , s2R , nflux )
!
!          end if
!
!          !=============================================================================================================!
!          !  Boundary post treatment :
!          !    - Feedback control of bathy_cell in ghost cells to properly control the Qin imposed
!          !    - Calculation of nflux sum for each inflow
!          !=============================================================================================================!
!
!          if ( mesh%edge(ie)%boundary ) then
!             call boundary_post( nflux(1) , iR , mesh )
!          endif
!
!          !=============================================================================================================!
!          !  Flux rotation and summation (as antisymmetric part to save time computation)
!          !=============================================================================================================!
!
!          lflux(1)  =                           nflux(1)
!          lflux(2)  =  mesh%edge(ie)%normal%x * nflux(2)  -  mesh%edge(ie)%normal%y * nflux(3)
!          lflux(3)  =  mesh%edge(ie)%normal%y * nflux(2)  +  mesh%edge(ie)%normal%x * nflux(3)
!
!          lflux(1:3)  =  lflux(1:3)  *  mesh%edge(ie)%length
!
!          tflux( 1 , iL )  =  tflux( 1 , iL )  +  lflux(1)
!          tflux( 2 , iL )  =  tflux( 2 , iL )  +  lflux(2)  -  mesh%edge(ie)%normal%x * s2L * mesh%edge(ie)%length
!          tflux( 3 , iL )  =  tflux( 3 , iL )  +  lflux(3)  -  mesh%edge(ie)%normal%y * s2L * mesh%edge(ie)%length
!
!          if ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain ) then
!
!             tflux( 1 , iR )  =  tflux( 1 , iR )  -  lflux(1)
!             tflux( 2 , iR )  =  tflux( 2 , iR )  -  lflux(2)  -  mesh%edge(ie)%normal%x * s2R * mesh%edge(ie)%length
!             tflux( 3 , iR )  =  tflux( 3 , iR )  -  lflux(3)  -  mesh%edge(ie)%normal%y * s2R * mesh%edge(ie)%length
!
!          end if
!
!          if ( mesh%edge(ie)%boundary ) then
!             if ( mesh%edgeb(mesh%edge(ie)%lim)%typlim == 'internal_2D' ) then
!
!                tflux( 1 , iR )  =  tflux( 1 , iR )  -  lflux(1)
!                tflux( 2 , iR )  =  tflux( 2 , iR )  -  lflux(2)  -  mesh%edge(ie)%normal%x * s2R * mesh%edge(ie)%length
!                tflux( 3 , iR )  =  tflux( 3 , iR )  -  lflux(3)  -  mesh%edge(ie)%normal%y * s2R * mesh%edge(ie)%length
!             endif
!          endif
!
!       end if
!
!    end do
!
!    !===================================================================================================================!
!    !  Cumulative rain Calculation
!    !===================================================================================================================!
!
!    do k=1,bc%nb_rn
!
!       bc%rain(k)%cumul = bc%rain(k)%cumul + dt*bc%rain(k)%qin
!
!    end do
!
!    !===================================================================================================================!
!    !  Euler Time Step
!    !===================================================================================================================!
!
!    do i = 1,mesh%nc
!
!       h  =  dof%h(i)
!       u  =  dof%u(i)
!       v  =  dof%v(i)
!
!       dof%h(i)  =  max( 0._rp , h  -  dt / SPorosity%Phi( SPorosity%land(i) ) * tflux(1,i) * mesh%cell(i)%invsurf )
!
!       ! Add rain source term
!       if (bc_rain == 1) then
!
!          k = bc%rain_land(i)!mesh%cell(i)%rain !Get rain group for current cell
!
!          if (k > 0) then !If the cell does have a rain value attributed
!
!             if ( bc_infil == 2 ) then ! If SCS-type infiltration is selected and this cell does have an infiltration value attributed
!
!                S = 25.4_rp * ( 1000._rp / abs(infil%SCS( infil%land( i ) )%CN) - 10._rp ) / 1000._rp
!
!                if ( bc%rain(k)%cumul > abs(infil%SCS( infil%land(i) )%lambdacn) * S ) then
!                   Fn1 = S * abs( infil%SCS( infil%land( i ) )%lambdacn ) + &
!                         S * ( bc%rain( k )%cumul -      abs(infil%SCS( infil%land( i ) )%lambdacn)   * S ) / &
!                             ( bc%rain( k )%cumul + (1 - abs(infil%SCS( infil%land( i ) )%lambdacn) ) * S )
!
!                else
!
!                   Fn1 = dof%infil(i) + dt*bc%rain( k )%qin
!
!                endif
!
!                dof%h( i     ) = dof%h( i ) + dt * bc%rain( k )%qin - Fn1 + dof%infil( i ) ! Output SCS-modified rain
!                dof%infil( i ) = Fn1
!
!             else !Unmodified rain
!
!                dof%h( i     ) = dof%h( i ) + dt * bc%rain( k )%qin
!
!             endif
!
!          endif
!
!       endif
!
!       if ( bc_infil == 1 ) then !If Green-Ampt infiltration is selected
!
!          if (infil%land(i) .ne. 0) then !If the current cell does have an infiltration value attributed
!
!
!             !!!FOR INVERSE MODEL, MOVE LATER
!             !         infil%GA( infil%land( i ) )%Ks = abs(infil%GA( infil%land( i ) )%Ks)
!             !         infil%GA( infil%land( i ) )%DeltaTheta = abs(infil%GA( infil%land( i ) )%DeltaTheta)
!             !         infil%GA( infil%land( i ) )%PsiF = abs(infil%GA( infil%land( i ) )%PsiF)
!             !         infil%h_infil_max( infil%land(i)) ) = abs(infil%h_infil_max( infil%land(i)) )
!             !!!
!
!             aFn1 = dof%infil(i) + dt * infil%GA( infil%land( i ) )%Ks * ( 1._rp - infil%GA( infil%land( i ) )%DeltaTheta )
!
!             bFn1 = infil%GA( infil%land( i ) )%Ks * dt * infil%GA( infil%land(i) )%DeltaTheta * &
!                            ( dof%infil(i) + dof%h(i) + infil%GA( infil%land(i) )%PsiF )
!
!             Fn1 = ( aFn1 + sqrt( aFn1**2._rp + 4._rp * bFn1 ) ) / 2._rp
!
! 		      h_infil = dof%h(i) + dof%infil(i) - Fn1
!
! 		      if (h_infil  <  0._rp ) then
!
!                Fn1 = dof%h( i ) + dof%infil( i )
! 			      h_infil = 0._rp
!
!             endif
!
!
!             if (dof%infil(i) < infil%h_infil_max( infil%land(i)) ) then!1.2_rp) then !TEST MAX INFIL
!
!                dof%h( i ) = h_infil  !Replace the local variable h_infil
!                dof%infil( i ) = Fn1
!
!             endif
!
!          endif
!
!       endif
!
!
!       !================================================================================================================!
!       !   Positivity cut-off
!       !================================================================================================================!
!
!       if ( dof%h(i) <= heps ) then
!
!          dof%u(i)  =  0._rp
!          dof%v(i)  =  0._rp
!
!       else
!
!          dof%u(i)  =  (  h * u  -  dt / SPorosity%Phi( SPorosity%land(i) ) * ( tflux(2,i) * mesh%cell(i)%invsurf )  )  /  dof%h(i)
!          dof%v(i)  =  (  h * v  -  dt / SPorosity%Phi( SPorosity%land(i) ) * ( tflux(3,i) * mesh%cell(i)%invsurf )  )  /  dof%h(i)
!
!          !=============================================================================================================!
!          !   Semi-Implicit Treatment of Friction Source Term (Manning/Strickler Formula)
!          !=============================================================================================================!
!
!          if ( friction == 1 ) then
!
!             vel  =  sqrt( dof%u( i )**2 + dof%v( i )**2 )
!
!             sfl  =  dof%h( i )**d2p3 + sqrt( dof%h(  i)**d4p3 + 4._rp * dt * g * &
!                   ( manning( land( i ) ) * dof%h( i )**manning_beta( land( i ) ))**2 * vel )
!
!             sfl  =  2._rp * dof%h( i )**d2p3 / sfl
!
!          else if ( friction == 2 ) then
!
!             sfl  =  one - dt * manning( land( i ) )
!
!          else
!
!             sfl  =  1._rp
!
!          end if
!
!          dof%u( i )  =  dof%u( i ) * sfl
!          dof%v( i )  =  dof%v( i ) * sfl
!
!       end if
!
!    end do
!
!    !===================================================================================================================!
!    !  Calling MPI and filling ghost cells
!    !===================================================================================================================!
!
!    call com_dof( dof , mesh )
!
!    call com_var_r( bathy_cell , mesh )                   ! Required MPI Communication due to inverse variable dependency
!
! !    do i = 1, size(bc%sum_mass_flux)
! ! !     write(*,*) proc, tc, "BEFORE MPI SUM bc%sum_mass_flux", bc%sum_mass_flux
! !    call mpi_sum_r ( bc%sum_mass_flux(i) )
! ! !    write(*,*) proc, tc, "AFTER MPI SUM bc%sum_mass_flux", bc%sum_mass_flux
! !    enddo
!
!
! END SUBROUTINE euler_time_step_first_b1_porosity
