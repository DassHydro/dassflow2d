!======================================================================================================================!
!
!                    DassFlow Version 2.0
!
!======================================================================================================================!
!
!  Copyright University of Toulouse-INSA & CNRS (France)
!
!  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
!  DassFlow is a computational software whose purpose is to simulate geophysical free surface flows,
!  designed for variational sensitivities and data assimilation (4D-var). Inverse capabilities are
!  based on the adjoint code generation by a source-to-source algorithmic differentiation (Tapenade software used).
!
!  DassFlow software includes few mostly independent "modules" with common architectures and structures:
!    - Shallow Module (Shallow Water Model, Finite Volume Method), i.e. the present code.
!    - 3D Module (Full Stokes Model, Finite Element Method, Mobile Gometries, ALE).
!  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
!
!  Many people have contributed to the DassFlow development from the initial version to the latest ones.
!  Current main developer:
!               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
!  with scientific and/or programming contributions of:
!               R. Madec   (Mathematics Institute of Toulouse IMT).
!               K. Larnier (Fluid Mechanics Institute of Toulouse IMFT).
!               J. Monnier (INSA & Mathematics Institute of Toulouse IMT).
!               J.-P. Vila (INSA & Mathematics Institute of Toulouse IMT).
!  and former other developers (M. Honnorat and J. Marin).
!
!  Scientific Contact : jerome.monnier@insa-toulouse.fr
!  Technical  Contact : frederic.couderc@math.univ-toulouse.fr
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


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Perform Euler Time Step dedicated to Shallow-Water Equations
!
!  Semi-Implicit Scheme with Picard fixed-point and/or Newton Methods to solve nonlinear system for water depth
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE euler_time_step_imp_muscl( dof , mesh )

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_time_screen
   USE m_model
   USE m_linear_solver
   USE m_numeric

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in   )  ::  mesh

   type( unk ), intent(inout)  ::  dof

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  iL , iR                                 ! Left and Right Cells indexes to Edge

   real(rp)  ::  hL(2) , uL(3) , vL(3) , zL , dL            ! Left  States in edge cell normal coordinates
   real(rp)  ::  hR(2) , uR(3) , vR(3) , zR , dR            ! Right States in edge cell normal coordinates

   real(rp)  ::  uM                                         ! Edge normal velocity

   real(rp)                        ::  u , v                ! Temporal velocities
   real(rp), dimension( mesh%nc )  ::  h                    ! Temporal depth

   TYPE up_vel

      real(rp)  ::  p , m                                   ! Structure defining positive and negative part for a real

   END TYPE up_vel

   type( up_vel ), dimension( mesh%ne )  ::  v1 , v2        ! Upwinded Edge normal velocities ( parts 1 & 2 )
   type( up_vel ), dimension( mesh%ne )  ::  vtot           ! Upwinded Edge normal velocities ( total )

   real(rp), dimension( mesh%ne )  ::  gam_e                ! Coefficient to stabilize semi-implicit scheme

   real(rp)                            ::  lflux(2)         ! Finite Volume edge flux in (x,y) coordinates
   real(rp), dimension( 2 , mesh%nc )  ::  tflux            ! Finite Volume total flux for each cell

   real(rp)  ::  norm_err                                   ! Stopping criteria for Picard iterations

   integer(ip)  ::  it_sol

   real(rp)  ::  lambda

   real(rp)  ::  aL , aR

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

#ifdef USE_MUMPS

   !===================================================================================================================!
   !  Storing h^n
   !===================================================================================================================!

   h(:)  =  dof%h(1:mesh%nc)

   !===================================================================================================================!
   !  Filling Ghost Cells with appropriate boundary condition
   !===================================================================================================================!

   call fill_bc( dof , mesh )

   !===================================================================================================================!
   !  Calculating no limited gradients by Least Square Method
   !===================================================================================================================!

   call Least_Square_Slope( dof%grad_h , dof%h , mesh )
   call Least_Square_Slope( dof%grad_u , dof%u , mesh )
   call Least_Square_Slope( dof%grad_v , dof%v , mesh )

!   call Slope_Limiter_Barth( dof%grad_h , dof%h , mesh )
!   call Slope_Limiter_Barth( dof%grad_u , dof%u , mesh )
!   call Slope_Limiter_Barth( dof%grad_v , dof%v , mesh )

!   dof%grad_h(:)%x  =  0._rp
!   dof%grad_h(:)%y  =  0._rp

   dof%grad_u(:)%x  =  0._rp
   dof%grad_u(:)%y  =  0._rp

   dof%grad_v(:)%x  =  0._rp
   dof%grad_v(:)%y  =  0._rp

   !===================================================================================================================!
   !  Calculating part 1 of upwinded edge velocity and gam_e coefficient (constant in Picard iterations)
   !===================================================================================================================!

   do ie = 1,mesh%ne

      iL  =  mesh%edge(ie)%cell(1)
      iR  =  mesh%edge(ie)%cell(2)

      hL(1)  =  dof%h( iL )
      hR(1)  =  dof%h( iR )

      zL     =  bathy_cell( iL )
      zR     =  bathy_cell( iR )

      uL(1)  =  dof%u( iL )
      uR(1)  =  dof%u( iR )

      vL(1)  =  dof%v( iL )
      vR(1)  =  dof%v( iR )

      hL(2)  =  hL(1)  +  ( dof%grad_h( iL )%x * mesh%edge(ie)%v_edge_cell(1)%x + &
                            dof%grad_h( iL )%y * mesh%edge(ie)%v_edge_cell(1)%y )

      uL(2)  =  uL(1)  +  ( dof%grad_u( iL )%x * mesh%edge(ie)%v_edge_cell(1)%x + &
                            dof%grad_u( iL )%y * mesh%edge(ie)%v_edge_cell(1)%y )

      vL(2)  =  vL(1)  +  ( dof%grad_v( iL )%x * mesh%edge(ie)%v_edge_cell(1)%x + &
                            dof%grad_v( iL )%y * mesh%edge(ie)%v_edge_cell(1)%y )

      hL(2)  =  MP_property( hL(2) , hL(1) , hR(1) )
      uL(2)  =  MP_property( uL(2) , uL(1) , uR(1) )
      vL(2)  =  MP_property( vL(2) , vL(1) , vR(1) )

      uL(3)  =  mesh%edge(ie)%normal%x * uL(2) + mesh%edge(ie)%normal%y * vL(2)
      vL(3)  =  mesh%edge(ie)%normal%x * vL(2) - mesh%edge(ie)%normal%y * uL(2)

      dL  =  mesh%cell( iL )%peri * mesh%cell( iL )%invsurf

      if ( mesh%edge(ie)%boundary ) then

         call calc_boundary_state( mesh , hL(2) , zL , uL(3) , vL(3) , &
                                          hR(2) , zR , uR(3) , vR(3) )

         dR  =  0._rp

      else

         hR(2)  =  hR(1)  +  ( dof%grad_h( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                               dof%grad_h( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )

         uR(2)  =  uR(1)  +  ( dof%grad_u( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                               dof%grad_u( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )

         vR(2)  =  vR(1)  +  ( dof%grad_v( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                               dof%grad_v( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )

         hR(2)  =  MP_property( hR(2) , hR(1) , hL(1) )
         uR(2)  =  MP_property( uR(2) , uR(1) , uL(1) )
         vR(2)  =  MP_property( vR(2) , vR(1) , vL(1) )

         uR(3)  =  mesh%edge(ie)%normal%x * uR(2) + mesh%edge(ie)%normal%y * vR(2)
         vR(3)  =  mesh%edge(ie)%normal%x * vR(2) - mesh%edge(ie)%normal%y * uR(2)

         dR  =  mesh%cell( iR )%peri * mesh%cell( iR )%invsurf

      end if

      uM  =  0.5_rp * ( uL(3) + uR(3) )

      v1(ie)%m  =  min( 0._rp , uM )
      v1(ie)%p  =  max( 0._rp , uM )

      gam_e(ie)  =  0.01_rp * dt * max( dL / hL(2) , dR / hR(2) )

   end do

   !===================================================================================================================!
   !  Peformaing Newton iterations to solve nonlinear system for h^(n+1)
   !===================================================================================================================!

   norm_err = hugem

   it_sol = 0

   do while ( it_sol <= 100 .and. norm_err > 10._rp * zerom )

      do i = 1,mesh%nc

         mumps_par%A  (i)  =  1._rp
         mumps_par%RHS(i)  =  h(i) - dof%h(i)

      end do

      k = mesh%nc

      do ie = 1,mesh%ne

         aL  =  0._rp
         aR  =  0._rp

         iL  =  mesh%edge(ie)%cell(1)
         iR  =  mesh%edge(ie)%cell(2)

         hL(1)  =  dof%h( iL )
         hR(1)  =  dof%h( iR )

         zL     =  bathy_cell( iL )
         zR     =  bathy_cell( iR )

         hL(2)  =  hL(1)  +  ( dof%grad_h( iL )%x * mesh%edge(ie)%v_edge_cell(1)%x + &
                               dof%grad_h( iL )%y * mesh%edge(ie)%v_edge_cell(1)%y )

         hL(2)  =  MP_property( hL(2) , hL(1) , hR(1) )

         if ( hL(2) /= hL(1) ) aL = ( hL(2) - hL(1) ) / ( hR(1) - hL(1) )

         if ( mesh%edge(ie)%boundary ) then

            call calc_boundary_state( mesh , hL(2) , zL , uL(3) , vL(3) , &
                                             hR(2) , zR , uR(3) , vR(3) )

         else

            hR(2)  =  hR(1)  +  ( dof%grad_h( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                                  dof%grad_h( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )

            hR(2)  =  MP_property( hR(2) , hR(1) , hL(1) )

            if ( hR(2) /= hR(1) ) aR = ( hR(2) - hR(1) ) / ( hL(1) - hR(1) )

         end if

         uM  =  gam_e(ie) * 0.5_rp * g * ( hR(2)**2 - hL(2)**2 )

         v2(ie)%m  =  min( 0._rp , uM )
         v2(ie)%p  =  max( 0._rp , uM )

         lambda  =  dt * mesh%edge(ie)%length

         vtot(ie)%p  =  ( v1(ie)%p - v2(ie)%m ) * lambda
         vtot(ie)%m  =  ( v1(ie)%m - v2(ie)%p ) * lambda

         lflux(1)  =  vtot(ie)%p * hL(2) + vtot(ie)%m * hR(2)

         if ( .not. mesh%edge(ie)%boundary ) then

            mumps_par%RHS( iL )  =  mumps_par%RHS( iL )  -  lflux(1) * mesh%cell( iL )%invsurf
            mumps_par%RHS( iR )  =  mumps_par%RHS( iR )  +  lflux(1) * mesh%cell( iR )%invsurf

            k = k + 1

            if ( hL(2) >= hR(2) ) then

               mumps_par%A( iL )  =  mumps_par%A( iL )  +  &

                                     (   vtot(ie)%p + lambda * gam_e(ie) * g * hL(1) * hL(1) ) * mesh%cell( iL )%invsurf

               mumps_par%A( k  )  =  (   vtot(ie)%m - lambda * gam_e(ie) * g * hR(1) * hL(1) ) * mesh%cell( iL )%invsurf

            else

               mumps_par%A( iL )  =  mumps_par%A( iL )  +  &

                                     (   vtot(ie)%p + lambda * gam_e(ie) * g * hL(1) * hR(1) ) * mesh%cell( iL )%invsurf

               mumps_par%A( k  )  =  (   vtot(ie)%m - lambda * gam_e(ie) * g * hR(1) * hR(1) ) * mesh%cell( iL )%invsurf

            end if

            k = k + 1

            if ( hR(2) >= hL(2) ) then

               mumps_par%A( iR )  =  mumps_par%A( iR )  +  &

                                     ( - vtot(ie)%m + lambda * gam_e(ie) * g * hR(1) * hR(1) ) * mesh%cell( iR )%invsurf

               mumps_par%A( k  )  =  ( - vtot(ie)%p - lambda * gam_e(ie) * g * hL(1) * hR(1) ) * mesh%cell( iR )%invsurf

            else

               mumps_par%A( iR )  =  mumps_par%A( iR )  +  &

                                     ( - vtot(ie)%m + lambda * gam_e(ie) * g * hR(1) * hL(1) ) * mesh%cell( iR )%invsurf

               mumps_par%A( k  )  =  ( - vtot(ie)%p - lambda * gam_e(ie) * g * hL(1) * hL(1) ) * mesh%cell( iR )%invsurf

            end if

         else

            mumps_par%RHS( iL )  =  mumps_par%RHS( iL )  -  lflux(1) * mesh%cell( iL )%invsurf

            if      ( hL(2) > hR(2) ) then

               mumps_par%A( iL )  =  mumps_par%A( iL )  +  &

                                     (   vtot(ie)%p + lambda * gam_e(ie) * g * hL(1) * hL(1) ) * mesh%cell( iL )%invsurf

            else if ( hR(2) > hL(2) ) then

               mumps_par%A( iL )  =  mumps_par%A( iL )  +  &

                                     (   vtot(ie)%p + lambda * gam_e(ie) * g * hL(1) * hR(1) ) * mesh%cell( iL )%invsurf

            end if

         end if

      end do

      norm_err  =  maxval( abs( mumps_par%RHS(:) ) )

      call DMUMPS( mumps_par )

      dof%h(1:mesh%nc)  =  dof%h(1:mesh%nc)  +  mumps_par%RHS(:)

      norm_err  =  min( norm_err , maxval( abs( mumps_par%RHS(:) ) ) )

      call fill_bc( dof , mesh )

!      call Least_Square_Slope( dof%grad_h , dof%h , mesh )

      it_sol = it_sol + 1

!      write(6,'(I3,ES15.7)') it_sol , norm_err

   end do

   write(6,'(I3,ES15.7)') it_sol , norm_err

   !===================================================================================================================!
   !  Calculating explicit fluxes for u^(n+1) and v^(n+1)
   !===================================================================================================================!

   tflux(:,:)  =  0._rp

   do ie = 1,mesh%ne

      iL  =  mesh%edge(ie)%cell(1)
      iR  =  mesh%edge(ie)%cell(2)

      hL(1)  =  dof%h( iL )
      hR(1)  =  dof%h( iR )

      zL     =  bathy_cell( iL )
      zR     =  bathy_cell( iR )

      uL(1)  =  dof%u( iL )
      uR(1)  =  dof%u( iR )

      vL(1)  =  dof%v( iL )
      vR(1)  =  dof%v( iR )

      hL(2)  =  hL(1)  +  ( dof%grad_h( iL )%x * mesh%edge(ie)%v_edge_cell(1)%x + &
                            dof%grad_h( iL )%y * mesh%edge(ie)%v_edge_cell(1)%y )

      hL(2)  =  MP_property( hL(2) , hL(1) , hR(1) )

      uL(2)  =  uL(1)  +  ( dof%grad_u( iL )%x * mesh%edge(ie)%v_edge_cell(1)%x + &
                            dof%grad_u( iL )%y * mesh%edge(ie)%v_edge_cell(1)%y )

      vL(2)  =  vL(1)  +  ( dof%grad_v( iL )%x * mesh%edge(ie)%v_edge_cell(1)%x + &
                            dof%grad_v( iL )%y * mesh%edge(ie)%v_edge_cell(1)%y )

      uL(2)  =  MP_property( uL(2) , uL(1) , uR(1) )
      vL(2)  =  MP_property( vL(2) , vL(1) , vR(1) )

      uL(3)  =  mesh%edge(ie)%normal%x * uL(2) + mesh%edge(ie)%normal%y * vL(2)
      vL(3)  =  mesh%edge(ie)%normal%x * vL(2) - mesh%edge(ie)%normal%y * uL(2)

      if ( mesh%edge(ie)%boundary ) then

         call calc_boundary_state( mesh , hL(2) , zL , uL(3) , vL(3) , &
                                          hR(2) , zR , uR(3) , vR(3) )

      else

         hR(2)  =  hR(1)  +  ( dof%grad_h( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                               dof%grad_h( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )

         hR(2)  =  MP_property( hR(2) , hL(1) , hR(1) )

         uR(2)  =  uR(1)  +  ( dof%grad_u( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                               dof%grad_u( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )

         vR(2)  =  vR(1)  +  ( dof%grad_v( iR )%x * mesh%edge(ie)%v_edge_cell(2)%x + &
                               dof%grad_v( iR )%y * mesh%edge(ie)%v_edge_cell(2)%y )

         uR(2)  =  MP_property( uR(2) , uR(1) , uL(1) )
         vR(2)  =  MP_property( vR(2) , vR(1) , vL(1) )

         uR(3)  =  mesh%edge(ie)%normal%x * uR(2) + mesh%edge(ie)%normal%y * vR(2)
         vR(3)  =  mesh%edge(ie)%normal%x * vR(2) - mesh%edge(ie)%normal%y * uR(2)

      end if

      uM  =  gam_e(ie) * 0.5_rp * g * ( hR(2)**2 - hL(2)**2 )

      v2(ie)%m  =  min( 0._rp , uM )
      v2(ie)%p  =  max( 0._rp , uM )

      vtot(ie)%p  =  v1(ie)%p - v2(ie)%m
      vtot(ie)%m  =  v1(ie)%m - v2(ie)%p

      lflux(1) = ( vtot(ie)%p * hL(2) * uL(2) + vtot(ie)%m * hR(2) * uR(2) + &
                   d1p4 * g * ( hL(2)**2 + hR(2)**2 ) * mesh%edge(ie)%normal%x ) * mesh%edge(ie)%length

      lflux(2) = ( vtot(ie)%p * hL(2) * vL(2) + vtot(ie)%m * hR(2) * vR(2) + &
                   d1p4 * g * ( hL(2)**2 + hR(2)**2 ) * mesh%edge(ie)%normal%y ) * mesh%edge(ie)%length

      tflux( 1 , iL )  =  tflux( 1 , iL )  +  lflux(1)
      tflux( 2 , iL )  =  tflux( 2 , iL )  +  lflux(2)

      if ( .not. mesh%edge(ie)%boundary ) then

         tflux( 1 , iR )  =  tflux( 1 , iR )  -  lflux(1)
         tflux( 2 , iR )  =  tflux( 2 , iR )  -  lflux(2)

      end if

   end do

   !===================================================================================================================!
   !  Updating u^(n+1) and v^(n+1)
   !===================================================================================================================!

   do i = 1,mesh%nc

      u  =  dof%u(i)
      v  =  dof%v(i)

      !================================================================================================================!
      !   Positivity cut-off
      !================================================================================================================!

      if ( dof%h(i) <= heps ) then

         dof%u(i)  =  0._rp
         dof%v(i)  =  0._rp

      else

         dof%u(i)  =  (  h(i) * u  -  dt * tflux(1,i) * mesh%cell(i)%invsurf ) / dof%h(i)
         dof%v(i)  =  (  h(i) * v  -  dt * tflux(2,i) * mesh%cell(i)%invsurf ) / dof%h(i)

      end if

   end do

#endif

END SUBROUTINE euler_time_step_imp_muscl
