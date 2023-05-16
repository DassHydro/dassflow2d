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


SUBROUTINE euler_time_step_imp_first( dof , mesh )

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_time_screen
   USE m_model
   USE m_linear_solver

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   TYPE( msh ), intent(in   )  ::  mesh
   TYPE( unk ), intent(inout)  ::  dof

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   integer(ip)  ::  iL , iR                                 ! Left and Right Cells indexes to Edge

   real(rp)  ::  hL , uL , vL , zL , dL                     ! Left  States
   real(rp)  ::  hR , uR , vR , zR , dR                     ! Right States

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

!======================================================================================================================!
!  Begin Subroutine
!======================================================================================================================!

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
   !  Calculating part 1 of upwinded edge velocity and gam_e coefficient (constant in Picard iterations)
   !===================================================================================================================!

   do ie = 1,mesh%ne

      iL  =  mesh%edge(ie)%cell(1)
      iR  =  mesh%edge(ie)%cell(2)

      uM  =  0.5_rp * ( mesh%edge(ie)%normal%x * ( dof%u( iL ) + dof%u( iR ) ) + &
                        mesh%edge(ie)%normal%y * ( dof%v( iL ) + dof%v( iR ) ) )

      v1(ie)%m  =  min( 0._rp , uM )
      v1(ie)%p  =  max( 0._rp , uM )

      hL  =  dof%h( iL )
      hR  =  dof%h( iR )

      dL  =  mesh%cell( iL )%peri * mesh%cell( iL )%invsurf

      if ( mesh%edge(ie)%boundary ) then

         dR  =  0._rp

      else

         dR  =  mesh%cell( iR )%peri * mesh%cell( iR )%invsurf

      end if

      gam_e(ie)  =  0.05_rp * dt * max( dL / hL , dR / hR )

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

         iL  =  mesh%edge(ie)%cell(1)
         iR  =  mesh%edge(ie)%cell(2)

         hL  =  dof%h( iL )
         hR  =  dof%h( iR )

         uM  =  gam_e(ie) * 0.5_rp * g * ( hR**2 - hL**2 )

         v2(ie)%m  =  min( 0._rp , uM )
         v2(ie)%p  =  max( 0._rp , uM )

         lambda  =  dt * mesh%edge(ie)%length

         vtot(ie)%p  =  ( v1(ie)%p - v2(ie)%m ) * lambda
         vtot(ie)%m  =  ( v1(ie)%m - v2(ie)%p ) * lambda

         lflux(1)  =  vtot(ie)%p * hL + vtot(ie)%m * hR

         if ( .not. mesh%edge(ie)%boundary ) then

            mumps_par%RHS( iL )  =  mumps_par%RHS( iL )  -  lflux(1) * mesh%cell( iL )%invsurf
            mumps_par%RHS( iR )  =  mumps_par%RHS( iR )  +  lflux(1) * mesh%cell( iR )%invsurf

            k = k + 1

            if ( hL >= hR ) then

               mumps_par%A( iL )  =  mumps_par%A( iL )  +  &

                                     (   vtot(ie)%p + lambda * gam_e(ie) * g * hL * hL ) * mesh%cell( iL )%invsurf

               mumps_par%A( k  )  =  (   vtot(ie)%m - lambda * gam_e(ie) * g * hR * hL ) * mesh%cell( iL )%invsurf

            else

               mumps_par%A( iL )  =  mumps_par%A( iL )  +  &

                                     (   vtot(ie)%p + lambda * gam_e(ie) * g * hL * hR ) * mesh%cell( iL )%invsurf

               mumps_par%A( k  )  =  (   vtot(ie)%m - lambda * gam_e(ie) * g * hR * hR ) * mesh%cell( iL )%invsurf

            end if

            k = k + 1

            if ( hR >= hL ) then

               mumps_par%A( iR )  =  mumps_par%A( iR )  +  &

                                     ( - vtot(ie)%m + lambda * gam_e(ie) * g * hR * hR ) * mesh%cell( iR )%invsurf

               mumps_par%A( k  )  =  ( - vtot(ie)%p - lambda * gam_e(ie) * g * hL * hR ) * mesh%cell( iR )%invsurf

            else

               mumps_par%A( iR )  =  mumps_par%A( iR )  +  &

                                     ( - vtot(ie)%m + lambda * gam_e(ie) * g * hR * hL ) * mesh%cell( iR )%invsurf

               mumps_par%A( k  )  =  ( - vtot(ie)%p - lambda * gam_e(ie) * g * hL * hL ) * mesh%cell( iR )%invsurf

            end if

         else

            mumps_par%RHS( iL )  =  mumps_par%RHS( iL )  -  lflux(1) * mesh%cell( iL )%invsurf

            if      ( hL > hR ) then

               mumps_par%A( iL )  =  mumps_par%A( iL )  +  &

                                     (   vtot(ie)%p + lambda * gam_e(ie) * g * hL * hL ) * mesh%cell( iL )%invsurf

            else if ( hR > hL ) then

               mumps_par%A( iL )  =  mumps_par%A( iL )  +  &

                                     (   vtot(ie)%p + lambda * gam_e(ie) * g * hL * hR ) * mesh%cell( iL )%invsurf

            end if

         end if

      end do

      norm_err  =  maxval( abs( mumps_par%RHS(:) ) )

      call DMUMPS( mumps_par )

      dof%h(1:mesh%nc)  =  dof%h(1:mesh%nc)  +  mumps_par%RHS(:)

      norm_err  =  min( norm_err , maxval( abs( mumps_par%RHS(:) ) ) )

      call fill_bc( dof , mesh )

      it_sol = it_sol + 1

!      write(6,'(I3,ES15.7)') it_sol , norm_err

   end do

!    write(6,'(I3,ES15.7)') it_sol , norm_err

   !===================================================================================================================!
   !  Calculating explicit fluxes for u^(n+1) and v^(n+1)
   !===================================================================================================================!

   tflux(:,:)  =  0._rp

   do ie = 1,mesh%ne

      iL  =  mesh%edge(ie)%cell(1)
      iR  =  mesh%edge(ie)%cell(2)

      hL  =  dof%h( iL )
      hR  =  dof%h( iR )

      uM  =  gam_e(ie) * 0.5_rp * g * ( hR**2 - hL**2 )

      v2(ie)%m  =  min( 0._rp , uM )
      v2(ie)%p  =  max( 0._rp , uM )

      vtot(ie)%p  =  v1(ie)%p - v2(ie)%m
      vtot(ie)%m  =  v1(ie)%m - v2(ie)%p

      lflux(1) = ( vtot(ie)%p * hL * dof%u( iL ) + vtot(ie)%m * hR * dof%u( iR ) + &
                   d1p4 * g * ( hL**2 + hR**2 ) * mesh%edge(ie)%normal%x ) * mesh%edge(ie)%length

      lflux(2) = ( vtot(ie)%p * hL * dof%v( iL ) + vtot(ie)%m * hR * dof%v( iR ) + &
                   d1p4 * g * ( hL**2 + hR**2 ) * mesh%edge(ie)%normal%y ) * mesh%edge(ie)%length

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

   return

END SUBROUTINE euler_time_step_imp_first
