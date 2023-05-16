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
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE euler_time_step_first_b2i( dof , mesh )

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_time_screen                                                                                              !NOADJ
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in   )  ::  mesh

   type( unk ), intent(inout)  ::  dof

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  iL , iR                              ! Left and Right cells indexes to edge

   real(rp)  ::  hL , uL , vL , zL                       ! Left  State in edge cell normal coordinates
   real(rp)  ::  hR , uR , vR , zR                       ! Right State in edge cell normal coordinates

   real(rp), dimension( sw_nb )  ::  nflux               ! Finite Volume normal edge flux
   real(rp), dimension( sw_nb )  ::  lflux               ! Finite Volume edge flux in (x,y) coordinates

   real(rp), dimension( sw_nb , mesh%nc )  ::  tflux     ! Finite Volume total flux for each cell

   real(rp)  ::  h , u , v                               ! Temporal primitive variables

   real(rp)  ::  vel                                     ! Velocity norm

   real(rp)  ::  sfl                                     ! Manning

   real(rp), dimension(2)  ::  bathy_balanced            ! Balanced Bathymetry

   type( vec2d ), dimension(mesh%nc)  ::  h_z_eq         ! Balanced Bathymetry

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   tflux(:,:)  =  0._rp

   h_z_eq(:)%x  =  0._rp
   h_z_eq(:)%y  =  0._rp

   do ie = 1,mesh%ne

      !================================================================================================================!
      !  Calculate Left and Right States
      !================================================================================================================!

      iL  =  mesh%edge(ie)%cell(1)
      iR  =  mesh%edge(ie)%cell(2)

      zL  =  bathy_cell( iL )
      zR  =  bathy_cell( iR )

      hL  =  dof%h( iL )

      uL  =  mesh%edge(ie)%normal%x * dof%u( iL ) + &
             mesh%edge(ie)%normal%y * dof%v( iL )

      vL  =  mesh%edge(ie)%normal%x * dof%v( iL ) - &
             mesh%edge(ie)%normal%y * dof%u( iL )

      if ( mesh%edge(ie)%boundary ) then

         call calc_boundary_state( mesh , hL , zL , uL , vL , &
                                          hR , zR , uR , vR )

      else

         hR  =  dof%h( iR )

         uR  =  mesh%edge(ie)%normal%x * dof%u( iR ) + &
                mesh%edge(ie)%normal%y * dof%v( iR )

         vR  =  mesh%edge(ie)%normal%x * dof%v( iR ) - &
                mesh%edge(ie)%normal%y * dof%u( iR )

      end if

      if ( hR > heps .or. hL > heps ) then

         !=============================================================================================================!
         !  Wet-Dry Front Treatment
         !=============================================================================================================!

         if ( zL + hL < zR .and. hR <= heps ) then

            uL  =  zero

            zR  =  zL + hL

         end if

         if ( zR + hR < zL .and. hL <= heps ) then

            uR  =  zero

            zL  =  zR + hR

         end if

         !=============================================================================================================!
         !  Calling the balanced HLLC Solver dedicated to Shallow-Water Equations
         !=============================================================================================================!

         call sw_balanced_hllc_src_out( hL , uL , vL , zL , &
                                        hR , uR , vR , zR , nflux , bathy_balanced )

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
         lflux(2)  =  mesh%edge(ie)%normal%x * nflux(2) - mesh%edge(ie)%normal%y * nflux(3)
         lflux(3)  =  mesh%edge(ie)%normal%y * nflux(2) + mesh%edge(ie)%normal%x * nflux(3)

         lflux(1:3)  =  lflux(1:3)  *  mesh%edge(ie)%length

         tflux( 1 , iL )  =  tflux( 1 , iL )  +  lflux(1)
         tflux( 2 , iL )  =  tflux( 2 , iL )  +  lflux(2)
         tflux( 3 , iL )  =  tflux( 3 , iL )  +  lflux(3)

         h_z_eq(iL)%x  =  h_z_eq(iL)%x + ( bathy_balanced(1) + bathy_balanced(2) ) * &
                                           mesh%edge(ie)%normal%x * &
                                           mesh%edge(ie)%length * &
                                           mesh%cell(iL)%invsurf

         h_z_eq(iL)%y  =  h_z_eq(iL)%y + ( bathy_balanced(1) + bathy_balanced(2) ) * &
                                           mesh%edge(ie)%normal%y * &
                                           mesh%edge(ie)%length * &
                                           mesh%cell(iL)%invsurf

         if ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain ) then

            tflux( 1 , iR )  =  tflux( 1 , iR )  -  lflux(1)
            tflux( 2 , iR )  =  tflux( 2 , iR )  -  lflux(2)
            tflux( 3 , iR )  =  tflux( 3 , iR )  -  lflux(3)

            h_z_eq(iR)%x  =  h_z_eq(iR)%x + ( bathy_balanced(1) - bathy_balanced(2) ) * &
                                              mesh%edge(ie)%normal%x * &
                                              mesh%edge(ie)%length * &
                                              mesh%cell(iR)%invsurf

            h_z_eq(iR)%y  =  h_z_eq(iR)%y + ( bathy_balanced(1) - bathy_balanced(2) ) * &
                                              mesh%edge(ie)%normal%y * &
                                              mesh%edge(ie)%length * &
                                              mesh%cell(iR)%invsurf

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

!      dof%h(i)  =  max( 0._rp , h  -  dt * tflux(1,i) * mesh%cell(i)%invsurf )

      dof%h(i)  =  h  -  dt * tflux(1,i) * mesh%cell(i)%invsurf

      if ( dof%h(i) < 0._rp ) call spread_mass_added( dof , mesh )

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

         if      ( friction == 1 ) then

            vel  =  sqrt( dof%u(i)**2 + dof%v(i)**2 )

            sfl  =  dof%h(i)**d2p3 + sqrt( dof%h(i)**d4p3 + 4._rp * dt * g * manning( land(i) )**2 * vel )

            sfl  =  2._rp * dof%h(i)**d2p3 / sfl

         else if ( friction == 2 ) then

            sfl  =  one - dt * manning( land(i) )

         else

            sfl  =  1._rp

         end if

         dof%u(i)  =  ( dof%u(i) - dt * h_z_eq(i)%x / dof%h(i) ) * sfl
         dof%v(i)  =  ( dof%v(i) - dt * h_z_eq(i)%y / dof%h(i) ) * sfl

      end if

   end do

   !===================================================================================================================!
   !  Calling MPI and filling ghost cells
   !===================================================================================================================!

   call com_dof( dof , mesh )

END SUBROUTINE euler_time_step_first_b2i
