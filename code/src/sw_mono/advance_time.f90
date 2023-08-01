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
!> \file advance_time.f90
!! \brief This file includes advance_time routine.
!! \details The file includes only advance_time routine (see doc advance_time routine ).

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Advancement of Time Variables / Eventual Adaptative Time Step Calculation
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!>  Advancement of time variables / Eventual adaptative time step calculation
!!
!! \details This subroutine compute the advancement of time variables and compute time step \f$dt\f$ from CFL value
!! if it is necessary.
!! \param[in]    dof Unknowns of the model.
!! \param[in]    mesh Mesh of the problem.
!! \return set up end_time_loop= .TRUE. if it is the case
!! either return new time (tc) and time step value (dt)
!! \details various cases are treated :
!!  - adapt_dt == 1 : adaptative time step first version
!!      - if fix_time_step_serie == 2, timestep are forced from a file
!!      - if fix_time_step_serie == 1 .and. proc == 0   we WRITE time step generated in a file
!!  - adapt_dt == 2 : adaptative time step second version : first do the same treatment as '==1' option, then compare with dt calculated from edge and keep smalest
!!  - adapt_dt == else : constant time step
SUBROUTINE advance_time( dof , mesh )

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_model
   USE m_time_screen                                                                                              !NOADJ

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   TYPE( unk ), intent(in)  ::  dof
   TYPE( msh ), intent(in)  ::  mesh

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   real(rp)  ::  c    !> wave celerity ? (sqrt(g*h) )
   real(rp)  ::  vel  !> flow velocity (sqrt(u2+v2) )
   real(rp)  ::  dist !> distance calculated as if it was a squared mesh  dist = 2*surface/perimeter
   real(rp)  ::  s    !> la réponse D
   real(rp) :: dt_min !> smalest timestep , temporary allocation for loop on mesh cells
   real(rp) :: test   ! store velocity value in this variable in some cases

   integer(ip)  ::  iL , iR

   real(rp)  ::  hL , hR , hM , uM , dL , dR

   integer(ip) :: imin

!======================================================================================================================!
!  Begin Subroutine
!======================================================================================================================!

   nt = nt + 1

   if ( nt == max_nt_for_direct ) end_time_loop = .true.

!======================================================================================================================!
!  Adaptive Time Step
!======================================================================================================================!

   if ( adapt_dt == 1 ) then

      dt = dt

!<NOADJ
      if ( fix_time_step_serie == 2 ) then

         read(80) tc , dt , end_time_loop

         return

      end if

      dt  =  hugem    ! hugem is  Machine overflow - precision limits numbers

      do i = 1,mesh%nc

         if ( dof%h(i) > heps ) then

            dist = 2._rp * mesh%cell(i)%surf / mesh%cell(i)%peri

            c = sqrt( g * dof%h(i) )

            vel = sqrt( dof%u(i)**2 + dof%v(i)**2 )

            dt_min = min( dt , dist / ( vel + c ) )
            if (dt_min <= dt) then
				dt = dt_min
				imin = i
				test = vel
            end if

         end if

      end do



      call mpi_min_r( dt )

      dt = cfl * dt


	if(dt >10) dt = 10

      if ( tc + dt < ts ) then

         tc = tc + dt

      else

         dt = ts - tc ; tc  = ts

         end_time_loop = .true.

      end if

      call write_scalar_in_time( dt , 'time_step' )
      if ( fix_time_step_serie == 1 .and. proc == 0 ) write(80) tc , dt , end_time_loop

!>NOADJ

   else if ( adapt_dt == 2 ) then

      dt = dt

!<NOADJ

      dt  =  hugem

      do i = 1,mesh%nc

         if ( dof%h(i) > heps ) then

            dist = 2._rp * mesh%cell(i)%surf / mesh%cell(i)%peri

            c = sqrt( g * dof%h(i) )

            vel = sqrt( dof%u(i)**2 + dof%v(i)**2 )

            dt = min( dt , dist / ( vel + c ) )

         end if

      end do

      dt = cfl * dt

     call write_scalar_in_time( dt , 'time_step_exp' )

      dt  =  hugem

      do ie = 1,mesh%ne

         if ( .not. mesh%edge(ie)%boundary ) then

            iL  =  mesh%edge(ie)%cell(1)
            iR  =  mesh%edge(ie)%cell(2)

            hL  =  dof%h( iL )
            hR  =  dof%h( iR )

            dL  =  mesh%cell( iL )%peri * mesh%cell( iL )%invsurf
            dR  =  mesh%cell( iR )%peri * mesh%cell( iR )%invsurf

            uM  =  0.5_rp * ( mesh%edge(ie)%normal%x * ( dof%u( iL ) + dof%u( iR ) ) + &
                              mesh%edge(ie)%normal%y * ( dof%v( iL ) + dof%v( iR ) ) )

            uM  =  abs( uM )

            hM  =  max( dL * hR / hL , dR * hL / hR )

            dt = min ( dt , 1._rp / max( zerom , uM * hM ) )

         end if

      end do

      dt = cfl * dt

      if ( tc + dt < ts ) then

         tc = tc + dt

      else

         dt = ts - tc ; tc  = ts

         end_time_loop = .true.

      end if
      call write_scalar_in_time( dt , 'time_step_imp' )

!>NOADJ

!======================================================================================================================!
!  Constant Time Step
!======================================================================================================================!

   else

      tc = real(nt,rp) * dt

      if ( tc + dt + zerom > ts ) end_time_loop = .true.

   end if

!======================================================================================================================!
!  Output Time Informations
!======================================================================================================================!
   call Print_Screen( 'dt' )   !NOADJ

END SUBROUTINE advance_time
