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
!> \file friction_imex.f90
!! \brief This file includes friction_imex routine.
!! \details The file includes only friction_imex routine (see doc friction_imex routine).

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> brief friction calculation for imex case 
SUBROUTINE friction_imex( dof_imp , dof_exp , rk1 , rk2 , mesh )

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_time_screen
   USE m_model

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   TYPE( msh ), intent(in   )  ::  mesh
   TYPE( unk ), intent(inout)  ::  dof_imp
   TYPE( unk ), intent(in   )  ::  dof_exp

   real(rp), intent(in)  ::  rk1 , rk2

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   real(rp)  ::  vel                                     ! Velocity norm

   real(rp)  ::  Ks

   real(rp)  ::  imp_part , exp_part_u , exp_part_v

!======================================================================================================================!
!  Begin Subroutine
!======================================================================================================================!

    do i = 1,mesh%nc

      !================================================================================================================!
      !   Positivity cut-off
      !================================================================================================================!

      if ( dof_imp%h(i) <= heps ) then

         dof_imp%u(i)  =  0._rp
         dof_imp%v(i)  =  0._rp

      else

         !=============================================================================================================!
         !   Semi-Implicit Treatment of Friction Source Term (Manning/Strickler Formula)
         !=============================================================================================================!

         if      ( friction == 1 ) then

            Ks  =  g * manning( land(i) )**2

            vel  =  sqrt( dof_exp%u(i)**2 + dof_exp%v(i)**2 )

            if ( dof_exp%h(i) <= heps ) then

               exp_part_u  =  0._rp
               exp_part_v  =  0._rp

            else

               exp_part_u  =  rk2 * dt * Ks * dof_exp%u(i) * vel / dof_exp%h(i)**d1p3
               exp_part_v  =  rk2 * dt * Ks * dof_exp%v(i) * vel / dof_exp%h(i)**d1p3

            end if

            exp_part_u  =  dof_imp%u(i) - exp_part_u / dof_imp%h(i)
            exp_part_v  =  dof_imp%v(i) - exp_part_v / dof_imp%h(i)

            vel  =  sqrt( exp_part_u**2 + exp_part_v**2 )

            imp_part  =  dof_imp%h(i)**d2p3 + sqrt( dof_imp%h(i)**d4p3 + 4._rp * rk1 * dt * Ks * vel )

            imp_part  =  2._rp * dof_imp%h(i)**d2p3 / imp_part

         else if ( friction == 2 ) then

            Ks  =  manning( land(i) )

            exp_part_u  =  rk2 * dt * Ks * dof_exp%h(i) * dof_exp%u(i)
            exp_part_v  =  rk2 * dt * Ks * dof_exp%h(i) * dof_exp%v(i)

            exp_part_u  =  dof_imp%u(i) - exp_part_u / dof_imp%h(i)
            exp_part_v  =  dof_imp%v(i) - exp_part_v / dof_imp%h(i)

            imp_part  =  one / ( one + rk1 * dt * Ks )

         end if

         dof_imp%u(i)  =  exp_part_u * imp_part
         dof_imp%v(i)  =  exp_part_v * imp_part

      end if

   end do

   !===================================================================================================================!
   !  Calling MPI and filling ghost cells
   !===================================================================================================================!

   call com_dof( dof_imp , mesh )

END SUBROUTINE friction_imex
