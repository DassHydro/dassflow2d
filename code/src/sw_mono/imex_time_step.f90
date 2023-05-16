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
!> \file imex_time_step.f90
!! \brief This file includes imex_time_step routine.
!! \details The file includes only imex_time_step routine (see doc imex_time_step routine).



!> \brief  Perform IMEX Time Step (dedicated to Shallow-Water Equations ?)
!! \return dof updated after this new timestep
SUBROUTINE imex_time_step( dof , mesh )

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
   TYPE( unk ), intent(inout)  ::  dof

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   TYPE( unk )  ::  dof_1 , dof_2 , dof_3 , dof_4

   integer(ip)  ::  icopy

   real(rp)  ::  fric2 , fric4

!======================================================================================================================!
!  Begin Subroutine
!======================================================================================================================!

   call alloc_dof( dof_1 , mesh )
   call alloc_dof( dof_2 , mesh )
   call alloc_dof( dof_3 , mesh )
   call alloc_dof( dof_4 , mesh )

   #ifdef USE_MPI
      icopy  =  mesh%nc + part_neighbs( proc , np )
   #else
      icopy  =  mesh%nc
   #endif

!======================================================================================================================!

   dof_1%h( 1 : icopy )  =  dof%h( 1 : icopy )
   dof_1%u( 1 : icopy )  =  dof%u( 1 : icopy )
   dof_1%v( 1 : icopy )  =  dof%v( 1 : icopy )

   dof_2%h( 1 : icopy )  =  dof%h( 1 : icopy )
   dof_3%h( 1 : icopy )  =  dof%h( 1 : icopy )

!======================================================================================================================!

   call friction_imex( dof_1 , dof , demi , zero , mesh )

!======================================================================================================================!

   dof_2%u( 1 : icopy )  =  2._rp * dof%u( 1 : icopy ) - dof_1%u( 1 : icopy )
   dof_2%v( 1 : icopy )  =  2._rp * dof%v( 1 : icopy ) - dof_1%v( 1 : icopy )

!======================================================================================================================!

   call friction_imex( dof_2 , dof , demi , zero , mesh )

!======================================================================================================================!

   dof_3%u( 1 : icopy )  =  dof_2%u( 1 : icopy )
   dof_3%v( 1 : icopy )  =  dof_2%v( 1 : icopy )

!======================================================================================================================!


    call muscl_aud_flux_n        ( dof_3 , dof , mesh )


!======================================================================================================================!

   dof_4%h( 1 : icopy )  =  dof_3%h( 1 : icopy )

   dof_4%u( 1 : icopy )  =  dof_1%u( 1 : icopy ) + dof_2%u( 1 : icopy ) + dof_3%u( 1 : icopy ) - 2._rp * dof%u( 1 : icopy )
   dof_4%v( 1 : icopy )  =  dof_1%v( 1 : icopy ) + dof_2%v( 1 : icopy ) + dof_3%v( 1 : icopy ) - 2._rp * dof%v( 1 : icopy )

!======================================================================================================================!

   call friction_imex( dof_4 , dof , demi , zero , mesh )

!======================================================================================================================!

   dof_1%h( 1 : icopy )  =  dof_4%h( 1 : icopy )
   dof_1%u( 1 : icopy )  =  dof_4%u( 1 : icopy )
   dof_1%v( 1 : icopy )  =  dof_4%v( 1 : icopy )

!======================================================================================================================!

	! remove the condition on the spatial scheme choice as only muscl_b1_b is used

   call muscl_aud_flux_n        ( dof_4 , dof , mesh )


!======================================================================================================================!

   dof%u( 1 : icopy )  =  0.5_rp * ( dof_4%h( 1 : icopy ) * dof_4%u( 1 : icopy ) - &
                                     dof_3%h( 1 : icopy ) * dof_3%u( 1 : icopy ) ) + &
                                     dof_1%h( 1 : icopy ) * dof_1%u( 1 : icopy )

   dof%v( 1 : icopy )  =  0.5_rp * ( dof_4%h( 1 : icopy ) * dof_4%v( 1 : icopy ) - &
                                     dof_3%h( 1 : icopy ) * dof_3%v( 1 : icopy ) ) + &
                                     dof_1%h( 1 : icopy ) * dof_1%v( 1 : icopy )

   dof%h( 1 : icopy )  =  0.5_rp * ( dof_3%h( 1 : icopy ) + dof_4%h( 1 : icopy ) )

!======================================================================================================================!

   where( dof%h( 1 : icopy ) <= heps )

      dof%u( 1 : icopy )  =  0._rp
      dof%v( 1 : icopy )  =  0._rp

   elsewhere

      dof%u( 1 : icopy )  =  dof%u( 1 : icopy ) / dof%h( 1 : icopy )
      dof%v( 1 : icopy )  =  dof%v( 1 : icopy ) / dof%h( 1 : icopy )

   end where

!======================================================================================================================!

   call dealloc_dof( dof_1 )
   call dealloc_dof( dof_2 )
   call dealloc_dof( dof_3 )
   call dealloc_dof( dof_4 )

END SUBROUTINE imex_time_step
