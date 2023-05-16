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


SUBROUTINE rk_time_step( dof , mesh )

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

   TYPE( unk )  ::  dof_temp

   integer(ip)  ::  icopy

!======================================================================================================================!
!  Begin Subroutine
!======================================================================================================================!

   call alloc_dof( dof_temp , mesh )

   #ifdef USE_MPI
      icopy  =  mesh%nc + part_neighbs( proc , np )
   #else
      icopy  =  mesh%nc
   #endif

   dof_temp%h( 1 : icopy )  =  dof%h( 1 : icopy )
   dof_temp%u( 1 : icopy )  =  dof%u( 1 : icopy )
   dof_temp%v( 1 : icopy )  =  dof%v( 1 : icopy )

   select case( spatial_scheme )

      case( 'first' )

         call euler_time_step_first( dof , mesh )
         call euler_time_step_first( dof , mesh )

      case( 'first_b1' )

         call euler_time_step_first_b1( dof , mesh )
         call euler_time_step_first_b1( dof , mesh )

      case( 'first_b2e' )

         call euler_time_step_first_b2e( dof , mesh )
         call euler_time_step_first_b2e( dof , mesh )

      case( 'first_b2i' )

         call euler_time_step_first_b2i( dof , mesh )
         call euler_time_step_first_b2i( dof , mesh )

      case( 'first_imp' )

         call euler_time_step_imp_first( dof , mesh )
         call euler_time_step_imp_first( dof , mesh )

      case( 'muscl' )

         call euler_time_step_muscl( dof , mesh )
         call euler_time_step_muscl( dof , mesh )

      case( 'muscl_b' )

         call euler_time_step_muscl( dof , mesh )
         call euler_time_step_muscl( dof , mesh )

      case( 'muscl_b1' )

         call euler_time_step_muscl_b1( dof , mesh )
         call euler_time_step_muscl_b1( dof , mesh )

      case( 'muscl_b1_b' )

         call euler_time_step_muscl_b1( dof , mesh )
         call euler_time_step_muscl_b1( dof , mesh )

      case( 'muscl_b2' )

         call euler_time_step_muscl_b2( dof , mesh )
         call euler_time_step_muscl_b2( dof , mesh )

      case( 'muscl_b3' )

         call euler_time_step_muscl_b3( dof , mesh )
         call euler_time_step_muscl_b3( dof , mesh )

      case( 'muscl_imp' )

         call euler_time_step_imp_muscl( dof , mesh )
         call euler_time_step_imp_muscl( dof , mesh )

      case default

         call Stopping_Program_Sub( 'Unknow spatial scheme' )

   end select

   dof%u( 1 : icopy )  =  0.5_rp * ( dof_temp%h( 1 : icopy ) * dof_temp%u( 1 : icopy ) + &
                                     dof%h     ( 1 : icopy ) * dof%u     ( 1 : icopy ) )

   dof%v( 1 : icopy )  =  0.5_rp * ( dof_temp%h( 1 : icopy ) * dof_temp%v( 1 : icopy ) + &
                                     dof%h     ( 1 : icopy ) * dof%v     ( 1 : icopy ) )

   dof%h( 1 : icopy )  =  0.5_rp * ( dof_temp%h( 1 : icopy ) + &
                                     dof%h     ( 1 : icopy ) )

   where( dof%h( 1 : icopy ) <= heps )

      dof%u( 1 : icopy )  =  0._rp
      dof%v( 1 : icopy )  =  0._rp

   elsewhere

      dof%u( 1 : icopy )  =  dof%u( 1 : icopy ) / dof%h( 1 : icopy )
      dof%v( 1 : icopy )  =  dof%v( 1 : icopy ) / dof%h( 1 : icopy )

   end where

   call dealloc_dof( dof_temp )

   return

END SUBROUTINE rk_time_step
