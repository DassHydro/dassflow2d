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
!> \file m_linear_solver.f90
!! \brief This file includes m_linear_solver module.
!! \details The file includes only m_linear_solver module (see doc m_linear_solver module).

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Module controlling Internal and External Solvers
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> Module m_linear_solver.
!!
!! \details It is in this module that linear solver (mumps is defined and DGTSV) are defined.
!! no further exploration has been made
MODULE m_linear_solver

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_linear_algebra
   USE m_model

   implicit none

   #ifdef USE_MUMPS

      include 'dmumps_struc.h'

      type( dmumps_struc )  ::  mumps_par

   #elif USE_AGMG

      integer(ip)  ::  agmg_iter

   #endif


CONTAINS


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Main Subroutine to Initialize Solver(s)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE Init_Linear_Solver(mesh)

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type(msh), intent(in)  ::  mesh

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      integer(ip)  ::  iL , iR

      !================================================================================================================!
      !  Begin Subroutine
      !================================================================================================================!
      
!       print *, "Init_Linear_Solver" ; read(*,*)

      #if defined USE_MUMPS

         mumps_par%JOB   = - 1
         mumps_par%SYM   =   0
         mumps_par%PAR   =   1

         call DMUMPS( mumps_par )

      #endif

      #if defined USE_SW_MONO && USE_MUMPS

         mumps_par%N   =  mesh%nc
         mumps_par%NZ  =  mesh%nc

         do i = 1,mesh%nc

            do j = 1,mesh%cell(i)%nbed

               if ( mesh%cell(i)%cell(j) >= 1 .and. &
                    mesh%cell(i)%cell(j) <= mesh%nc ) mumps_par%NZ = mumps_par%NZ + 1

            end do

         end do

         allocate( mumps_par%IRN ( mumps_par%NZ ) )
         allocate( mumps_par%JCN ( mumps_par%NZ ) )

         k = 0

         do i = 1,mesh%nc

            k = k + 1

            mumps_par%IRN(k) = i
            mumps_par%JCN(k) = i

         end do

         do ie = 1,mesh%ne

            iL  =  mesh%edge(ie)%cell(1)
            iR  =  mesh%edge(ie)%cell(2)

            if ( .not. mesh%edge(ie)%boundary ) then

               k = k + 1

               mumps_par%IRN(k) = iL
               mumps_par%JCN(k) = iR

               k = k + 1

               mumps_par%IRN(k) = iR
               mumps_par%JCN(k) = iL

            end if

         end do

         allocate( mumps_par%A  ( mumps_par%NZ ) )
         allocate( mumps_par%RHS( mumps_par%N  ) )
         
         mumps_par%JOB  =  1

         call DMUMPS( mumps_par )

         mumps_par%JOB  =  5

!          mumps_par%JOB  =  6
          mumps_par%ICNTL(1) = 0
          mumps_par%ICNTL(2) = 0
          mumps_par%ICNTL(3) = 0
         mumps_par%ICNTL(4) = 0
         
         print *, "Init_Linear_Solver :: SW_MONO", mumps_par%N, loc(mumps_par%RHS), mumps_par%ICNTL(4)  ; read(*,*)

      #endif

      #if defined USE_NS_MULTIFLUID

         call Init_Jacobian( mesh , 'large_bandwith' )

         #ifdef USE_MUMPS

            call Init_Mumps_struct( jacobian )

            mumps_par%ICNTL(1) = 0
            mumps_par%ICNTL(2) = 0
            mumps_par%ICNTL(3) = 0
            mumps_par%ICNTL(7) = 3

            mumps_par%SYM  =  1
            mumps_par%JOB  =  1

            call DMUMPS( mumps_par )

            mumps_par%JOB  =  5

         #elif USE_AGMG

            call matrix_COO_to_CSR( jacobian )

         #endif

      #endif

   END SUBROUTINE Init_Linear_Solver


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Main Subroutine to Initialize Solver(s)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   #ifdef USE_MUMPS

      SUBROUTINE Init_MUMPS_struct( mat )

         implicit none

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         type( sys_lin_sparse ), intent(in)  ::  mat

         !=============================================================================================================!
         !  Local Variables
         !=============================================================================================================!

         integer(ip)  ::  elt

         !=============================================================================================================!
         !  Begin Subroutine
         !=============================================================================================================!

         mumps_par%N   =  mat%n
         mumps_par%NZ  =  mat%nz

         allocate( mumps_par%IRN ( mumps_par%NZ ) )
         allocate( mumps_par%JCN ( mumps_par%NZ ) )

         do elt = 1,mat%nz

            mumps_par%IRN( elt )  =  mat%ia( elt )
            mumps_par%JCN( elt )  =  mat%ja( elt )

         end do

         allocate( mumps_par%A  ( mumps_par%NZ ) )
         allocate( mumps_par%RHS( mumps_par%N  ) )

      END SUBROUTINE Init_MUMPS_struct


      SUBROUTINE Fill_MUMPS_struct( mat )

         implicit none

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         type( sys_lin_sparse ), intent(in)  ::  mat

         !=============================================================================================================!
         !  Begin Subroutine
         !=============================================================================================================!

         mumps_par%A  ( 1:mat%nz )  =  mat%a  ( 1:mat%nz )
         mumps_par%RHS( 1:mat%n  )  =  mat%rhs( 1:mat%n  )

      END SUBROUTINE Fill_MUMPS_struct

   #endif


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Convert a Sparse Matrix in COO format into CSR format
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE matrix_COO_to_CSR( mat )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( sys_lin_sparse ), intent(inout)  ::  mat

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      type( sys_lin_sparse )  ::  copy

      integer(ip)  ::  elt , row , ct( mat%nz )

      !================================================================================================================!
      !  Begin Subroutine
      !================================================================================================================!

      allocate( copy%ia ( mat%nz ) ) ; copy%ia  =  mat%ia
      allocate( copy%ja ( mat%nz ) ) ; copy%ja  =  mat%ja
!      allocate( copy%a  ( mat%nz ) ) ; copy%a   =  mat%a

!      copy%a = (/ ( i , i=1,mat%nz ) /)

      call reallocate_i( mat%ia , mat%n + 1 )

      mat%ia(:) = 1

      do elt = 1,mat%nz

         row = copy%ia( elt )

         mat%ia( row + 1 ) = mat%ia( row + 1 ) + 1

      end do

      do elt = 1,mat%n

         mat%ia( elt + 1 ) = mat%ia( elt + 1 ) + mat%ia( elt ) - 1

      end do

      ct(:) = 0

      do elt = 1,mat%nz

         row = mat%ia( copy%ia( elt ) )

         mat%ja( row + ct( row ) ) = copy%ja( elt )

!         mat%a( row + ct( row ) ) = copy%a( elt )

         mat%swap( row + ct( row ) ) = elt

         ct( row ) = ct( row ) + 1

      end do

!      write(6,'(200I3)') int(copy%a)
!      write(6,'(200I3)') int(copy%a( mat%swap(1:mat%nz) ) )
!      write(6,'(200I3)') int(mat%a)

!      STOP

      deallocate( copy%ia )
      deallocate( copy%ja )
!      deallocate( copy%a  )

   END SUBROUTINE matrix_COO_to_CSR


END MODULE m_linear_solver
