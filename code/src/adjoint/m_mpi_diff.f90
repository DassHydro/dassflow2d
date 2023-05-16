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


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Module m_mpi
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


MODULE m_mpi_diff

   USE m_mpi

   #ifdef USE_SW_MONO
      USE m_model
   #endif

   implicit none


CONTAINS


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE com_var_i_diff( var , var_diff , mesh )

      USE m_mesh

      implicit none

      type(msh), intent(in)  ::  mesh

      integer(ip), dimension( mesh%nc + mesh%ncb ), intent(inout)  ::  var
      integer(ip), dimension( mesh%nc + mesh%ncb ), intent(inout)  ::  var_diff

      #ifdef USE_MPI

         call com_var_i( var      , mesh )
         call com_var_i( var_diff , mesh )

      #endif

   END SUBROUTINE com_var_i_diff


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE com_var_r_diff( var , var_diff , mesh )

      USE m_mesh

      implicit none

      type(msh), intent(in)  ::  mesh

      real(rp), dimension( mesh%nc + mesh%ncb ), intent(inout)  ::  var
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(inout)  ::  var_diff

      #ifdef USE_MPI

         call com_var_r( var      , mesh )
         call com_var_r( var_diff , mesh )

      #endif

   END SUBROUTINE com_var_r_diff


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE mpi_sum_r_diff( val , val_diff )

      implicit none

      real(rp), intent(inout)  ::  val
      real(rp), intent(inout)  ::  val_diff

      #ifdef USE_MPI

         call mpi_sum_r( val      )
         call mpi_sum_r( val_diff )

      #endif

   END SUBROUTINE mpi_sum_r_diff


   SUBROUTINE mpi_sum_i_diff( val , val_diff )

      implicit none

      integer(ip), intent(inout)  ::  val
      integer(ip), intent(inout)  ::  val_diff

      #ifdef USE_MPI

         call mpi_sum_i( val      )
         call mpi_sum_i( val_diff )

      #endif

   END SUBROUTINE mpi_sum_i_diff


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE mpi_max_r_diff( val , val_diff )

      implicit none

      real(rp), intent(inout)  ::  val
      real(rp), intent(inout)  ::  val_diff

      real(rp), dimension(np)  ::  temp , temp_diff

      integer(ip)  ::  index_m

      #ifdef USE_MPI

         if ( np > 1 ) then

            call MPI_ALLGATHER( val,       1 , realtype , temp      , 1 , realtype , mpi_comm_world , code )
            call MPI_ALLGATHER( val_diff , 1 , realtype , temp_diff , 1 , realtype , mpi_comm_world , code )

         else

            temp     (1)  =  val
            temp_diff(1)  =  val_diff

         end if

         index_m  =  maxloc( temp , 1 )

         val_diff  =  temp_diff( index_m )

         call mpi_max_r( val )

      #endif

   END SUBROUTINE mpi_max_r_diff


   SUBROUTINE mpi_max_i_diff( val , val_diff )

      implicit none

      integer(ip), intent(inout)  ::  val
      integer(ip), intent(inout)  ::  val_diff

      integer(ip), dimension(np)  ::  temp , temp_diff

      integer(ip)  ::  index_m

      #ifdef USE_MPI

         if ( np > 1 ) then

            call MPI_ALLGATHER( val,       1 , inttype , temp      , 1 , inttype , mpi_comm_world , code )
            call MPI_ALLGATHER( val_diff , 1 , inttype , temp_diff , 1 , inttype , mpi_comm_world , code )

         else

            temp     (1)  =  val
            temp_diff(1)  =  val_diff

         end if

         index_m  =  maxloc( temp , 1 )

         val_diff  =  temp_diff( index_m )

         call mpi_max_i( val )

      #endif

   END SUBROUTINE mpi_max_i_diff


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE mpi_min_r_diff( val , val_diff )

      implicit none

      real(rp), intent(inout)  ::  val
      real(rp), intent(inout)  ::  val_diff

      real(rp), dimension(np)  ::  temp , temp_diff

      integer(ip)  ::  index_m

      #ifdef USE_MPI

         if ( np > 1 ) then

            call MPI_ALLGATHER( val,       1 , realtype , temp      , 1 , realtype , mpi_comm_world , code )
            call MPI_ALLGATHER( val_diff , 1 , realtype , temp_diff , 1 , realtype , mpi_comm_world , code )

         else

            temp     (1)  =  val
            temp_diff(1)  =  val_diff

         end if

         index_m  =  minloc( temp , 1 )

         val_diff  =  temp_diff( index_m )

         call mpi_min_r( val )

      #endif

   END SUBROUTINE mpi_min_r_diff


   SUBROUTINE mpi_min_i_diff( val , val_diff )

      implicit none

      integer(ip), intent(inout)  ::  val
      integer(ip), intent(inout)  ::  val_diff

      integer(ip), dimension(np)  ::  temp , temp_diff

      integer(ip)  ::  index_m

      #ifdef USE_MPI

         if ( np > 1 ) then

            call MPI_ALLGATHER( val,       1 , inttype , temp      , 1 , inttype , mpi_comm_world , code )
            call MPI_ALLGATHER( val_diff , 1 , inttype , temp_diff , 1 , inttype , mpi_comm_world , code )

         else

            temp     (1)  =  val
            temp_diff(1)  =  val_diff

         end if

         index_m  =  minloc( temp , 1 )

         val_diff  =  temp_diff( index_m )

         call mpi_min_i( val )

      #endif

   END SUBROUTINE mpi_min_i_diff


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Calling MPI and communicating dof to fill ghost cells
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE com_dof_diff( dof , dof_diff , mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( msh ), intent(in   )  ::  mesh

      type( unk ), intent(inout)  ::  dof
      type( unk ), intent(inout)  ::  dof_diff

      !================================================================================================================!
      !
      !================================================================================================================!

      #ifdef USE_SW_MONO

         call com_var_r_diff( dof%h(:) , dof_diff%h(:) , mesh )
         call com_var_r_diff( dof%u(:) , dof_diff%u(:) , mesh )
         call com_var_r_diff( dof%v(:) , dof_diff%v(:) , mesh )

      #endif

   END SUBROUTINE com_dof_diff


END MODULE m_mpi_diff
