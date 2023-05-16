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
!  Module m_mpi backward
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


MODULE m_mpi_back

   USE m_mpi
   USE m_time_screen

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


   SUBROUTINE com_var_i_back( var , var_back , mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( msh ), intent(in)  ::  mesh

      integer(ip), dimension( mesh%nc + mesh%ncb ), intent(inout)  ::  var
      integer(ip), dimension( mesh%nc + mesh%ncb ), intent(inout)  ::  var_back

      !================================================================================================================!
      !
      !================================================================================================================!

      #ifdef USE_MPI

         integer(ip), dimension(                   2*np )  ::  req
         integer(ip), dimension( MPI_status_size , 2*np )  ::  status_mpi

         integer(ip)  ::  ireq

         !=============================================================================================================!
         !  temporary vector that will contain every proc's neighbour's contribution
         !=============================================================================================================!

         integer(ip), dimension( mesh%nc + mesh%ncb , nneighb )  ::  temp_var

         integer(ip)  ::  count_v

         call Time_Init_Part(80)

         temp_var = 0

         ireq = 0

         do k = 0,np-1

            if ( part_neighbs( proc , k ) == 0 ) cycle

            ireq = ireq + 1

            !==========================================================================================================!
            !  RECV in direct become SEND in adjoint
            !==========================================================================================================!

            call MPI_ISEND( var_back( part_size( proc ) + 1 + sum( part_neighbs( proc , 0:k-1 ) , dim = 1 ) ) , &
                            part_neighbs( proc , k ) , inttype , k , k , mpi_comm_world , req(ireq) , code )

         end do

         count_v = 0

         do k = 0,np-1

            if ( part_neighbs( proc , k ) == 0 ) cycle

            ireq = ireq + 1

            count_v = count_v + 1

            !==========================================================================================================!
            !  SEND in direct become RECV in adjoint
            !==========================================================================================================!

            call MPI_IRECV( temp_var( 1 , count_v ) , 1 , type_com(k) , k , proc , mpi_comm_world , req(ireq) , code )

         end do

         call MPI_WAITALL( ireq , req , status_mpi , code )

         var_back(1:mesh%nc)  =  var_back(1:mesh%nc)  +  sum( temp_var(1:mesh%nc,:) , dim = 2 )

         !=============================================================================================================!
         !  Set to 0 the adjoint variables sent
         !=============================================================================================================!

         var_back( mesh%nc + 1 : mesh%nc + sum( part_neighbs( proc , 0:np-1 ) , dim = 1 ) )  =  0_ip

         call Time_End_Part(80)

      #endif

   END SUBROUTINE com_var_i_back


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE com_var_r_back( var , var_back , mesh )

      USE m_mesh

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type(msh), intent(in)  ::  mesh

      real(rp), dimension( mesh%nc + mesh%ncb ), intent(inout)  ::  var
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(inout)  ::  var_back

      !================================================================================================================!
      !
      !================================================================================================================!

      #ifdef USE_MPI

         integer(ip), dimension(                   2*np )  ::  req
         integer(ip), dimension( MPI_status_size , 2*np )  ::  status_mpi

         integer(ip)  ::  ireq

         !=============================================================================================================!
         !  temporary vector that will contain every proc's neighbour's contribution
         !=============================================================================================================!

         real(rp), dimension( mesh%nc + mesh%ncb , nneighb )  ::  temp_var

         integer(ip)  ::  count_v

         call Time_Init_Part(80)

         temp_var = 0._rp

         ireq = 0

         do k = 0,np-1

            if ( part_neighbs( proc , k ) == 0 ) cycle

            ireq = ireq + 1

            !==========================================================================================================!
            !  RECV in direct become SEND in adjoint
            !==========================================================================================================!

            call MPI_ISEND( var_back( part_size( proc ) + 1 + sum( part_neighbs( proc , 0:k-1 ) , dim = 1 ) ) , &
                            part_neighbs( proc , k ) , realtype , k , k , mpi_comm_world , req(ireq) , code )

         end do

         count_v = 0

         do k = 0,np-1

            if ( part_neighbs( proc , k ) == 0 ) cycle

            ireq = ireq + 1

            count_v = count_v + 1

            !==========================================================================================================!
            !  SEND in direct become RECV in adjoint
            !==========================================================================================================!

            call MPI_IRECV( temp_var( 1 , count_v ) , 1 , type_com(k) , k , proc , mpi_comm_world , req(ireq) , code )

         end do

         call MPI_WAITALL( ireq , req , status_mpi , code )

         var_back(1:mesh%nc)  =  var_back(1:mesh%nc)  +  sum( temp_var(1:mesh%nc,:) , dim = 2 )

         !=============================================================================================================!
         !  Set to 0 the adjoint variables sent
         !=============================================================================================================!

         var_back( mesh%nc + 1 : mesh%nc + sum( part_neighbs( proc , 0:np-1 ) , dim = 1 ) )  =  0._rp

         call Time_End_Part(80)

      #endif

   END SUBROUTINE com_var_r_back


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE mpi_sum_r_back( val , val_back )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      real(rp), intent(inout)  ::  val
      real(rp), intent(inout)  ::  val_back

      !================================================================================================================!
      !
      !================================================================================================================!

      #ifdef USE_MPI

         call mpi_sum_r( val_back )

      #endif

   END SUBROUTINE mpi_sum_r_back


   SUBROUTINE mpi_sum_i_back( val , val_back )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      integer(ip), intent(inout)  ::  val
      integer(ip), intent(inout)  ::  val_back

      !================================================================================================================!
      !
      !================================================================================================================!

      #ifdef USE_MPI

         call mpi_sum_i( val_back )

      #endif

   END SUBROUTINE mpi_sum_i_back


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE mpi_max_r_back( val , val_back )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      real(rp), intent(inout)  ::  val
      real(rp), intent(inout)  ::  val_back

      !================================================================================================================!
      !
      !================================================================================================================!

      #ifdef USE_MPI

         real(rp), dimension(np)  ::  temp

         integer(ip)  ::  index_m

        if ( np > 1 ) then

           call MPI_ALLGATHER( val , 1 , realtype , temp , 1 , realtype , mpi_comm_world , code )

        else

           temp(1)  =  val

        end if

        index_m  =  maxloc( temp , 1 )

        call mpi_sum_r( val_back )

        if ( proc /= index_m - 1 )  val_back = 0._rp

      #endif

   END SUBROUTINE mpi_max_r_back


   SUBROUTINE mpi_max_i_back( val , val_back )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      integer(ip), intent(inout)  ::  val
      integer(ip), intent(inout)  ::  val_back

      !================================================================================================================!
      !
      !================================================================================================================!

      #ifdef USE_MPI

         real(rp), dimension(np)  ::  temp

         integer(ip)  ::  index_m

        if ( np > 1 ) then

           call MPI_ALLGATHER( val , 1 , inttype , temp , 1 , inttype , mpi_comm_world , code )

        else

           temp(1)  =  val

        end if

        index_m  =  maxloc( temp , 1 )

        call mpi_sum_i( val_back )

        if ( proc /= index_m - 1 )  val_back = 0

      #endif

   END SUBROUTINE mpi_max_i_back


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE mpi_min_r_back( val , val_back )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      real(rp), intent(inout)  ::  val
      real(rp), intent(inout)  ::  val_back

      !================================================================================================================!
      !
      !================================================================================================================!

      #ifdef USE_MPI

         real(rp), dimension(np)  ::  temp

         integer(ip)  ::  index_m

        if ( np > 1 ) then

           call MPI_ALLGATHER( val , 1 , realtype , temp , 1 , realtype , mpi_comm_world , code )

        else

           temp(1)  =  val

        end if

        index_m  =  minloc( temp , 1 )

        call mpi_sum_r( val_back )

        if ( proc /= index_m - 1 )  val_back = 0._rp

      #endif

   END SUBROUTINE mpi_min_r_back


   SUBROUTINE mpi_min_i_back( val , val_back )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      integer(ip), intent(inout)  ::  val
      integer(ip), intent(inout)  ::  val_back

      !================================================================================================================!
      !
      !================================================================================================================!

      #ifdef USE_MPI

         integer(ip), dimension(np)  ::  temp

         integer(ip)  ::  index_m

        if ( np > 1 ) then

           call MPI_ALLGATHER( val , 1 , realtype , temp , 1 , realtype , mpi_comm_world , code )

        else

           temp(1)  =  val

        end if

        index_m  =  minloc( temp , 1 )

        call mpi_sum_i( val_back )

        if ( proc /= index_m - 1 )  val_back = 0

      #endif

   END SUBROUTINE mpi_min_i_back


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Calling MPI and communicating dof to fill ghost cells
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE com_dof_back( dof , dof_back , mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( msh ), intent(in   )  ::  mesh

      type( unk ), intent(inout)  ::  dof
      type( unk ), intent(inout)  ::  dof_back

      !================================================================================================================!
      !
      !================================================================================================================!

      #ifdef USE_SW_MONO

         call com_var_r_back( dof%h(:) , dof_back%h(:) , mesh )
         call com_var_r_back( dof%u(:) , dof_back%u(:) , mesh )
         call com_var_r_back( dof%v(:) , dof_back%v(:) , mesh )

      #endif

   END SUBROUTINE com_dof_back


END MODULE m_mpi_back
