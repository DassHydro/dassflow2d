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
!> \file m_time_screen.f90
!! \brief This file includes m_time_screen module.
!! \details The file includes only m_time_screen module (see doc m_time_screen module).

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Module controlling Time and Screen Outputs
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> Module m_time_screen.
!!
!! \details  Module controlling Time and Screen Outputs
MODULE m_time_screen
   USE m_common
   USE m_mesh
   USE m_mpi

   implicit none

   integer, parameter  ::  max_time_index = 100

   real(rp), dimension( max_time_index )  ::  t_ini
   real(rp), dimension( max_time_index )  ::  t_end
   real(rp), dimension( max_time_index )  ::  time


CONTAINS


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Initialization of the ième time clock
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief  Initialization of the ième time clock
!> \details return time value of the ieme time clock ??? 
!! force time value to zero ?
   SUBROUTINE Time_Init( i )

      implicit none

      intrinsic cpu_time

      integer(ip), intent(in)  ::  i

      time( i ) = 0._rp

      #ifdef USE_MPI
         t_ini( i ) = MPI_WTIME()
      #else
         call CPU_TIME( t_ini( i ) )
      #endif

   END SUBROUTINE Time_Init


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Stopping of the ième time clock
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief   Stopping of the ième time clock
!> \details return the total time  between initial time and end time
   SUBROUTINE Time_End( i )

      implicit none

      intrinsic cpu_time

      integer(ip), intent(in)  ::  i

      #ifdef USE_MPI
         t_end( i ) = MPI_WTIME()
      #else
         call CPU_TIME( t_end( i ) )
      #endif

      time( i )  =  t_end( i )  -  t_ini( i )

   END SUBROUTINE Time_End


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Initialization of the ième time clock
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief  Initialization of the ième time clock
!> \details return time value of the ieme time clock ???
!! DO NOT force time value to zero ?
   SUBROUTINE Time_Init_Part( i )

      implicit none

      intrinsic cpu_time

      integer(ip), intent(in)  ::  i

      #ifdef USE_MPI
         t_ini( i ) = MPI_WTIME()
      #else
         call CPU_TIME( t_ini( i ) )
      #endif

   END SUBROUTINE Time_Init_Part


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Stopping of the ième time clock
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief   Stopping of the ième time clock
!> \details return the total time  between initial time and end time
!!      time( i )  =  time( i )  +  t_end( i )  -  t_ini( i )
   SUBROUTINE Time_End_Part( i )

      implicit none

      intrinsic cpu_time

      integer(ip), intent(in)  ::  i

      #ifdef USE_MPI
         t_end( i ) = MPI_WTIME()
      #else
         call CPU_TIME( t_end( i ) )
      #endif

      time( i )  =  time( i )  +  t_end( i )  -  t_ini( i )

   END SUBROUTINE Time_End_Part


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Output to screen time clock(s)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief  Output to screen time clock(s)
!> \details following tags are treated
!!  - USE_MPI
!!  - USE_NS_MULTIFLUID && USE_MUMPS
!! following options are treated
!!  - proc == 0
!!  -  spatial_scheme(1:5) == 'muscl'
!!  - proc == 0
   SUBROUTINE Time_Screen( mesh )

      implicit none

      TYPE( msh ), intent(in)  ::  mesh

      integer(ip)  ::  mesh_total_size

      #ifdef USE_MPI
         call MPI_ALLREDUCE( mesh%nc , mesh_total_size , 1 , inttype , MPI_SUM , MPI_COMM_WORLD , code )
      #else
         mesh_total_size = mesh%nc
      #endif

      if ( proc == 0 ) then

         write(6,'(A)'        )
         write(6,'(A)'        ) '********************************************************************************'
         write(6,'(A,F15.2)'  ) ' Time of simulation   =   ' , time(1)
         write(6,'(A)'        ) '********************************************************************************'
         write(6,'(A,F15.2,A)') ' Performance          =   ' , time(1) * np / mesh_total_size / nt * 1.d6 , &
                                    ' microsecond / dx / dt / proc'
         write(6,'(A)'        ) '********************************************************************************'

         if ( spatial_scheme(1:5) == 'muscl' ) then

            write(6,'(A,F5.2)'   ) ' MUSCL reconstruction in percent  =  ' , 100._rp * time(2) / time(1)
            write(6,'(A,F5.2)'   ) ' MUSCL limitation     in percent  =  ' , 100._rp * time(3) / time(1)
            write(6,'(A,F5.2)'   ) ' Fluxes computation   in percent  =  ' , 100._rp * time(4) / time(1)
            write(6,'(A,F5.2)'   ) ' Variables update     in percent  =  ' , 100._rp * time(5) / time(1)
            write(6,'(A)'        ) '********************************************************************************'

         end if ! end if   ( spatial_scheme(1:5) == 'muscl' ) 

         #if defined USE_NS_MULTIFLUID && USE_MUMPS

            time(3) = time(3) - time(2)

            write(6,'(A,F5.2)') ' Mumps in percent        =   ' , 100._rp * time(2) / time(1)
            write(6,'(A,F5.2)') ' Transfert to Mumps      =   ' , 100._rp * time(4) / time(1)
            write(6,'(A,F5.2)') ' Assembling in percent   =   ' , 100._rp * time(3) / time(1)
            write(6,'(A)'     ) '********************************************************************************'

         #endif

         #if defined USE_MPI

            write(6,'(A,F5.2)') ' MPI communications   in percent  =  ' , 100._rp * time(80) / time(1)
            write(6,'(A)'     ) '********************************************************************************'

         #endif

      end if ! end if proc==0

   END SUBROUTINE Time_Screen


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Subroutine dedicated to Screen Output Control
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief  Subroutine dedicated to Screen Output Control)
!> \details following tags are treated
!! following options are treated
!!  - proc == 0
!!  - verbose >= 0
!!  - test_dt_just_after( 0.01_rp * ts )
!! following tags are treated
!!  - USE_SW_MONO
!!  - USE_NS_MULTIFLUID
!! following case can be proposed
!!  - start_mesh
!!  - end_mesh
!!  - start_direct
!!  - end_direct
!!  - start_testadj
!!  - end_testadj
!!  - start_grad_cost
!!  - end_grad_cost
!!  - start_minimize
!!  - end_minimize
!!  - result  (if verbose >=0)
!!  - dt    (  if verbose >=0 && test_dt_just_after( 0.01_rp * ts ))
!!  - number_limits
!!  - norms_h
!!  - norms_u
!!  - norms_q
!!  - tangent_gradient_test
!!  - backward_gradient_test 
!!  - backward_scalar_product_test
!!  - cost
!!  - cost_diff
!!  - cost_back
   SUBROUTINE Print_Screen( screen_case , var )

      implicit none

      intrinsic range , precision

      character(len=*), intent(in)  ::  screen_case

      real(rp), optional, intent(in)  ::  var


      if ( proc == 0 ) then

         select case( screen_case )

            case( 'start_mesh' )

               write(6,'(A)')
               write(6,'(A80)') '================================================================================'
               write(6,'(A80)') '*  Mesh Loading                                                                *'
               write(6,'(A80)') '================================================================================'

            case( 'end_mesh' )

               write(6,'(A80)') '================================================================================'
               write(6,'(A80)') '*  Mesh Loading is OK                                                          *'
               write(6,'(A80)') '================================================================================'

            case( 'start_direct' )

               write(6,'(A)')

              #if defined USE_SW_MONO

                  write(6,'(A80)') '================================================================================'
                  write(6,'(A80)') '*  Running Shallow-Water Model in direct mode                                  *'
                  write(6,'(A80)') '================================================================================'

               #elif defined USE_NS_MULTIFLUID

                  write(6,'(A80)') '================================================================================'
                  write(6,'(A80)') '*  Running Multifluid Navier-Stokes Model in direct mode                       *'
                  write(6,'(A80)') '================================================================================'

               #endif

            case( 'end_direct' )

              #if defined USE_SW_MONO

                  write(6,'(A80)') '================================================================================'
                  write(6,'(A80)') '*  End of run of the Shallow-Water Model in direct mode                        *'
                  write(6,'(A80)') '================================================================================'

               #elif defined USE_NS_MULTIFLUID

                  write(6,'(A80)') '================================================================================'
                  write(6,'(A80)') '*  End of run of the Multifluid Navier-Stokes Model in direct mode             *'
                  write(6,'(A80)') '================================================================================'

               #endif

            case( 'start_testadj' )

               write(6,'(A)')
               write(6,'(A80)') '================================================================================'
               write(6,'(A80)') '*  Testing the generated Adjoint to Shallow-Water Model                        *'
               write(6,'(A80)') '================================================================================'

            case( 'end_testadj' )

               write(6,'(A80)') '================================================================================'
               write(6,'(A80)') '*  End of test of the generated Adjoint to Shallow-Water Model                 *'
               write(6,'(A80)') '================================================================================'

            case( 'start_grad_cost' )

               write(6,'(A)')
               write(6,'(A80)') '================================================================================'
               write(6,'(A80)') '*  Running the Shallow-Water Model Adjoint to calculate a cost-gradient        *'
               write(6,'(A80)') '================================================================================'

            case( 'end_grad_cost' )

               write(6,'(A80)') '================================================================================'
               write(6,'(A80)') '*  End of run of the Shallow-Water Model Adjoint                               *'
               write(6,'(A80)') '================================================================================'

            case( 'start_minimize' )

               write(6,'(A)')
               write(6,'(A80)') '================================================================================'
               write(6,'(A80)') '*  Running the Shallow-Water loop to minimize the cost function                *'
               write(6,'(A80)') '================================================================================'

            case( 'end_minimize' )

               write(6,'(A80)') '================================================================================'
               write(6,'(A80)') '*  End of run of the minimization Shallow-Water Model loop                     *'
               write(6,'(A80)') '================================================================================'

            case( 'result' )

               if ( verbose >= 0 ) then

                  write(6,'(A80)') '================================================================================'
                  write(6,'(A80)') '*  Writing Result File                                                         *'
                  write(6,'(A80)') '================================================================================'

               end if

            case( 'dt' )

               if ( test_dt_just_after( 0.01_rp * ts ) .and. verbose >= 0 ) then

                  write(6,'("nt = ",I8," t = ",ES12.5," / ",ES12.5," ( ",F5.1," % ) , dt = ",ES13.6)') &

                  nt , tc , ts , 100._rp * tc / ts , dt

               end if

            case( 'number_limits' )

               write(6,'(A)')
               write(6,'(A80)') '********************************************************************************'
               write(6,'(A,I15  )') ' range     = ' , range(1._rp)
               write(6,'(A,I15  )') ' precision = ' , precision(1._rp)
               write(6,'(A,ES15.7)') ' zerom     = ' , zerom
               write(6,'(A,ES15.7)') ' pinfm     = ' , pinfm
               write(6,'(A,ES15.7)') ' minfm     = ' , minfm
               write(6,'(A,ES15.7)') ' hugem     = ' , hugem
               write(6,'(A,ES15.7)') ' tinym     = ' , tinym
               write(6,'(A80)') '********************************************************************************'

            case( 'norms_h' )

               write(6,'(A80)') '********************************************************************************'
               write(6,'(A,ES15.7)') ' norm_h_inf  =  ' , norm_inf(1)
               write(6,'(A,ES15.7)') ' norm_h_L1   =  ' , norm_L1 (1)
               write(6,'(A,ES15.7)') ' norm_h_L2   =  ' , norm_L2 (1)
               write(6,'(A80)') '********************************************************************************'

            case( 'norms_u' )

               write(6,'(A80)') '********************************************************************************'
               write(6,'(A,ES15.7)') ' norm_u_inf  =  ' , norm_inf(2)
               write(6,'(A,ES15.7)') ' norm_u_L1   =  ' , norm_L1 (2)
               write(6,'(A,ES15.7)') ' norm_u_L2   =  ' , norm_L2 (2)
               write(6,'(A80)') '********************************************************************************'

            case( 'norms_q' )

               write(6,'(A80)') '********************************************************************************'
               write(6,'(A,ES15.7)') ' norm_q_inf  =  ' , norm_inf(3)
               write(6,'(A,ES15.7)') ' norm_q_L1   =  ' , norm_L1 (3)
               write(6,'(A,ES15.7)') ' norm_q_L2   =  ' , norm_L2 (3)
               write(6,'(A80)') '********************************************************************************'

            case( 'tangent_gradient_test' )

               write(6,'(A)')
               write(6,'(A80)') '================================================================================'
               write(6,'(A80)') '*  Tangent Gradient Test                                                       *'
               write(6,'(A80)') '================================================================================'

            case( 'backward_gradient_test' )

               write(6,'(A)')
               write(6,'(A80)') '================================================================================'
               write(6,'(A80)') '*  Backward Gradient Test                                                      *'
               write(6,'(A80)') '================================================================================'

            case( 'backward_scalar_product_test' )

               write(6,'(A)')
               write(6,'(A80)') '================================================================================'
               write(6,'(A80)') '*  Backward Scalar Product Test                                                *'
               write(6,'(A80)') '================================================================================'

            case( 'cost' )

               write(6,'(A80)') '================================================================================'
               write(6,'(A,ES24.15)') ' cost function       =  ' , var
               write(6,'(A80)') '================================================================================'

            case( 'cost_diff' )

               write(6,'(A80)') '================================================================================'
               write(6,'(A,ES24.15)') ' diff cost function  =  ' , var
               write(6,'(A80)') '================================================================================'

            case( 'cost_back' )

               write(6,'(A80)') '================================================================================'
               write(6,'(A,ES24.15)') ' back cost function  =  ' , var
               write(6,'(A80)') '================================================================================'

         end select

      end if

      call mpi_wait_all

   END SUBROUTINE Print_Screen


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Stopping Program in case of bad parameter or numeric convergence
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief  PRINT WHY Stopping Program in case of bad parameter or numeric convergence
!> \details print error message for different error cases
!!  - unknow_mesh : say unknown mesh
!!  - default  : say stop dassflow run
!! f90wrap_abort is called here to message wrapping error
   SUBROUTINE Stopping_Program( screen_case )

      implicit none

      character(len=*), intent(in)  ::  screen_case

      if ( proc == 0 ) then

         select case( screen_case )

            case( 'unknow_mesh' )

               write(6,*)
               write(6,'(A80)') '================================================================================'
               write(6,'(A80)') '*  Unknow mesh type in input.txt'
               write(6,'(A80)') '================================================================================'

            case default

         end select

      end if

      call mpi_wait_all

      if ( proc == 0 ) then

         write(6,*)
         write(6,'(A)') '================================================================================'
         write(6,'(A)') '*'
         write(6,'(A)') '*  STOPPING DASSFLOW RUN'
         write(6,'(A)') '*'
         write(6,'(A)') '================================================================================'

      end if

		call f90wrap_abort(screen_case)
!~       call End_MPI ; STOP

   END SUBROUTINE Stopping_Program


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Write at screen a progress bar
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief  Write at screen a progress bar
   SUBROUTINE waitbar( iter , niter )

	   implicit none

	   integer(ip)  ::  niter , iter , ch , per

      per = 100 * iter / niter

      call mpi_min_i( per )

      if ( proc == 0  .and. per /= 100 * ( iter - 1 ) / niter ) then

      	write( 6 , '(                  256a1)' , advance='no' ) ( char(8) , ch = 1 , per/2 + 9 )
         write( 6 , '(2x,1i3,1a1,2x,1a1,256a1)' , advance='no' ) per , '%' , '|' , ( '=' , ch = 1,per/2 )
         close( 6 )
         open ( 6 )

      end if

      if ( proc == 0 .and. iter == niter ) write(6,'(a)') '|'

   END SUBROUTINE waitbar

END MODULE m_time_screen
