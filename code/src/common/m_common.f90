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
!> \file m_common.f90
!! \brief This file includes m_common module.
!! \details The file includes only m_common module (see doc m_common module).

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Module m_common
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> Module m_common
!!
!! \details This module includes common variable used in all software (loop indexes, constant, arguments, ...) and some
!! function for the allocation/reallocation of memories.
MODULE m_common

   implicit none

!======================================================================================================================!
!  A set of variables precision (space in memory management)
!======================================================================================================================!

   integer, parameter  ::  ip = 4                              !< Fix the integer numbers machine precision
   integer, parameter  ::  rp = 8                              !<  Fix the real    numbers machine precision

!   integer, parameter  ::  selected_int_kind ( 18       )     ! Fix the integer numbers machine precision
!   integer, parameter  ::  selected_real_kind( 15 , 307 )     ! Fix the real    numbers machine precision

   integer, parameter  ::  lchar = 128                         !< Fix maximun of the string size (for name file for example).

!======================================================================================================================!
!  Local loop index (dangerous to manipulate !)
!======================================================================================================================!

   integer(ip)  ::  i , ie , iK , iKe , ib                  !<  Local loop index  ???
   integer(ip)  ::  j , je , jK , jKe , jb                  !<  Local loop index  ???
   integer(ip)  ::  k , ke , kK , kKe , kb                  !<  Local loop index  ???

!======================================================================================================================!
!  Input variables (to enter in input.txt)
!======================================================================================================================!

   character(len=lchar)  ::  mesh_type                        !< calling cartesian mesh or reader type
   character(len=lchar)  ::  mesh_name                        !< mesh file name if not basic

   character(len=lchar)  ::  bc_N                             !< Type of boundary condition at North mesh boundary if mesh = 'basic'
   character(len=lchar)  ::  bc_S                             !< Type of boundary condition at South mesh boundary if mesh = 'basic'
   character(len=lchar)  ::  bc_W                             !< Type of boundary condition at West  mesh boundary if mesh = 'basic'
   character(len=lchar)  ::  bc_E                             !< Type of boundary condition at East  mesh boundary if mesh = 'basic'

   integer(ip) :: bc_rain                                      ! rain condition
   integer(ip) :: bc_infil                                     ! infiltration condition
   
   real(rp)     ::  lx                                        !< Lenght of computational domain x horizontal direction if mesh = 'basic'
   real(rp)     ::  ly                                        !< Lenght of computational domain y vertical   direction if mesh = 'basic'

   integer(ip)  ::  nx                                        !< Number of nodes in x horizontal direction if mesh = 'basic'
   integer(ip)  ::  ny                                        !< Number of nodes in y vertical   direction if mesh = 'basic'

   real(rp)     ::  ts                                        !< User defined simulation time (total simulation time)

   integer(ip)  ::  adapt_dt                                  !< Activate Adaptative time step
   real(rp)     ::  dt                                        !< Time step if not adaptative
   real(rp)     ::  cfl                                       !< CFL value if adaptative
   
   logical  ::  do_warmup                                      !> Toggle Warmup run for GR4 module

   real(rp)     ::  dtw                                       !<  Time step to Output Result Files
   real(rp)     ::  dtp                                       !<  Time step to Output Post Variables
   real(rp)     ::  dta                                       !<  Time Step to Generate BC (for Data Assimilation)

   integer(ip)  ::  w_tecplot                                 !< write Tecplot Output File
   integer(ip)  ::  w_vtk                                     !< write VTK Output File
   integer(ip)  ::  w_gnuplot                                 !< write Gnuplot Output File
   integer(ip)  ::  w_bin                                     !< write Binary Output File

   integer(ip)  ::  w_exact                                   !< write Exact Solution Output File (calculated from m_user_data)
   integer(ip)  ::  w_norm                                    !< Error Norms Calculation (need exact solution)

   integer(ip)  ::  w_obs                                     !< Gen Observation Output File
   integer(ip)  ::  use_obs                                   !< Use Observations in cost function definition
   integer(ip)  ::  use_Qobs                                   ! Use Observations in cost function definition (model outflow)
   integer(ip)  ::  use_Qobs_gr4                               ! Use Observations in cost function definition (gr4 basins outflow)
   
   character(len=lchar)  ::  spatial_scheme                   !< Name of Spatial  Discretization Scheme ('muscl' only)
   character(len=lchar)  ::  temp_scheme                      !< Name of Temporal Discretization Scheme ('imex' or 'euler')

   character(len=lchar), dimension(:), allocatable  ::  args  !< Arguments passed on the command line

   integer(ip)  ::  max_nt_for_direct                         !< Maximum iterations to perform the direct model
   integer(ip)  ::  max_nt_for_adjoint                        !< Maximum iterations to perform the direct model in view
                                                              !< to bound the memory of the adjoint model

   integer(ip)  ::  length_real                               !< Real record length in a direct access unformatted file

!======================================================================================================================!
!  Working variables
!======================================================================================================================!

   integer(ip)  ::  nt , nt0                                  !< Global Time Iteration Index

   real(rp)  ::  tc , tc0                                     !< Simulation Time

   logical  ::  end_time_loop                                 !< Time Loop Stopping Criterion Logical

   real(rp)  ::  dx                                           !< Space step in x horizontal direction
   real(rp)  ::  dy                                           !< Space step in y vertical   direction

   character(len=lchar)  ::  is_file_open(10000)               !< Helping to create automaticly basic output files

   integer(ip)  ::  file_open_counter = 0                     !< Counter of opened basic output files

   logical  ::  file_exist(10)                                !< Logical to test a file existence

   character(len=1028)  ::  buffer                            !< Temporal buffer to read formated files

   logical  ::  logic_test                                    !< Time Loop Stopping Criterion Logical

   real(rp)  ::  norm_inf(3)                                  !< To store Linf error norm
   real(rp)  ::  norm_L1 (3)                                  !< To store L1   error norm
   real(rp)  ::  norm_L2 (3)                                  !< To store L2   error norm

   integer(ip)  ::  verbose                                   !< Level of verbosity

   integer(ip)  ::  restart_min                               !< Maximum number of gradients calculation called

   real(rp)  ::  eps_min                                      !< Precision stopping criterion for the gradient norm

!======================================================================================================================!
!  Machine precision limits numbers
!======================================================================================================================!

   real(rp)  ::  zerom                                         !< Machine zero  - precision limits numbers
   real(rp)  ::  pinfm , minfm                                 !< Machine infinities - precision limits numbers
   real(rp)  ::  hugem , tinym                                 !< Machine overflow/undeflow - precision limits numbers

!======================================================================================================================!
!  Constants ( preferring multiplying rather than division to accelerate program )
!======================================================================================================================!

   real(rp), parameter  ::  zero   =   0._rp                    !<  Constant : zero
   real(rp), parameter  ::  one    =   1._rp                    !<  Constant : 1
   real(rp), parameter  ::  two    =   2._rp                    !<  Constant : 2
   real(rp), parameter  ::  demi   =   0.5_rp                   !<  Constant : 1/2
   real(rp), parameter  ::  d1p4   =   0.25_rp                  !<  Constant : 1/4

   real(rp), parameter  ::  d1p3   =   1._rp / 3._rp            !<  Constant : 1/3
   real(rp), parameter  ::  d2p3   =   2._rp / 3._rp            !<  Constant : 2/3
   real(rp), parameter  ::  d4p3   =   4._rp / 3._rp            !<  Constant : 4/3
   real(rp), parameter  ::  d5p3   =   5._rp / 3._rp            !<  Constant : 5/3
   real(rp), parameter  ::  d7p3   =   7._rp / 3._rp            !<  Constant : 7/3
   real(rp), parameter  ::  d8p3   =   8._rp / 3._rp            !<  Constant : 8/3
   real(rp), parameter  ::  d10p3  =  10._rp / 3._rp            !<  Constant : 10/3

   real(rp), parameter  ::  d3p2   =   3._rp / 2._rp            !<  Constant : 3/2
   real(rp), parameter  ::  d3p5   =   3._rp / 5._rp            !<  Constant : 3/5
   real(rp), parameter  ::  d3p8   =   3._rp / 8._rp            !<  Constant : 3/8
    
   real(rp), parameter  ::  pi  =  3.14159265358979_rp          !<  Constant : pi


!======================================================================================================================!
!  Structure defining a vector of weights
!======================================================================================================================!

   TYPE weights

      real(rp), dimension(:), allocatable  ::  weights          !< Weight of the observations. -->  IF VAR_chg other treatment ?

   END TYPE

!===================================================================================================================!
!  Recording Structures
!===================================================================================================================!
   
   TYPE station_obsQ

      integer(ip) :: nb
      
      real(rp)  :: weight    ! weight of observations

      integer(ip) :: ind_t   ! Index observation time

      integer(ip) :: nb_dt   ! Number of observation time
      
      real(rp), dimension(:), allocatable  ::  t , q

   END TYPE station_obsQ
      
!===================================================================================================================!
!  Recording Variables
!===================================================================================================================!

   type( station_obsQ ), dimension(:), allocatable ::  stationQ
   
CONTAINS


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Function testing best matching simulation times each given dt
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!
   !> Function testing best matching simulation times each given dt.
   !!
   !! \details If the dt_to_test is the more close to time step of the dt possible, this function return true,
   !! else return false. 
   !! \param[in] dt_to_test Time step to test.
   !! \returns Boolean returned
   logical FUNCTION test_dt_nearest( dt_to_test )

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      real(rp), intent(in)  ::  dt_to_test

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      real(rp)  ::  t_test

      !================================================================================================================!
      !  Begin
      !================================================================================================================!

      if ( dt_to_test < dt ) then

         test_dt_nearest  =  .true.

         return

      end if

      t_test  =  real( floor ( ( tc + 0.5_rp * dt ) / dt_to_test ) , rp ) * dt_to_test

      if ( abs( tc - t_test ) + zerom  <  0.5_rp * dt ) then

         test_dt_nearest  =  .true.

      else

         test_dt_nearest  =  .false.

      end if

   END FUNCTION test_dt_nearest

   !> Function testing best matching simulation times each given dt.
   !!
   !! \details If the dt_to_test is the dt just after the actual dt, this function return true, else return false. 
   !! \param[in] dt_to_test Time step to test.
   !! \returns Boolean returned
   logical FUNCTION test_dt_just_after( dt_to_test )

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      real(rp), intent(in)  ::  dt_to_test

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      real(rp)  ::  t_test

      !================================================================================================================!
      !  Begin
      !================================================================================================================!

      if ( dt_to_test < dt ) then

         test_dt_just_after  =  .true.

         return

      end if

      t_test  =  real( floor ( ( tc + zerom ) / dt_to_test ) , rp ) * dt_to_test

      if ( abs( tc - t_test ) + zerom  <  dt ) then

         test_dt_just_after  =  .true.

      else

         test_dt_just_after  =  .false.

      end if

   END FUNCTION test_dt_just_after


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Reading eventual arguments passed on the command line when the program was invoked
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !> Reading eventual arguments passed on the command line when the program was invoked.
   !> \detail NO ADJ IN THE 1D version 
   SUBROUTINE reading_args

      implicit none

      integer(ip)  ::  iarg , nb_args

      nb_args = command_argument_count()

      if ( nb_args > 0 ) then

         allocate( args( 0 : nb_args ) )

         do iarg = 0,nb_args

            call get_command_argument( number = iarg , VALUE  = args( iarg ) )

         end do

      else

         allocate( args( 1 ) )

         args(1) = ""

      end if

   END SUBROUTINE reading_args


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Increasing or Decreasing Allocated Array in Memory
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   !> Subroutine which allocate or reaallocate space to a variable (real case).
   !!
   !! \details This function allocate more space to an original real variable. If this space is already allocate to this 
   !! variable, space is added. If space is not allocate, the function allocates space to the time for this variable (of
   !! lenght new).
   !! \param[inout] var Original variable.
   !! \param[in] new Size to allocate.
   SUBROUTINE alloc_or_realloc_r( var , new )

      implicit none

      real(rp), dimension(:), allocatable, intent(inout)  ::  var

      integer(ip), intent(in)  ::  new

      if ( .not. allocated( var ) ) then

         allocate( var( 1 : new ) )

      else

         call reallocate_r( var , new )

      end if

   END SUBROUTINE alloc_or_realloc_r

   !> Subroutine which allocate or reaallocate space to a variable (real case).
   !!
   !! \details This function allocate more space to an original real variable. If this space is already allocate to this 
   !! variable, space is added only if new is greater than lenght of the original variable. If space is not allocate, 
   !! the function allocates space to the time for this variable (of lenght new).
   !! \param[inout] var Original variable.
   !! \param[in] new Size to allocate.
   SUBROUTINE alloc_or_larger_r( var , new )

      implicit none

      real(rp), dimension(:), allocatable, intent(inout)  ::  var

      integer(ip), intent(in)  ::  new

      integer(ip)  ::  old

      old = size(var)

      if ( .not. allocated( var ) ) then

         allocate( var( 1 : new ) )

      else if ( new > old ) then

         call reallocate_r( var , new )

      end if

   END SUBROUTINE alloc_or_larger_r

   !> Subroutine which allocate or reaallocate space to a variable (real case).
   !!
   !! \details This function allocate more space to an original real variable. 
   !!    - If the allocate space length is equal to new: do nothing.
   !!    - else if the allocate space length is smaller than the original variable: change structuration of the
   !! original variable (original to the end).
   !!    - else if the allocate space length is smaller than the original variable: change structuration of the
   !! original variable (original to the beginning).
   !! \param[inout] var Original variable.
   !! \param[in] new Size to allocate.
   SUBROUTINE reallocate_r( var , new )

      implicit none

      real(rp), dimension(:), allocatable, intent(inout)  ::  var

      integer(ip), intent(in)  ::  new

      integer(ip)  ::  old

      real(rp), dimension(:), allocatable  ::  temp

      intrinsic move_alloc

      old = size(var)

      if ( new == old ) then

         return

      else if ( new  < old ) then

         allocate( temp( new ) )

         temp( 1 : new ) = var( 1 : new )

         call move_alloc( temp , var )

      else

         allocate( temp( new ) )

         temp( 1 : old ) = var( 1 : old )

         call move_alloc( temp , var )

      end if

      if ( allocated( temp ) ) deallocate( temp )

   END SUBROUTINE reallocate_r

   !> Subroutine which allocate or reaallocate space to a variable (integer case).
   !!
   !! \details This function allocate more space to an original integer variable. If this space is already allocate to this 
   !! variable, space is added. If space is not allocate, the function allocates space to the time for this variable (of
   !! lenght new).
   !! \param[inout] var Original variable.
   !! \param[in] new Size to allocate.
   SUBROUTINE alloc_or_realloc_i( var , new )

      implicit none

      integer(ip), dimension(:), allocatable, intent(inout)  ::  var

      integer(ip), intent(in)  ::  new

      if ( .not. allocated( var ) ) then

         allocate( var( 1 : new ) )

      else

         call reallocate_i( var , new )

      end if

   END SUBROUTINE alloc_or_realloc_i

   !> Subroutine which allocate or reaallocate space to a variable (integer case).
   !!
   !! \details This function allocate more space to an original integer variable. If this space is already allocate to this 
   !! variable, space is added only if new is greater than lenght of the original variable. If space is not allocate, 
   !! the function allocates space to the time for this variable (of lenght new).
   !! \param[inout] var Original variable.
   !! \param[in] new Size to allocate.
   SUBROUTINE alloc_or_larger_i( var , new )

      implicit none

      integer(ip), dimension(:), allocatable, intent(inout)  ::  var

      integer(ip), intent(in)  ::  new

      integer(ip)  ::  old

      old = size(var)

      if ( .not. allocated( var ) ) then

         allocate( var( 1 : new ) )

      else if ( new > old ) then

         call reallocate_i( var , new )

      end if

   END SUBROUTINE alloc_or_larger_i

   !> Subroutine which allocate or reaallocate space to a variable (integer case).
   !!
   !! \details This function allocate more space to an original integer variable. 
   !!    - If the allocate space length is equal to new: do nothing.
   !!    - else if the allocate space length is smaller than the original variable: change structuration of the
   !! original variable (original to the end).
   !!    - else if the allocate space length is smaller than the original variable: change structuration of the
   !! original variable (original to the beginning).
   !! \param[inout] var Original variable.
   !! \param[in] new Size to allocate.
   SUBROUTINE reallocate_i( var , new )

      implicit none

      integer(ip), dimension(:), allocatable, intent(inout)  ::  var

      integer(ip), intent(in)  ::  new

      integer(ip)  ::  old

      integer(ip), dimension(:), allocatable  ::  temp

      intrinsic move_alloc

      old = size(var)

      if      ( new == old ) then

         return

      else if ( new  < old ) then

         allocate( temp( new ) )

         temp( 1 : new ) = var( 1 : new )

         call move_alloc( temp , var )

      else

         allocate( temp( new ) )

         temp( 1 : old ) = var( 1 : old )

         call move_alloc( temp , var )

      end if

      if ( allocated( temp ) ) deallocate( temp )

   END SUBROUTINE reallocate_i


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Deallocation of one variable
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !> Deallocation of one variable (real case).
   !!
   !! \details This function deallocate space of a real variable.
   !! \param[inout] var Original variable.
   SUBROUTINE dealloc_r( var )

      implicit none

      real(rp), dimension(:), allocatable, intent(inout)  ::  var

      if ( allocated( var ) ) deallocate( var )

   END SUBROUTINE dealloc_r


   SUBROUTINE dealloc_i( var )

      implicit none

      integer(ip), dimension(:), allocatable, intent(inout)  ::  var

      if ( allocated( var ) ) deallocate( var )

   END SUBROUTINE dealloc_i


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Machine precision limits calculation
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !> Machine precision limits calculation.
   !!
   !! \details Compute the value of the machine infinities, machine zeo and machine overflow/underflow.
   SUBROUTINE Machine_Number_Limits

      implicit none

      intrinsic huge , tiny

      !================================================================================================================!
      !  Machine zero
      !================================================================================================================!

      zerom  =  one

      do while ( one + zerom > one )

         zerom  =  0.9_rp * zerom

      end do

      !================================================================================================================!
      !  Machine +inf
      !================================================================================================================!

      pinfm  =  one

      do while ( pinfm + one > pinfm )

         pinfm  =  1.1_rp * pinfm

      end do

      !================================================================================================================!
      !  Machine -inf
      !================================================================================================================!

      minfm  =  - one

      do while ( minfm - one < minfm )

         minfm = 1.1_rp * minfm

      end do

      !================================================================================================================!
      !  Greatest Machine Real Value
      !================================================================================================================!

      hugem = huge( 0._rp )

      !================================================================================================================!
      !  Lowest Machine Real Value
      !================================================================================================================!

      tinym = tiny( 0._rp )

   END SUBROUTINE Machine_Number_Limits


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Swapping Subroutines
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !> Swapping Subroutines (real case).
   !!
   !! \details Swap two real variables.
   !! \param[inout] a real variable to swap.
   !! \param[inout] b real variable to swap.
   SUBROUTINE swap_r( a , b )

      implicit none

      real(rp), intent(inout)  ::  a , b

      real(rp)  ::  temp

      temp = a ; a = b ; b = temp

   END SUBROUTINE swap_r

   !> Swapping Subroutines (integer case).
   !!
   !! \details Swap two integer variables.
   !! \param[inout] a integer variable to swap.
   !! \param[inout] b integer variable to swap.
   SUBROUTINE swap_i( a , b )

      implicit none

      integer(ip), intent(inout)  ::  a , b

      integer(ip)  ::  temp

      temp = a ; a = b ; b = temp

   END SUBROUTINE swap_i

   !> Swapping Subroutines (real case).
   !!
   !! \details Swap two real vectors.
   !! \param[inout] vec real vector variable to swap.
   !! \param[inout] swap real vector variable to swap.
   SUBROUTINE swap_vec_r( vec , swap )

      implicit none

      integer(ip), dimension(:), intent(in   )  ::  swap
      real(rp)   , dimension(:), intent(inout)  ::  vec

      real(rp), dimension(:), allocatable :: temp

      allocate( temp( size(vec) ) ) ; temp = vec

      do k = 1,size(swap)

         if ( swap(k) == 0 .or. swap(k) == k ) cycle

         vec( k ) = temp( swap(k) )

      end do

      deallocate( temp )

   END SUBROUTINE swap_vec_r

   !> Swapping Subroutines (integer case).
   !!
   !! \details Swap two integer vectors.
   !! \param[inout] vec integer vector variable to swap.
   !! \param[inout] swap integer vector variable to swap.
   SUBROUTINE swap_vec_i( vec , swap )

      implicit none

      integer(ip), dimension(:), intent(in   )  ::  swap
      integer(ip), dimension(:), intent(inout)  ::  vec

      integer(ip), dimension(:), allocatable  ::  temp

      allocate( temp( size(vec) ) ) ; temp = vec

      do k = 1,size(swap)

         if ( swap(k) == 0 .or. swap(k) == k ) cycle

         vec( k ) = temp( swap(k) )

      end do

      deallocate( temp )

   END SUBROUTINE swap_vec_i


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Basic Functions and Subroutines
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !> Divise function (continuity in 0).
   !!
   !! \details This function does the division between two real, if the denominator is equal to 0 the result is 0.
   !! \param[in] a numerator.
   !! \param[in] b denominator.
   !! \returns results of a/b
	real(rp) FUNCTION div_by_except_0( a , b )

      implicit none

		real(rp), intent(in)  ::  a , b

      if ( abs( b ) > zerom ) then

         div_by_except_0  =  a  /  b

      else

         div_by_except_0  =  0._rp

      end if

	END FUNCTION div_by_except_0


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Ascending sorts - Author : John Burkardt (GNU LGPL license)
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !> Ascending sorts.
   !!
   !! \details In lexicographic order, the statement "X < Y", applied to two real vectors X and Y of length M,
   !! means that there is some index I, with 1 <= I <= M, with the property that X(J) = Y(J) for J < I, and 
   !! X(I) < Y(I).In other words, the first time they differ, X is smaller.
   !! \author John Burkardt (GNU LGPL license)
   !! \param[in] m length of x.
   !! \param[in] n lenght of y.
   !! \param[inout] a maxtrix x,y.
   SUBROUTINE i4col_sort_a( m , n , a )

      implicit none

      integer(ip), intent(in   )  ::  m , n
      integer(ip), intent(inout)  ::  a(m,n)

      integer(ip)  ::  i , j , indx , isgn

      if ( m <= 0 .or. n <= 1 ) return

      !================================================================================================================!
      !  Initialize.
      !================================================================================================================!

      i = 0 ; j = 0 ; indx = 0 ; isgn = 0

      !================================================================================================================!
      !  Call the external heap sorter.
      !================================================================================================================!

      do

         call sort_heap_external( n , indx , i , j , isgn )

         !=============================================================================================================!
         !  Interchange the I and J objects.
         !=============================================================================================================!

         if ( 0 < indx ) then

            call i4col_swap( m , n , a , i , j )

            !=============================================================================================================!
            !  Compare the I and J objects.
            !=============================================================================================================!

         else if ( indx < 0 ) then

            call i4col_compare( m , n , a , i , j , isgn )

         else if ( indx == 0 ) then

            exit

         end if

      end do

   END SUBROUTINE i4col_sort_a


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Externally sorts a list of items into ascending order - Author : John Burkardt (GNU LGPL license)
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integer (ip)s, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   !> Externally sorts a list of items into ascending order.
   !!
   !! \details The actual list of data is not passed to the routine.  Hence this routine may be used to sort integer 
   !! (ip)s, reals, numbers, names, dates, shoe sizes, and so on.  After each call, the routine asks the user to compare
   !! or interchange two items, until a special return value signals that the sorting is completed.
   !! \author John Burkardt (GNU LGPL license)
   !! \param[in]  n integer.
   !! \param[inout] indx integer.
   !! \param[out] i integer.
   !! \param[out] j integer.
   !! \param[in] isgn integer.
   SUBROUTINE sort_heap_external( n , indx , i , j , isgn )

      implicit none

      integer(ip), intent(in   )  ::  isgn , n
      integer(ip), intent(inout)  ::  indx
      integer(ip), intent(  out)  ::  i , j

      integer(ip), save  ::  i_save = 0
      integer(ip), save  ::  j_save = 0
      integer(ip), save  ::  k      = 0
      integer(ip), save  ::  k1     = 0
      integer(ip), save  ::  n1     = 0

      !================================================================================================================!
      !  INDX = 0: This is the first call.
      !================================================================================================================!

      if ( indx == 0 ) then

         i_save  =  0
         j_save  =  0
         k       =  n / 2
         k1      =  k
         n1      =  n

      !================================================================================================================!
      !  INDX < 0: The user is returning the results of a comparison.
      !================================================================================================================!

      else if ( indx < 0 ) then

         if ( indx == -2 ) then

            if ( isgn < 0 ) i_save = i_save + 1

            j_save  =   k1
            k1      =   i_save
            indx    = - 1
            i       =   i_save
            j       =   j_save

            return

         end if

         if ( 0 < isgn ) then

            indx  =  2
            i     =  i_save
            j     =  j_save

            return

         end if

         if ( k <= 1 ) then

            if ( n1 == 1 ) then

               i_save = 0
               j_save = 0
               indx = 0

            else

               i_save = n1
               n1 = n1 - 1
               j_save = 1
               indx = 1

            end if

            i = i_save
            j = j_save

            return

         end if

         k   =  k - 1
         k1  =  k

      !================================================================================================================!
      !  0 < INDX, the user was asked to make an interchange.
      !================================================================================================================!

      else if ( indx == 1 ) then

         k1  =  k

      end if

      do

         i_save = 2 * k1

         if ( i_save == n1 ) then

            j_save  =   k1
            k1      =   i_save
            indx    = - 1
            i       =   i_save
            j       =   j_save

            return

         else if ( i_save <= n1 ) then

            j_save  =   i_save + 1
            indx    = - 2
            i       =   i_save
            j       =   j_save

            return

         end if

         if ( k <= 1 ) exit

         k   =  k - 1
         k1  =  k

      end do

      if ( n1 == 1 ) then

         i_save  =  0
         j_save  =  0
         indx    =  0
         i       =  i_save
         j       =  j_save

      else

         i_save  =  n1
         n1      =  n1 - 1
         j_save  =  1
         indx    =  1
         i       =  i_save
         j       =  j_save

      end if

   END SUBROUTINE sort_heap_external


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Swaps columns I and J - Author : John Burkardt (GNU LGPL license)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !> Swaps columns I and J.
   !!
   !! \author John Burkardt (GNU LGPL license)
   !! \param[in]  m integer.
   !! \param[in] n integer.
   !! \param[inout] a integer.
   !! \param[in] i integer.
   !! \param[in] j integer.
   SUBROUTINE i4col_swap( m , n , a , i , j )

      implicit none

      integer(ip), intent(in   )  ::  m , n , i , j
      integer(ip), intent(inout)  ::  a(m,n)

      integer(ip)  ::  col(m)

      if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) STOP 'Problem with i4col_swap'

      if ( i == j ) return

      col(1:m  )  =  a  (1:m,i)
      a  (1:m,i)  =  a  (1:m,j)
      a  (1:m,j)  =  col(1:m  )

   END SUBROUTINE i4col_swap


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Compares columns I and J - Author : John Burkardt (GNU LGPL license)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !> Compares columns I and J.
   !!
   !! \author John Burkardt (GNU LGPL license)
   !! \param[in]  m integer.
   !! \param[in] n integer.
   !! \param[in] a integer.
   !! \param[in] i integer.
   !! \param[in] j integer.
   !! \param[out] isgn integer.
   SUBROUTINE i4col_compare( m , n , a , i , j , isgn )

      implicit none

      integer(ip), intent(in )  ::   m , n , a(m,n) , i , j
      integer(ip), intent(out)  ::  isgn

      integer(ip)  ::  k

      if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) STOP 'Problem with i4col_compare'

      isgn = 0

      if ( i == j ) return

      k = 1

      do while ( k <= m )

         if      ( a(k,i) < a(k,j) ) then

            isgn = - 1

            return

         else if ( a(k,j) < a(k,i) ) then

            isgn =   1

            return

         end if

         k = k + 1

      end do

   END SUBROUTINE i4col_compare


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Complete a Filename with a extension
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !>  Complete a filename with a extension.
   !!
   !! \details From a file name and a extension ('tecplot', 'gnuplot', ...) complete file name.
   !! \param[in]  file_name Name of the file.
   !! \param[in]  typ Type of the file.
   !! \return The complete file name.
   FUNCTION file_name_ext( file_name, typ ) RESULT( file_name_res )

      character(len=*), intent(in)  ::  typ , file_name

      character(len=lchar)  ::  file_name_res

      select case( trim( typ ) )

         case( 'tecplot' )

            file_name_res = trim(file_name)//'.plt'

         case( 'gnuplot' )

            file_name_res = trim(file_name)//'.dat'

         case( 'bin' )

            file_name_res = trim(file_name)//'.bin'

         case default

            file_name_res = trim(file_name)//'.txt'

      end select

   END FUNCTION


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Count the Number of Lines in a File
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !>  Count the number of lines in a file.
   !!
   !! \param[in]  file_name Name of the file.
   !! \return Number of the lines.
   FUNCTION count_lines( file_name ) RESULT( nb_lines )

      character(len=*), intent(in)  ::  file_name

      integer(ip)  ::  nb_lines , icod

      open(100,file=trim(file_name),status='old',form='formatted')

      nb_lines = 0

      do

         read(100,*,iostat=icod)

         if ( icod >= 0 ) then

            nb_lines = nb_lines + 1

         else

            exit

         end if

      end do

      close(100)

   END FUNCTION


END MODULE m_common
