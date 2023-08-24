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
!  Module using Tapenade generated Output Files in /tap directory
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


MODULE m_adjoint

   USE m_common
   USE m_linear_algebra
   USE m_mesh
   USE m_mpi
   USE m_random
   USE m_time_screen

   USE m_tap_vars
   USE m_obs

#ifdef USE_SW_MONO
      USE m_model
#endif

#ifdef USE_HYDRO
!       USE m_gr4
#endif

   implicit none

   type(unk), target  ::  dof0_diff
   type(unk), target  ::  dof0_back

   type(unk), target  ::  dof_diff
   type(unk), target  ::  dof_back

   real(rp), dimension(:), allocatable  ::  control             !              Control Vector
   real(rp), dimension(:), allocatable  ::  control_back        ! Adjoint      Control Vector
   real(rp), dimension(:), allocatable  ::  control_diff        ! Perturbation Control Vector
   real(rp), dimension(:), allocatable  ::  control_perturb     ! Perturbated  Control Vector
   real(rp), dimension(:), allocatable  ::  control_lbound      !< Lower bounds vector.
   real(rp), dimension(:), allocatable  ::  control_ubound      !< Upper bounds vector.

   real(rp)  ::  cost                                           !              Cost function
   real(rp)  ::  cost_back                                      ! Adjoint      Cost function
   real(rp)  ::  cost_diff                                      ! Perturbation Cost function
   real(rp)  ::  cost_perturb                                   ! Perturbated  Cost function



   integer(ip)  ::  ic                                          ! Control vector index

   integer(ip)  ::  ite_min                                     ! Iteration of Minimization Procedure

   integer(ip)  ::  nb_vars_in_control                          ! Number of variables in the control vector

   integer(ip), dimension(100)  ::  dim_vars_in_control         ! Dimension of each variable in the control vector

   integer(ip) :: dim_all 										! Sum of Dimensions of all variables in the control vector

   type( unk )  ::  dof0_copy

#ifdef USE_SW_MONO

      real(rp), dimension(:), allocatable  ::  bathy_cell_copy

#endif

#ifdef USE_LBFGSB3
   !===================================================================================================================!
   !  Input bounds for LBFGSB3 minimizer
   !===================================================================================================================!

    real(rp)     ::  x1_lbound                      !< Lower bound for GR4 hydrological parameters using LBFGSB-3.0 minimizer
    real(rp)     ::  x1_ubound                      !< Upper bound for GR4 hydrological parameters using LBFGSB-3.0 minimizer
    real(rp)     ::  x2_lbound
    real(rp)     ::  x2_ubound
    real(rp)     ::  x3_lbound
    real(rp)     ::  x3_ubound
    real(rp)     ::  x4_lbound
    real(rp)     ::  x4_ubound

    real(rp)     ::  Ks_lbound
    real(rp)     ::  Ks_ubound
    real(rp)     ::  PsiF_lbound
    real(rp)     ::  PsiF_ubound
    real(rp)     ::  DeltaTheta_lbound
    real(rp)     ::  DeltaTheta_ubound

    integer(ip)  ::  n_bathyb
    integer(ip)  ::  n_manningb

    TYPE spatial_bounds

        real(rp), dimension(:), allocatable ::  manning_lbounds, manning_beta_lbounds
        real(rp), dimension(:), allocatable ::  manning_ubounds, manning_beta_ubounds
        real(rp), dimension(:), allocatable ::  bathy_lbounds
        real(rp), dimension(:), allocatable ::  bathy_ubounds
        real(rp), dimension(1) ::  global_bathy_lbounds, global_bathy_ubounds

    END TYPE spatial_bounds

    type( spatial_bounds ) ::  spatial_b
#endif

CONTAINS


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Compute Direct Model considering :
!     - input  : control vector ( control )
!     - output : cost function  ( cost    )
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE model_direct( mesh , dof0 , dof )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type(msh), intent(in   )  ::  mesh
      type(unk), intent(inout)  ::  dof0
      type(unk), intent(inout)  ::  dof

      !================================================================================================================!
      !  Read control vectror
      !================================================================================================================!

      call read_control( dof0 , mesh )

      !================================================================================================================!
      !  Calling direct model
      !================================================================================================================!

      #ifdef USE_SW_MONO
         call run_model( mesh , dof0 , dof , cost )
      #endif

   END SUBROUTINE model_direct


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Compute Direct Model considering :
!     - input  : perturbated control vector ( control_perturb )
!     - output : perturbated cost function  ( cost_perturb    )
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE model_direct_perturb( mesh , dof0 , dof )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type(msh), intent(in   )  ::  mesh
      type(unk), intent(inout)  ::  dof0
      type(unk), intent(inout)  ::  dof

      !================================================================================================================!
      !  Read perturbated control vector
      !================================================================================================================!

      call read_control_perturb( dof0 , mesh )

      !================================================================================================================!
      !  Calling direct model
      !================================================================================================================!

      call Time_Init(1_ip)

      #ifdef USE_SW_MONO
         call run_model( mesh , dof0 , dof , cost_perturb )
      #endif

      call Time_End(1_ip)

   END SUBROUTINE model_direct_perturb


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Compute Tangent Linear Model considering :
!     - input  : control vector / perturbation control vector ( control / control_diff )
!     - output : cost function  / perturbation cost           ( cost    / cost_diff    )
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE linear_tangent_model( mesh , dof0 , dof )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( msh ), intent(in   )  ::  mesh
      type( unk ), intent(inout)  ::  dof0
      type( unk ), intent(inout)  ::  dof

      !================================================================================================================!
      !  Read control vector
      !================================================================================================================!

      call read_control( dof0 , mesh )

      !================================================================================================================!
      !  Read perturbation control vector
      !================================================================================================================!

      call read_control_diff( dof0_diff , mesh )

      !================================================================================================================!
      !  Calling linear tangent model
      !================================================================================================================!

      call Time_Init(1_ip)

      #ifdef USE_SW_MONO

         call run_model_diff( mesh , dof0 , dof0_diff , dof , dof_diff , cost , cost_diff )

      #endif

      call Time_End(1_ip)

   END SUBROUTINE linear_tangent_model


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Compute Adjoint Model considering :
!     - input  :         control vector / adjoint cost function ( control      / cost_back )
!     - output : adjoint control vector /         cost function ( control_back / cost      )
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


	SUBROUTINE adjoint_model( mesh , dof0 , dof )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type(msh), intent(in   )  ::  mesh
      type(unk), intent(inout)  ::  dof0
      type(unk), intent(inout)  ::  dof

		integer(ip) :: i_loc
      !================================================================================================================!
      !  Read control vector
      !================================================================================================================!

      call read_control( dof0 , mesh )

      !================================================================================================================!
      !  Calling adjoint model
      !================================================================================================================!

      call Time_Init(1_ip)

      if ( proc /= 0 ) cost_back = 0._rp

      #ifdef USE_SW_MONO

         dof_back%h(:) = 0._rp
         dof_back%u(:) = 0._rp
         dof_back%v(:) = 0._rp

         dof_back%infil(:) = 0._rp

         infil_back%GA(:)%Ks           = 0._rp
         infil_back%GA(:)%PsiF         = 0._rp
         infil_back%GA(:)%DeltaTheta   = 0._rp
         infil_back%SCS(:)%lambdacn       = 0._rp
         infil_back%SCS(:)%CN           = 0._rp

         bathy_cell_back(:) = 0._rp

         XSshape_back(1)%xleft = 0._rp
         XSshape_back(1)%xcenter = 0._rp
         XSshape_back(1)%xright = 0._rp
         XSshape_back(1)%s = 0._rp
         XSshape_back(1)%hmax = 0._rp
         XSshape_back(1)%topz = 0._rp

         manning_back(:) = 0._rp
         manning_beta_back(:) = 0._rp

         slope_y_back = 0._rp
         slope_x_back = 0._rp

         call run_model_back( mesh , dof0 , dof0_back , dof , dof_back , cost , cost_back )

      #endif

      call Time_End(1_ip)

     !================================================================================================================!
      !  Filling control_back vector (cost gradient vector)
      !================================================================================================================!

      call write_control_back( dof0_back , mesh )
      
   END SUBROUTINE adjoint_model


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Calc cost function gradient using Adjoint Model
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE calc_grad_cost( mesh , dof0 , dof )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( msh ), intent(in   )  ::  mesh
      type( unk ), intent(inout)  ::  dof0
      type( unk ), intent(inout)  ::  dof

      !================================================================================================================!
      !  Initialization
      !================================================================================================================!
        !call dealloc_back_vars()   ! might become usefull to call it for python, hope to user
		call alloc_back_vars(dof0_back, dof_back, mesh )

      call write_control( dof0 , mesh )

      verbose  =  -1

      cost_back  =  1._rp

      !================================================================================================================!
      !  Calc cost function gradient using Adjoint Model
      !================================================================================================================!

      call adjoint_model( mesh , dof0 , dof )

      !================================================================================================================!
      !  Output in proper Files the Gradient of the Cost Function
      !================================================================================================================!

      call output_control_back( mesh )

   END SUBROUTINE calc_grad_cost


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Testing Adjoint :
!
!        0 : all
!        1 : two direct run
!        2 : gradient test on tangent linear
!        3 : gradient test on adjoint
!        4 : scalar product
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE test_adjoint( mesh , dof0 , dof , arg )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( msh ), intent(in   )  ::  mesh
      type( unk ), intent(inout)  ::  dof0
      type( unk ), intent(inout)  ::  dof

      character(len=*), intent(in)  ::  arg

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      integer(ip)  ::  test_type

      real(rp)  ::  eps , scal1 , scal2

      !================================================================================================================!
      !  Initialization
      !================================================================================================================!

      call alloc_diff_vars( dof0_diff , dof_diff , mesh )
      call alloc_back_vars(dof0_back, dof_back, mesh )

      !================================================================================================================!
      !  Gradient Direction
      !================================================================================================================!

      call write_control     ( dof0 , mesh )
      call write_control_diff( dof0 , mesh )

      allocate( control_perturb ( size( control ) ) )
      allocate( control_back    ( size( control ) ) )

      !================================================================================================================!
      !  Direct Test (2 run that should return the same p_Y result
      !================================================================================================================!

      read(arg,*) test_type

      verbose  =  -1

      if ( test_type == 0 .or. test_type == 1 ) then

         if ( proc == 0 ) write(6,*)

         call model_direct( mesh , dof0 , dof )

         call Print_Screen( 'cost' , cost )

         call model_direct( mesh , dof0 , dof )

         call Print_Screen( 'cost' , cost )

         if ( test_type == 1 ) call End_MPI

      end if

      !================================================================================================================!
      !  Tangent Gradient Test
      !================================================================================================================!

      if ( test_type == 0 .or. test_type == 2 ) then

         call Print_Screen( 'tangent_gradient_test' )

         if ( proc == 0 ) open(80,file='time_step_serie_for_adjoint.bin',form='unformatted',status='replace')

         fix_time_step_serie = 1

         call linear_tangent_model( mesh , dof0 , dof )

         if ( proc == 0 ) then

            close(80)

            write(6,*)

         end if

         call Print_Screen( 'cost'      , cost      )
         call Print_Screen( 'cost_diff' , cost_diff )

         eps = 1._rp

         if ( proc == 0 ) then

            write(6,*)
            write(6,'(A)') '--------------------------------------------------------------------------------'
            write(6,'(A,TR21,A,TR6,A)') ' eps' , 'cost_perturb_diff' , 'relative error'
            write(6,'(A)') '--------------------------------------------------------------------------------'
            write(6,*)

         end if

         do while( eps >= 1.d-8 )

            control_perturb(:)  =  control(:)  +  eps * control_diff(:)

            open(80,file='time_step_serie_for_adjoint.bin',form='unformatted',status='old')

            fix_time_step_serie = 2

            call model_direct_perturb( mesh , dof0 , dof )

            if ( proc == 0 ) write(6,'(3ES23.15)') eps , &
                                                   ( cost_perturb - cost ) / eps, &
                                                   abs( one -  ( cost_perturb - cost ) / ( eps * cost_diff ) )

            close(80)

            eps = 0.1_rp * eps

         end do

         if ( test_type == 2 ) call End_MPI

      end if

      !================================================================================================================!
      !  Backward Gradient Test
      !================================================================================================================!

      if ( test_type == 0 .or. test_type == 3 ) then

         call Print_Screen( 'backward_gradient_test' )

         if ( proc == 0 ) open(80,file='time_step_serie_for_adjoint.bin',form='unformatted',status='replace')

         fix_time_step_serie = 1

         cost_back  =  one

         call adjoint_model( mesh , dof0 , dof )

         cost_diff  =  sum( control_diff(:) * control_back(:) )

         if ( proc == 0 ) then

            close(80)

            write(6,*)

         end if

         call Print_Screen( 'cost'      , cost      )
         call Print_Screen( 'cost_diff' , cost_diff )

         eps = 1._rp

         if ( proc == 0 ) then

            write(6,*)
            write(6,'(A)') '--------------------------------------------------------------------------------'
            write(6,'(A,TR21,A,TR6,A)') ' eps' , 'cost_perturb_diff' , 'relative error'
            write(6,'(A)') '--------------------------------------------------------------------------------'
            write(6,*)

         end if

         do while( eps >= 1.d-8 )

            control_perturb(:)  =  control(:)  +  eps * control_diff(:)

            open(80,file='time_step_serie_for_adjoint.bin',form='unformatted',status='old')

            fix_time_step_serie = 2

            call model_direct_perturb( mesh , dof0 , dof )

            if ( proc == 0 ) write(6,'(3ES23.15)') eps , &
                                                   ( cost_perturb - cost ) / eps, &
                                                   abs( one -  ( cost_perturb - cost ) / ( eps * cost_diff ) )

            close(80)

            eps = 0.1_rp * eps

         end do

         if ( test_type == 3 ) call End_MPI

      end if

      !================================================================================================================!
      !  Backward Scalar Product Test
      !================================================================================================================!

      if ( test_type == 0 .or. test_type == 4 ) then

         fix_time_step_serie = 0

         call Print_Screen( 'backward_scalar_product_test' )

         call linear_tangent_model( mesh , dof0 , dof )

         call Print_Screen( 'cost'      , cost      )
         call Print_Screen( 'cost_diff' , cost_diff )

         cost_back  =  cost_diff

         scal1  =  cost_back * cost_diff

         call adjoint_model( mesh , dof0 , dof )

         call Print_Screen( 'cost'      , cost      )
         call Print_Screen( 'cost_back' , cost_back )

         scal2 = sum( control_diff(:) * control_back(:) )

         if ( proc == 0 ) write(6,'(A,2ES22.15)') ' scal1 = ' , scal1
         if ( proc == 0 ) write(6,'(A,2ES22.15)') ' scal2 = ' , scal2

         if ( proc == 0 ) write(6,'(A,2ES22.15)') ' relative error = '  , abs( ( scal2 - scal1 ) / scal1 )

         if ( test_type == 4 ) call End_MPI

      end if

   END SUBROUTINE test_adjoint


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Filling the Model Input Control Vector
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE write_control( dof0 , mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type(msh), intent(in   )  ::  mesh
      type(unk), intent(inout)  ::  dof0
      real(rp) :: x1min,x1max,x2min,x2max,x3min,x3max,x4min,x4max
      !================================================================================================================!
      !  Filling desired Control Vector (with c_ in input.txt)
      !================================================================================================================!

      nb_vars_in_control = 0

      ic = 1

      #ifdef USE_SW_MONO

         if ( c_shape_s == 1 ) then
            call var_2_control( XSshape(:)%s    , size(XSshape)   , 0 )
        endif
        if ( c_hmax == 1 ) then
            call var_2_control( XSshape(:)%hmax    , size(XSshape)   , 0 )
        endif
        if ( c_xcenter == 1 ) then
            call var_2_control( XSshape(:)%xcenter    , size(XSshape)   , 0 )
         endif

         if ( c_manning == 1 ) call var_2_control( manning    , nland   , manning_data_glob )
         if ( c_manning_beta == 1 ) call var_2_control( manning_beta    , nland   , manning_data_glob )
         if ( c_bathy   == 1 ) call var_2_control( bathy_cell, mesh%nc , 0                 )

         if ( c_slope_y == 1 ) call var_2_control( slope_y , size(slope_y) , 0                 )
         if ( c_slope_x == 1 ) call var_2_control( slope_x , size(slope_x) , 0                 )

         if ( c_ic      == 1 ) call var_2_control( dof0%h     , mesh%nc , 0                 )
         if ( c_ic      == 1 ) call var_2_control( dof0%u     , mesh%nc , 0                 )
         if ( c_ic      == 1 ) call var_2_control( dof0%v     , mesh%nc , 0                 )

#ifdef USE_HYDRO
        if ( c_gr4params == 1 ) then

          x1min = 225.0_rp
          x1max = 625.0_rp
          x2min = -20.0_rp
          x2max = -3_rp
          x3min = 60.0_rp
          x3max = 160.0_rp
          x4min = 0.1_rp
          x4max = 0.2_rp

      do k = 1,bc%nb_gr4in

       write(*,*) ''
       write(*,*) 'old gr4 params'
       write(*,*) bc%gr4( k )%params(:)

           !Normalize here for now
            bc%gr4( k )%params(1) = (bc%gr4( k )%params(1)   - x1min) / (x1max - x1min)
            bc%gr4( k )%params(2) = (bc%gr4( k )%params(2)   - x2min) / (x2max - x2min)
            bc%gr4( k )%params(3) = (bc%gr4( k )%params(3)   - x3min) / (x3max - x3min)
            bc%gr4( k )%params(4) = (bc%gr4( k )%params(4)   - x4min) / (x4max - x4min)

       write(*,*) ''
       write(*,*) 'old gr4 params (normalized)'
       write(*,*) bc%gr4( k )%params(:)


               call var_2_control( bc%gr4( k )%params(:) , size(bc%gr4( k )%params(:)) , 1 )

          end do

         endif
#endif



         if ( c_Ks == 1         ) call var_2_control( infil%GA(:)%Ks , infil%nland , 0 )
         if ( c_PsiF == 1       ) call var_2_control( infil%GA(:)%PsiF , infil%nland , 0 )
         if ( c_DeltaTheta == 1 ) call var_2_control( infil%GA(:)%DeltaTheta, infil%nland , 0 )
         if ( c_lambda == 1     ) call var_2_control( infil%SCS(:)%lambdacn , infil%nland , 0 )
         if ( c_CN == 1         ) call var_2_control( infil%SCS(:)%CN, infil%nland , 0 )


         if ( c_hydrograph == 1 ) then
            do k = 1,bc%nb_in
               call var_2_control( bc%hyd( k )%q(:) , size( bc%hyd( k )%q(:) ) , 1 )
            end do
         end if

         if      ( c_ratcurve == 1 ) then
            do k = 1,bc%nb_out
               call var_2_control( bc%rat( k )%q(:) , size( bc%rat( k )%q(:) ) , 1 )
            end do
         else if ( c_ratcurve == 2 ) then

            do k = 1,bc%nb_out
               call var_2_control( bc%rat( k )%pow(:) , 2 , 1 )
            end do

         end if

         if ( c_rain == 1 ) then
            do k = 1,bc%nb_rn
               call var_2_control( bc%rain( k )%q(:) , size( bc%rain( k )%q(:) ) , 1 )
            end do
         end if

         if(allocated(bathy_cell_copy)) deallocate(bathy_cell_copy)
         allocate( bathy_cell_copy( size( bathy_cell ) ) )

         bathy_cell_copy(:)  =  bathy_cell(:)

         call alloc_dof( dof0_copy , mesh )

         dof0_copy  =  dof0

      #endif

   CONTAINS


      SUBROUTINE var_2_control( var , n , data_glob )

         implicit none

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         integer(ip), intent(in)  ::  n , data_glob

         real(rp), dimension(n), intent(in)  ::  var

         !=============================================================================================================!
         !
         !=============================================================================================================!

         integer(ip)  ::  n_k

         integer(ip)  :: k_loc

         !=============================================================================================================!
         !
         !=============================================================================================================!

         nb_vars_in_control  =  nb_vars_in_control  +  1



         do k_loc = 0,np-1

            if ( k_loc == 0 ) then

               call alloc_or_larger_r( control , ic + n - 1 )

               control( ic : ic + n - 1 )  =  var( 1 : n )

               ic = ic + n

               dim_vars_in_control( nb_vars_in_control )  =  n

            else if ( data_glob == 0 ) then

               call mpi_send_recv_scal_i( n , n_k , k_loc , 0 )

               if ( proc == 0 ) call alloc_or_larger_r( control , ic + n_k - 1 )

               call mpi_send_recv_array_r( var( 1 : n ) , control( ic : ic + n_k - 1 ) , k_loc , 0 )

               if ( proc == 0 ) ic = ic + n_k

               if ( proc == 0 ) dim_vars_in_control( nb_vars_in_control ) = &
                                dim_vars_in_control( nb_vars_in_control ) + n_k

            end if

            call mpi_wait_all

         end do

      END SUBROUTINE var_2_control


   END SUBROUTINE write_control


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Filling the Model Input Perturbation Control Vector
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE write_control_diff( dof0 , mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( unk ), intent(in)  ::  dof0
      type( msh ), intent(in)  ::  mesh

      !================================================================================================================!
      !  Filling input control vector ( control )
      !================================================================================================================!
		write(*,*) "in write_control_diff( dof0 , mesh )"
      ic = 1

      #ifdef USE_SW_MONO

        if ( c_shape_s == 1 ) then
            call var_2_control_diff( eps_bathy * XSshape(:)%s    , size(XSshape)   , 0 )
        endif
        if ( c_hmax == 1 ) then
            call var_2_control_diff( eps_bathy * XSshape(:)%hmax    , size(XSshape)   , 0 )
        endif
        if ( c_xcenter == 1 ) then
            call var_2_control_diff( eps_bathy * XSshape(:)%xcenter    , size(XSshape)   , 0 )
         endif

         if ( c_manning == 1 ) call var_2_control_diff( eps_manning * manning    , nland   , manning_data_glob )
         if ( c_manning_beta == 1 ) call var_2_control_diff( eps_manning * manning_beta    , nland   , manning_data_glob ) ! should have manning_BETA_data_glob ??
         if ( c_bathy   == 1 ) call var_2_control_diff( eps_bathy   * bathy_cell , mesh%nc , 0                 )

         if ( c_slope_y == 1 ) call var_2_control_diff( slope_y , size(slope_y) , 0                 )
         if ( c_slope_x == 1 ) call var_2_control_diff( slope_x , size(slope_x) , 0                 )

         if ( c_ic      == 1 ) call var_2_control_diff( eps_ic      * dof0%h     , mesh%nc , 0                 )
         if ( c_ic      == 1 ) call var_2_control_diff( eps_ic      * dof0%u     , mesh%nc , 0                 )
         if ( c_ic      == 1 ) call var_2_control_diff( eps_ic      * dof0%v     , mesh%nc , 0                 )

#ifdef USE_HYDRO
         if ( c_gr4params == 1 ) then

          do k = 1,bc%nb_gr4in

               call var_2_control_diff( eps_gr4params * bc%gr4( k )%params(:) , size(bc%gr4( k )%params(:)) , 1 )

          end do

         endif
#endif


         if ( c_Ks == 1 ) call var_2_control_diff( eps_Ks * infil%GA(:)%Ks , infil%nland , 0 )
         if ( c_PsiF == 1 ) call var_2_control_diff( eps_PsiF * infil%GA(:)%PsiF , infil%nland , 0 )
         if ( c_DeltaTheta == 1 ) call var_2_control_diff( eps_DeltaTheta * infil%GA(:)%DeltaTheta, infil%nland , 0 )
         if ( c_lambda == 1 ) call var_2_control_diff( infil%SCS(:)%lambdacn , infil%nland , 0 )
         if ( c_CN == 1 ) call var_2_control_diff(infil%SCS(:)%CN, infil%nland , 0 )

         if ( c_hydrograph == 1 ) then
            do k = 1,bc%nb_in
               call var_2_control_diff( eps_hydrograph * bc%hyd( k )%q(:) , size( bc%hyd( k )%q(:) ) , 1 )
            end do
         end if

         if      ( c_ratcurve == 1 ) then
            do k = 1,bc%nb_out
               call var_2_control_diff( eps_ratcurve * bc%rat( k )%q(:) , size( bc%rat( k )%q(:) ) , 1 )
            end do

         else if ( c_ratcurve == 2 ) then
            do k = 1,bc%nb_out
               call var_2_control_diff( eps_ratcurve * bc%rat( k )%pow(:) , 2 , 1 )
            end do

         end if

         if ( c_rain == 1 ) then
            do k = 1,bc%nb_rn
               call var_2_control_diff( eps_rain * bc%rain( k )%q(:) , size( bc%rain( k )%q(:) ) , 1 )
            end do
         end if

      #endif


   CONTAINS


      SUBROUTINE var_2_control_diff( var , n , data_glob )

         implicit none

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         integer(ip), intent(in)  ::  n , data_glob

         real(rp), dimension(n), intent(in)  ::  var

         !=============================================================================================================!
         !
         !=============================================================================================================!

         integer(ip)  ::  n_k,k_loc

         real(rp), dimension(n)  ::  rn

         !=============================================================================================================!
         !
         !=============================================================================================================!

         call init_random_seed

         do k_loc = 0,np-1

            if ( k_loc == 0 ) then

               call alloc_or_larger_r( control_diff , ic + n - 1 )

               call random_number( rn )

               control_diff( ic : ic + n - 1 )  =  ( -1._rp + 2._rp * rn ( 1 : n ) ) * var( 1 : n )

               ic = ic + n

            else if ( data_glob == 0 ) then

               call mpi_send_recv_scal_i( n , n_k , k_loc , 0 )

               if ( proc == 0 ) call alloc_or_larger_r( control_diff , ic + n_k - 1 )

               call random_number( rn )

               call mpi_send_recv_array_r( ( -1._rp + 2._rp * rn ( 1 : n ) ) * var( 1 : n ) , &
                                           control_diff( ic : ic + n_k - 1 ) , k_loc , 0 )

               if ( proc == 0 ) ic = ic + n_k

            end if

            call mpi_wait_all

         end do

      END SUBROUTINE var_2_control_diff


   END SUBROUTINE write_control_diff


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Filling the Model Output Adjoint Control Vector, the first called, control_back is allocated little by little
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE write_control_back( dof0 , mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type(unk), intent(in)  ::  dof0
      type(msh), intent(in)  ::  mesh

      !================================================================================================================!
      !  Filling input control vector ( control )
      !================================================================================================================!

      ic = 1

      #ifdef USE_SW_MONO

!          if ( manning_data_glob == 1 ) then

            do i = 1,nland
               call mpi_sum_r( manning_back(i) )
               call mpi_sum_r( manning_beta_back(i) )
            end do

!          end if

#ifdef USE_HYDRO
         do k = 1,size(bc_back%gr4)

            do i = 1,4
               call mpi_sum_r( bc_back%gr4( k )%params(i) )
            end do

         end do
#endif

         do k = 1,bc%nb_in

            do i = 1,size( bc_back%hyd( k )%q(:) )
!             write(*,*) "k, i , bc_back%hyd( k )%q(:)", k, i , bc_back%hyd( k )%q(:)
               call mpi_sum_r( bc_back%hyd( k )%q(i) )
            end do

         end do

         do k = 1,bc%nb_out

            do i = 1,size( bc_back%rat( k )%q(:) )
               call mpi_sum_r( bc_back%rat( k )%q(i) )
            end do

            do i = 1,size( bc_back%rat( k )%pow(:) )
               call mpi_sum_r( bc_back%rat( k )%pow(i) )
            end do

         end do

         do k = 1,bc%nb_rn

            do i = 1,size( bc_back%rain( k )%q(:) )
               call mpi_sum_r( bc_back%rain( k )%q(i) )
            end do

         end do


        if ( c_shape_s == 1 ) then
            call var_2_control_back( XSshape_back(:)%s    ,  size(XSshape_back)   , 0 )
        endif
        if ( c_hmax == 1 ) then
            call var_2_control_back( XSshape_back(:)%hmax    ,  size(XSshape_back)   , 0 )
        endif
        if ( c_xcenter == 1 ) then
            call var_2_control_back( XSshape_back(:)%xcenter    ,  size(XSshape_back)   , 0 )
         endif

         if ( c_manning == 1 ) call var_2_control_back( manning_back    , nland   , manning_data_glob )
         if ( c_manning_beta == 1 ) call var_2_control_back( manning_beta_back    , nland   , manning_data_glob )
         if ( c_bathy   == 1 ) call var_2_control_back( bathy_cell_back, mesh%nc , 0                 )

         if ( c_slope_y == 1 ) call var_2_control_back( slope_y , size(slope_y) , 0                 )
         if ( c_slope_x == 1 ) call var_2_control_back( slope_x , size(slope_x) , 0                 )

         if ( c_ic      == 1 ) call var_2_control_back( dof0_back%h     , mesh%nc , 0                 )
         if ( c_ic      == 1 ) call var_2_control_back( dof0_back%u     , mesh%nc , 0                 )
         if ( c_ic      == 1 ) call var_2_control_back( dof0_back%v     , mesh%nc , 0                 )

#ifdef USE_HYDRO
         if ( c_gr4params == 1 ) then

          do k = 1,bc%nb_gr4in

               call var_2_control_back( bc_back%gr4( k )%params(:) , size(bc_back%gr4( k )%params(:)) , 1 )

          end do

         endif
#endif

         if ( c_Ks              == 1 ) call var_2_control_back( infil_back%GA%Ks   , infil%nland      , manning_data_glob )
         if ( c_PsiF            == 1 ) call var_2_control_back( infil_back%GA%PsiF , infil%nland      , 0                 )
         if ( c_DeltaTheta      == 1 ) call var_2_control_back( infil_back%GA%DeltaTheta   , infil%nland , 0      )
         if ( c_lambda          == 1 ) call var_2_control_back( infil_back%SCS%lambdacn, infil%nland    , 0                 )
         if ( c_CN              == 1 ) call var_2_control_back( infil_back%SCS%CN    , infil%nland    , 0                 )

         if ( c_hydrograph == 1 ) then
            do k = 1,bc%nb_in
               call var_2_control_back( bc_back%hyd( k )%q(:) , size( bc_back%hyd( k )%q(:) ) , 1 )
            end do
         end if

         if ( c_ratcurve == 1 ) then
            do k = 1,bc%nb_out
               call var_2_control_back( bc_back%rat( k )%q(:) , size( bc_back%rat( k )%q(:) ) , 1 )
            end do

         else if ( c_ratcurve == 2 ) then
            do k = 1,bc%nb_out
               call var_2_control_back( bc_back%rat( k )%pow(:) , 2 , 1 )
            end do

         end if

         if ( c_rain == 1 ) then
            do k = 1,bc%nb_rn
               call var_2_control_back( bc_back%rain( k )%q(:) , size( bc_back%rain( k )%q(:) ) , 1 )
            end do
         end if

      #endif


   CONTAINS


      SUBROUTINE var_2_control_back( var , n , data_glob )

         implicit none

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         integer(ip), intent(in)  ::  n , data_glob

         real(rp), dimension(n), intent(in)  ::  var

         !=============================================================================================================!
         !
         !=============================================================================================================!

         integer(ip)  ::  n_k=0,k_loc

         !=============================================================================================================!
         !
         !=============================================================================================================!

         do k_loc = 0,np+1

            if ( k_loc == 0 ) then

               call alloc_or_larger_r( control_back , ic + n - 1 )

               control_back( ic : ic + n - 1 )  =  var( 1 : n )

               ic = ic + n

            else if ( data_glob == 0 ) then

               call mpi_send_recv_scal_i( n , n_k , k_loc , 0 )

               if ( proc == 0 ) call alloc_or_larger_r( control_back , ic + n_k - 1 )

               call mpi_send_recv_array_r( var( 1 : n ) , control_back( ic : ic + n_k - 1 ) , k_loc , 0 )

               if ( proc == 0 ) ic = ic + n_k

            end if

            call mpi_wait_all

         end do

      END SUBROUTINE var_2_control_back


   END SUBROUTINE write_control_back



!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Filling the Model Boundaries Control Vector
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !>  Filling the model input adjoint control vector.
   !!
   !! \details Add to the control vector k variable. In the end of this function the vector k is created.
   !! \param[in]    mesh Mesh of the model.
   !! \param[inout] dof0 Initial conditions.
   SUBROUTINE write_control_bounds( dof0 , mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( unk ), intent(in)  ::  dof0
      type( msh ), intent(in)  ::  mesh

      !================================================================================================================!
      !  Filling input control vector ( control )
      !================================================================================================================!

      ic = 1

#ifdef USE_LBFGSB3

x1_lbound = 0_rp
x1_ubound = 1_rp
x2_lbound = 0_rp
x2_ubound = 1_rp
x3_lbound = 0_rp
x3_ubound = 1_rp
x4_lbound = 0_rp
x4_ubound = 1_rp


      if ( c_gr4params == 1 ) then
        do k = 1,size(bc%gr4)
         call var_2_control_scal_bounds( x1_lbound, x1_ubound, 1 )
         call var_2_control_scal_bounds( x2_lbound, x2_ubound, 1 )
         call var_2_control_scal_bounds( x3_lbound, x3_ubound, 1 )
         call var_2_control_scal_bounds( x4_lbound, x4_ubound, 1 )
        enddo
      endif

#endif


   CONTAINS

      !>  Filling the model input control bounds vectors.
      !!
      !! \details This subroutine adds to the control bounds vectors the bounds for a variable array of dimension n.
      !! \param[in] lb Lower bound to append to the control lower-bounds vectors.
      !! \param[in] ub Upper bound to append to the control upper-bounds vectors.
      !! \param[in] n Number of components in the control vector.
      !! \param[in] data_glob Variable not used.
      SUBROUTINE var_2_control_bounds( lb, ub, n, data_glob )

         implicit none

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         integer(ip), intent(in)  ::  n , data_glob

         real(rp), intent(in)  ::  lb
         real(rp), intent(in)  ::  ub

         call alloc_or_larger_r( control_lbound , ic + n - 1 )
         call alloc_or_larger_r( control_ubound , ic + n - 1 )

         control_lbound( ic : ic + n - 1 )  =  lb
         control_ubound( ic : ic + n - 1 )  =  ub

         ic = ic + n


      END SUBROUTINE var_2_control_bounds


      !BELOW: 1D routine, not adapted for 2D

      !>  Filling the model input control bounds vectors.
      !!
      !! \details This subroutine adds to the control bounds vectors the bounds for the bathymetry.
      !! \param[in] lb Lower bound to append to the control lower-bounds vectors.
      !! \param[in] ub Upper bound to append to the control upper-bounds vectors.
      !! \param[in] var Array of initial values of bathymetry.
      !! \param[in] delta Admissible delta around initial values of bathymetry.
      !! \param[in] n Size of the control vector.
      !! \param[in] data_glob Variable not used.
!       SUBROUTINE bathy_2_control_bounds( lb, ub, var, delta, n, data_glob )
!
!          implicit none
!
!          !=============================================================================================================!
!          !  Interface Variables
!          !=============================================================================================================!
!
!          integer(ip), intent(in)  ::  n , data_glob
!
!          real(rp), intent(in)  ::  lb
!          real(rp), intent(in)  ::  ub
!          real(rp), intent(in)  ::  delta
!          real(rp), dimension(n), intent(in)  ::  var
!
!
!          call alloc_or_larger_r( control_lbound , ic + nb_bathy_control_pts - 1 )
!          call alloc_or_larger_r( control_ubound , ic + nb_bathy_control_pts - 1 )
!          if (delta > 1e-6) then
!             control_lbound(ic:ic+nb_bathy_control_pts-1)  =  var(bathy_first:bathy_last) - 0.5 * delta
!             control_ubound(ic:ic+nb_bathy_control_pts-1)  =  var(bathy_first:bathy_last) + 0.5 * delta
!          else
!             control_lbound(ic:ic+nb_bathy_control_pts-1)  =  lb
!             control_lbound(ic:ic+nb_bathy_control_pts-1)  =  ub
!          end if
!          ic = ic + nb_bathy_control_pts
!
!       END SUBROUTINE bathy_2_control_bounds


      !>  Append boundaries for a scalar component of the control vector to the array of control boundaries
      !!
      !! \param[in] lb Lower bound for the scalar component of the control vector
      !! \param[in] ub Upper bound for the scalar component of the control vector
      !! \param[in] data_glob Unused variable.
      SUBROUTINE var_2_control_scal_bounds( lb, ub, data_glob )

         implicit none

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         integer(ip), intent(in)  ::  data_glob

         real(rp), intent(in)  ::  lb
         real(rp), intent(in)  ::  ub

         call alloc_or_larger_r( control_lbound , ic )
         call alloc_or_larger_r( control_ubound , ic )

         control_lbound( ic  )  =  lb
         control_ubound( ic  )  =  ub

         ic = ic + 1

      END SUBROUTINE var_2_control_scal_bounds


   END SUBROUTINE write_control_bounds

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Reading the Model Input Control Vector
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE read_control( dof0 , mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type(unk), intent(inout)  ::  dof0
      type(msh), intent(in   )  ::  mesh

      character(len=lchar)  ::  file_name
      real(rp) :: x1min,x1max,x2min,x2max,x3min,x3max,x4min,x4max
      real(rp) :: Ksmin, Ksmax, PsiFmin, PsiFmax, DeltaThetamin, DeltaThetamax
      !================================================================================================================!
      !  Filling input control vector ( control )
      !================================================================================================================!

      ic = 1

      #ifdef USE_SW_MONO

         bathy_cell(:) = bathy_cell_copy(:)

        if ( c_shape_s == 1 ) then
            call control_2_var( XSshape(:)%s    , size(XSshape)   , 0 )
        endif
        if ( c_hmax == 1 ) then
            call control_2_var( XSshape(:)%hmax    , size(XSshape)   , 0 )
        endif
        if ( c_xcenter == 1 ) then
            call control_2_var( XSshape(:)%xcenter    , size(XSshape)   , 0 )
        endif

         if ( c_manning == 1 ) call control_2_var( manning    , nland   , manning_data_glob )
         if ( c_manning_beta == 1 ) call control_2_var( manning_beta    , nland   , manning_data_glob )
         if ( c_bathy   == 1 ) call control_2_var( bathy_cell, mesh%nc , 0                 )

         if ( c_slope_y == 1 ) call control_2_var( slope_y , size(slope_y) , 0                 )
         if ( c_slope_x == 1 ) call control_2_var( slope_x , size(slope_x) , 0                 )

         do ie = 1,mesh%neb

            i = mesh%edge( mesh%edgeb(ie)%ind )%cell(2)
            j = mesh%edge( mesh%edgeb(ie)%ind )%cell(1)

            if ( mesh%edgeb(ie)%typlim == 'wall' ) bathy_cell(i)  =  bathy_cell(j)

         end do

         if ( c_ic      == 1 ) call control_2_var( dof0%h     , mesh%nc , 0 )
         if ( c_ic      == 1 ) call control_2_var( dof0%u     , mesh%nc , 0 )
         if ( c_ic      == 1 ) call control_2_var( dof0%v     , mesh%nc , 0 )

#ifdef USE_HYDRO
         if ( c_gr4params == 1 ) then

          x1min = 225.0_rp
          x1max = 625.0_rp
          x2min = -20.0_rp
          x2max = -3_rp
          x3min = 60.0_rp
          x3max = 160.0_rp
          x4min = 0.1_rp
          x4max = 0.2_rp

          do k = 1,bc%nb_gr4in

               call control_2_var( bc%gr4( k )%params(:) , size(bc%gr4( k )%params(:)) , 1 )


               write(file_name,'(A,I3.3,A,I3.3)') 'min/gr4ite_' , k , '-' , k

                inquire( file = file_name , exist = file_exist(1) )


               if ( .not. file_exist(1) ) then
                open(10,file=file_name,position='append',form='formatted')
               else
                open(10,file=file_name,status='old',position='append',form='formatted')
               endif
                  write(10,*) bc%gr4( k )%params(:)
               close(10)

           !Denormalize here for now
           bc%gr4( k )%params(1) = (bc%gr4( k )%params(1) ) * (x1max - x1min) + x1min
           bc%gr4( k )%params(2) = (bc%gr4( k )%params(2) ) * (x2max - x2min) + x2min
           bc%gr4( k )%params(3) = (bc%gr4( k )%params(3) ) * (x3max - x3min) + x3min
           bc%gr4( k )%params(4) = (bc%gr4( k )%params(4) ) * (x4max - x4min) + x4min

          end do

         end if
#endif

         if ( c_Ks          == 1 ) call control_2_var( infil%GA(:)%Ks , infil%nland , 0 )
         if ( c_PsiF        == 1 ) call control_2_var( infil%GA(:)%PsiF , infil%nland , 0 )
         if ( c_DeltaTheta  == 1 ) call control_2_var( infil%GA(:)%DeltaTheta, infil%nland , 0 )
         if ( c_lambda      == 1 ) call control_2_var( infil%SCS(:)%lambdacn , infil%nland , 0 )
         if ( c_CN          == 1 ) call control_2_var( infil%SCS(:)%CN, infil%nland , 0 )

         if ( c_hydrograph == 1 ) then
            do k = 1,bc%nb_in
               call control_2_var( bc%hyd( k )%q(:) , size( bc%hyd( k )%q(:) ) , 1 )
            end do
         end if

         if      ( c_ratcurve == 1 ) then
            do k = 1,bc%nb_out
               call control_2_var( bc%rat( k )%q(:) , size( bc%rat( k )%q(:) ) , 1 )
            end do

         else if ( c_ratcurve == 2 ) then
            do k = 1,bc%nb_out
               call control_2_var( bc%rat( k )%pow(:) , 2 , 1 )
            end do

         end if

         if ( c_rain == 1 ) then
            do k = 1,bc%nb_rn
               call control_2_var( bc%rain( k )%q(:) , size( bc%rain( k )%q(:) ) , 1 )
            end do
         end if

      #endif


   CONTAINS


      SUBROUTINE control_2_var( var , n , data_glob )

         implicit none

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         integer(ip), intent(in)  ::  n , data_glob

         real(rp), dimension(n), intent(out)  ::  var

         !=============================================================================================================!
         !
         !=============================================================================================================!

         integer(ip)  ::  n_k=0,k_loc

         !=============================================================================================================!
         !
         !=============================================================================================================!

         do k_loc = 0,np-1

            if ( k_loc == 0 ) then

               var( 1 : n )  =  control( ic : ic + n - 1 )

               if ( data_glob == 0 ) ic = ic + n

            else

               call mpi_send_recv_scal_i( n , n_k , k_loc , 0 )

               call mpi_send_recv_array_r( control( ic : ic + n_k - 1 ) , var( 1 : n ) , 0 , k_loc )

               if ( data_glob == 0 ) ic = ic + n_k

            end if

            call mpi_wait_all

         end do

         if ( data_glob == 1 ) ic = ic + n

      END SUBROUTINE control_2_var


   END SUBROUTINE read_control


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Reading the Model Input Perturbation Control Vector
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE read_control_diff( dof0 , mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( unk ), intent(inout)  ::  dof0
      type( msh ), intent(in   )  ::  mesh

      !================================================================================================================!
      !  Filling input control vector ( control )
      !================================================================================================================!

      ic = 1

      #ifdef USE_SW_MONO

         manning_diff   (:)  =  0._rp
         manning_beta_diff(:) = 0._rp
         bathy_cell_diff(:)  =  0._rp
         dof0_diff%h    (:)  =  0._rp
         dof0_diff%u    (:)  =  0._rp
         dof0_diff%v    (:)  =  0._rp

         slope_y_diff = 0._rp
         slope_x_diff = 0._rp

         infil_diff%GA(:)%Ks           = 0._rp
         infil_diff%GA(:)%PsiF         = 0._rp
         infil_diff%GA(:)%DeltaTheta   = 0._rp
         infil_diff%SCS(:)%lambdacn       = 0._rp
         infil_diff%SCS(:)%CN           = 0._rp

#ifdef USE_HYDRO
         do k = 1,bc%nb_gr4in
            bc_diff%gr4( k )%params(:)  =  0._rp
         end do
#endif

#ifdef USE_HYDRO
         do k = 1,bc%nb_gr4in
            bc_diff%gr4( k )%params(:)  =  0._rp
         end do
#endif

         do k = 1,bc%nb_in
            bc_diff%hyd( k )%q(:)  =  0._rp
         end do

         do k = 1,bc%nb_out
            bc_diff%rat( k )%q  (:)  =  0._rp
            bc_diff%rat( k )%pow(:)  =  0._rp
         end do

         do k = 1,bc%nb_rn
            bc_diff%rain( k )%q(:)  =  0._rp
         end do


        if ( c_shape_s == 1 ) then
            call control_diff_2_var( XSshape_diff(:)%s    ,size(XSshape_diff)   , 0 )
        endif
        if ( c_hmax == 1 ) then
            call control_diff_2_var( XSshape_diff(:)%hmax    ,size(XSshape_diff)   , 0 )
        endif
        if ( c_xcenter == 1 ) then
            call control_diff_2_var( XSshape_diff(:)%xcenter    ,size(XSshape_diff)   , 0 )
        endif

         if ( c_manning == 1 ) call control_diff_2_var( manning_diff    , nland   , manning_data_glob )
         if ( c_manning_beta == 1 ) call control_diff_2_var( manning_beta_diff    , nland   , manning_data_glob )
         if ( c_bathy   == 1 ) call control_diff_2_var( bathy_cell_diff , mesh%nc , 0                 )

         if ( c_slope_y == 1 ) call control_diff_2_var( slope_y , size(slope_y) , 0                 )
         if ( c_slope_x == 1 ) call control_diff_2_var( slope_x , size(slope_x) , 0                 )

         if ( c_ic      == 1 ) call control_diff_2_var( dof0_diff%h     , mesh%nc , 0                 )
         if ( c_ic      == 1 ) call control_diff_2_var( dof0_diff%u     , mesh%nc , 0                 )
         if ( c_ic      == 1 ) call control_diff_2_var( dof0_diff%v     , mesh%nc , 0                 )

#ifdef USE_HYDRO
         if ( c_gr4params == 1 ) then

            do k = 1,bc%nb_gr4in
               call control_diff_2_var( bc_diff%gr4( k )%params(:) , size(bc_diff%gr4( k )%params(:)) , 1 )
            end do

         end if
#endif

         if ( c_Ks          == 1 ) call control_diff_2_var( infil%GA(:)%Ks , infil%nland , 0 )
         if ( c_PsiF        == 1 ) call control_diff_2_var( infil%GA(:)%PsiF , infil%nland , 0 )
         if ( c_DeltaTheta  == 1 ) call control_diff_2_var( infil%GA(:)%DeltaTheta, infil%nland , 0 )
         if ( c_lambda      == 1 ) call control_diff_2_var( infil%SCS(:)%lambdacn , infil%nland , 0 )
         if ( c_CN          == 1 ) call control_diff_2_var(infil%SCS(:)%CN, infil%nland , 0 )

         if ( c_hydrograph == 1 ) then
            do k = 1,bc%nb_in
               call control_diff_2_var( bc_diff%hyd( k )%q(:) , size( bc_diff%hyd( k )%q(:) ) , 1 )
            end do
         end if

         if      ( c_ratcurve == 1 ) then
            do k = 1,bc%nb_out
               call control_diff_2_var( bc_diff%rat( k )%q(:) , size( bc_diff%rat( k )%q(:) ) , 1 )
            end do

         else if ( c_ratcurve == 2 ) then
            do k = 1,bc%nb_out
               call control_diff_2_var( bc_diff%rat( k )%pow(:) , 2 , 1 )
            end do

         end if

         if ( c_rain == 1 ) then
            do k = 1,bc%nb_rn
               call control_diff_2_var( bc_diff%rain( k )%q(:) , size( bc_diff%rain( k )%q(:) ) , 1 )
            end do
         end if

      #endif


   CONTAINS


      SUBROUTINE control_diff_2_var( var , n , data_glob )

         implicit none

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         integer(ip), intent(in)  ::  n , data_glob

         real(rp), dimension(n), intent(out)  ::  var

         !=============================================================================================================!
         !
         !=============================================================================================================!

         integer(ip)  ::  n_k=0,k_loc

         !=============================================================================================================!
         !
         !=============================================================================================================!

         do k_loc = 0,np-1

            if ( k_loc == 0 ) then

               var( 1 : n )  =  control_diff( ic : ic + n - 1 )

               if ( data_glob == 0 ) ic = ic + n

            else

               call mpi_send_recv_scal_i( n , n_k , k_loc , 0 )

               call mpi_send_recv_array_r( control_diff( ic : ic + n_k - 1 ) , var( 1 : n ) , 0 , k_loc )

               if ( data_glob == 0 ) ic = ic + n_k

            end if

            call mpi_wait_all

         end do

         if ( data_glob == 1 ) ic = ic + n

      END SUBROUTINE control_diff_2_var


   END SUBROUTINE read_control_diff


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Reading the Model Input Perturbated Control Vector
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE read_control_perturb( dof0 , mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( unk ), intent(inout)  ::  dof0
      type( msh ), intent(in   )  ::  mesh

      !================================================================================================================!
      !  Filling input control vector ( control )
      !================================================================================================================!

      ic = 1

      #ifdef USE_SW_MONO

         bathy_cell(:) = bathy_cell_copy(:)

         dof0 = dof0_copy


        if ( c_shape_s == 1 ) then
            call control_perturb_2_var( XSshape(:)%s    , size(XSshape)   , 0 )
        endif
        if ( c_hmax == 1 ) then
            call control_perturb_2_var( XSshape(:)%hmax    , size(XSshape)   , 0 )
        endif
        if ( c_xcenter == 1 ) then
            call control_perturb_2_var( XSshape(:)%xcenter    , size(XSshape)   , 0 )
        endif

         if ( c_manning == 1 ) call control_perturb_2_var( manning    , nland   , manning_data_glob )
         if ( c_manning_beta == 1 ) call control_perturb_2_var( manning_beta    , nland   , manning_data_glob )
         if ( c_bathy   == 1 ) call control_perturb_2_var( bathy_cell , mesh%nc , 0                 )

         if ( c_slope_y == 1 ) call control_perturb_2_var( slope_y , size(slope_y) , 0                 )
         if ( c_slope_x == 1 ) call control_perturb_2_var( slope_x , size(slope_x) , 0                 )

         if ( c_ic      == 1 ) call control_perturb_2_var( dof0%h     , mesh%nc , 0                 )
         if ( c_ic      == 1 ) call control_perturb_2_var( dof0%u     , mesh%nc , 0                 )
         if ( c_ic      == 1 ) call control_perturb_2_var( dof0%v     , mesh%nc , 0                 )

         if ( c_Ks          == 1 ) call control_perturb_2_var( infil%GA(:)%Ks , infil%nland , 0 )
         if ( c_PsiF        == 1 ) call control_perturb_2_var( infil%GA(:)%PsiF , infil%nland , 0 )
         if ( c_DeltaTheta  == 1 ) call control_perturb_2_var(infil%GA(:)%DeltaTheta, infil%nland , 0 )
         if ( c_lambda      == 1 ) call control_perturb_2_var( infil%SCS(:)%lambdacn , infil%nland , 0 )
         if ( c_CN          == 1 ) call control_perturb_2_var(infil%SCS(:)%CN, infil%nland , 0 )

#ifdef USE_HYDRO
         if ( c_gr4params == 1 ) then

            do k = 1,bc%nb_gr4in
               call control_perturb_2_var( bc%gr4( k )%params(:) , size(bc%gr4( k )%params(:)) , 1 )
            end do

         end if
#endif

         if ( c_hydrograph == 1 ) then
            do k = 1,bc%nb_in
               call control_perturb_2_var( bc%hyd( k )%q(:) , size( bc%hyd( k )%q(:) ) , 1 )
            end do
         end if

         if      ( c_ratcurve == 1 ) then
            do k = 1,bc%nb_out
               call control_perturb_2_var( bc%rat( k )%q(:) , size( bc%rat( k )%q(:) ) , 1 )
            end do

         else if ( c_ratcurve == 2 ) then
            do k = 1,bc%nb_out
               call control_perturb_2_var( bc%rat( k )%pow(:) , 2 , 1 )
            end do

         end if

         if ( c_rain == 1 ) then
            do k = 1,bc%nb_rn
               call control_perturb_2_var( bc%rain( k )%q(:) , size( bc%rain( k )%q(:) ) , 1 )
            end do
         end if

      #endif


   CONTAINS


      SUBROUTINE control_perturb_2_var( var , n , data_glob )

         implicit none

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         integer(ip), intent(in)  ::  n , data_glob

         real(rp), dimension(n), intent(out)  ::  var

         !=============================================================================================================!
         !
         !=============================================================================================================!

         integer(ip)  ::  n_k=0,k_loc

         !=============================================================================================================!
         !
         !=============================================================================================================!
		write(*,*) "in control_perturb_2_var( var , n , data_glob )"

         do k_loc = 0,np-1

            if ( k_loc == 0 ) then

               var( 1 : n )  =  control_perturb( ic : ic + n - 1 )

               if ( data_glob == 0 ) ic = ic + n

            else

               call mpi_send_recv_scal_i( n , n_k , k_loc , 0 )

               call mpi_send_recv_array_r( control_perturb( ic : ic + n_k - 1 ) , var( 1 : n ) , 0 , k_loc )

               if ( data_glob == 0 ) ic = ic + n_k

            end if

            call mpi_wait_all

         end do

         if ( data_glob == 1 ) ic = ic + n

      END SUBROUTINE control_perturb_2_var


   END SUBROUTINE read_control_perturb


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Output in proper Files the Control Vector
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE output_control( dof0 , mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type(msh), intent(in)  ::  mesh

      type(unk), intent(inout)  ::  dof0

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      character(len=lchar)  ::  file_name

      !================================================================================================================!
      !  Filling input control vector ( control )
      !================================================================================================================!

      call system('mkdir -p min')

      call read_control( dof0 , mesh )

      #ifdef USE_SW_MONO

#ifdef USE_HYDRO
          if ( proc == 0 .and. c_gr4params == 1 ) then
             do k = 1,bc%nb_gr4in

               write(file_name,'(A,I3.3,A,I3.3)') 'min/gr4params_' , k , '.' , ite_min
               open(10,file=file_name,status='replace',form='formatted')

               do i = 1,4
                  write(10,*) bc%gr4( k )%params(i)
               end do

               close(10)

            end do
         end if
#endif

         if ( proc == 0 .and. c_hydrograph == 1 ) then

            do k = 1,bc%nb_in

               write(file_name,'(A,I3.3,A,I3.3)') 'min/hydrograph_' , k , '.' , ite_min

               open(10,file=file_name,status='replace',form='formatted')

               do i = 1,size( bc%hyd( k )%t(:) )
                  write(10,*) bc%hyd( k )%t(i) , bc%hyd( k )%q(i)
               end do

               close(10)

            end do

         end if

        if ( proc == 0 .and. c_rain == 1 ) then

            do k = 1,bc%nb_rn

               write(file_name,'(A,I3.3,A,I3.3)') 'min/rain' , k , '.' , ite_min

               open(10,file=file_name,status='replace',form='formatted')

               do i = 1,size( bc%rain( k )%t(:) )
                  write(10,*) bc%rain( k )%t(i) , bc%rain( k )%q(i)
               end do

               close(10)

            end do

         end if

         if ( proc == 0 .and. c_manning == 1 ) then

            write(file_name,'(A,I3.3)') 'min/manning.' , ite_min

            open(10,file=file_name,status='replace',form='formatted')

            do i = 1,nland

               write(10,*) i , manning(i)

            end do

            close(10)



!            if ( manning_data_glob == 1 ) then

!               do i = 1,nland

!                  write(buffer,'(A,I4.4,A)') 'min/land_' , i , '.txt'

!                  call write_pscalar( (/ real(ite_min,8) , manning(i) /) , &
!                                      (/ 'Ite'   , 'Man'      /) , buffer , 2 )

!               end do

!            end if

         end if

        if ( proc == 0 .and. c_manning_beta == 1 ) then

            write(file_name,'(A,I3.3)') 'min/manning_beta.' , ite_min

            open(10,file=file_name,status='replace',form='formatted')

            do i = 1,nland

               write(10,*) i , manning_beta(i)

            end do

            close(10)

         end if

         if ( proc == 0 .and. c_bathy == 1 ) then

            write(file_name,'(A,I3.3)') 'min/bathy.' , ite_min

            open(10,file=file_name,status='replace',form='formatted')

!             call write_scalar_field( bathy_cell , mesh , file_name )
              do i = 1,mesh%nc
                write(10,*) i , bathy_cell(i)
              enddo

            close(10)
         end if


         if ( proc == 0 .and. c_Ks == 1 ) then

            write(file_name,'(A,I3.3)') 'min/Ks.' , ite_min

            open(10,file=file_name,status='replace',form='formatted')

            do i=1,infil%nland

               write(10,*) i , infil%GA(i)%Ks

            enddo

         endif

         if ( proc == 0 .and. c_PsiF == 1 ) then

            write(file_name,'(A,I3.3)') 'min/PsiF.' , ite_min

            open(10,file=file_name,status='replace',form='formatted')

            do i=1,infil%nland

               write(10,*) i , infil%GA(i)%PsiF

            enddo

         endif

         if ( proc == 0 .and. c_DeltaTheta == 1  ) then

            write(file_name,'(A,I3.3)') 'min/DeltaTheta.' , ite_min

            open(10,file=file_name,status='replace',form='formatted')

            do i=1,infil%nland

               write(10,*) i , infil%GA(i)%DeltaTheta

            enddo

         endif

         if ( proc == 0 .and. c_lambda == 1  ) then

            write(file_name,'(A,I3.3)') 'min/lambdacn.' , ite_min

            open(10,file=file_name,status='replace',form='formatted')

            do i=1,infil%nland

               write(10,*) i , infil%SCS(i)%lambdacn

            enddo

         endif

         if ( proc == 0 .and. c_CN == 1  ) then

            write(file_name,'(A,I3.3)') 'min/CN.' , ite_min

            open(10,file=file_name,status='replace',form='formatted')

            do i=1,infil%nland

               write(10,*) i , infil%SCS(i)%CN

            enddo

         endif

         if ( proc == 0 .and. ((c_shape_s == 1) .or. (c_hmax == 1) .or. (c_xcenter == 1))) then

            write(file_name,'(A,I3.3)') 'min/xshape.' , ite_min

            open(10,file=file_name,status='replace',form='formatted')

!             call write_scalar_field( bathy_cell , mesh , file_name )
              do i = 1,size(XSshape)
                write(10,*) i , XSshape(i)%s , XSshape(i)%xcenter , XSshape(i)%hmax
              enddo

            close(10)
         end if

         if ( proc == 0 .and. c_slope_y == 1 ) then

            write(file_name,'(A,I3.3)') 'min/slope_y.' , ite_min

            open(10,file=file_name,status='replace',form='formatted')

              do i = 1,size(slope_y)
                write(10,*) i , slope_y(i)
              enddo

            close(10)
         end if

         if ( proc == 0 .and. c_slope_x == 1 ) then

            write(file_name,'(A,I3.3)') 'min/slope_x.' , ite_min

            open(10,file=file_name,status='replace',form='formatted')

              do i = 1,size(slope_x)
                write(10,*) i , slope_x(i)
              enddo

            close(10)
         end if

      #endif

   END SUBROUTINE output_control


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Output in proper Files the Back Control Vector (the gradient of the cost function)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE output_control_back( mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type(msh), intent(in)  ::  mesh

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      character(len=lchar)  ::  file_name

      real(rp), dimension( mesh%nc + mesh%ncb )  ::  temp

      !================================================================================================================!
      !  Filling input control vector ( control )
      !================================================================================================================!

      call system('mkdir -p grad')

      #ifdef USE_SW_MONO

#ifdef USE_HYDRO
         if ( proc == 0 .and. c_gr4params == 1) then
            do k = 1,bc%nb_gr4in

               write(file_name,'(A,I3.3,A)') 'grad/gr4params' , k , '_grad'
               open(10,file=file_name,status='replace',form='formatted')

               do i = 1,4
                  write(10,*) bc_back%gr4( k )%params(i)
               end do

               close(10)

            end do
         end if
#endif

         if ( proc == 0 ) then

            do k = 1,bc%nb_in

               write(file_name,'(A,I3.3,A)') 'grad/hydrograph' , k , '_grad'

               open(10,file=file_name,status='replace',form='formatted')

               do i = 1,size( bc%hyd( k )%t(:) )
                  write(10,*) bc%hyd( k )%t(i) , bc_back%hyd( k )%q(i) 
               end do

               close(10)

            end do

            do k = 1,bc%nb_rn

               write(file_name,'(A,I3.3,A)') 'grad/rain' , k , '_grad'

               open(10,file=file_name,status='replace',form='formatted')

               do i = 1,size( bc%rain( k )%t(:) )
                  write(10,*) bc%rain( k )%t(i) , bc_back%rain( k )%q(i)
               end do

               close(10)

            end do

         end if

!          call write_scalar_field( bathy_cell_back , mesh , 'grad/bathy_grad' )



               write(file_name,'(A,I3.3,A)') 'grad/bathy_grad'

               open(10,file=file_name,status='replace',form='formatted')

               do i = 1,mesh%nc
                  write(10,*) i, bathy_cell_back(i)
               end do

               close(10)


         do i = 1,mesh%nc
            temp(i)  =  manning_back( land(i) )
         end do

         call write_scalar_field( temp , mesh , 'grad/manning_grad' )

         do i = 1,mesh%nc
            temp(i)  =  manning_beta_back( land(i) )
         end do

         call write_scalar_field( temp , mesh , 'grad/manning_beta_grad' )

      #endif

   END SUBROUTINE output_control_back


!################# !################# !################# !################# !#################
!################# !################# !################# !################# !#################
! temporaryly, add minimize cost here
!################# !################# !################# !################# !#################
!################# !################# !################# !################# !#################






!      implicit none

!      !================================================================================================================!
!      !  Interface Variables
!      !================================================================================================================!

!      type( msh ), intent(in   )  ::  mesh
!      type( unk ), intent(inout)  ::  dof0
!      type( unk ), intent(inout)  ::  dof

!      !================================================================================================================!
!      !  Local Variables
!      !================================================================================================================!

!      external  simul_rc , euclid , ctonbe , ctcabe

!	  integer(ip) :: ite_line_search
!      integer(ip)  ::  n            ! Control vector size
!      real(rp)     ::  dxmin        ! Resolution for x in Linf norm
!      real(rp)     ::  dJ_expect    ! Estimation of the expected decrease in cost during the first iteration
!      real(rp)     ::  epsg         ! Stopping criterion
!      character*3  ::  normtype     ! specifies the norm used to test optimality
!      integer(ip)  ::  impres       ! Print option for m1qn3
!      integer(ip)  ::  io           ! File output label
!      integer(ip)  ::  imode(3)     ! Input mode
!      integer(ip)  ::  omode        ! Output mode
!      integer(ip)  ::  niter        ! Maximum number of minimization iterations
!      integer(ip)  ::  nsim         ! Maximum number of simulator calls
!      integer(ip)  ::  iz(5)        ! Adress of working array for m1qn3
!      integer(ip)  ::  ndz          ! Dimension of dz
!      integer(ip)  ::  reverse      ! Specifies direct or reverse mode
!      integer(ip)  ::  indic        ! m1qn3 indicates it needs cost/gradient

!      integer(ip)  ::  izs(1)       !
!      real(rp)     ::  rzs(1)       ! Working areas not needed in reverse mode
!      real(rp)     ::  dzs(1)       !

!      integer(ip)  :: max_ite_line_search
!      real(rp), dimension(:), allocatable  ::  dz     ! Adress of working array for m1qn3

!      real(rp)  ::  norm_gradJ

!      real(rp)  ::  cost_ini , norm_grad_cost_ini(100) , norm_grad_cost(100)

!!~       real(rp), dimension(:), allocatable :: control_temp

!      !================================================================================================================!
!      !  Initialization
!      !================================================================================================================!
!		write(*,*) "in minimize cost"
!      call system('mkdir -p min')

!      call alloc_back_vars(dof0_back, dof_back, mesh ) ! allocate _back variables (both dof0_back and other control vector components, to be call only on time

!      call write_control( dof0 , mesh )

!!~       allocate(control_temp(size(control)))

!      !================================================================================================================!
!      !  m1qn3 arguments
!      !================================================================================================================!

!      n          =   size( control )
!      dxmin      =   1.d-5
!      epsg       =   eps_min
!      normtype   =   'two'
!      impres     =   5
!      io         =   40
!      imode(1)   =   0
!      imode(2)   =   0
!      imode(3)   =   1
!      ndz        =   4 * n + 100 * ( 2 * n + 1 )
!      reverse    =   1
!      indic      =   4
!      max_ite_line_search = 15

!      allocate( dz( ndz ) ) ; dz(:) = 0._rp

!      nsim  =  100

!      if ( restart_min == 0 ) then

!         niter  =  nsim

!      else

!         niter  =  restart_min

!      end if

!      !================================================================================================================!
!      !  Initialization / Restarting
!      !================================================================================================================!

!      inquire( file = 'min/min_cost.txt' , exist = file_exist(1) )

!      call mpi_wait_all

!      if ( .not. file_exist(1) ) then

!         ite_min         = 0
!         ite_line_search = 0

!         open(30,file='min/min_cost.txt'    ,status='replace',form = 'formatted')
!         open(40,file='min/m1qn3_output.txt',status='replace',form = 'formatted')

!         write(30,'(A)') '# ite cost norm_grad_cost'

!      else

!         open(30,file='min/min_cost.txt'    ,status='old',position = 'append',form='formatted')
!         open(40,file='min/m1qn3_output.txt',status='old',position = 'append',form='formatted')

!         allocate( control_back( size(control) ) )

!         call read_restart_m1qn3_v2

!         imode(2)   =   1
!         epsg       =   eps_min / epsg
!         reverse    =   1

!      end if



!# --> reverse
!# --> indic
!      !================================================================================================================!
!      !  Minimization ... M1QN3 loop in reverse mode
!      !================================================================================================================!

!      do while( reverse >= 0 .and. indic /= 0 )
!         !write(*,*) 'indic ' , indic, 'reverse', reverse, ' omode ', omode

!         !=============================================================================================================!
!         !  Case indic = 4 -> M1QN3 needs new values of the cost and its gradient to compute
!         !  either the step descent (line-search) or the new iterate
!         !=============================================================================================================!

!         if ( indic == 4 ) then

!            cost_back = one

!            verbose = -1
!write(*,*) '--control before adjoint', control
!            call adjoint_model( mesh , dof0 , dof )
!write(*,*) '--control after adjoint' , control
!            ite_line_search = ite_line_search + 1

!         end if

!         !=============================================================================================================!
!         !  As M1QN3 is sequential, only master thread perform the minimization
!         !  The control vector is managed in the module m_adjoint
!         !=============================================================================================================!

!         if ( proc == 0 ) then

!            !==========================================================================================================!
!            !  Store the initial cost and its gradient norm
!            !==========================================================================================================!

!            if ( ite_min == 0 .and. indic == 4 ) then

!               cost_ini  =  cost

!               i = 1

!               do k = 1,nb_vars_in_control

!                  norm_grad_cost_ini(k) = sqrt( sum( control_back( i : i - 1 + dim_vars_in_control(k) ) * &
!                                                     control_back( i : i - 1 + dim_vars_in_control(k) ) ) )

!                  i = i + dim_vars_in_control(k)

!               end do

!            end if

!            !==========================================================================================================!
!            !  Normalization of the control vector
!            !==========================================================================================================!

!            i = 1

!            do k = 1,nb_vars_in_control

!               control( i : i - 1 + dim_vars_in_control(k) ) = &
!               control( i : i - 1 + dim_vars_in_control(k) ) * norm_grad_cost_ini(k) / cost_ini

!               i = i + dim_vars_in_control(k)

!            end do

!            !==========================================================================================================!
!            !  Normalization of the cost, its gradient and the control vector
!            !==========================================================================================================!

!            if ( indic == 4 ) then

!               cost = cost / cost_ini

!               i = 1

!               do k = 1,nb_vars_in_control

!                  control_back( i : i - 1 + dim_vars_in_control(k) ) = &
!                  control_back( i : i - 1 + dim_vars_in_control(k) ) / norm_grad_cost_ini(k)

!                  norm_grad_cost(k) = sqrt( sum( control_back( i : i - 1 + dim_vars_in_control(k) ) * &
!                                                 control_back( i : i - 1 + dim_vars_in_control(k) ) ) )

!                  i = i + dim_vars_in_control(k)

!               end do

!               norm_gradJ  =  sqrt( sum( control_back(:) * control_back(:) ) )

!               write(6,'("ite ",I3," , J =",ES13.6," , |grad J| =",ES13.6, " , time =",ES13.6)') &
!               ite_min , cost , norm_gradJ , time(1)

!            end if

!            !==========================================================================================================!
!            !  dJ_expect is df1 in the M1QN3 documentation
!            !==========================================================================================================!

!            dJ_expect = 0.5_rp * cost

!            !==========================================================================================================!
!            !  Call of M1QN3
!            !==========================================================================================================!

!!~ 			control_temp(:) = control(:)
!            call m1qn3( simul_rc     , &
!                        euclid       , &
!                        ctonbe       , &
!                        ctcabe       , &
!                        n            , &
!                        control      , &
!                        cost         , &
!                        control_back , &
!                        dxmin        , &
!                        dJ_expect    , &
!                        epsg         , &
!                        normtype     , &
!                        impres       , &
!                        io           , &
!                        imode        , &
!                        omode        , &
!                        niter        , &
!                        nsim         , &
!                        iz           , &
!                        dz           , &
!                        ndz          , &
!                        reverse      , &
!                        indic        , &
!                        izs          , &
!                        rzs          , &
!                        dzs )

!!~ 			if (all(control(:) == control_temp(:) )) print *, "control is not changed here l344"
!            !==========================================================================================================!
!            !  Back normalization of the control vector in order to output it
!            !==========================================================================================================!

!            i = 1

!            do k = 1,nb_vars_in_control

!               control( i : i - 1 + dim_vars_in_control(k) ) = &
!               control( i : i - 1 + dim_vars_in_control(k) ) / norm_grad_cost_ini(k) * cost_ini

!               i = i + dim_vars_in_control(k)

!            end do

!         end if

!         call mpi_wait_all

!         call mpi_bcast_i( reverse , 0 )
!         call mpi_bcast_i( indic   , 0 )
!         call mpi_bcast_i( omode   , 0 )

!         call write_restart_m1qn3_v2

!         !=============================================================================================================!
!         !  Case indic = 1 -> M1QN3 has finished an iterate ( imode(3) = 1 )
!         !=============================================================================================================!

!         if      ( indic == 1 ) then

!            call output_control( dof0 , mesh )

!            if ( proc == 0 ) write(30,'(I4,102ES15.7)') ite_min , cost , norm_grad_cost(1:nb_vars_in_control)

!            ite_min = ite_min + 1

!            ite_line_search = 0

!         else if ( indic == 4 .and. ite_line_search > max_ite_line_search) then

!            call Stopping_Program_Sub( 'Stopping Minimization, too much iterations in line-search' )

!         end if

!      end do

!      close(30)
!      close(40)

!      call mpi_wait_all

!      deallocate( dz )

!      if ( omode == 1 .or. &
!           omode == 3 .or. &
!           omode == 4 .or. &
!           omode == 5 .or. &
!           omode == 6 ) then

!         verbose = 1

!         call model_direct( mesh , dof0 , dof )

!         call Print_Screen( 'cost' , cost )

!      end if

!      if      ( omode == 1 ) then
!         call Stopping_Program_Sub( 'Successful Minimization, the test on the gradient norm (eps_min) is satisfied' )

!      else if ( omode == 2 ) then
!         call Stopping_Program_Sub( 'Stopping Minimization, input parmeter not well initialized (omode=2 in the M1QN3 documentation)' )

!      else if ( omode == 3 ) then
!         call Stopping_Program_Sub( 'Stopping Minimization, the line-search is blocked (omode=3 in the M1QN3 documentation)' )

!      else if ( omode == 4 ) then
!         call Stopping_Program_Sub( 'Stopping Minimization, maximal number of iterations reached (omode=4 in the M1QN3 documentation)' )

!      else if ( omode == 5 ) then
!         call Stopping_Program_Sub( 'Stopping Minimization, maximal number of simulations reached (omode=5 in the M1QN3 documentation)' )

!      else if ( omode == 6 ) then
!         call Stopping_Program_Sub( 'Stopping Minimization, stop on dxmin during the line-search (omode=6 in the M1QN3 documentation)' )

!      else

!         call Stopping_Program_Sub( 'Stopping Minimization, unknown reason !?)' )

!      end if

!CONTAINS


!	SUBROUTINE write_restart_m1qn3_v2

!	 implicit none

!	 if ( proc == 0 ) then

!		open(50,file='min/restart_min.bin',status='replace',form='unformatted')

!		write(50) ite_min , ite_line_search , reverse , indic , epsg , iz , dz, &
!				  control , cost , control_back , cost_ini , norm_grad_cost_ini

!		close(50)

!	 end if

!	END SUBROUTINE write_restart_m1qn3_v2


!	SUBROUTINE read_restart_m1qn3_v2

!	 implicit none

!	 open(50,file='min/restart_min.bin',status='old',form='unformatted')

!	 read(50) ite_min , ite_line_search , reverse , indic , epsg , iz , dz, &
!			  control , cost , control_back , cost_ini , norm_grad_cost_ini

!	 close(50)

!	END SUBROUTINE read_restart_m1qn3_v2

! END SUBROUTINE minimize_cost_v2


END MODULE m_adjoint
