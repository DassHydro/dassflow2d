MODULE m_adjoint
   USE m_common
   USE m_linear_algebra
   USE m_mesh
   USE m_random
   USE m_time_screen
   USE m_tap_vars
   USE m_obs
   USE m_model
   implicit none
   type(unk), target :: dof0_diff
   type(unk), target :: dof0_back
   type(unk), target :: dof_diff
   type(unk), target :: dof_back
   real(rp), dimension(:), allocatable :: control ! Control Vector
   real(rp), dimension(:), allocatable :: control_back ! Adjoint Control Vector
   real(rp), dimension(:), allocatable :: control_diff ! Perturbation Control Vector
   real(rp), dimension(:), allocatable :: control_perturb ! Perturbated Control Vector
   real(rp), dimension(:), allocatable :: control_lbound !< Lower bounds vector.
   real(rp), dimension(:), allocatable :: control_ubound !< Upper bounds vector.
   real(rp) :: cost ! Cost function
   real(rp) :: cost_back ! Adjoint Cost function
   real(rp) :: cost_diff ! Perturbation Cost function
   real(rp) :: cost_perturb ! Perturbated Cost function
   integer(ip) :: ic ! Control vector index
   integer(ip) :: ite_min ! Iteration of Minimization Procedure
   integer(ip) :: nb_vars_in_control ! Number of variables in the control vector
   integer(ip), dimension(100) :: dim_vars_in_control ! Dimension of each variable in the control vector
   integer(ip) :: dim_all ! Sum of Dimensions of all variables in the control vector
   type( unk ) :: dof0_copy
  real(rp), dimension(:), allocatable :: bathy_cell_copy
CONTAINS
   SUBROUTINE model_direct( mesh , dof0 , dof )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(in ) :: mesh
      type(unk), intent(inout) :: dof0
      type(unk), intent(inout) :: dof
      !================================================================================================================!
      ! Read control vectror
      !================================================================================================================!
      call read_control( dof0 , mesh )
      !================================================================================================================!
      ! Calling direct model
      !================================================================================================================!
         call run_model( mesh , dof0 , dof , cost )
   END SUBROUTINE model_direct
   SUBROUTINE model_direct_perturb( mesh , dof0 , dof )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(in ) :: mesh
      type(unk), intent(inout) :: dof0
      type(unk), intent(inout) :: dof
      !================================================================================================================!
      ! Read perturbated control vector
      !================================================================================================================!
      call read_control_perturb( dof0 , mesh )
      !================================================================================================================!
      ! Calling direct model
      !================================================================================================================!
      call Time_Init(1_ip)
         call run_model( mesh , dof0 , dof , cost_perturb )
      call Time_End(1_ip)
   END SUBROUTINE model_direct_perturb
   SUBROUTINE linear_tangent_model( mesh , dof0 , dof )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in ) :: mesh
      type( unk ), intent(inout) :: dof0
      type( unk ), intent(inout) :: dof
      !================================================================================================================!
      ! Read control vector
      !================================================================================================================!
      call read_control( dof0 , mesh )
      !================================================================================================================!
      ! Read perturbation control vector
      !================================================================================================================!
      call read_control_diff( dof0_diff , mesh )
      !================================================================================================================!
      ! Calling linear tangent model
      !================================================================================================================!
      call Time_Init(1_ip)
      call run_model_diff( mesh , dof0 , dof0_diff , dof , dof_diff , cost , cost_diff )
      call Time_End(1_ip)
   END SUBROUTINE linear_tangent_model
 SUBROUTINE adjoint_model( mesh , dof0 , dof )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(in ) :: mesh
      type(unk), intent(inout) :: dof0
      type(unk), intent(inout) :: dof
  integer(ip) :: i_loc
      !================================================================================================================!
      ! Read control vector
      !================================================================================================================!
      call read_control( dof0 , mesh )
write(*,*) "control read in ajoint model", control
      !================================================================================================================!
      ! Calling adjoint model
      !================================================================================================================!
      call Time_Init(1_ip)
      if ( proc /= 0 ) cost_back = 0._rp
         dof_back%h(:) = 0._rp
         dof_back%u(:) = 0._rp
         dof_back%v(:) = 0._rp
         dof_back%infil(:) = 0._rp
         infil_back%GA(:)%Ks = 0._rp
         infil_back%GA(:)%PsiF = 0._rp
         infil_back%GA(:)%DeltaTheta = 0._rp
         infil_back%SCS(:)%lambda = 0._rp
         infil_back%SCS(:)%CN = 0._rp
         bathy_cell_back(:) = 0._rp
         manning_back(:) = 0._rp
         manning_beta_back(:) = 0._rp
write(*,*) "control before run_model_back", control
write(*,*) "manning", manning(:)
         call run_model_back( mesh , dof0 , dof0_back , dof , dof_back , cost , cost_back )
write(*,*) "control after run_model_back", control
write(*,*) "manning", manning(:)
        !write(*,*) "DONE call run_model_back"
      call Time_End(1_ip)
     !================================================================================================================!
      ! Filling control_back vector (cost gradient vector)
      !================================================================================================================!
      call write_control_back( dof0_back , mesh )
   END SUBROUTINE adjoint_model
   SUBROUTINE calc_grad_cost( mesh , dof0 , dof )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in ) :: mesh
      type( unk ), intent(inout) :: dof0
      type( unk ), intent(inout) :: dof
      !================================================================================================================!
      ! Initialization
      !================================================================================================================!
        !call dealloc_back_vars() ! might become usefull to call it for python, hope to user
  call alloc_back_vars(dof0_back, dof_back, mesh )
      call write_control( dof0 , mesh )
      verbose = -1
      cost_back = 1._rp
      !================================================================================================================!
      ! Calc cost function gradient using Adjoint Model
      !================================================================================================================!
      call adjoint_model( mesh , dof0 , dof )
      !================================================================================================================!
      ! Output in proper Files the Gradient of the Cost Function
      !================================================================================================================!
      call output_control_back( mesh )
   END SUBROUTINE calc_grad_cost
   SUBROUTINE test_adjoint( mesh , dof0 , dof , arg )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in ) :: mesh
      type( unk ), intent(inout) :: dof0
      type( unk ), intent(inout) :: dof
      character(len=*), intent(in) :: arg
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      integer(ip) :: test_type
      real(rp) :: eps , scal1 , scal2
      !================================================================================================================!
      ! Initialization
      !================================================================================================================!
      call alloc_diff_vars( dof0_diff , dof_diff , mesh )
      call alloc_back_vars(dof0_back, dof_back, mesh )
      !================================================================================================================!
      ! Gradient Direction
      !================================================================================================================!
      call write_control ( dof0 , mesh )
      call write_control_diff( dof0 , mesh )
      allocate( control_perturb ( size( control ) ) )
      allocate( control_back ( size( control ) ) )
      !================================================================================================================!
      ! Direct Test (2 run that should return the same p_Y result
      !================================================================================================================!
      read(arg,*) test_type
      verbose = -1
      if ( test_type == 0 .or. test_type == 1 ) then
         if ( proc == 0 ) write(6,*)
         call model_direct( mesh , dof0 , dof )
         call Print_Screen( 'cost' , cost )
         call model_direct( mesh , dof0 , dof )
         call Print_Screen( 'cost' , cost )
      end if
      !================================================================================================================!
      ! Tangent Gradient Test
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
         call Print_Screen( 'cost' , cost )
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
            control_perturb(:) = control(:) + eps * control_diff(:)
            open(80,file='time_step_serie_for_adjoint.bin',form='unformatted',status='old')
            fix_time_step_serie = 2
            call model_direct_perturb( mesh , dof0 , dof )
            if ( proc == 0 ) write(6,'(3ES23.15)') eps , &
                                                   ( cost_perturb - cost ) / eps, &
                                                   abs( one - ( cost_perturb - cost ) / ( eps * cost_diff ) )
            close(80)
            eps = 0.1_rp * eps
         end do
      end if
      !================================================================================================================!
      ! Backward Gradient Test
      !================================================================================================================!
      if ( test_type == 0 .or. test_type == 3 ) then
         call Print_Screen( 'backward_gradient_test' )
         if ( proc == 0 ) open(80,file='time_step_serie_for_adjoint.bin',form='unformatted',status='replace')
         fix_time_step_serie = 1
         cost_back = one
         call adjoint_model( mesh , dof0 , dof )
         cost_diff = sum( control_diff(:) * control_back(:) )
         if ( proc == 0 ) then
            close(80)
            write(6,*)
         end if
         call Print_Screen( 'cost' , cost )
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
            control_perturb(:) = control(:) + eps * control_diff(:)
            open(80,file='time_step_serie_for_adjoint.bin',form='unformatted',status='old')
            fix_time_step_serie = 2
            call model_direct_perturb( mesh , dof0 , dof )
            if ( proc == 0 ) write(6,'(3ES23.15)') eps , &
                                                   ( cost_perturb - cost ) / eps, &
                                                   abs( one - ( cost_perturb - cost ) / ( eps * cost_diff ) )
            close(80)
            eps = 0.1_rp * eps
         end do
      end if
      !================================================================================================================!
      ! Backward Scalar Product Test
      !================================================================================================================!
      if ( test_type == 0 .or. test_type == 4 ) then
         fix_time_step_serie = 0
         call Print_Screen( 'backward_scalar_product_test' )
         call linear_tangent_model( mesh , dof0 , dof )
         call Print_Screen( 'cost' , cost )
         call Print_Screen( 'cost_diff' , cost_diff )
         cost_back = cost_diff
         scal1 = cost_back * cost_diff
         call adjoint_model( mesh , dof0 , dof )
         call Print_Screen( 'cost' , cost )
         call Print_Screen( 'cost_back' , cost_back )
         scal2 = sum( control_diff(:) * control_back(:) )
         if ( proc == 0 ) write(6,'(A,2ES22.15)') ' scal1 = ' , scal1
         if ( proc == 0 ) write(6,'(A,2ES22.15)') ' scal2 = ' , scal2
         if ( proc == 0 ) write(6,'(A,2ES22.15)') ' relative error = ' , abs( ( scal2 - scal1 ) / scal1 )
      end if
   END SUBROUTINE test_adjoint
   SUBROUTINE write_control( dof0 , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(in ) :: mesh
      type(unk), intent(inout) :: dof0
      real(rp) :: x1min,x1max,x2min,x2max,x3min,x3max,x4min,x4max
      !================================================================================================================!
      ! Filling desired Control Vector (with c_ in input.txt)
      !================================================================================================================!
      nb_vars_in_control = 0
      ic = 1
         if ( c_manning == 1 ) call var_2_control( manning , nland , manning_data_glob )
         if ( c_manning_beta == 1 ) call var_2_control( manning_beta , nland , manning_data_glob )
         if ( c_bathy == 1 ) call var_2_control( bathy_cell , mesh%nc , 0 )
         if ( c_ic == 1 ) call var_2_control( dof0%h , mesh%nc , 0 )
         if ( c_ic == 1 ) call var_2_control( dof0%u , mesh%nc , 0 )
         if ( c_ic == 1 ) call var_2_control( dof0%v , mesh%nc , 0 )
         if ( c_Ks == 1 ) call var_2_control( infil%GA(:)%Ks , infil%nland , 0 )
         if ( c_PsiF == 1 ) call var_2_control( infil%GA(:)%PsiF , infil%nland , 0 )
         if ( c_DeltaTheta == 1 ) call var_2_control( infil%GA(:)%DeltaTheta, infil%nland , 0 )
         if ( c_lambda == 1 ) call var_2_control( infil%SCS(:)%lambda , infil%nland , 0 )
         if ( c_CN == 1 ) call var_2_control( infil%SCS(:)%CN, infil%nland , 0 )
         if ( c_hydrograph == 1 ) then
            do k = 1,bc%nb_in
               call var_2_control( bc%hyd( k )%q(:) , size( bc%hyd( k )%q(:) ) , 1 )
            end do
         end if
         if ( c_ratcurve == 1 ) then
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
         bathy_cell_copy(:) = bathy_cell(:)
         call alloc_dof( dof0_copy , mesh )
         dof0_copy = dof0
   CONTAINS
      SUBROUTINE var_2_control( var , n , data_glob )
         implicit none
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         integer(ip), intent(in) :: n , data_glob
         real(rp), dimension(n), intent(in) :: var
         !=============================================================================================================!
         !
         !=============================================================================================================!
         integer(ip) :: n_k
         integer(ip) :: k_loc
         !=============================================================================================================!
         !
         !=============================================================================================================!
         nb_vars_in_control = nb_vars_in_control + 1
         do k_loc = 0,np-1
            if ( k_loc == 0 ) then
               call alloc_or_larger_r( control , ic + n - 1 )
               control( ic : ic + n - 1 ) = var( 1 : n )
               ic = ic + n
               dim_vars_in_control( nb_vars_in_control ) = n
            else if ( data_glob == 0 ) then
               if ( proc == 0 ) call alloc_or_larger_r( control , ic + n_k - 1 )
               if ( proc == 0 ) ic = ic + n_k
               if ( proc == 0 ) dim_vars_in_control( nb_vars_in_control ) = &
                                dim_vars_in_control( nb_vars_in_control ) + n_k
            end if
         end do
      END SUBROUTINE var_2_control
   END SUBROUTINE write_control
   SUBROUTINE write_control_diff( dof0 , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( unk ), intent(in) :: dof0
      type( msh ), intent(in) :: mesh
      !================================================================================================================!
      ! Filling input control vector ( control )
      !================================================================================================================!
  write(*,*) "in write_control_diff( dof0 , mesh )"
      ic = 1
         if ( c_manning == 1 ) call var_2_control_diff( eps_manning * manning , nland , manning_data_glob )
         if ( c_manning_beta == 1 ) call var_2_control_diff( eps_manning * manning_beta , nland , manning_data_glob ) ! should have manning_BETA_data_glob ??
         if ( c_bathy == 1 ) call var_2_control_diff( eps_bathy * bathy_cell , mesh%nc , 0 )
         if ( c_ic == 1 ) call var_2_control_diff( eps_ic * dof0%h , mesh%nc , 0 )
         if ( c_ic == 1 ) call var_2_control_diff( eps_ic * dof0%u , mesh%nc , 0 )
         if ( c_ic == 1 ) call var_2_control_diff( eps_ic * dof0%v , mesh%nc , 0 )
         if ( c_Ks == 1 ) call var_2_control_diff( eps_Ks * infil%GA(:)%Ks , infil%nland , 0 )
         if ( c_PsiF == 1 ) call var_2_control_diff( eps_PsiF * infil%GA(:)%PsiF , infil%nland , 0 )
         if ( c_DeltaTheta == 1 ) call var_2_control_diff( eps_DeltaTheta * infil%GA(:)%DeltaTheta, infil%nland , 0 )
         if ( c_lambda == 1 ) call var_2_control_diff( infil%SCS(:)%lambda , infil%nland , 0 )
         if ( c_CN == 1 ) call var_2_control_diff(infil%SCS(:)%CN, infil%nland , 0 )
         if ( c_hydrograph == 1 ) then
            do k = 1,bc%nb_in
               call var_2_control_diff( eps_hydrograph * bc%hyd( k )%q(:) , size( bc%hyd( k )%q(:) ) , 1 )
            end do
         end if
         if ( c_ratcurve == 1 ) then
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
   CONTAINS
      SUBROUTINE var_2_control_diff( var , n , data_glob )
         implicit none
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         integer(ip), intent(in) :: n , data_glob
         real(rp), dimension(n), intent(in) :: var
         !=============================================================================================================!
         !
         !=============================================================================================================!
         integer(ip) :: n_k,k_loc
         real(rp), dimension(n) :: rn
         !=============================================================================================================!
         !
         !=============================================================================================================!
         call init_random_seed
         do k_loc = 0,np-1
            if ( k_loc == 0 ) then
               call alloc_or_larger_r( control_diff , ic + n - 1 )
               call random_number( rn )
               control_diff( ic : ic + n - 1 ) = ( -1._rp + 2._rp * rn ( 1 : n ) ) * var( 1 : n )
               ic = ic + n
            else if ( data_glob == 0 ) then
               if ( proc == 0 ) call alloc_or_larger_r( control_diff , ic + n_k - 1 )
               call random_number( rn )
               if ( proc == 0 ) ic = ic + n_k
            end if
         end do
      END SUBROUTINE var_2_control_diff
   END SUBROUTINE write_control_diff
   SUBROUTINE write_control_back( dof0 , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(unk), intent(in) :: dof0
      type(msh), intent(in) :: mesh
      !================================================================================================================!
      ! Filling input control vector ( control )
      !================================================================================================================!
      ic = 1
         if ( c_manning == 1 ) call var_2_control_back( manning_back , nland , manning_data_glob )
         if ( c_manning_beta == 1 ) call var_2_control_back( manning_beta_back , nland , manning_data_glob )
         if ( c_bathy == 1 ) call var_2_control_back( bathy_cell_back , mesh%nc , 0 )
         if ( c_ic == 1 ) call var_2_control_back( dof0_back%h , mesh%nc , 0 )
         if ( c_ic == 1 ) call var_2_control_back( dof0_back%u , mesh%nc , 0 )
         if ( c_ic == 1 ) call var_2_control_back( dof0_back%v , mesh%nc , 0 )
         if ( c_Ks == 1 ) call var_2_control_back( infil_back%GA%Ks , infil%nland , manning_data_glob )
         if ( c_PsiF == 1 ) call var_2_control_back( infil_back%GA%PsiF , infil%nland , 0 )
         if ( c_DeltaTheta == 1 ) call var_2_control_back( infil_back%GA%DeltaTheta , infil%nland , 0 )
         if ( c_lambda == 1 ) call var_2_control_back( infil_back%SCS%lambda, infil%nland , 0 )
         if ( c_CN == 1 ) call var_2_control_back( infil_back%SCS%CN , infil%nland , 0 )
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
   CONTAINS
      SUBROUTINE var_2_control_back( var , n , data_glob )
         implicit none
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         integer(ip), intent(in) :: n , data_glob
         real(rp), dimension(n), intent(in) :: var
         !=============================================================================================================!
         !
         !=============================================================================================================!
         integer(ip) :: n_k=0,k_loc
         !=============================================================================================================!
         !
         !=============================================================================================================!
         do k_loc = 0,np-1
            if ( k_loc == 0 ) then
               call alloc_or_larger_r( control_back , ic + n - 1 )
               control_back( ic : ic + n - 1 ) = var( 1 : n )
               ic = ic + n
            else if ( data_glob == 0 ) then
               if ( proc == 0 ) call alloc_or_larger_r( control_back , ic + n_k - 1 )
               if ( proc == 0 ) ic = ic + n_k
            end if
         end do
      END SUBROUTINE var_2_control_back
   END SUBROUTINE write_control_back
   !> Filling the model input adjoint control vector.
   !!
   !! \details Add to the control vector k variable. In the end of this function the vector k is created.
   !! \param[in] mesh Mesh of the model.
   !! \param[inout] dof0 Initial conditions.
   SUBROUTINE write_control_bounds( dof0 , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( unk ), intent(in) :: dof0
      type( msh ), intent(in) :: mesh
      !================================================================================================================!
      ! Filling input control vector ( control )
      !================================================================================================================!
      ic = 1
   CONTAINS
      !> Filling the model input control bounds vectors.
      !!
      !! \details This subroutine adds to the control bounds vectors the bounds for a variable array of dimension n.
      !! \param[in] lb Lower bound to append to the control lower-bounds vectors.
      !! \param[in] ub Upper bound to append to the control upper-bounds vectors.
      !! \param[in] n Number of components in the control vector.
      !! \param[in] data_glob Variable not used.
      SUBROUTINE var_2_control_bounds( lb, ub, n, data_glob )
         implicit none
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         integer(ip), intent(in) :: n , data_glob
         real(rp), intent(in) :: lb
         real(rp), intent(in) :: ub
         call alloc_or_larger_r( control_lbound , ic + n - 1 )
         call alloc_or_larger_r( control_ubound , ic + n - 1 )
         control_lbound( ic : ic + n - 1 ) = lb
         control_ubound( ic : ic + n - 1 ) = ub
         ic = ic + n
      END SUBROUTINE var_2_control_bounds
      !BELOW: 1D routine, not adapted for 2D
      !> Filling the model input control bounds vectors.
      !!
      !! \details This subroutine adds to the control bounds vectors the bounds for the bathymetry.
      !! \param[in] lb Lower bound to append to the control lower-bounds vectors.
      !! \param[in] ub Upper bound to append to the control upper-bounds vectors.
      !! \param[in] var Array of initial values of bathymetry.
      !! \param[in] delta Admissible delta around initial values of bathymetry.
      !! \param[in] n Size of the control vector.
      !! \param[in] data_glob Variable not used.
      !> Append boundaries for a scalar component of the control vector to the array of control boundaries
      !!
      !! \param[in] lb Lower bound for the scalar component of the control vector
      !! \param[in] ub Upper bound for the scalar component of the control vector
      !! \param[in] data_glob Unused variable.
      SUBROUTINE var_2_control_scal_bounds( lb, ub, data_glob )
         implicit none
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         integer(ip), intent(in) :: data_glob
         real(rp), intent(in) :: lb
         real(rp), intent(in) :: ub
         call alloc_or_larger_r( control_lbound , ic )
         call alloc_or_larger_r( control_ubound , ic )
         control_lbound( ic ) = lb
         control_ubound( ic ) = ub
         ic = ic + 1
      END SUBROUTINE var_2_control_scal_bounds
   END SUBROUTINE write_control_bounds
   SUBROUTINE read_control( dof0 , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(unk), intent(inout) :: dof0
      type(msh), intent(in ) :: mesh
      character(len=lchar) :: file_name
      real(rp) :: x1min,x1max,x2min,x2max,x3min,x3max,x4min,x4max
      !================================================================================================================!
      ! Filling input control vector ( control )
      !================================================================================================================!
      ic = 1
         bathy_cell(:) = bathy_cell_copy(:)
         if ( c_manning == 1 ) call control_2_var( manning , nland , manning_data_glob )
         if ( c_manning_beta == 1 ) call control_2_var( manning_beta , nland , manning_data_glob )
         if ( c_bathy == 1 ) call control_2_var( bathy_cell , mesh%nc , 0 )
         do ie = 1,mesh%neb
            i = mesh%edge( mesh%edgeb(ie)%ind )%cell(2)
            j = mesh%edge( mesh%edgeb(ie)%ind )%cell(1)
            if ( mesh%edgeb(ie)%typlim == 'wall' ) bathy_cell(i) = bathy_cell(j)
         end do
         if ( c_ic == 1 ) call control_2_var( dof0%h , mesh%nc , 0 )
         if ( c_ic == 1 ) call control_2_var( dof0%u , mesh%nc , 0 )
         if ( c_ic == 1 ) call control_2_var( dof0%v , mesh%nc , 0 )
         if ( c_Ks == 1 ) call control_2_var( infil%GA(:)%Ks , infil%nland , 0 )
         if ( c_PsiF == 1 ) call control_2_var( infil%GA(:)%PsiF , infil%nland , 0 )
         if ( c_DeltaTheta == 1 ) call control_2_var( infil%GA(:)%DeltaTheta, infil%nland , 0 )
         if ( c_lambda == 1 ) call control_2_var( infil%SCS(:)%lambda , infil%nland , 0 )
         if ( c_CN == 1 ) call control_2_var( infil%SCS(:)%CN, infil%nland , 0 )
         if ( c_hydrograph == 1 ) then
            do k = 1,bc%nb_in
               call control_2_var( bc%hyd( k )%q(:) , size( bc%hyd( k )%q(:) ) , 1 )
            end do
         end if
         if ( c_ratcurve == 1 ) then
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
write(*,*) "manning value after  read control", manning
   CONTAINS
      SUBROUTINE control_2_var( var , n , data_glob )
         implicit none
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         integer(ip), intent(in) :: n , data_glob
         real(rp), dimension(n), intent(out) :: var
         !=============================================================================================================!
         !
         !=============================================================================================================!
         integer(ip) :: n_k=0,k_loc
         !=============================================================================================================!
         !
         !=============================================================================================================!
         do k_loc = 0,np-1
            if ( k_loc == 0 ) then
               var( 1 : n ) = control( ic : ic + n - 1 )
               if ( data_glob == 0 ) ic = ic + n
            else
               if ( data_glob == 0 ) ic = ic + n_k
            end if
         end do
         if ( data_glob == 1 ) ic = ic + n
      END SUBROUTINE control_2_var
   END SUBROUTINE read_control
   SUBROUTINE read_control_diff( dof0 , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( unk ), intent(inout) :: dof0
      type( msh ), intent(in ) :: mesh
      !================================================================================================================!
      ! Filling input control vector ( control )
      !================================================================================================================!
      ic = 1
         manning_diff (:) = 0._rp
         manning_beta_diff(:) = 0._rp
         bathy_cell_diff(:) = 0._rp
         dof0_diff%h (:) = 0._rp
         dof0_diff%u (:) = 0._rp
         dof0_diff%v (:) = 0._rp
         infil_diff%GA(:)%Ks = 0._rp
         infil_diff%GA(:)%PsiF = 0._rp
         infil_diff%GA(:)%DeltaTheta = 0._rp
         infil_diff%SCS(:)%lambda = 0._rp
         infil_diff%SCS(:)%CN = 0._rp
         do k = 1,bc%nb_in
            bc_diff%hyd( k )%q(:) = 0._rp
         end do
         do k = 1,bc%nb_out
            bc_diff%rat( k )%q (:) = 0._rp
            bc_diff%rat( k )%pow(:) = 0._rp
         end do
         do k = 1,bc%nb_rn
            bc_diff%rain( k )%q(:) = 0._rp
         end do
         if ( c_manning == 1 ) call control_diff_2_var( manning_diff , nland , manning_data_glob )
         if ( c_manning_beta == 1 ) call control_diff_2_var( manning_beta_diff , nland , manning_data_glob )
         if ( c_bathy == 1 ) call control_diff_2_var( bathy_cell_diff , mesh%nc , 0 )
         if ( c_ic == 1 ) call control_diff_2_var( dof0_diff%h , mesh%nc , 0 )
         if ( c_ic == 1 ) call control_diff_2_var( dof0_diff%u , mesh%nc , 0 )
         if ( c_ic == 1 ) call control_diff_2_var( dof0_diff%v , mesh%nc , 0 )
         if ( c_Ks == 1 ) call control_diff_2_var( infil%GA(:)%Ks , infil%nland , 0 )
         if ( c_PsiF == 1 ) call control_diff_2_var( infil%GA(:)%PsiF , infil%nland , 0 )
         if ( c_DeltaTheta == 1 ) call control_diff_2_var( infil%GA(:)%DeltaTheta, infil%nland , 0 )
         if ( c_lambda == 1 ) call control_diff_2_var( infil%SCS(:)%lambda , infil%nland , 0 )
         if ( c_CN == 1 ) call control_diff_2_var(infil%SCS(:)%CN, infil%nland , 0 )
         if ( c_hydrograph == 1 ) then
            do k = 1,bc%nb_in
               call control_diff_2_var( bc_diff%hyd( k )%q(:) , size( bc_diff%hyd( k )%q(:) ) , 1 )
            end do
         end if
         if ( c_ratcurve == 1 ) then
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
   CONTAINS
      SUBROUTINE control_diff_2_var( var , n , data_glob )
         implicit none
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         integer(ip), intent(in) :: n , data_glob
         real(rp), dimension(n), intent(out) :: var
         !=============================================================================================================!
         !
         !=============================================================================================================!
         integer(ip) :: n_k=0,k_loc
         !=============================================================================================================!
         !
         !=============================================================================================================!
         do k_loc = 0,np-1
            if ( k_loc == 0 ) then
               var( 1 : n ) = control_diff( ic : ic + n - 1 )
               if ( data_glob == 0 ) ic = ic + n
            else
               if ( data_glob == 0 ) ic = ic + n_k
            end if
         end do
         if ( data_glob == 1 ) ic = ic + n
      END SUBROUTINE control_diff_2_var
   END SUBROUTINE read_control_diff
   SUBROUTINE read_control_perturb( dof0 , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( unk ), intent(inout) :: dof0
      type( msh ), intent(in ) :: mesh
      !================================================================================================================!
      ! Filling input control vector ( control )
      !================================================================================================================!
      ic = 1
         bathy_cell(:) = bathy_cell_copy(:)
         dof0 = dof0_copy
         if ( c_manning == 1 ) call control_perturb_2_var( manning , nland , manning_data_glob )
         if ( c_manning_beta == 1 ) call control_perturb_2_var( manning_beta , nland , manning_data_glob )
         if ( c_bathy == 1 ) call control_perturb_2_var( bathy_cell , mesh%nc , 0 )
         if ( c_ic == 1 ) call control_perturb_2_var( dof0%h , mesh%nc , 0 )
         if ( c_ic == 1 ) call control_perturb_2_var( dof0%u , mesh%nc , 0 )
         if ( c_ic == 1 ) call control_perturb_2_var( dof0%v , mesh%nc , 0 )
         if ( c_Ks == 1 ) call control_perturb_2_var( infil%GA(:)%Ks , infil%nland , 0 )
         if ( c_PsiF == 1 ) call control_perturb_2_var( infil%GA(:)%PsiF , infil%nland , 0 )
         if ( c_DeltaTheta == 1 ) call control_perturb_2_var(infil%GA(:)%DeltaTheta, infil%nland , 0 )
         if ( c_lambda == 1 ) call control_perturb_2_var( infil%SCS(:)%lambda , infil%nland , 0 )
         if ( c_CN == 1 ) call control_perturb_2_var(infil%SCS(:)%CN, infil%nland , 0 )
         if ( c_hydrograph == 1 ) then
            do k = 1,bc%nb_in
               call control_perturb_2_var( bc%hyd( k )%q(:) , size( bc%hyd( k )%q(:) ) , 1 )
            end do
         end if
         if ( c_ratcurve == 1 ) then
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
   CONTAINS
      SUBROUTINE control_perturb_2_var( var , n , data_glob )
         implicit none
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         integer(ip), intent(in) :: n , data_glob
         real(rp), dimension(n), intent(out) :: var
         !=============================================================================================================!
         !
         !=============================================================================================================!
         integer(ip) :: n_k=0,k_loc
         !=============================================================================================================!
         !
         !=============================================================================================================!
  write(*,*) "in control_perturb_2_var( var , n , data_glob )"
         do k_loc = 0,np-1
            if ( k_loc == 0 ) then
               var( 1 : n ) = control_perturb( ic : ic + n - 1 )
               if ( data_glob == 0 ) ic = ic + n
            else
               if ( data_glob == 0 ) ic = ic + n_k
            end if
         end do
         if ( data_glob == 1 ) ic = ic + n
      END SUBROUTINE control_perturb_2_var
   END SUBROUTINE read_control_perturb
   SUBROUTINE output_control( dof0 , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(in) :: mesh
      type(unk), intent(inout) :: dof0
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      character(len=lchar) :: file_name
      !================================================================================================================!
      ! Filling input control vector ( control )
      !================================================================================================================!
      call system('mkdir -p min')
      call read_control( dof0 , mesh )
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
         end if
        if ( proc == 0 .and. c_manning_beta == 1 ) then
            write(file_name,'(A,I3.3)') 'min/manning_beta.' , ite_min
            open(10,file=file_name,status='replace',form='formatted')
            do i = 1,nland
               write(10,*) i , manning_beta(i)
            end do
            close(10)
         end if
         if ( c_bathy == 1 ) then
            write(file_name,'(A,I3.3)') 'min/bathy.' , ite_min
          ! call write_scalar_field( bathy_cell , mesh , file_name )
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
         if ( proc == 0 .and. c_DeltaTheta == 1 ) then
            write(file_name,'(A,I3.3)') 'min/DeltaTheta.' , ite_min
            open(10,file=file_name,status='replace',form='formatted')
            do i=1,infil%nland
               write(10,*) i , infil%GA(i)%DeltaTheta
            enddo
         endif
         if ( proc == 0 .and. c_lambda == 1 ) then
            write(file_name,'(A,I3.3)') 'min/lambda.' , ite_min
            open(10,file=file_name,status='replace',form='formatted')
            do i=1,infil%nland
               write(10,*) i , infil%SCS(i)%lambda
            enddo
         endif
         if ( proc == 0 .and. c_CN == 1 ) then
            write(file_name,'(A,I3.3)') 'min/CN.' , ite_min
            open(10,file=file_name,status='replace',form='formatted')
            do i=1,infil%nland
               write(10,*) i , infil%SCS(i)%CN
            enddo
         endif
   END SUBROUTINE output_control
   SUBROUTINE output_control_back( mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(in) :: mesh
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      character(len=lchar) :: file_name
      real(rp), dimension( mesh%nc + mesh%ncb ) :: temp
      !================================================================================================================!
      ! Filling input control vector ( control )
      !================================================================================================================!
      call system('mkdir -p grad')
         if ( proc == 0 ) then
            do k = 1,bc%nb_in
               write(file_name,'(A,I3.3,A)') 'grad/hydrograph' , k , '_grad'
               open(10,file=file_name,status='replace',form='formatted')
               do i = 1,size( bc%hyd( k )%t(:) )
                  write(10,*) bc%hyd( k )%t(i) , bc_back%hyd( k )%q(i) !,bc_diff%hyd( k )%q(i)
                  !write(*,*) bc%hyd( k )%t(i) , bc_back%hyd( k )%q(i)
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
      ! call write_scalar_field( bathy_cell_back , mesh , 'grad/bathy_grad' )
         do i = 1,mesh%nc
            temp(i) = manning_back( land(i) )
         end do
      ! call write_scalar_field( temp , mesh , 'grad/manning_grad' )
         do i = 1,mesh%nc
            temp(i) = manning_beta_back( land(i) )
         end do
     ! call write_scalar_field( temp , mesh , 'grad/manning_beta_grad' )
   END SUBROUTINE output_control_back
END MODULE m_adjoint
