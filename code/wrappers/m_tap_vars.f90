MODULE m_tap_vars
   USE m_common
   USE m_obs
   USE m_model
   real(rp) :: dt_diff , tc_diff
   real(rp) :: dt_back , tc_back
   type(bcs), target :: bc_diff
   type(bcs), target :: bc_back
   type(infiltration_data), target :: infil_diff
   type(infiltration_data), target :: infil_back
   type( innovation_obs ), dimension(:), allocatable, target :: innovation_diff
   type( innovation_obs ), dimension(:), allocatable, target :: innovation_back
    real(rp), dimension(:), allocatable :: manning_diff , manning_beta_diff, bathy_cell_diff
    real(rp), dimension(:), allocatable :: manning_back ,manning_beta_back, bathy_cell_back
CONTAINS
   SUBROUTINE alloc_diff_vars( dof0_diff , dof_diff , mesh )
      type( unk ), intent(inout) :: dof0_diff , dof_diff
      type( msh ), intent(in) :: mesh
      call alloc_dof( dof0_diff , mesh )
      call alloc_dof( dof_diff , mesh )
      allocate( bc_diff%inflow ( size( bc%inflow ) ) )
      allocate( bc_diff%outflow( size( bc%outflow ) ) )
      bc_diff%inflow (:) = 0._rp
      bc_diff%outflow(:) = 0._rp
      allocate( bc_diff%hyd( bc%nb_in ) )
      do i = 1,bc%nb_in
         allocate( bc_diff%hyd( i )%t( size( bc%hyd( i )%t ) ) )
         allocate( bc_diff%hyd( i )%q( size( bc%hyd( i )%q ) ) )
         bc_diff%hyd( i )%t(:) = 0._rp
         bc_diff%hyd( i )%q(:) = 0._rp
      end do
      allocate( bc_diff%rat( bc%nb_out ) )
      do i = 1,bc%nb_out
         allocate( bc_diff%rat( i )%h( size( bc%rat( i )%h ) ) )
         allocate( bc_diff%rat( i )%q( size( bc%rat( i )%q ) ) )
         bc_diff%rat( i )%h(:) = 0._rp
         bc_diff%rat( i )%q(:) = 0._rp
      end do
      do i = 1,bc%nb_rn
         allocate( bc_diff%rain( i )%t( size( bc%rain( i )%t ) ) )
         allocate( bc_diff%rain( i )%q( size( bc%rain( i )%q ) ) )
         bc_diff%rain( i )%t(:) = 0._rp
         bc_diff%rain( i )%q(:) = 0._rp
      end do
      allocate( innovation_diff ( size( innovation ) ) )
      do iobs = 1,size( innovation )
         allocate( innovation_diff ( iobs )%diff( size( innovation ( iobs )%diff ) ) )
         innovation_diff ( iobs )%diff(:) = 0._rp
      end do
      dt_diff = 0._rp
     allocate( manning_diff ( size( manning ) ) )
     allocate( manning_beta_diff ( size( manning_beta ) ) )
     allocate( bathy_cell_diff ( size( bathy_cell ) ) )
     manning_diff(:) = 0._rp
     manning_beta_diff(:)= 0._rp
     bathy_cell_diff(:) = 0._rp
     allocate( bc_diff%sum_mass_flux( bc%nb ) )
   END SUBROUTINE alloc_diff_vars
   SUBROUTINE alloc_back_vars(dof0_back, dof_back, mesh )
  type(unk) , intent(inout) :: dof0_back, dof_back
      type(msh), intent(in) :: mesh
      call alloc_dof( dof0_back , mesh )
      call alloc_dof( dof_back , mesh )
      !===============================================
      ! Classical hydraulic BCs
      !===============================================
      allocate( bc_back%inflow ( size( bc%inflow ) ) )
      allocate( bc_back%outflow( size( bc%outflow ) ) )
      bc_back%inflow (:) = 0._rp
      bc_back%outflow(:) = 0._rp
      allocate( bc_back%hyd( bc%nb_in ) )
      do i = 1,bc%nb_in
         allocate( bc_back%hyd( i )%t( size( bc%hyd( i )%t ) ) )
         allocate( bc_back%hyd( i )%q( size( bc%hyd( i )%q ) ) )
         bc_back%hyd( i )%t(:) = 0._rp
         bc_back%hyd( i )%q(:) = 0._rp
      end do
      allocate( bc_back%rat( bc%nb_out ) )
      do i = 1,bc%nb_out
        write(*,*) "rating curve must always be allocated while doing minimization (even if not used)"
         allocate( bc_back%rat( i )%h( size( bc%rat( i )%h ) ) )
         allocate( bc_back%rat( i )%q( size( bc%rat( i )%q ) ) )
         bc_back%rat( i )%h(:) = 0._rp
         bc_back%rat( i )%q(:) = 0._rp
      end do
      !===============================================
      ! Rain BCs
      !===============================================
      if (bc_rain == 1) then
        allocate( bc_back%rain( bc%nb_rn ) )
        do i = 1,bc%nb_rn
                ! dirty
            allocate( bc_back%rain( i )%t( size( bc%rain( i )%t ) ) )
            allocate( bc_back%rain( i )%q( size( bc%rain( i )%q ) ) )
            bc_back%rain( i )%t(:) = 0._rp
            bc_back%rain( i )%q(:) = 0._rp
        enddo
      else
            allocate( bc_back%rain( 1 ) )
            allocate( bc_back%rain( 1 )%t( 1 ) )
            allocate( bc_back%rain( 1 )%q( 1 ) )
            bc_back%rain( 1 )%t(1) = 0._rp
            bc_back%rain( 1 )%q(1) = 0._rp
      endif
      !===============================================
      ! Infiltration params
      !===============================================
      allocate( infil_back%GA ( size( infil%GA ) ) )
      allocate( infil_back%SCS( size( infil%SCS ) ) )
      infil_back%GA(:)%Ks = 0._rp
      infil_back%GA(:)%PsiF = 0._rp
      infil_back%GA(:)%DeltaTheta = 0._rp
      infil_back%SCS(:)%lambda = 0._rp
      infil_back%SCS(:)%CN = 0._rp
      allocate( innovation_back ( size( innovation ) ) )
      do iobs = 1,size( innovation )
         allocate( innovation_back ( iobs )%diff( size( innovation( iobs )%diff ) ) )
         innovation_back ( iobs )%diff(:) = 0._rp
      end do
      dt_back = 0._rp
         allocate( manning_back ( size( manning ) ) )
         allocate( manning_beta_back ( size( manning_beta ) ) )
         allocate( bathy_cell_back( size( bathy_cell ) ) )
         manning_back(:) = 0._rp
         manning_beta_back(:) = 0._rp
         bathy_cell_back(:) = 0._rp
         allocate( bc_back%sum_mass_flux( bc%nb ) )
   END SUBROUTINE alloc_back_vars
SUBROUTINE dealloc_back_vars()
      !type(unk) , intent(inout) :: dof0_back, dof_back
      !if(allocated(dof0_back )) deallocate(dof0_back)
      !if(allocated(dof_back )) deallocate(dof_back)
      if(allocated(bc_back%inflow )) deallocate(bc_back%inflow )
      if(allocated(bc_back%outflow )) deallocate(bc_back%outflow )
      if(allocated(bc_back%rat )) deallocate(bc_back%rat )
      if(allocated(bc_back%hyd )) deallocate(bc_back%hyd )
      if(allocated( bc_back%rain)) deallocate( bc_back%rain)
      if(allocated(infil_back%GA )) deallocate(infil_back%GA )
      if(allocated(infil_back%SCS )) deallocate(infil_back%SCS)
      if(allocated(innovation_back)) deallocate(innovation_back)
      if(allocated(manning_back)) deallocate(manning_back)
      if(allocated(manning_beta_back)) deallocate(manning_beta_back)
      if(allocated(bathy_cell_back)) deallocate(bathy_cell_back)
      if(allocated(bc_back%sum_mass_flux)) deallocate(bc_back%sum_mass_flux)
      !if(allocated(control )) deallocate(control ) !to check this does not generate issues lilian
      !if(allocated(control_back )) deallocate(control_back ) !to check this does not generate issues lilian
      !if(allocated(control_diff )) deallocate(control_diff ) !to check this does not generate issues lilian
END SUBROUTINE dealloc_back_vars
END MODULE m_tap_vars
