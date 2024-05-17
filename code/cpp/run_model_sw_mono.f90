SUBROUTINE run_model( mesh , dof0 , dof , cost )
   USE m_common
   USE m_linear_algebra
   USE m_mesh
   USE m_mpi
   USE m_time_screen !NOADJ
   USE m_numeric
   USE m_model
   USE m_obs
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(inout) :: mesh
   type( unk ), intent(inout) :: dof0
   type( unk ), intent(inout) :: dof
   real(rp), intent(out) :: cost
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   integer(ip) :: sub_nt, iR, iL
    real(rp) :: bathy_temp
    real(rp) :: Ks_ubound, Ks_lbound, z, lambda, shift
   !===================================================================================================================!
   ! Define parameterized bathymetry
   !===================================================================================================================!
   if (use_xsshp == 1) then
      if (xsshp_along_y == 1) then
        do ie = 1,mesh%nc
            bathy_temp = bathy_cell(ie)
            if (mesh%cell(ie)%grav%x <= XSshape(1)%xcenter) then
            bathy_cell( ie ) = XSshape(1)%topz + &
            min( 0._rp, - abs(XSshape(1)%hmax) *&
            (1 - ( abs((mesh%cell(ie)%grav%x - XSshape(1)%xleft) -&
            (XSshape(1)%xcenter-XSshape(1)%xleft))/&
            ((XSshape(1)%xcenter-XSshape(1)%xleft)) )**XSshape(1)%s)) &
            + bathy_temp &
            + slope_y(1) * mesh%cell(ie)%grav%y &
            + slope_x(1) * mesh%cell(ie)%grav%x
            else
            bathy_cell( ie ) = XSshape(1)%topz + &
            min( 0._rp, - abs(XSshape(1)%hmax) *&
            (1 - ( abs((mesh%cell(ie)%grav%x - XSshape(1)%xleft) -&
            (XSshape(1)%xcenter-XSshape(1)%xleft))/&
            ((XSshape(1)%xright-XSshape(1)%xcenter)) )**XSshape(1)%s)) &
            + bathy_temp &
            + slope_y(1) * mesh%cell(ie)%grav%y &
            + slope_x(1) * mesh%cell(ie)%grav%x
            endif
        enddo
     elseif (xsshp_along_x == 1) then
        do ie = 1,mesh%nc
            bathy_temp = bathy_cell(ie)
            if (mesh%cell(ie)%grav%y <= XSshape(1)%xcenter) then
            bathy_cell( ie ) = XSshape(1)%topz + &
            min( 0._rp, - abs(XSshape(1)%hmax) *&
            (1 - ( abs((mesh%cell(ie)%grav%y - XSshape(1)%xleft) -&
            (XSshape(1)%xcenter-XSshape(1)%xleft))/&
            ((XSshape(1)%xcenter-XSshape(1)%xleft)) )**XSshape(1)%s)) &
            + bathy_temp &
            + slope_y(1) * mesh%cell(ie)%grav%y &
            + slope_x(1) * mesh%cell(ie)%grav%x
            else
            bathy_cell( ie ) = XSshape(1)%topz + &
            min( 0._rp, - abs(XSshape(1)%hmax) *&
            (1 - ( abs((mesh%cell(ie)%grav%y - XSshape(1)%xleft) -&
            (XSshape(1)%xcenter-XSshape(1)%xleft))/&
            ((XSshape(1)%xright-XSshape(1)%xcenter)) )**XSshape(1)%s)) &
            + bathy_temp &
            + slope_y(1) * mesh%cell(ie)%grav%y &
            + slope_x(1) * mesh%cell(ie)%grav%x
            endif
        enddo
      endif
   endif
      !===================================================================================================================!
   ! Define parameterized bathymetry
   !===================================================================================================================!
   if (use_ptf == 1) then
    if (bc_infil == 1) then
   !===================================================================================================================!
   ! PTF from From Imhoff et al. 2020
   !===================================================================================================================!
   !===================================================================================================================!
   ! PTF from From Cooper et al. 2021
   !===================================================================================================================!
      do i = 1,size(phys_desc%soil)
        !infil%GA(i)%DeltaTheta = PTF( phys_desc%ptf_land(i) )%Kappa(1) &
        ! - PTF( phys_desc%ptf_land(i) )%Kappa(2) * phys_desc%soil(i)%clay &
        ! - PTF( phys_desc%ptf_land(i) )%Kappa(3) * phys_desc%soil(i)%sand !&
        !infil%GA(i)%PsiF = 0.01_rp * 10._rp ** (&
        ! PTF( phys_desc%ptf_land(i) )%Kappa(4) &
        ! - PTF( phys_desc%ptf_land(i) )%Kappa(5) * phys_desc%soil(i)%clay &
        ! - PTF( phys_desc%ptf_land(i) )%Kappa(6) * phys_desc%soil(i)%sand )
        infil%GA(i)%Ks = 25.4_rp / 3600._rp * 10._rp ** (&
                                - PTF( phys_desc%ptf_land(i) )%Kappa(7) &
                                - PTF( phys_desc%ptf_land(i) )%Kappa(8) * phys_desc%soil(i)%clay &
                                + PTF( phys_desc%ptf_land(i) )%Kappa(9) * phys_desc%soil(i)%sand ) &
                                / 1000._rp !kg/m2 to m
      ! Sigmoid parameterization
      enddo
    elseif (bc_infil == 2) then
        call Stopping_Program_Sub( "PTF for SCS-CN method not defined yet" )
    endif
   endif
   !===================================================================================================================!
   ! ZSinit.txt
   !===================================================================================================================!
    inquire( file = 'zs_init.txt' , exist = file_exist(1) )
    if ( file_exist(1) ) then
        do i = 1, mesh%nc
            dof0%h(i) = max(0._rp, dof0%h(i) - bathy_cell(i) &
                     + slope_y(1) * mesh%cell(i)%grav%y &
                     + slope_x(1) * mesh%cell(i)%grav%x)
        enddo
    endif
    !======================================================================================================================!
   ! Run the GR4 scheme for all hourly time steps
   !======================================================================================================================!
   !===================================================================================================================!
   ! Model Loop Time Initialization
   !===================================================================================================================!
   dof%h = dof0%h
   dof%u = dof0%u
   dof%v = dof0%v
   tc = tc0
   nt = nt0
    if ( use_obs == 1 ) then
        if ( use_Zobs == 1 ) then
            innovation(:)%ind_t = 1_ip
        end if
        if ( use_UVobs == 1 ) then
            innovUV(:)%ind_t = 1_ip
        end if
        if ( (use_Qobs == 1) .or. (use_Qobs_gr4 == 1) ) then
            innovQ(:)%ind_t = 1_ip
        endif
    endif
   if ( use_Qobs == 1 .or. use_Qobs_gr4 == 1) then
      innovQ(:)%ind_t = 1_ip
   endif
   cost = 0._rp
   !===================================================================================================================!
   ! Writing Initial Condition Output File
   !===================================================================================================================!
   call write_results( dof0 , mesh ) !NOADJ
   !===================================================================================================================!
   ! Initializing post treatment variables
   !===================================================================================================================!
   call sw_pre_treatment( dof0 , mesh ) !NOADJ
   !===================================================================================================================!
   ! SW Model Loop Time
   !===================================================================================================================!
   end_time_loop = .false.
   do while( .not. end_time_loop)
      call sub_run_model
   end do
   !===================================================================================================================!
   ! Cost Function Calculation using Innovation Vector
   !===================================================================================================================!
write(*,*) "call calc_cost_function( cost , mesh )"
   call calc_cost_function( cost , mesh )
CONTAINS
   SUBROUTINE sub_run_model
      implicit none
      sub_nt = 0
      do while ( .not. end_time_loop .and. sub_nt < max_nt_for_adjoint )
         !=============================================================================================================!
         ! Boundary Conditions
         !=============================================================================================================!
         call set_bc( dof , mesh )
         !=============================================================================================================!
         ! Time Control
         !=============================================================================================================!
         call advance_time( dof , mesh )
         !=============================================================================================================!
         ! Reinitialization of post variables
         !=============================================================================================================!
         bc%sum_mass_flux(:) = 0._rp
         !=============================================================================================================!
         ! Time Stepping Performing ( Euler + first_b1, IMEX + )
         !=============================================================================================================!
 select case( temp_scheme )
    case( 'euler' )
       select case( spatial_scheme )
      case( 'first_b1' )
       call euler_time_step_first_b1( dof , mesh )
      case default
       call Stopping_Program_Sub( 'Unknow spatial scheme' )
     end select
                              !<NOADJ
    case( 'imex' )
       select case( spatial_scheme )
      case( 'muscl_b1_b' )
       call imex_time_step( dof , mesh )
      case default
       call Stopping_Program_Sub( 'If temporal scheme is imex, spatial scheme MUST BE muscl_b1_b ' )
     end select
    case( 'low_froude' )
       call low_froude_time_step( dof , mesh )
                              !>NOADJ
    case default
       call Stopping_Program_Sub( 'Unknow temporal scheme' )
 end select
         !============================================================================!
         ! Post-processing
         !=============================================================================================================!
         call sw_post_treatment( dof , mesh )
         !=============================================================================================================!
         ! Filling Innovation Vector
         !=============================================================================================================!
    if ( use_obs == 1 ) then
         if ( use_Zobs == 1) then
            call calc_innovation( dof,mesh )
        endif
        if ( use_UVobs == 1 ) then
            call calc_innovUV( dof,mesh )
        endif
         if ( use_Qobs == 1 ) then
            call calc_innovQ( dof,mesh )
        endif
         if ( .not. use_obs == 1 ) then
            call update_cost_function( dof , cost )
         endif
         sub_nt = sub_nt + 1
    endif
         !=============================================================================================================!
         ! Writing Output Result File, now with innovUV
         !=============================================================================================================!
         call write_results( dof , mesh ) !NOADJ
   end do
   END SUBROUTINE sub_run_model
END SUBROUTINE run_model
