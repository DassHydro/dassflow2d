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
!> \file call_run_model.f90
!! \brief This file includes call_run_model module.
!! \details The file includes only call_run_model module (see doc call_run_model module).

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Module call_run_model : This file is for wrapper part.
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!>  Module call_run_model : This file is for wrapper part.
!! \details This module replace the previous main.f90, this is for wrapper part
MODULE call_model

   USE m_common
   USE m_linear_algebra
   USE m_mesh
   USE m_model
   USE m_linear_solver
   USE m_numeric
   USE m_obs

   #ifdef USE_ADJ
   USE m_adjoint
   USE m_tap_vars
   USE m_minimization
   #endif

   implicit none





!>  Type Model : contain global variables (not initialised in modules) of the Dassflow-2D
   TYPE Model

!>  Type Model : contain global variables (not initialised in modules) of the Dassflow-2D
!>   blablablal
	type(msh)  ::  mesh             !> mesh of the model

   	type(unk)  ::  dof0            !> initial conditions of the unknown of the probelm

   	type(unk)  ::  dof             !> unkown of the probelm

   	real(rp)   ::  cost            !> cost value


 	type(Input_Param) :: param     !> not used yet
	type(friction_data) :: my_friction !> friction parameterisation
	type(infiltration_data) :: my_infiltration !> infiltration parameterisation
!     type(greenampt)  :: my_greenampt

	type(input_data)  :: my_phys_desc !> Physical descriptors


	type(param_model) ::  my_param_model         !> Bathymetry at mesh cells gravity center

    type(bcs) :: my_bc

   END TYPE



CONTAINS

	!> \brief initialise the model
    !! \return mdl : Model type
   subroutine Model_initialise(mdl)

	implicit none

	type(Model), intent(out) :: mdl

	mdl%cost = 0._rp

   end subroutine Model_initialise

    !> \brief destroy the model
    !! \details destroy mesh, dof and dof0 (do not destroy cost)
   subroutine Model_finalise(mdl)

	implicit none

	type(Model), intent(inout) :: mdl

	! deallocate wrapped (interfaced variables)
	call unk_finalise(mdl%dof0)
	call unk_finalise(mdl%dof)
	call msh_finalise(mdl%mesh)

    call infiltration_finalise(mdl%my_infiltration)
    call phys_desc_finalise(mdl%my_phys_desc)
    call bc_finalise(mdl%my_bc)
	! deallocate fortran

   call dealloc_model()

   call dealloc_m_numeric()

   if(allocated(innovation)) deallocate(innovation)
   if(allocated(innovUV)) deallocate(innovUV)
   if(allocated(innovW))	deallocate(innovW)
   if(allocated(innovQ))	deallocate(innovQ)

   end subroutine Model_finalise


   subroutine init_solver(mdl)
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   !>FORTRAN DOCUMENTATION
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   	   !>init_solver
	   !>
	   !> initialize solver
	   !>
	   !> -------------------
	   !> details
	   !> -------------------
	   !>
     !> calls Init_Linear_Solver(model) from  m_numeric.f90

	implicit none

	!<input/ output variables
	type(Model), intent(inout) :: mdl

	!< ======================================================================================================================!
	!<   Machine precision limits numbers
	!< ======================================================================================================================!

    call Machine_Number_Limits

	call Print_Screen( 'number_limits' )

	!<======================================================================================================================!
	!<  Initialization of the linear solver (MUMPS, AGMG for the moment)
	!<======================================================================================================================!

	call Init_Linear_Solver(mdl%mesh)

   end subroutine init_solver

   subroutine init_friction(mdl)
	   !>> init_friction solver
	   !>call f90 routine friction_initialise
	   !>blablas
	   !>blablasS
	   !>blablassss
		implicit none

		! input/ output variables
		type(Model), intent(inout) :: mdl

		!<======================================================================================================================!
		!<  Model Variables Initialization (Using m_user_data.f90 file provided in the /bin directory)
		!<======================================================================================================================!

 		call friction_initialise(mdl%my_friction, mdl%mesh) ! can be called in the python

! 		call my_param_model_initialise(mdl%my_param_model, mdl%mesh) ! can be called in the python
! 		call my_bc_initialise(mdl%my_bc, bc, mdl%mesh,mdl%dof0)

   end subroutine init_friction


   subroutine init_infiltration(mdl)
	   !>> init_infiltration solver
	   !>call f90 routine infiltration_initialise,

		implicit none

		! input/ output variables
		type(Model), intent(inout) :: mdl

		!<======================================================================================================================!
		!<  Model Variables Initialization (Using m_user_data.f90 file provided in the /bin directory)
		!<======================================================================================================================!

 		call infiltration_initialise(mdl%my_infiltration, mdl%mesh) ! can be called in the python

   end subroutine init_infiltration


   subroutine init_phys_desc(mdl)
	   !>> init_phys_desc solver
	   !>call f90 routine phys_desc_initialise

		implicit none

		! input/ output variables
		type(Model), intent(inout) :: mdl

		!<======================================================================================================================!
		!<  Model Variables Initialization (Using m_user_data.f90 file provided in the /bin directory)
		!<======================================================================================================================!

        call phys_desc_initialise(mdl%my_phys_desc, mdl%mesh)

   end subroutine init_phys_desc


   subroutine init_bc(mdl)
	   !>> init_bc solver
	   !>call f90 routine bc_initialise

		implicit none

		! input/ output variables
		type(Model), intent(inout) :: mdl

		!<======================================================================================================================!
		!<  Model Variables Initialization (Using m_user_data.f90 file provided in the /bin directory)
		!<======================================================================================================================!

        call bc_initialise(mdl%my_bc, mdl%mesh)

   end subroutine init_bc


   subroutine init_fortran(mdl)
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   !>FORTRAN DOCUMENTATION
   !>++++++++++++++++++++++++++++++++++++++++++++++++
	   !>init_fortran
	   !>
	   !>Initialize fortran derived types
	   !>- bc
	   !>- land use
	   !>- eventual bc value (bc%hyd, bc%rat, bc%zpresc, bc%hpresc, etc...)
	   !>(manning, manning_beta, bathy_cell are ALREADY defined with mesh)
	   !>
	   !>-------------------
	   !>details
	   !>-------------------
	   !>
	   !>call not interfaced  f90 subroutines
	   !>-  Initial()      from initialization.f90
	   !>-  Init_Schemes() from m_numeric.f90
		implicit none
		! input/ output variables
		type(Model), intent(inout) :: mdl

		call Initial(mdl%dof0, mdl%mesh, mdl%my_friction, mdl%my_infiltration, &
		mdl%my_param_model, mdl%my_phys_desc, mdl%my_bc)

		!>======================================================================================================================!
		!>  Fill arrays usefull to save time computation for some schemes (MUSCL, Diamond Scheme, etc ... )
		!>======================================================================================================================!

		call Init_Schemes(mdl%mesh)   !  in m_numeric.f90, which is not interaced

!~ 		call write_mesh(mdl%mesh)
!~ 		call write_land_uses()
!~ 		call write_hydrograph()
!~ 		call write_bc()

   end subroutine init_fortran


!    subroutine Init_all
!    !>> Init_all solver
!    !> call f90 routine friction_initialise, infiltration_initialise,
!
! 	implicit none
!
!
! ! 	call init_solver(mdl)
! ! 	call init_friction(mdl)
! ! 	call init_fortran(mdl)
!
! 	end subroutine  Init_all


#ifdef USE_ADJ
    !> \brief  initialize "back" variables for adjoint model
    !! \details  initialize for adjoint model
    !! - allocate _back variables (both dof0_back and other control vector components, to be call only on time)
    !! - write control (at initialization)
 subroutine init_back(mdl)
 	type(Model), intent(inout) :: mdl

 	call alloc_back_vars(dof0_back, dof_back, mdl%mesh ) ! allocate _back variables (both dof0_back and other control vector components, to be call only on time)

     call write_control( mdl%dof0 , mdl%mesh )

     dim_all = sum(dim_vars_in_control)

 end subroutine init_back
#endif



   !======================================================================================================================!
   ! main routine
   !======================================================================================================================!

   subroutine run(mdl, arg)
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   !>FORTRAN DOCUMENTATION
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   !>run
   !>-------------------
   !>details
   !>-------------------
   !!run the model depending on arg value, the different cases are treated
   !!- 'direct' : direct simulation
   !!- 'testadj' : for debug ??
   !!- 'grad' : for returning grad_cost for sensibility analisys
   !!- 'min' : perform a minimization
   !!- else : f90wrap_abort

	implicit none

	! input/ output variables

	type(Model), intent(inout) :: mdl                  !> model instance
	character(len=*), optional, intent(in) :: arg      !> arguments precising the type of simulation we want
	character(len=32) :: arg_value                     !> arg variable treated for fortran reading

	! local variables

   integer(ip)  ::  num_points_display

	arg_value(:)=' '
	if (present(arg)) then
		arg_value = arg
	else
		arg_value(1:6) = 'direct'
	endif

	!<======================================================================================================================!
	!<  time loop
	!<======================================================================================================================!
	select case( trim(arg_value) )
		case('direct')
			call Time_Init(1_ip)

			call Print_Screen( 'start_direct' )

			call run_model(mdl%mesh, mdl%dof0, mdl%dof, mdl%cost)

            call Print_Screen( 'end_direct' )

			call Time_End(1_ip)

			call Time_Screen( mdl%mesh )
		#ifdef USE_ADJ
		case('testadj')

			call Print_Screen( 'start_testadj' )

            call test_adjoint( mdl%mesh , mdl%dof0 , mdl%dof , "0" )

            call Print_Screen( 'end_testadj' )

        case('grad')
			call Print_Screen( 'start_grad_cost' )

            call calc_grad_cost( mdl%mesh , mdl%dof0 , mdl%dof )

            call Print_Screen( 'end_grad_cost' )

        case('min')

			call Print_Screen( 'start_minimize' )

            call minimize_cost( mdl%mesh , mdl%dof0 , mdl%dof )

            call Print_Screen( 'end_minimize' )
		#endif

		case default
			call f90wrap_abort('unknown case')
	end select

   end subroutine run

    !> \brief
	subroutine clean_model(mdl)
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   !>FORTRAN DOCUMENTATION
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   !>run
   !>-------------------
   !>details
   !>-------------------
   !> Dealocate the model variable
   !> so that it will be safe to run another simulation
   !> WARNING: if not done, allocation fault will happen

		implicit none

		type(Model), intent(inout) :: mdl

		! deallocate the
		call Model_finalise(mdl)

	end subroutine

	!======================================================================================================================!
	! some routines for outputs (gnuplot, vtk, tecplot)
    !======================================================================================================================!

	!> \brief write txt file for gnuplote
	subroutine output_gnu(mdl, dof, filename)
		implicit none

		type(Model), intent(in) :: mdl
		type(unk), intent(in) :: dof
		character(len=128), intent(in) :: filename

		call v_gnuplot(dof, mdl%mesh, trim(filename))

	end subroutine output_gnu

		!> \brief write vtk file for paraview
	subroutine output_vtk(mdl, dof, filename)
		implicit none

		type(Model), intent(in) :: mdl
		type(unk), intent(in) :: dof
		character(len=128), intent(in) :: filename

		call v_vtk(dof, mdl%mesh, trim(filename))

	end subroutine output_vtk
    !> \brief write tec file for texplot ?
	subroutine output_tec(mdl, dof, filename)
		implicit none

		type(Model), intent(in) :: mdl
		type(unk), intent(in) :: dof
		character(len=128), intent(in) :: filename

		call v_tecplot(dof, mdl%mesh, trim(filename))

	end subroutine output_tec





!==================================================================================================================!
! initialise and finalize  friction
!==================================================================================================================!

    SUBROUTINE friction_initialise(my_friction, mesh)
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   !>FORTRAN DOCUMENTATION
   !>++++++++++++++++++++++++++++++++++++++++++++++++
    !> define my_friction fortran type (linked to mdl%my_friction)
    !> allocate my_friction%xxx and set values from text file

    implicit none
    !  Global Variables
    type(msh), intent(in)  ::  mesh
    type(friction_data), intent(inout) :: my_friction
    !  Local Variables
    real(rp)  ::  tmp                                     !
    integer(ip) :: mesh_total_cells
    
   inquire( file = 'land_uses.txt' , exist = file_exist(1) )
   if ( file_exist(1) ) write(*,*) "WARNING: you are trying to allocate manning from land_uses.txt and from init_friction"
    
   if(mesh_type=="dassflow") then

      mesh_total_cells = mesh%nc
      call mpi_sum_i( mesh_total_cells )

!===================================================================================================================!
! set the correspondance between cell and land value
!===================================================================================================================!

      allocate( my_friction%land( mesh_total_cells ) )

      open(10,file=mesh%file_name ,status='old')
   ! skip rows to access wished info
      read(10,*) ! comment line (# Simple channel mesh)
      read(10,*) ! nc, nn, ne line
      read(10,*)  ! comment line (# node)
      do i = 1,mesh%nn
         read(10,*) ! node coordinates
      end do
      read(10,*) ! comment line (# cells)

      do i = 1,mesh_total_cells

         read(10,*)   k , &
                   tmp , & !node1 but useless here
                   tmp , & !node2 but useless here
                   tmp , & !node3 but useless here
                   tmp , & !node4 but useless here
                   my_friction%land(k) , &
                   tmp   ! bathymetry value but useless here

      enddo

   close(10)

      !my_friction%nland = !maxval(my_friction%land(:))
      allocate( my_friction%manning( my_friction%nland ) )
      allocate( my_friction%manning_beta( my_friction%nland ) )

      my_friction%manning(:) = 0._rp
      my_friction%manning_beta(:) = 0._rp

   else
        write(*,*) "only 'dassflow' mesh_type formats is implemented"
   endif


  END SUBROUTINE friction_initialise

  SUBROUTINE friction_finalise(my_friction)
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   !>FORTRAN DOCUMENTATION
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   !> deallocate my_friction%xxx derived type
      implicit none
      !  Global Variables
      type(friction_data), intent(inout)  ::  my_friction

      if (allocated(my_friction%land))         deallocate(my_friction%land)
      if (allocated(my_friction%manning_beta)) deallocate(my_friction%manning_beta)
      if (allocated(my_friction%manning))      deallocate(my_friction%manning)

  END SUBROUTINE friction_finalise


SUBROUTINE infiltration_initialise(my_infiltration, mesh)
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   !>FORTRAN DOCUMENTATION
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   ! initialise friction (allocate and set values from file)
   implicit none
   !  Global Variables
   type(msh), intent(in)  ::  mesh
   type(infiltration_data), intent(inout) :: my_infiltration
   !  Local Variables
   integer(ip) :: mesh_total_cells

   inquire( file = 'land_uses_GA.txt' , exist = file_exist(1) )
   if ( file_exist(1) ) write(*,*) "WARNING: you are trying to allocate infiltration from land_uses_GA.txt and from init_infiltration"
   inquire( file = 'land_uses_SCS.txt' , exist = file_exist(1) )
   if ( file_exist(1) ) write(*,*) "WARNING: you are trying to allocate infiltration from land_uses_SCS.txt and from init_infiltration"
   
   if (mesh_type=='dassflow') then

      mesh_total_cells = mesh%nc
      call mpi_sum_i( mesh_total_cells )

   !===================================================================================================================!
   ! set the correspondance between cell and land value
   !===================================================================================================================!

      if (bc_infil == 1) then

        allocate( my_infiltration%land( mesh_total_cells ) )
        #allocate( my_infiltration%coord( 4, my_infiltration%nland )   )
        allocate( my_infiltration%GA( my_infiltration%nland ) )
        allocate( my_infiltration%SCS( 1 ) )

      elseif ( bc_infil == 2 ) then

        allocate( my_infiltration%land( mesh_total_cells ) )
        allocate( my_infiltration%coord( 4, my_infiltration%nland )   )
        allocate( my_infiltration%SCS( my_infiltration%nland ) )
        allocate( my_infiltration%GA( 1 ) )

      else

        allocate( my_infiltration%land ( 1 ) )
        allocate( my_infiltration%coord( 4, 1 )   )
        allocate( my_infiltration%SCS  ( 1 ) )
        allocate( my_infiltration%GA   ( 1 ) )

      endif


      my_infiltration%GA( : )%ks = -1._rp
      my_infiltration%GA( : )%psif = -1._rp
      my_infiltration%GA( : )%deltatheta = -1._rp
      my_infiltration%SCS( : )%lambdacn = -1._rp
      my_infiltration%SCS( : )%CN = -1._rp

    endif

!     do i = 1, mesh%nc
!       do j  = 1,my_infiltration%nland
!         if (mesh%cell(i)%grav%x < my_infiltration%coord(1,j) .and. mesh%cell(i)%grav%x > my_infiltration%coord(2,j) .and.  mesh%cell(i)%grav%y < my_infiltration%coord(3,j) .and. mesh%cell(i)%grav%y > my_infiltration%coord(4,j)) then
!             my_infiltration%land(i) = j
!       enddo
!     enddo

  END SUBROUTINE infiltration_initialise

  SUBROUTINE infiltration_finalise(my_infiltration)

      implicit none
      !  Global Variables
      type(infiltration_data), intent(inout)  ::  my_infiltration

      if (allocated(my_infiltration%land)) deallocate(my_infiltration%land)
      if (allocated(my_infiltration%GA))   deallocate(my_infiltration%GA)
      if (allocated(my_infiltration%SCS))  deallocate(my_infiltration%SCS)
      if (allocated(my_infiltration%coord)) deallocate(my_infiltration%coord)

  END SUBROUTINE infiltration_finalise


  SUBROUTINE phys_desc_initialise(my_phys_desc, mesh)

      implicit none
      !  Global Variables
      type(input_data), intent(inout)  ::  my_phys_desc
      type(msh), intent(in)  ::  mesh
      !  Local Variables
      integer(ip) :: mesh_total_cells

      mesh_total_cells = mesh%nc
      call mpi_sum_i( mesh_total_cells )

      allocate(my_phys_desc%soil_land(mesh_total_cells))!mesh%nc))
      allocate(my_phys_desc%soil(my_phys_desc%soil_nland))!my_phys_desc%soil_nland))
      my_phys_desc%soil(:)%clay = 0._rp
      my_phys_desc%soil(:)%silt = 0._rp
      my_phys_desc%soil(:)%sand = 0._rp

      allocate(my_phys_desc%surf_land(mesh_total_cells))!mesh%nc))
      allocate(my_phys_desc%surf(my_phys_desc%surf_nland))!my_phys_desc%surf_nland))
      my_phys_desc%surf(:)%imperm = 0._rp
      my_phys_desc%surf(:)%Dmax = 0._rp

      allocate(my_phys_desc%struct_land(mesh_total_cells))!mesh%nc))
      allocate(my_phys_desc%structures(my_phys_desc%struct_nland))!my_phys_desc%struct_nland))
      my_phys_desc%structures(:)%C1 = 0._rp
      my_phys_desc%structures(:)%C2 = 0._rp
      my_phys_desc%structures(:)%C3 = 0._rp
      my_phys_desc%structures(:)%true_x = 0._rp
      my_phys_desc%structures(:)%true_y = 0._rp

  END SUBROUTINE phys_desc_initialise

  SUBROUTINE phys_desc_finalise(my_phys_desc)

      implicit none
      !  Global Variables
      type(input_data), intent(inout)  ::  my_phys_desc

      if (allocated(my_phys_desc%soil)) deallocate(my_phys_desc%soil)
      if (allocated(my_phys_desc%soil_land)) deallocate(my_phys_desc%soil_land)

      if (allocated(my_phys_desc%surf))   deallocate(my_phys_desc%surf)
      if (allocated(my_phys_desc%surf_land))   deallocate(my_phys_desc%surf_land)

      if (allocated(my_phys_desc%structures))  deallocate(my_phys_desc%structures)
      if (allocated(my_phys_desc%struct_land))  deallocate(my_phys_desc%struct_land)

  END SUBROUTINE phys_desc_finalise


  SUBROUTINE bc_initialise(my_bc, mesh)

      implicit none
      !  Global Variables
      type(bcs), intent(inout)  ::  my_bc
      type(msh), intent(inout)  ::  mesh
      !  Local Variables
      integer(ip) :: mesh_total_cells

      
      inquire( file = 'rain.txt' , exist = file_exist(1) )
       if ( file_exist(1) ) write(*,*) "WARNING: you are trying to allocate rain from rain.txt and from init_bc"
   
      mesh_total_cells = mesh%nc
      call mpi_sum_i( mesh_total_cells )


      allocate(my_bc%rain(my_bc%nb_rn))
      allocate(my_bc%rain_land(mesh_total_cells))
      my_bc%rain_land(:) = 0._ip

      do i = 1, my_bc%nb_rn

        allocate(my_bc%rain(i)%t(my_bc%nb_rn_t))
        allocate(my_bc%rain(i)%q(my_bc%nb_rn_t))
      enddo


  END SUBROUTINE bc_initialise

  SUBROUTINE bc_finalise(my_bc)

    implicit none
    !   Global Variables
    type(bcs), intent(inout)  ::  my_bc

    if (allocated(my_bc%rain)) deallocate(my_bc%rain)
    if (allocated(my_bc%rain_land)) deallocate(my_bc%rain_land)

  END SUBROUTINE bc_finalise

 	subroutine func(mdl,ctrl_in, cost_func, grad_func)
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   !>FORTRAN DOCUMENTATION
   !>++++++++++++++++++++++++++++++++++++++++++++++++
   !> func
   !>compute the cost function and the gradient of cost function
   !>-------------------
   !>details
   !>-------------------
    !> set control = ctrl_in, and compute the cost function and the gradient of cost function
	!> dim_all is defined in m_adjoint, is calculated once the control vector is allocated (write_control)
	!> \param[inout] mdl : Model, the model instance
    !> \param[in] ctrl_in : real(rp), , the control vector to define as new control vector
    !> \param[out] cost_func : real(rp),dimension(1), the curent value of cost funtion (cost in fortran, cost_func in python)
    !> \param[out] grad_func : real(rp),dimension(control),the curent value of gradient of cost funtion  (control_back in fortran, grad_func in python)

 		implicit none
 		type(Model), intent(inout) :: mdl
 		real(rp), dimension(dim_all), intent(in) :: ctrl_in

 		real(rp), intent(out) :: cost_func
 		real(rp), dimension(dim_all), intent(out) :: grad_func

 		! set the input control vector

 		control(:) = ctrl_in(:)

 		cost_back = one
 		verbose = -1

 		! call adjoint model to compute cost function and gradient of cost function
 		call adjoint_model(mdl%mesh, mdl%dof0, mdl%dof)

 		! save the current value of cost function and the gradient as outputs
 		cost_func = cost
 		grad_func(:) = control_back(:)

 	end subroutine func


!=============================================!
! additionnal write functions
!=============================================!

! > write_mesh()
!! very similar to read_dass_mesh, but write the msh object in the 'inhouse dassflow' convention
!! use  mesh and :
! 	-land
! 	-bathy_cell
!	-bc
SUBROUTINE write_mesh( mesh)

   USE m_common
   USE m_mesh
   USE m_model
   USE m_time_screen
   USE m_mpi

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in)  ::  mesh

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!
! To check size of inlet and outlet bcs manualy
   integer(ip)  ::  nb_bc_in, nb_cell_in
   integer(ip)  ::  nb_bc_out, nb_cell_out

   ! logical, used to check if the boundary cell is a boundary cell a the given bc
   logical(ip)  ::  is_boundary_cell

    !>>>>>>>>> boundary_edge_index is the corresponding index to global_edge_index.
    ! boundary_edge_index store the corresponding id in the object mesh%edge
    integer(ip)  :: boundary_edge_index
    ! global_edge_index store the corresponding id in the object mesh%edgeb
    integer(ip)  :: global_edge_index
    !>>>>>>>>> to store the index of a cell identified
    integer(ip)  :: cell_index
    ! >>>>>>>> store the id of the boundary edge related to cell
    integer(ip) :: edge_incell_index

    character(51)  ::  tmp_mesh_name ! automaticaly_generated_mesh =
   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

    tmp_mesh_name = "res/generated_files/automaticaly_generated_mesh.txt"

    call system('mkdir -p res/generated_files')



   open(10,file=tmp_mesh_name,status='replace',form='formatted')

   !===================================================================================================================!
   !  Writing header
   !===================================================================================================================!

   write(10,fmt='(A)') "# Generated mesh with subroutine write_mesh() ||| number of nodes |   number of cells | mesh scale == 0 always"
   write(10,'(I10, I10, A)') mesh%nn , mesh%nc , '   0'
   !===================================================================================================================!
   !  Writing Nodes
   !===================================================================================================================!

   write(10,fmt='(A)') "# Nodes||| id node, x coord, y coord, bathy (x,y) "

   do i = 1,mesh%nn
      write (10,*)  i , mesh%node(i)%coord%x , mesh%node(i)%coord%y , 0
   end do

   !===================================================================================================================!
   !  Writing Nodes connectivity
   !===================================================================================================================!

     write (10,*) '# Cells|||  id of cell, id node 1, id node 2, id  node 3,  id node 4 , land_type, bathymetry'

   do k = 1,mesh%nc

    write (10,*)   k , &
                   mesh%cell(k)%node(1) , &
                   mesh%cell(k)%node(2) , &
                   mesh%cell(k)%node(3) , &
                   mesh%cell(k)%node(4) , &
                   land(k) , &
                   bathy_cell(k)


   end do


   !===================================================================================================================!
   !  Prepare necessary inlet and outlet informations
   !===================================================================================================================!


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     count the  number of boundaries defined as inlet or  outlet
nb_bc_in = 0
nb_bc_out = 0

do k = 1,bc%nb
  if (bc%typ(k,1 ) == 'discharg1'  .or. bc%typ(k,1 ) == 'discharg2' ) then
    nb_bc_in = nb_bc_in +1
  else if(bc%typ(k,1 )(1:8) =='ratcurve' .or. bc%typ(k,1 )(1:6) == 'transm' .or. bc%typ(k,1)(1:6) =='zspresc' .or.  bc%typ(k,1 )(1:6) == 'hpresc') then
    nb_bc_out = nb_bc_out +1

  end if
end do


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    count the  number of cells defined on a boundary of type  inlet or  outlet
nb_cell_in = 0
nb_cell_out = 0

! A) LOOP ON ALL BOUNDARY EDGES  in mesh object
    do boundary_edge_index = 1,mesh%neb

    ! get corresponding index in object mesh%edge
        global_edge_index =  mesh%edgeb(boundary_edge_index)%ind
!         write(*,*) mesh%edgeb( boundary_edge_index  )%typlim(1:9)

    ! check if this boundary typlim is one of discharg1 or dischard2
        if (mesh%edgeb( boundary_edge_index  )%typlim(1:9) == 'discharg1'  .or. mesh%edgeb( boundary_edge_index   )%typlim(1:9) == 'discharg2' ) then

        nb_cell_in = nb_cell_in + 1

        else if(mesh%edgeb( boundary_edge_index  )%typlim(1:8) ==  'ratcurve' .or. mesh%edgeb( boundary_edge_index   )%typlim(1:6) == 'transm' &
          .or. mesh%edgeb( boundary_edge_index   )%typlim(1:7) =='zspresc' .or.  mesh%edgeb( boundary_edge_index   )%typlim(1:6) == 'hpresc') then

        nb_cell_out = nb_cell_out +1

        end if
    end do


   !!===================================================================================================================!
   !! Write inlet INLETS
   !!===================================================================================================================!
write(10,*) '# Boundaries'

write(10,*) 'INLET', nb_cell_in, nb_bc_in



! A) Loop on all the defined boundaries in bc object
do k = 1,bc%nb
! A) LOOP ON ALL BOUNDARY EDGES  in mesh object
    do boundary_edge_index = 1,mesh%neb

    ! get corresponding index in object mesh%edge
        global_edge_index =  mesh%edgeb(boundary_edge_index)%ind

    ! check if this boundary typlim is one of discharg1 or dischard2
        if (mesh%edgeb( boundary_edge_index  )%typlim == 'discharg1'  .or. mesh%edgeb( boundary_edge_index   )%typlim == 'discharg2' ) then

    ! check if this edge of of correct grogp
        if (mesh%edgeb( boundary_edge_index  )%group == bc%grpf(k) ) then
        cell_index =  mesh%edge(global_edge_index)%cell(1)
                do i = 1, size(mesh%cell(cell_index)%edge)
                if(mesh%cell(cell_index)%edge(i) == global_edge_index) then
                    edge_incell_index = i
                endif
            end do
    ! write gathered data                                      ! mesh%edge(global_edge_index)%cell(2) is the id of the ghost cell
        write(10,*) cell_index, edge_incell_index, bc%grpf(k) ,  bathy_cell(mesh%edge(global_edge_index)%cell(2)), bc%grpf(k)

        end if
        end if
    end do


end do

   !!===================================================================================================================!
   !! Write Outlet
   !!===================================================================================================================!

write(10,*) 'OUTLET', nb_cell_out, nb_bc_out


! A) Loop on all the defined boundaries in bc object
do k = 1,bc%nb

! A) LOOP ON ALL BOUNDARY EDGES  in mesh object
    do boundary_edge_index = 1,mesh%neb

    ! get corresponding index in object mesh%edge
        global_edge_index =  mesh%edgeb(boundary_edge_index)%ind
       !     write(*,*) mesh%edgeb( boundary_edge_index  )%group
       !     write(*,*) mesh%edgeb( boundary_edge_index  )%typlim

    ! check if this boundary typlim is one of discharg1 or dischard2
        if (  mesh%edgeb( boundary_edge_index  )%typlim(1:8) ==  'ratcurve' &
        .or. mesh%edgeb( boundary_edge_index  )%typlim(1:6) == 'transm' &
        .or. mesh%edgeb( boundary_edge_index   )%typlim(1:7) == 'zspresc' &
        .or. mesh%edgeb( boundary_edge_index   )%typlim(1:6) == 'hpresc' ) then
    ! check if this edge of of correct grogp
        if (mesh%edgeb( boundary_edge_index  )%group == bc%grpf(k) ) then
        cell_index =  mesh%edge(global_edge_index)%cell(1)

    !       identify id of edge in mesh%cell%edge corresponding to global_edge index
            do i = 1, size(mesh%cell(cell_index)%edge)
                if(mesh%cell(cell_index)%edge(i) == global_edge_index) then
                    edge_incell_index = i
                endif
            end do
    ! write gathered data                                                                    ! mesh%edge(global_edge_index)%cell(2) is the id of the ghost cell
        write(10,*) cell_index, edge_incell_index, bc%grpf(k) ,  bathy_cell(mesh%edge(global_edge_index)%cell(2)), bc%grpf(k)

        end if
        end if
    end do


end do
close(10)

END SUBROUTINE write_mesh



SUBROUTINE write_land_uses( )
   USE m_common ! necessary ?
   USE m_model
   implicit none


    open(10,file="res/generated_files/land_uses.txt",status='replace',form='formatted')

    write(10,*) "!=================================================================================================!"
    write(10,*) "!	Number	of	Land	Uses				                                                        "
    write(10,*) "!=================================================================================================!"
    write(10,fmt = '(I10)') nland
    write(10,*) "!=================================================================================================!"
    write(10,*) "!	List	of	Land	Uses, for each patch, the values defined in the file are : 				                                                        "
    write(10,*) "!====================== id patch, manning, manning beta ==========================================!"

    !define values for each paches
      do i = 1,nland
      write(10,fmt = "(I10, 2ES15.7)") i , manning(i), manning_beta(i)
      end do
    close(10)


END SUBROUTINE write_land_uses

   !===================================================================================================================!
   !===================================================================================================================!
   ! # routines to ckeck INITIALization
   !===================================================================================================================!
   !===================================================================================================================!
SUBROUTINE write_bc()
   USE m_model

   implicit none
   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!
    ! all information needed is in object bc.
    ! the type bcs is defined in module m_model
    ! bcs is initialized in either
    !   - read_dass_mesh() --> call read_bc_file()
    !   - Initial()
    ! most subroutines modifying and using bc information are located in in src/sw_mono/boundary.f90 file

   !===================================================================================================================!
   !  Local variables
   !===================================================================================================================!
    character(len=50) ::  filename
    integer :: bc_nb ! the number of inflow + number of outflow

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   !calculate bc nb
   bc_nb = 0
      do i = 1,bc%nb
        if(bc%typ(i,1)(1:8)=='discharg' .or. bc%typ(i,1)(1:8)  ==  'ratcurve' &
        .or. bc%typ(i,1)(1:7) == 'zspresc' .or.bc%typ(i,1)(1:6) == 'hpresc' &
        .or. bc%typ(i,1)(1:6) == 'transm') then
            bc_nb = bc_nb +1
        end if
      end do

   ! if directory where we wrtie files not exist, create if
   call system('mkdir -p res/generated_files/')

   ! define name of the file
   filename='res/generated_files/bc.txt'

   ! write information in bc.txt located in res/bc/
   open(20,file=filename,status='replace',form='formatted')

   write(20,*) '!===============================================!'
   write(20,*) '!  Number of boundary conditions'
   write(20,*) '!===============================================!'

   write(20,*) bc_nb

   write(20,*) '!===============================================!'
   write(20,*) '!  List of boundary conditions (all unspecified boundaries are set as "wall" boundary condition type ) '
   write(20,*) '!===============================================!'

!    allocate( bc%typ ( bc%nb , 2 ) )
!    allocate( bc%grpf( bc%nb     ) )
!    bc%typ (:,:)  =  ''


   do i = 1,bc%nb
   ! discharge1 adn sicharge2 type treatment
   if(bc%typ(i,1)(1:8)=='discharg' ) then
   !  ================= WRITE()'s 'fmt' argument explanation :   =================
   !  type   I = integer, A = character
   !  the "sufix" values (3,4,9) represent the number of space allocated to the write function (ie, the first "sufix" values are going to be printed")
   !  the  "prefix" values indicates teh number of repetition of this write format
   !                                  I3 ,  A4,   A9      ,   2A4
    write(20,fmt='(I3,A4, A9, 2A4)') bc%grpf(i), '    ',  bc%typ(i,1),'    ', 'file'
   else if(bc%typ(i,1)(1:8)  ==  'ratcurve') then
       write(20,fmt='(I3,A4, A9, 2A4)') bc%grpf(i), '    ',  bc%typ(i,1),'    ', 'file'
   else if(bc%typ(i,1)(1:6) == 'transm') then
       write(20,fmt='(I3,A4, A9, 2A4)') bc%grpf(i), '    ',  bc%typ(i,1),'    ', '    '
   else if(bc%typ(i,1)(1:7) == 'zspresc') then
       write(20,fmt='(I3,A4, A9, 2A4)') bc%grpf(i), '    ',  bc%typ(i,1),'    ', '    '
   else if(bc%typ(i,1)(1:6) == 'hpresc') then
       write(20,fmt='(I3,A4, A9, 2A4)') bc%grpf(i), '    ',  bc%typ(i,1),'    ', '    '
   else
   ! wall conditions dont need to be specified
       write(20,fmt='(A5, I3, A4, A9)') '!    ', bc%grpf(i), '    ',  bc%typ(i,1)
   end if
   end do

   write(20,*) '!===============================================!'
   write(20,*) '!   NOTES                                       !'
   write(20,*) '!===============================================!'
   write(20,*) '!  File generated using subroutine write_bc()    '
   write(20,*) '!======= integer values of bc'
   write(20,*) '!bc%nb=', bc%nb
   write(20,*) '!bc%nb_in=', bc%nb_in
   write(20,*) '!bc%nb_out=', bc%nb_out
   write(20,fmt='(100A13)') '!bc%typ(:,1)=', bc%typ(:,1)
   write(20,fmt='(100A13)') '!bc%typ(:,2)=', bc%typ(:,2)
   write(20,fmt='(A13,100I10)') '!bc%grpf(:)=', bc%grpf(:)
   close(20)


END SUBROUTINE write_bc


!> very similar to gen_hydrograph
!! the difference resides in the fact that in this function, the data already available in bc is used to generate the hydrograph (instead of using m_used_data routines)
!! this function serve the purpose to check what is read by the model at any time
SUBROUTINE write_hydrograph()
   USE m_model

   implicit none
   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!
    ! all information needed is in object bc.
    ! the type bcs is defined in module m_model
    ! bcs is initialized in either
    !   - read_dass_mesh() --> call read_bc_file()
    !   - Initial()
    ! most subroutines modifying and using bc information are located in in src/sw_mono/boundary.f90 file

   !===================================================================================================================!
   !  Local variables
   !===================================================================================================================!
    character(len=50) ::  filename
    integer :: nb_hydrograph ! numberof hydrograf
    integer :: size_hyd ! size of the hydrograph (number of rows)
    nb_hydrograph = size(bc%hyd)

    ! define name of the file
   filename='res/generated_files/hydrograph.txt'

      if ( proc  == 0 ) then
   ! write information in hydrograph.txt located in res/bc/hydrograph.txt
   open(10,file=filename,status='replace',form='formatted')


         write(10,'(A)') '!===============================================!'
         write(10,'(A)') '!  Number of hydrograph'
         write(10,'(A)') '!===============================================!'

        write(10,'(I1)') nb_hydrograph

        do k = 1,nb_hydrograph

         size_hyd = size(bc%hyd( k )%t( : ))
         write(10,'(A)') '!===============================================!'
         write(10,'(A)') '!  Hydrograph'
         write(10,'(A)') '!===============================================!'

         write(10,'(I10)') size_hyd

         do i = 1,size_hyd
           write(10,'(2ES15.7)') bc%hyd( k )%t( i ) , bc%hyd( k )%q( i )
         end do ! END  do i = 1,size_hyd
     end do ! END do k = 1,nb_hydrograph


   write(10,*) '!===============================================!'
   write(10,*) '!   NOTES                                       !'
   write(10,*) '!===============================================!'
   write(10,*) '!  File generated using subroutine write_hydrograph()    '

         close(10)

      end if !  if ( proc  == 0 )

END SUBROUTINE write_hydrograph

! ---------------------------------------------------------------------- !
! ENABLE COPY FOR DERIVED TYPES
! ---------------------------------------------------------------------- !

!> Allow user to make a copy of the derived type bcs
!> @param i bcs derived type to copy
!> @param o bcs derived type copy

  SUBROUTINE BOUNDARIES_COPY(i, o)
    IMPLICIT NONE
    TYPE(bcs), INTENT(IN) :: i
    TYPE(bcs), INTENT(OUT) :: o
    o = i
  END SUBROUTINE BOUNDARIES_COPY



  !> Allow user to make a copy of the derived type dof
  !> @param i bcs derived type to copy
  !> @param o bcs derived type copy

  SUBROUTINE DOF_COPY(i, o)
    IMPLICIT NONE
    TYPE(unk), INTENT(IN) :: i
    TYPE(unk), INTENT(OUT) :: o
    o = i
  END SUBROUTINE DOF_COPY


  !> Allow user to make a copy of the derived type dof
  !> @param i bcs derived type to copy
  !> @param o bcs derived type copy

  SUBROUTINE MESH_COPY(i, o)
    IMPLICIT NONE
    TYPE(msh), INTENT(IN) :: i
    TYPE(msh), INTENT(OUT) :: o
    o = i
  END SUBROUTINE MESH_COPY

  ! ---------------------------------------------------------------------- !
  ! ENABLE REALLOCATION OF  DERIVED TYPES
  ! ---------------------------------------------------------------------- !

 ! reallocate m_model.manning and m_model.manning_beta global variables
  SUBROUTINE reallocate_manning(new_size)!, manning)     !size_in, size_out, manning_out )
        USE m_model
        IMPLICIT NONE
        integer, INTENT(in):: new_size
        !real(rp), dimension(:), INTENT(INOUT), allocatable :: manning

          if (allocated(manning)) then
            deallocate(manning)
          end if
          if (allocated(manning_beta)) then
            deallocate(manning_beta)
          end if

          allocate( manning( new_size ) )
          allocate( manning_beta( new_size ) )

  END SUBROUTINE reallocate_manning

END MODULE call_model
