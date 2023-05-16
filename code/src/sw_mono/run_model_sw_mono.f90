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
!  Main Routine to run a Shallow-Water simulation
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE run_model( mesh , dof0 , dof , cost )

   USE m_common
   USE m_linear_algebra
   USE m_mesh
   USE m_mpi
   USE m_time_screen                                                                                              !NOADJ
   USE m_numeric
   USE m_model
   USE m_obs

#ifdef USE_HYDRO
    USE m_gr4
#endif

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(inout)  ::  mesh
   type( unk ), intent(in   )  ::  dof0
   type( unk ), intent(inout)  ::  dof

   real(rp), intent(out)  ::  cost

   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  sub_nt


    !======================================================================================================================!
   !  Run the GR4 scheme for all hourly time steps
   !======================================================================================================================!
#ifdef USE_HYDRO
        do catchnb = 1, bc%nb_gr4in
!write(*,*) catchnb, 'params', bc%gr4(catchnb)%params
!write(*,*) 'state 1 and 13 before warmup', bc%gr4(catchnb)%state(1), bc%gr4(catchnb)%state(13)
!             if (do_warmup) then
              call GR4_main("warmup")
!             endif
!write(*,*) 'state 1 and 13 after  warmup', bc%gr4(catchnb)%state(1), bc%gr4(catchnb)%state(13)
            call GR4_main("launch")
!write(*,*) 'state 1 and 13 after     run', bc%gr4(catchnb)%state(1), bc%gr4(catchnb)%state(13)
        enddo
#endif
   !===================================================================================================================!
   !  Model Loop Time Initialization
   !===================================================================================================================!

   dof%h  =  dof0%h
   dof%u  =  dof0%u
   dof%v  =  dof0%v

   tc  =  tc0
   nt  =  nt0

   if ( use_obs == 1 ) then
      innovation (:)%ind_t    =  1_ip
      innovW(:)%ind_t  =  1_ip
   end if

   if ( use_Qobs == 1 .or. use_Qobs_gr4 == 1) then
      innovQ(:)%ind_t  =  1_ip
   endif

   cost  =  0._rp

   !===================================================================================================================!
   !  Writing Initial Condition Output File
   !===================================================================================================================!
!write(*,*) "+ do  call write_results( dof0 , mesh )      "
   call write_results( dof0 , mesh )                                              !NOADJ
!write(*,*) "+ done  call write_results( dof0 , mesh )      "
   !===================================================================================================================!
   !  Initializing post treatment variables
   !===================================================================================================================!
!write(*,*) "+ do  call sw_pre_treatment( dof0 , mesh  )      "
   call sw_pre_treatment( dof0 , mesh )                                           !NOADJ
!write(*,*) "+ done call sw_pre_treatment( dof0 , mesh )      "

   !===================================================================================================================!
   !  SW Model Loop Time
   !===================================================================================================================!

   end_time_loop = .false.

   do while( .not. end_time_loop)
!write(*,*) "+ do call sub_run_model"
!write(*,*) 'w_obs=', w_obs

      call sub_run_model
!write(*,*) "+ done call sub_run_model"

   end do


   !===================================================================================================================!
   !  Cost Function Calculation using Innovation Vector
   !===================================================================================================================!
   call calc_cost_function( cost , mesh )

CONTAINS


   SUBROUTINE sub_run_model

      implicit none

      sub_nt = 0

! write(*,*) "++ ENTER WHILE LOOP: while ( .not. end_time_loop .and. sub_nt < max_nt_for_adjoint )"
      do while ( .not. end_time_loop .and. sub_nt < max_nt_for_adjoint )

         !=============================================================================================================!
         !  Boundary Conditions
         !=============================================================================================================!

!write(*,*) "--- do  call set_bc( dof , mesh )"
         call set_bc( dof , mesh )
!write(*,*) "--- donr call set_bc( dof , mesh )"

         !=============================================================================================================!
         !  Time Control
         !=============================================================================================================!

!write(*,*) "--- do  call advance_time( dof , mesh )"
         call advance_time( dof , mesh )
!write(*,*) "--- done  call advance_time( dof , mesh )"

         !=============================================================================================================!
         !  Reinitialization of post variables
         !=============================================================================================================!

         bc%sum_mass_flux(:) = 0._rp

         !=============================================================================================================!
         !  Time Stepping Performing ( Euler + first_b1, IMEX +  )
         !=============================================================================================================!

! write(*,*) "--- do  call appropriate scheme"
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

!write(*,*) "--- done  call appropriate scheme"
         !============================================================================!
         !  Post-processing
         !=============================================================================================================!
!write(*,*) ">>> do      call sw_post_treatment( dof , mesh )"
         call sw_post_treatment( dof , mesh )
!write(*,*) ">>> done      call sw_post_treatment( dof , mesh )"
         !=============================================================================================================!
         !  Writing Output Result File
         !=============================================================================================================!
!write(*,*) ">>> do      call write_results( dof , mesh )          "
         call write_results( dof , mesh )                                                                         !NOADJ
!write(*,*) ">>> done      call write_results( dof , mesh ) "

         !=============================================================================================================!
         !  Filling Innovation Vector
         !=============================================================================================================!

         if ( use_obs == 1 ) then
!write(*,*) ">>> do  call calc_innovation( dof,mesh )         "
            call calc_innovation( dof,mesh )
           ! call calc_innovW( dof,mesh )
!write(*,*) ">>> done  call calc_innovation( dof,mesh )         "
        endif

#ifdef USE_HYDRO
         if ( use_Qobs_gr4 == 1 ) then
!write(*,*) ">>> do  call calc_innovQ_gr4( dof,mesh )      "
!             write(*,*) 'Using flow observations at GR4 outputs.'
            call calc_innovQ_gr4( dof,mesh )
!write(*,*) ">>> done  call calc_innovQ_gr4( dof,mesh )      "
        endif
#endif

         if ( use_Qobs == 1 ) then
!write(*,*) ">>> do  call calc_innovQ( dof,mesh )     "
            !write(*,*) 'Calling flow observations at hydraulic station.t'
            call calc_innovQ( dof,mesh )
!write(*,*) ">>> done   call calc_innovQ( dof,mesh )      "
        endif

         if ( .not. use_obs == 1 .and. .not. use_Qobs_gr4 == 1 .and. .not. use_Qobs == 1) then
!write(*,*) ">>> do  call update_cost_function( dof , cost )     "
            call update_cost_function( dof , cost )
!write(*,*) ">>> done  call update_cost_function( dof , cost )     "
         endif

         sub_nt = sub_nt + 1

!write(*,*) "tc=", tc
!write(*,*) "ts=", ts

      end do
!write(*,*) ">>>> OUT sub_run_model"
   END SUBROUTINE sub_run_model


END SUBROUTINE run_model
