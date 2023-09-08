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


MODULE m_tap_vars

   USE m_common
   USE m_obs

#ifdef USE_SW_MONO
      USE m_model
#endif

#ifdef USE_HYDRO
!       USE m_gr4
#endif

   real(rp) ::  dt_diff , tc_diff
   real(rp) ::  dt_back , tc_back

   type(bcs), target  ::  bc_diff
   type(bcs), target  ::  bc_back

   type(xsshp), dimension(:), allocatable, target  ::  XSshape_back
   type(xsshp), dimension(:), allocatable, target  ::  XSshape_diff

   type(infiltration_data), target  ::  infil_diff
   type(infiltration_data), target  ::  infil_back

   type(ptf_data), dimension(:), allocatable, target  ::  ptf_diff
   type(ptf_data), dimension(:), allocatable, target  ::  ptf_back
   
   type(innovation_obs), dimension(:), allocatable, target  ::  innovation_diff
   type(innovation_obs), dimension(:), allocatable, target  ::  innovation_back
   type(innovation_obs), dimension(:), allocatable, target  ::  innovW_diff
   type(innovation_obs), dimension(:), allocatable, target  ::  innovW_back
   type(innovation_obs), dimension(:), allocatable, target  ::  innovUV_diff
   type(innovation_obs), dimension(:), allocatable, target  ::  innovUV_back
   type(innovation_obs ), dimension(:), allocatable  ::  innovQ_diff
   type(innovation_obs ), dimension(:), allocatable  ::  innovQ_back


#ifdef USE_HYDRO
!    real(rp), dimension(13)  ::  State_diff , State_back
#endif

#ifdef USE_SW_MONO
      real(rp), dimension(:), allocatable  ::  manning_diff , manning_beta_diff, bathy_cell_diff
      real(rp), dimension(:), allocatable  ::  manning_back , manning_beta_back, bathy_cell_back
      real(rp), dimension(:), allocatable  ::  slope_y_back, slope_y_diff
      real(rp), dimension(:), allocatable  ::  slope_x_back, slope_x_diff
#endif


CONTAINS


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Allocation of needed variables to run the linear tangent model
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE alloc_diff_vars( dof0_diff , dof_diff , mesh )

      type( unk ), intent(inout)  ::  dof0_diff , dof_diff

      type( msh ), intent(in)  ::  mesh

      call alloc_dof( dof0_diff , mesh )
      call alloc_dof( dof_diff  , mesh )

      allocate( bc_diff%inflow ( size( bc%inflow  ) ) )
      allocate( bc_diff%outflow( size( bc%outflow ) ) )

      bc_diff%inflow (:)  =  0._rp
      bc_diff%outflow(:)  =  0._rp

#ifdef USE_HYDRO
      allocate( bc_diff%gr4( bc%nb_gr4in ) )

      do i = 1,bc%nb_gr4in

         allocate( bc_diff%gr4( i )%t( size( bc%gr4( i )%t ) ) )
         allocate( bc_diff%gr4( i )%Q( size( bc%gr4( i )%Q ) ) )

         bc_diff%gr4( i )%t(:)  =  0._rp
         bc_diff%gr4( i )%Q(:)  =  0._rp
         bc_diff%gr4( i )%params(:)  =  0._rp
         bc_diff%gr4( i )%state(:)  =  0._rp

      end do
#endif

      allocate( bc_diff%hyd( bc%nb_in ) )

      do i = 1,bc%nb_in

         allocate( bc_diff%hyd( i )%t( size( bc%hyd( i )%t ) ) )
         allocate( bc_diff%hyd( i )%q( size( bc%hyd( i )%q ) ) )

         bc_diff%hyd( i )%t(:)  =  0._rp
         bc_diff%hyd( i )%q(:)  =  0._rp

      end do


      allocate( bc_diff%rat( bc%nb_out ) )

      do i = 1,bc%nb_out
         allocate( bc_diff%rat( i )%h( size( bc%rat( i )%h ) ) )
         allocate( bc_diff%rat( i )%q( size( bc%rat( i )%q ) ) )

         bc_diff%rat( i )%h(:)  =  0._rp
         bc_diff%rat( i )%q(:)  =  0._rp

      end do

      do i = 1,bc%nb_rn

         allocate( bc_diff%rain( i )%t( size( bc%rain( i )%t ) ) )
         allocate( bc_diff%rain( i )%q( size( bc%rain( i )%q ) ) )

         bc_diff%rain( i )%t(:)  =  0._rp
         bc_diff%rain( i )%q(:)  =  0._rp

      end do

      allocate( innovation_diff ( size( innovation ) ) )
      allocate( innovW_diff( size( innovW ) ) )
      allocate( innovQ_diff( size( innovQ ) ) )

      do iobs = 1,size( innovation )
         allocate( innovation_diff ( iobs )%diff( size( innovation ( iobs )%diff ) ) )
         innovation_diff ( iobs )%diff(:)  =  0._rp
      end do

      do iobs = 1,size( innovW )
         allocate( innovW_diff( iobs )%diff( size( innovW( iobs )%diff ) ) )
         innovW_diff( iobs )%diff(:)  =  0._rp
      end do

      do iobs = 1,size( innovUV )
         allocate( innovUV_diff( iobs )%diff( size( innovUV( iobs )%diff ) ) )
         innovUV_diff( iobs )%diff(:)  =  0._rp
      end do

      do iobs = 1,size( innovQ )
         allocate( innovQ_diff( iobs )%diff( size( innovQ( iobs )%diff ) ) )
         innovQ_diff( iobs )%diff(:)  =  0._rp
      end do

      dt_diff = 0._rp
      #ifdef USE_SW_MONO

         allocate( XSshape_diff        ( size( XSshape    ) ) )

         allocate( manning_diff        ( size( manning    ) ) )
         allocate( manning_beta_diff   ( size( manning_beta    ) ) )
         allocate( bathy_cell_diff     ( size( bathy_cell ) ) )
         allocate( slope_y_diff (1_ip))
         allocate( slope_x_diff (1_ip))

         XSshape_diff(1)%xleft = 0._rp
         XSshape_diff(1)%xcenter = 0._rp
         XSshape_diff(1)%xright = 0._rp
         XSshape_diff(1)%s = 0._rp
         XSshape_diff(1)%hmax = 0._rp
         XSshape_diff(1)%topz = 0._rp

         manning_diff(:)     =  0._rp
         manning_beta_diff(:)=  0._rp
         bathy_cell_diff(:)  =  0._rp
         slope_y_diff(:) = 0._rp
         slope_x_diff(:) = 0._rp

         allocate( bc_diff%sum_mass_flux( bc%nb ) )

      #endif

   END SUBROUTINE alloc_diff_vars


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Allocation of needed variables to run the adjoint model
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE alloc_back_vars(dof0_back, dof_back, mesh )
	 type(unk) , intent(inout) :: dof0_back, dof_back
      type(msh), intent(in)  ::  mesh

      call alloc_dof( dof0_back , mesh )
      call alloc_dof( dof_back  , mesh )

      !===============================================
      ! Classical hydraulic BCs
      !===============================================
      allocate( bc_back%inflow ( size( bc%inflow  ) ) )
      allocate( bc_back%outflow( size( bc%outflow ) ) )

      bc_back%inflow (:)  =  0._rp
      bc_back%outflow(:)  =  0._rp

#ifdef USE_HYDRO
      allocate( bc_back%gr4( size(bc%gr4) ))

      do i = 1,size(bc%gr4)

         allocate( bc_back%gr4( i )%t( size( bc%gr4( i )%t ) ) )
         allocate( bc_back%gr4( i )%Q( size( bc%gr4( i )%Q ) ) )

         bc_back%gr4( i )%t(:)  =  0._rp
         bc_back%gr4( i )%Q(:)  =  0._rp
         bc_back%gr4( i )%params(:) =  0._rp
         bc_back%gr4( i )%state(:) =  0._rp

      end do
#endif

      allocate( bc_back%hyd( bc%nb_in ) )

      do i = 1,bc%nb_in

         allocate( bc_back%hyd( i )%t( size( bc%hyd( i )%t ) ) )
         allocate( bc_back%hyd( i )%q( size( bc%hyd( i )%q ) ) )

         bc_back%hyd( i )%t(:)  =  0._rp
         bc_back%hyd( i )%q(:)  =  0._rp

      end do


      allocate( bc_back%rat( bc%nb_out ) )

      do i = 1,bc%nb_out
        write(*,*) "rating curve must always be allocated while doing minimization (even if not used)"

         allocate( bc_back%rat( i )%h( size( bc%rat( i )%h ) ) )
         allocate( bc_back%rat( i )%q( size( bc%rat( i )%q ) ) )

         bc_back%rat( i )%h(:)  =  0._rp
         bc_back%rat( i )%q(:)  =  0._rp

      end do

      !===============================================
      ! Rain BCs
      !===============================================

      if (bc_rain == 1) then

        allocate( bc_back%rain( bc%nb_rn ) )

        do i = 1,bc%nb_rn
                ! dirty
!             if (.not. allocated(bc_back%rain( i )%t) ) then
            allocate( bc_back%rain( i )%t( size( bc%rain( i )%t ) ) )
!                 end if
!             if (.not. allocated(bc_back%rain( i )%q) ) then
            allocate( bc_back%rain( i )%q( size( bc%rain( i )%q ) ) )
!                 end if

            bc_back%rain( i )%t(:)  =  0._rp
            bc_back%rain( i )%q(:)  =  0._rp
        enddo

      else

            allocate( bc_back%rain( 1 ) )
            allocate( bc_back%rain( 1 )%t( 1 ) )
            allocate( bc_back%rain( 1 )%q( 1 ) )

            bc_back%rain( 1 )%t(1)  =  0._rp
            bc_back%rain( 1 )%q(1)  =  0._rp

      endif

      !===============================================
      ! Infiltration params
      !===============================================

      allocate( infil_back%GA ( size( infil%GA ) ) )
      allocate( infil_back%SCS( size( infil%SCS ) ) )

      infil_back%GA(:)%Ks           =  0._rp
      infil_back%GA(:)%PsiF         =  0._rp
      infil_back%GA(:)%DeltaTheta   =  0._rp
      infil_back%SCS(:)%lambdacn      =  0._rp
      infil_back%SCS(:)%CN          =  0._rp
      
      allocate( ptf_back( size( PTF ) ) )
      do i = 1, size( PTF )
        ptf_back(i)%kappa(:) = 0._rp
      enddo

      allocate( innovation_back ( size( innovation ) ) )
      allocate( innovW_back     ( size( innovW ) ) )
      allocate( innovUV_back     ( size( innovUV ) ) )
      allocate( innovQ_back( size( innovQ ) ) )

      do iobs = 1,size( innovation )
         allocate( innovation_back ( iobs )%diff( size( innovation( iobs )%diff ) ) )
         innovation_back ( iobs )%diff(:)  =  0._rp
      end do

      do iobs = 1,size( innovW )
         allocate( innovW_back( iobs )%diff( size( innovW( iobs )%diff ) ) )
         innovW_back( iobs )%diff(:)  =  0._rp
      end do

      do iobs = 1,size( innovUV )
         allocate( innovUV_back( iobs )%diff( size( innovUV( iobs )%diff ) ) )
         innovUV_back( iobs )%diff(:)  =  0._rp
      end do

      do iobs = 1,size( innovQ )
         allocate( innovQ_back( iobs )%diff( size( innovQ( iobs )%diff ) ) )
         innovQ_back( iobs )%diff(:)  =  0._rp
      end do


      dt_back = 0._rp

      #ifdef USE_SW_MONO

         allocate( XSshape_back   ( size( XSshape    ) ) )
         allocate( manning_back   ( size( manning    ) ) )
         allocate( manning_beta_back   ( size( manning_beta    ) ) )
         allocate( bathy_cell_back( size( bathy_cell ) ) )
         allocate( slope_y_back (size(slope_y)))
         allocate( slope_x_back (size(slope_x)))

         XSshape_back(:)%xleft = 0._rp
         XSshape_back(:)%xcenter = 0._rp
         XSshape_back(:)%xright = 0._rp
         XSshape_back(:)%s = 0._rp
         XSshape_back(:)%hmax = 0._rp
         XSshape_back(:)%topz = 0._rp

         manning_back(:)      =  0._rp
         manning_beta_back(:) =  0._rp
         bathy_cell_back(:)   =  0._rp
         slope_y_back(:)   =  0._rp
         slope_x_back(:)   =  0._rp

         allocate( bc_back%sum_mass_flux( bc%nb ) )

      #endif

   END SUBROUTINE alloc_back_vars


SUBROUTINE dealloc_back_vars()
      !type(unk) , intent(inout) :: dof0_back, dof_back
      !if(allocated(dof0_back )) deallocate(dof0_back)
      !if(allocated(dof_back )) deallocate(dof_back)

      if(allocated(bc_back%inflow )) deallocate(bc_back%inflow )
      if(allocated(bc_back%outflow )) deallocate(bc_back%outflow )
      if(allocated(bc_back%rat )) deallocate(bc_back%rat )
      if(allocated(bc_back%hyd )) deallocate(bc_back%hyd )
      if(allocated(XSshape_back )) deallocate(XSshape_back )
      if(allocated(slope_x_back )) deallocate(slope_x_back )
      if(allocated(slope_y_back )) deallocate(slope_y_back )
      #ifdef USE_HYDRO
      if(allocated(bc_back%gr4 )) deallocate(bc_back%gr4 )
      #endif USE_HYDRO
      if(allocated( bc_back%rain)) deallocate( bc_back%rain)
      if(allocated(infil_back%GA )) deallocate(infil_back%GA )
      if(allocated(infil_back%SCS  )) deallocate(infil_back%SCS)
      if(allocated(ptf_back )) deallocate(ptf_back )

      if(allocated(innovation_back)) deallocate(innovation_back)
      if(allocated(innovUV_back)) deallocate(innovUV_back)
      if(allocated(innovW_back)) deallocate(innovW_back)
      if(allocated(innovQ_back)) deallocate(innovQ_back)

      #ifdef USE_SW_MONO
      if(allocated(manning_back))  deallocate(manning_back)
      if(allocated(manning_beta_back))  deallocate(manning_beta_back)
      if(allocated(bathy_cell_back))  deallocate(bathy_cell_back)
      if(allocated(bc_back%sum_mass_flux)) deallocate(bc_back%sum_mass_flux)
      #endif USE_SW_MONO


      !if(allocated(control )) deallocate(control )            !to check this does not generate issues lilian
      !if(allocated(control_back )) deallocate(control_back )  !to check this does not generate issues lilian
      !if(allocated(control_diff )) deallocate(control_diff )  !to check this does not generate issues lilian


END SUBROUTINE dealloc_back_vars


END MODULE m_tap_vars
