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
!> \file initialization.f90
!! \brief This file includes Initial routine.
!! \details The file includes only Initial routine (see doc Initial routine ).

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Initialization Subroutine specific to Shallow-Water Model
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!>  Initialization Subroutine specific to Shallow-Water Model
!!
!! \details
SUBROUTINE Initial( dof0 , mesh, my_friction, my_infiltration, my_param_model, my_phys_desc)

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_model
   USE m_numeric
   USE m_obs
   USE m_time_screen

#ifdef USE_HYDRO
!    USE m_gr4
#endif

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(inout)  ::  mesh
   type( unk ), intent(inout)  ::  dof0
   type( friction_data )    , intent(in   )  ::  my_friction
   type( infiltration_data ), intent(in   )  ::  my_infiltration
   type( param_model ), intent(in   )  ::  my_param_model
   type( input_data ), intent(in   )  ::  my_phys_desc

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   character(len=lchar)  ::  filename

   integer(ip)  ::  icod, i_loc, l

   integer(ip), dimension(:), allocatable  ::  index_bathy_min

   real(rp), dimension(:), allocatable  ::  bathy_min , bathy_min_glob,temp

   !===================================================================================================================!
   !  Bathymetry Initialization --> CLEAN done via text file 'dassflow format'
   !===================================================================================================================!

!     call my_bathy_2_fortran() !(my_param_model)
!
!     call my_friction_2_fortran(my_friction) ! propagate definition of friction from fortran to manning,
!
     if (bc_infil .ne. 0) call my_infiltration_2_fortran(my_infiltration)

!      call my_phys_desc_2_fortran(my_phys_desc)

#ifdef USE_MPI
         call swap_vec_i  ( land , swap_index( 1 : mesh%nc ) )
         call reallocate_i( land ,                 mesh%nc   )
         call swap_vec_r  ( bathy_cell , swap_index( 1 : mesh%nc + mesh%ncb ) )
#endif

   !===================================================================================================================!
   !  DOF Initialization --> dof pushed from python
   ! just keep parralel order com_dof
   !===================================================================================================================!

     inquire( file = 'ic.bin'      , exist = file_exist(1) )
     inquire( file = 'restart.bin' , exist = file_exist(2) )
     inquire( file = 'dof_init.txt' , exist = file_exist(3) )
     if      ( file_exist(1) ) then

        open(10,file='ic.bin',form='unformatted',status='old',access='direct',recl=3*length_real)

        read(10,rec=1) tc0

        tc0 = 0._rp

     else if ( file_exist(2) ) then

       open(10,file='restart.bin',form='unformatted',status='old',access='direct',recl=3*length_real)

        read(10,rec=1) tc0

        if ( abs( ts - tc0 ) < zerom ) call Stopping_Program_Sub( 'End of simulation time reached' )

     end if

     if ( file_exist(1) .or. file_exist(2)) then

        do i = 1,mesh%nc

  		!read(10,rec=1) dof0%h(i) , dof0%u(i) , dof0%v(i)
          read(10,rec=1+swap_index(i)) dof0%h(i) , dof0%u(i) , dof0%v(i)

        end do

        close(10)

     endif


! COMMENTED, need to replace zs0_user, u0_user and v0_user
!      else
!
!         do i = 1,mesh%nc
!
!            dof0%h(i)  = max( 0._rp , zs0_user( mesh%cell(i)%grav%x , mesh%cell(i)%grav%y ) - bathy_cell(i) ) ! h_ex(mesh%cell(i)%grav%x , mesh%cell(i)%grav%y)!
!
!            if ( dof0%h(i) > heps ) then
!
!               dof0%u(i)  =  u0_user( mesh%cell(i)%grav%x , mesh%cell(i)%grav%y )
!               dof0%v(i)  =  v0_user( mesh%cell(i)%grav%x , mesh%cell(i)%grav%y )
!            else
!               dof0%u(i)  =  0._rp
!               dof0%v(i)  =  0._rp
!
!            end if
!
!         end do
!      end if


       if(file_exist(3)) then

       open(10,file='dof_init.txt',status='old',form='formatted')
  		do i = 1,mesh%nc

  			read(10,*) dof0%h(i), dof0%u(i), dof0%v(i)

  		end do

       end if

   call com_dof( dof0 , mesh )

   !===================================================================================================================!
   !  Boundary Condition Initialization
   !===================================================================================================================!


!~    if ( mesh_type == 'basic' ) then

!~       bc%nb = 4

!~       allocate( bc%typ ( bc%nb , 2 ) )
!~       allocate( bc%grpf( bc%nb     ) )

!~       bc%typ(1,1)  =  bc_W
!~       bc%typ(2,1)  =  bc_E
!~       bc%typ(3,1)  =  bc_S
!~       bc%typ(4,1)  =  bc_N

!~       bc%typ(:,2)  =  ''

!~    end if

   allocate( bc%inflow ( 2 * mesh%neb ) )
   allocate( bc%outflow( 2 * mesh%neb ) )
   allocate( bc%sum_mass_flux( bc%nb ) )


#ifdef USE_HYDRO
   !===================================================================================================================!
   !  Reading input data for GR4 module
   !===================================================================================================================!

   inquire( file = 'GR4params.txt' , exist = file_exist(1) )

   call mpi_wait_all

   if ( file_exist(1) ) then

      open(10,file='GR4params.txt',status='old')

      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*) bc%nb_gr4in

      if (.not. allocated(bc%gr4)) allocate( bc%gr4( bc%nb_gr4in ))

      if (use_Qobs_gr4 == 1) then
         if (.not. allocated(stationQ)) allocate( stationQ(bc%nb_gr4in))
      else
        write(*,*) "Attempt to allocate already allocated variable stationQ when reading GR4params.txt"
      endif

      do i=1,bc%nb_gr4in

        read(10,*)
        read(10,*) bc%gr4( i )%params(:)
        read(10,*)
        read(10,*) bc%gr4( i )%state(:)
        read(10,*)
        read(10,*) bc%gr4( i )%cell_id, bc%gr4( i )%surf

      enddo

      close(10)

      do i=1,bc%nb_gr4in

       write(filename,'("GR4forcings_",1I1.1,".txt")') i

      open(10,file=filename,status='old',form='formatted')

            read(10,*)
            read(10,*) l

            allocate( bc%gr4( i )%t( l ))
            allocate( bc%gr4( i )%P( l ))
            allocate( bc%gr4( i )%E( l ))
            allocate( bc%gr4( i )%Q( l ))

        do j=1,l
!             write(*,*) i, j
            read(10,*) bc%gr4( i )%t( j ), bc%gr4( i )%P( j ), bc%gr4( i )%E( j )
!             write(*,*) bc%gr4( i )%t( j ), bc%gr4( i )%P( j ), bc%gr4( i )%E( j )
        enddo

      close(10)
       enddo

   else

!     allocate( bc%gr4( 1 ) )
!     allocate( bc%gr4( 1 )%P( 1 ) )
!     allocate( bc%gr4( 1 )%E( 1 ) )
!     allocate( bc%gr4( 1 )%t( 1 ) )
!     allocate( bc%gr4( 1 )%Q( 1 ) )
!
!     else

      call Stopping_Program_Sub( 'File GR4data.txt not provided ... Stopping.' )

   endif




   !===================================================================================================================!
   !  Initialize hydrological reservoir states (1 year warmup)
   !===================================================================================================================!
   inquire( file = 'GR4_PEwarmup.txt' , exist = file_exist(1) )

   call mpi_wait_all

   if ( file_exist(1) ) then

      open(10,file='GR4_PEwarmup.txt',status='old')

      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*) l

      if (l .ne. bc%nb_gr4in) write(*,*) 'There are ',l,' GR4 warmup catchements and ',&
     bc%nb_gr4in, ' GR4 launch catchments'

      do i=1,bc%nb_gr4in

        read(10,*)
        read(10,*) l !should be 8670 , i.e. 365 days, for warmup

        allocate( bc%gr4( i )%P0( l ))
        allocate( bc%gr4( i )%E0( l ))

        do j=1,l
         read(10,*) bc%gr4( i )%P0( l ), bc%gr4( i )%E0( l )
        enddo

     enddo

   else
        write(*,*)'No GR4_PEwarmup.txt file found.'
   endif


!    do catchnb = 1, nb_gr4in
!         call GR4_main("warmup") !1 year run using P0 and E0 to generate reservoir statest consistent
!    enddo
#endif

   !===================================================================================================================!
   !  Loading/Creating hydrograph File
   !===================================================================================================================!

   inquire( file = 'hydrograph.txt' , exist = file_exist(1) )

   call mpi_wait_all

   if ( file_exist(1) ) then

      open(10,file='hydrograph.txt',status='old')

      read(10,*)
      read(10,*)
      read(10,*)

      read(10,*) bc%nb_in

    allocate( bc%hyd( bc%nb_in ) )

      do i = 1,bc%nb_in
     ! write(*,*) "i=",i
        read(10,*)
        read(10,*)
        read(10,*)

        read(10,*) j
       ! write(*,*)"j = ", j



         allocate( bc%hyd( i )%t( j ) )
         allocate( bc%hyd( i )%q( j ) )


         do k = 1,j
			!write(*,*) k
            read(10,*) bc%hyd( i ) %t( k ) , bc%hyd( i )%q( k )

         end do

         if ( mesh_type == 'basic' ) then

            do k = 1,4

               if ( bc%typ(k,1)(1:8)  ==  'discharg' ) then

                  bc%typ ( k , 2 )  =  'file'

                  bc%grpf( k     )  =  i

                  exit

               end if

               if ( k == 4 .and. bc%typ(k,2) == '' ) &

               !call Stopping_Program_Sub( 'Hydrograph given and no discharg boundary is specified' )
                write(*,*) 'BEWARE : Hydrograph given and no discharg boundary is specified'

            end do

         end if

      end do


		close(10)

	else if ( any( bc%typ(:,1)(1:8)  ==  'discharg' ) ) then
				write(*,*) " <ERROR> no File hydrograph.txt provided and any( bc%typ(:,1)(1:8)  ==  'discharg' )"
!~ 	!#write(*,*) dta
!~       if ( dta > zerom ) then

!~          call gen_hydrograph
!~          write(*,*) 'File hydrograph.txt generated with provided dta and inflow_user'
!~          !call Stopping_Program_Sub( 'File hydrograph.txt generated with provided dta and inflow_user' )

!~       else
!~          call gen_hydrograph
!~          write(*,*) ' <WARNING> File hydrograph.txt generated with provided dta and inflow_user'
!~          !call Stopping_Program_Sub( 'File hydrograph.txt not provided ... (can be genetated with dta and inflow_user)' )

!~       end if

    end if

   !===================================================================================================================!
   !  Loading/Creating hpresc File
   !====================================================================================================
   inquire( file = 'hpresc.txt' , exist = file_exist(1) )

   call mpi_wait_all

   if ( file_exist(1) ) then

      open(10,file='hpresc.txt',status='old')

      read(10,*)
      read(10,*)
      read(10,*)

      read(10,*) bc%nb_out

    allocate( bc%hpresc( bc%nb_out ) )

      do i = 1,bc%nb_out
		! 3 comment lines of the file
        read(10,*)
        read(10,*)
        read(10,*)
		! size of hpresc
        read(10,*) j

         allocate( bc%hpresc( i )%t( j ) )
         allocate( bc%hpresc( i )%h( j ) )

		! read all available value
         do k = 1,j
            read(10,*) bc%hpresc( i ) %t( k ) , bc%hpresc( i )%h( k )
         end do
      end do
		close(10)

    end if


   !===================================================================================================================!
   !  Loading/Creating hpresc File
   !===================================================================================================================!

   inquire( file = 'zspresc.txt' , exist = file_exist(2) )

   call mpi_wait_all

   if ( file_exist(2) ) then

      open(10,file='zspresc.txt',status='old')

      read(10,*)
      read(10,*)
      read(10,*)

      read(10,*) bc%nb_out

    allocate( bc%zspresc( bc%nb_out ) )

      do i = 1,bc%nb_out
		! 3 comment lines of the file
        read(10,*)
        read(10,*)
        read(10,*)
		! size of zspresc
        read(10,*) j

         allocate( bc%zspresc( i )%t( j ) )
         allocate( bc%zspresc( i )%z( j ) )

		! read all available value
         do k = 1,j
            read(10,*) bc%zspresc( i ) %t( k ) , bc%zspresc( i )%z( k )
         end do
      end do
		close(10)

    end if

   !===================================================================================================================!
   !  Loading/Creating rain File
   !===================================================================================================================!

   if (bc_rain  == 1) then

   inquire( file = 'rain.txt' , exist = file_exist(1) )

   call mpi_wait_all

   if ( file_exist(1) ) then

      open(10,file='rain.txt',status='old')

      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*) bc%nb_rn
      read(10,*)
      read(10,*)
      read(10,*)

      allocate(bc%rain(bc%nb_rn))

      do i = 1,bc%nb_rn
         read(10,*) bc%rain(i)%x_min, bc%rain(i)%x_max, bc%rain(i)%y_min, bc%rain(i)%y_max
      enddo

      do i = 1,bc%nb_rn

         read(10,*)
         read(10,*)
         read(10,*)
         read(10,*) j

         allocate( bc%rain(i)%t(j) )
         allocate( bc%rain(i)%q(j) )

         do k = 1,j
            read(10,*) bc%rain(i)%t(k) , bc%rain(i)%q(k)
         enddo

      end do

      close(10)

    else
         call Stopping_Program_Sub( 'File rain.txt not provided ...')
    endif

   else

     allocate( bc%rain(1) )
     allocate( bc%rain(1)%q(1) )
     allocate( bc%rain(1)%t(1) )

     bc%rain(1)%q = 0._rp
     bc%rain(1)%t = 0._rp
     bc%rain(1)%cumul = 0._rp
     bc%rain(1)%qin = 0._rp

   end if
   !===================================================================================================================!
   !  Rain condition Initialization
   !===================================================================================================================!

   if (bc_rain == 1) then

      do i=1,mesh%nc

         do k=1,bc%nb_rn

            if ((mesh%cell(i)%grav%x > bc%rain(k)%x_min .and. mesh%cell(i)%grav%x < bc%rain(k)%x_max) .and. &
                  (mesh%cell(i)%grav%y > bc%rain(k)%y_min .and. mesh%cell(i)%grav%y < bc%rain(k)%y_max)) then

                  mesh%cell(i)%rain = k
                  exit
            else

                  mesh%cell(i)%rain = 0_ip

            endif

         enddo


      enddo

        open(10,file='rain_post.dat',status='replace',form='formatted')

        write(10,*) '# Gnuplot DataFile Version'
        write(10,*) '# x y rain'

        do i=1,mesh%nc
                    write(10,'(2ES15.8,i3)') mesh%cell(i)%grav%x    , &
                                        mesh%cell(i)%grav%y    , &
                                        mesh%cell(i)%rain

        end do

        close(10)

   else

        mesh%cell(:)%rain = 0_ip

   endif


   !===================================================================================================================!
   !  Loading rating curve
   !===================================================================================================================!

   inquire( file = 'rating_curve.txt' , exist = file_exist(1) )

  call mpi_wait_all

   if ( file_exist(1) ) then

		open(10,file='rating_curve.txt',status='old')

		read(10,*)
		read(10,*)
		read(10,*)

		read(10,*) bc%nb_out

		allocate( bc%rat( bc%nb_out ) )

		do i = 1,bc%nb_out

			read(10,*)
			read(10,*)
			read(10,*)

			read(10,*) j , bc%rat( i )%z_rat_ref


			allocate( bc%rat( i )%h( j ) )
			allocate( bc%rat( i )%q( j ) )

			do k = 1,j

				read(10,*) bc%rat( i )%h( k ) , bc%rat( i )%q( k )

			end do

			if ( mesh_type == 'basic' ) then

				do k = 1,4

					if ( bc%typ(k,1)(1:8)  ==  'ratcurve' ) then

						bc%typ ( k , 2 )  =  'file'

						bc%grpf( k     )  =  i

					exit

					end if

					if ( k == 4 .and. bc%typ(k,2) /= 'file' ) &

						!call Stopping_Program_Sub( 'Ratcurve given and no ratcurve boundary is specified' )
                write(*,*) 'BEWARE : Ratcurve given and no ratcurve boundary is specified'
				end do

			end if

		end do

		close(10)

	else if ( any( bc%typ(:,1)  ==  'ratcurve' ) ) then

		call Stopping_Program_Sub( 'File ratcurve.txt not provided ...')

	else ! for debug inverse
		allocate( bc%rat( 1 ) )
		allocate( bc%rat(1)%q( 1 ) )
		allocate( bc%rat(1)%h( 1 ) )

		bc%rat(1)%q( 1 ) = 1
		bc%rat(1)%h( 1 ) = 1
		bc%rat(1)%z_rat_ref = 1
		bc%rat(1)%zout = 1
		bc%rat(1)%c1 = 1
		bc%rat(1)%c2 = 1
		bc%rat(1)%pow(:) = 0
		bc%rat(1)%group = 0
	end if

   !===================================================================================================================!
   !  Setting minimal condition to initialize discharg bc
   !===================================================================================================================!

   allocate( index_bathy_min( bc%nb ) )

   allocate( bathy_min( bc%nb ) )

   allocate( bathy_min_glob( np ) )

   bathy_min(:)        =   hugem
   index_bathy_min(:)  = - 1

   do ie = 1,mesh%neb

      if ( mesh%edgeb(ie)%typlim(1:8) == 'discharg' ) then

         i = mesh%edgeb(ie)%group
         j = mesh%edge( mesh%edgeb(ie)%ind )%cell(1)

         if ( bathy_cell(j) < bathy_min(i) ) then

            bathy_min(i) = bathy_cell(j)

            index_bathy_min(i) = j

         end if

      end if

   end do

   do i = 1,bc%nb

      if ( bc%typ(i,1)(1:8) == 'discharg' ) then

         call mpi_allgather_r( bathy_min(i) , bathy_min_glob )

         if ( index_bathy_min(i) > 0 .and. proc == minloc( bathy_min_glob , 1 ) - 1 ) then

            if ( dof0%h( index_bathy_min(i) ) == 0._rp ) then ! security : h must be strictly positive to corectly impose discharge

               dof0%h( index_bathy_min(i) ) =  1.d-1

            end if

         end if

      end if

   end do

   deallocate( index_bathy_min , bathy_min , bathy_min_glob )

   !===================================================================================================================!
   !  Cell Bathymetry Gradient
   !===================================================================================================================!

   allocate( grad_z ( mesh%nc + mesh%ncb ) )
   allocate( grad_z2( mesh%nc + mesh%ncb ) )
   allocate( z_eq   ( mesh%nc + mesh%ncb ) )

   call FV_Cell_Grad ( grad_z  , bathy_cell , mesh )
   call FV_Cell_Grad2( grad_z2 , bathy_cell , mesh )

   !===================================================================================================================!
   !  Set heps to be at minimum to zero machine precision
   !===================================================================================================================!

   heps = heps + zerom

   !===================================================================================================================!
   !  Initilization of Stations and Sections to Output it Reading the obs.txt File
   !===================================================================================================================!

   nb_obs = 0


   inquire( file = 'obs.txt' , exist = file_exist(1) )

   if ( w_obs == 1 .and. file_exist(1) ) then
      !================================================================================================================!
      !  Opening Data File concerning Stations and Sections Recording at prescribed frequency
      !================================================================================================================!

      open(10,file='obs.txt',form='formatted',status='old')

      !================================================================================================================!
      !  Treating First Time Stations or Sections Data Record
      !================================================================================================================!

      buffer = '' ; icod = 1

      do while ( buffer(1:8) /= 'stations'          .and. &
                buffer(1:8) /= 'stations_with_grp' .and. &
                 buffer(1:8) /= 'sections'          .and. icod >= 0 )


		read(10,'(A)',iostat=icod) buffer

      end do
      if ( icod >= 0 ) call read_obs_file
		!================================================================================================================!
		!  Treating Second Time Stations or Sections Data Record
		!================================================================================================================!

		buffer = '' ; icod = 1

		do while ( buffer(1:8) /= 'stations'          .and. &
				 buffer(1:8) /= 'stations_with_grp' .and. &
				 buffer(1:8) /= 'sections'          .and. icod >= 0 )

		read(10,'(A)',iostat=icod) buffer

! 		write(*,*) buffer

		end do

		if ( icod >= 0 ) call read_obs_file
		!================================================================================================================!
		!  Treating Thrid Time Stations or Sections Data Record
		!================================================================================================================!

		buffer = '' ; icod = 1

		do while ( buffer(1:8) /= 'stations'          .and. &
				 buffer(1:8) /= 'stations_with_grp' .and. &
				 buffer(1:8) /= 'sections'          .and. icod >= 0 )

		 read(10,'(A)',iostat=icod) buffer

! 		write(*,*) buffer
		end do

		if ( icod >= 0 ) call read_obs_file

		!================================================================================================================!
		!  Closing Data File concerning Stations and Sections Recording at prescribed frequency
		!================================================================================================================!

		close(10)


    elseif ( w_obs == 1 .and. .not. file_exist(1) ) then

 		call Stopping_Program_Sub( 'No obs.txt is provided ') ! to add f90wrap abbort

   end if

!    ===================================================================================================================!
!     Reading param_obs.txt
!    ===================================================================================================================!
   inquire( file = 'param_obs.txt' , exist = file_exist(1) )
   if ( w_obs == 1 .and. file_exist(1) ) then
      open(10,file='param_obs.txt',form='formatted',status='old')
      read(10,'(A)',iostat=icod) buffer

      do iobs = 1,size( station )
         read(10,*) station( iobs )%length , station(iobs)%dt_offset, temp,temp,temp,temp
      end do
      close(10)

   endif


   !===================================================================================================================!
   !  Array dt obs
   !===================================================================================================================!

   if ( w_obs == 1 ) then

    do iobs = 1,size(station)

        station(iobs)%nb_dt=int((ts- station(iobs)%dt_offset)/station(iobs)%dt)+1_ip

        allocate( station(iobs)%dt_obs(station(iobs)%nb_dt  ) )

        station(iobs)%dt_obs(1)=station(iobs)%dt_offset

        do i =2,station(iobs)%nb_dt
          station(iobs)%dt_obs(i)=station(iobs)%dt_obs(i-1)+station(iobs)%dt ! Create array of observation time ( offset + repitivity satellite)
        end do

        station( iobs )%ind_t=1_ip

    end do

  end if

   !===================================================================================================================!
   !  Initialize obs
   !===================================================================================================================!

   if ( use_obs == 1 ) then

	call read_stations()

      if (allocated( innovation )) deallocate( innovation )
      if (allocated( innovW )) deallocate( innovW )
      if (allocated( innovQ )) deallocate( innovQ )

	  allocate( innovation ( size( station ) ) )
	  allocate( innovW( size( station ) ) )

      do iobs = 1,size( station )

     innovation ( iobs )%nb_dt  =  size( station( iobs )%t(:) )
		 innovW( iobs )%nb_dt  =  size( station( iobs )%t(:) )

         innovation ( iobs )%nb_dx  =  1
         innovW( iobs )%nb_dx  =  1

		if (allocated (innovation ( iobs )%diff)) then
				write(*,*) "already allocated innovation ( iobs )%diff"
		else
				allocate( innovation ( iobs )%diff( innovation( iobs )%nb_dt * &
                                            innovation( iobs )%nb_dx ) )
				allocate( innovW( iobs )%diff( innovW( iobs )%nb_dt * &
                                            innovW( iobs )%nb_dx ) )
		end if

         innovation ( iobs )%diff(:)  =  0._rp
         innovW( iobs )%diff(:)  =  0._rp

      end do

  else


   end if

!#ifdef USE_HYDRO !NOT CLEAN, ADD CORRECT FLAGS
   if ( use_Qobs == 1 .or. use_Qobs_gr4 == 1) then

      call read_stationsQ

      allocate( innovQ( size( stationQ ) ) )

      do iobs = 1,size( stationQ )

         innovQ( iobs )%nb_dt  =  size( stationQ( iobs )%t(:) )

         innovQ( iobs )%nb_dx  =  1

         allocate( innovQ( iobs )%diff( innovQ( iobs )%nb_dt ))

         innovQ( iobs )%diff(:)  =  0._rp

      end do

   else
      if (allocated( innovQ )) deallocate( innovQ )
      allocate( innovQ(1) )

   end if
!#endif

   time(:) = 0._rp

   call display_mesh_cell(mesh)

CONTAINS


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Reading Stations or Sections informations in obs.txt file and Creating associated structure
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE read_obs_file

      implicit none

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      type( vec2d )  ::  pt_temp(2)

      real(rp)  ::  length

      integer(ip)  ::  nb_rec , nb_pt , pt , search_cell_inc_point , etq , cell

      logical  ::  point_in_cell

      character(len=8)  ::  typ

      !================================================================================================================!
      !  Beginning Treatment of File 'obs.txt'
      !================================================================================================================!

      backspace 10

      read(10,*) buffer , nb_rec
      select case( trim(buffer) )

         !=============================================================================================================!
         !  Stations at One Point Case
         !=============================================================================================================!

		case( 'stations' )
!~ 			print *, buffer, nb_rec
            call alloc_or_realloc_station( station , nb_obs + nb_rec )

            do iobs = nb_obs+1 , nb_obs+nb_rec

               !=======================================================================================================!
               !  Allocate Point Coordinate (only one in this case)
               !=======================================================================================================!

               allocate( station( iobs )%pt(1) )

               !=======================================================================================================!
               !  Reading Stations Informations in 'obs.txt' File
               !=======================================================================================================!

               read(10,*) station( iobs )%pt(1)%coord%x , &
                          station( iobs )%pt(1)%coord%y , &
                          station( iobs )%dt            , &
                          station( iobs )%weight


               !=======================================================================================================!
               !  Searching Cells for Stations
               !=======================================================================================================!

               station( iobs )%pt(1)%cell = search_cell_inc_point( mesh , station( iobs )%pt(1)%coord )

               !=======================================================================================================!
               !  Writing Ouput Files to Visualize Stations
               !=======================================================================================================!

!                write(buffer,'(A,I4.4)') 'res/obs/pos_station_', iobs
                 if ( proc == 0 ) call system('mkdir -p res/obs')
!                call write_station_pos_in_file( buffer , station( iobs )%pt(1)%coord )

            end do

         !=============================================================================================================!
         !  Stations with a Group of Points Case
         !=============================================================================================================!

         case( 'stations_with_grp' )
            call alloc_or_realloc_station( station , nb_obs + nb_rec )

            do iobs = nb_obs+1 , nb_obs+nb_rec

               !=======================================================================================================!
               !  Reading Stations Informations in 'obs.txt' File
               !=======================================================================================================!

               read(10,*) etq , station( nb_obs + etq )%dt , station( nb_obs + etq )%weight


               write(buffer,'(A,I4.4,A)') 'station_grp_' , iobs - nb_obs , '.txt'

               nb_pt = count_lines( buffer ) - 1

               allocate( station( iobs )%pt( nb_pt ) )

               open(20,file=buffer,status='old',form='formatted')

               read(20,'(A)') typ

               do pt = 1,nb_pt

                  !====================================================================================================!
                  !  Reading Stations Informations in 'station_grp_xxxx.txt' File
                  !====================================================================================================!

                  if      ( trim(typ) == 'points' ) then
                     read(20,*) station( iobs )%pt( pt )%coord%x , &
                                station( iobs )%pt( pt )%coord%y

                     !=================================================================================================!
                     !  Searching Cells for Stations
                     !=================================================================================================!

                     station( iobs )%pt( pt )%cell = search_cell_inc_point( mesh , station( iobs )%pt( pt )%coord )

                  else if ( trim(typ) == 'indexes' ) then

                     read(20,*) cell

                     station( iobs )%pt( pt )%cell = cell
                     station( iobs )%pt( pt )%coord  =  mesh%cell( station( iobs )%pt( pt )%cell )%grav

                     if ( part( cell ) == proc ) then

                        station( iobs )%pt( pt )%cell = inv_swap_index( cell )

                        station( iobs )%pt( pt )%coord  =  mesh%cell( station( iobs )%pt( pt )%cell )%grav

                     else

                        station( iobs )%pt( pt )%cell = -1

                     end if

                     call mpi_bcast_r( station( iobs )%pt( pt )%coord%x , part( cell ) )
                     call mpi_bcast_r( station( iobs )%pt( pt )%coord%y , part( cell ) )

                  else

                     call Stopping_Program_Sub( 'One station_grp_xxxx.txt file has a wrong header (points or indexes)' )

                  end if

                  !====================================================================================================!
                  !  Writing Ouput Files to Visualize Stations
                  !====================================================================================================!

!                   write(buffer,'(A,I4.4,A,I4.4)') 'res/obs/pos_station_', iobs , '_' , pt

!                   call write_station_pos_in_file( buffer , station( iobs )%pt( pt )%coord )

               end do

               close(20)

            end do

         !=============================================================================================================!
         !  Sections Case
         !=============================================================================================================!

         case( 'sections' )

            allocate( section( nb_rec ) )

            do iobs = 1,nb_rec

               !====================================================================================================!
               !  Reading Sections Informations in 'obs.txt' File
               !====================================================================================================!

               read(10,*) pt_temp(1)%x , &
                          pt_temp(1)%y , &
                          pt_temp(2)%x , &
                          pt_temp(2)%y , &
                          nb_pt , &
                          section( iobs )%dt

               !====================================================================================================!
               !  Calculating Sections points and normal
               !====================================================================================================!

               allocate( section( iobs )%pt( nb_pt ) )

               section( iobs )%pt( 1     )%coord  =  pt_temp(1)
               section( iobs )%pt( nb_pt )%coord  =  pt_temp(2)

               do pt = 2,nb_pt-1

                  section( iobs )%pt( pt )%coord%x  =  section( iobs )%pt( 1     )%coord%x + real( pt - 1 , rp ) * &
                                                     ( section( iobs )%pt( nb_pt )%coord%x - &
                                                       section( iobs )%pt( 1     )%coord%x ) / real( nb_pt-1 , rp )

                  section( iobs )%pt( pt )%coord%y  =  section( iobs )%pt( 1     )%coord%y + real( pt - 1 , rp ) * &
                                                     ( section( iobs )%pt( nb_pt )%coord%y - &
                                                       section( iobs )%pt( 1     )%coord%y ) / real( nb_pt-1 , rp )

               end do

               length  =  sqrt( ( pt_temp(2)%x - pt_temp(1)%x )**2 + &
                                ( pt_temp(2)%y - pt_temp(1)%y )**2 )

               section( iobs )%normal%x = ( pt_temp(1)%y - pt_temp(2)%y ) / length
               section( iobs )%normal%y = ( pt_temp(2)%x - pt_temp(1)%x ) / length

               section( iobs )%dx  =  length / real( nb_pt - 1 , rp )

              !====================================================================================================!
               !  Searching Cells for all points in Sections
               !====================================================================================================!

               do pt = 1,nb_pt

                  section( iobs )%pt( pt )%cell = search_cell_inc_point( mesh , section( iobs )%pt( pt )%coord )

               end do

            end do

         end select

         nb_obs = nb_obs + nb_rec

   END SUBROUTINE read_obs_file


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Generating an hydrograph.txt file using provided inflow_user (in m_user_data.f90) and dta (in input.txt)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE write_station_pos_in_file( file_name , coord )

      character(len=*), intent(in)  ::  file_name

      type( vec2d ), intent(in)  ::  coord

      if ( proc == 0 ) then

         buffer = file_name_ext( file_name, 'tecplot' )

         call system('mkdir -p res/obs')

         open(100,file=buffer,status='replace',form='formatted')

         write(100,'(A)') 'TITLE = "DassFlow Station Position"'
         write(100,'(A)') 'VARIABLES = "x","y"'

         write(100,'(2ES15.7)') coord%x , coord%y

         close(100)

      end if

   END SUBROUTINE write_station_pos_in_file


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Generating an hydrograph.txt file using provided inflow_user (in m_user_data.f90) and dta (in input.txt)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


!~    SUBROUTINE gen_hydrograph

!~       implicit none

!~       if ( proc  == 0 ) then

!~          allocate( bc%hyd( 1 ) )

!~         j = 1

!~          tc = 0._rp

!~          do while ( tc < ts )

!~            j = j + 1

!~             tc = tc + dta

!~          end do

!~          allocate( bc%hyd( 1 )%t( j ) )
!~          allocate( bc%hyd( 1 )%q( j ) )

!~          i = 1

!~          tc = 0._rp

!~          do while ( tc < ts )

!~             bc%hyd( 1 )%t( i )  =  tc
!~             bc%hyd( 1 )%q( i )  =  inflow_user( tc , 0._rp , 0._rp )

!~             tc = tc + dta

!~             i = i + 1

!~          end do

!~          bc%hyd( 1 )%t( i )  =  ts
!~          bc%hyd( 1 )%q( i )  =  inflow_user( ts , 0._rp , 0._rp )

!~          open(10,file='hydrograph.txt',status='replace')

!~          write(10,'(A)') '!===============================================!'
!~          write(10,'(A)') '!  Number of hydrograph'
!~          write(10,'(A)') '!===============================================!'

!~         write(10,'(I1)') 1

!~          write(10,'(A)') '!===============================================!'
!~          write(10,'(A)') '!  Hydrograph'
!~          write(10,'(A)') '!===============================================!'

!~          write(10,'(I10)') j

!~          do i = 1,j

!~            write(10,'(2ES15.7)') bc%hyd( 1 )%t( i ) , bc%hyd( 1 )%q( i )

!~          end do

!~          close(10)

!~       end if

!~    END SUBROUTINE gen_hydrograph





! use variable my_friction (wrapped varible)
! to set up fortran variables (manning, manning_beta, land)

SUBROUTINE my_friction_2_fortran(my_friction)

implicit none

   type( friction_data ), intent(in   )  ::  my_friction

     nland = my_friction%nland
!       write(*,*) proc, nland, size( my_friction%land )
!       allocate( land( size( my_friction%land ) ) )
!       allocate( manning( nland ) )
!       allocate( manning_beta( nland ) )

      ! loop on all cells to define patch correspondance
      do i = 1,mesh%nc
          land( i )  =  my_friction%land( i )
      end do

      !define values for each patch
      do i = 1,nland
        manning(i) = my_friction%manning(i)
        manning_beta(i) = my_friction%manning_beta(i)
      end do

END SUBROUTINE my_friction_2_fortran
!
!
!
!
! use variable my_infiltration (wrapped varible)
! to set up fortran variables (infil)
SUBROUTINE my_infiltration_2_fortran(my_infiltration)

implicit none

   type( infiltration_data ), intent(in   )  ::  my_infiltration

      infil%nland = my_infiltration%nland

if( infil%nland  > 0) then ! if not zero

      allocate( infil%land( size(my_infiltration%land ) ) )
      if (allocated (my_infiltration%GA)) then
        allocate( infil%GA  ( infil%nland ) )
        do i = 1, size(my_infiltration%GA)
          infil%GA(i)  = my_infiltration%GA(i)
        enddo
      else
        allocate( infil%GA  ( 1 ) )
        infil%GA  ( 1 )%Ks = 0._rp
        infil%GA  ( 1 )%DeltaTheta = 0._rp
        infil%GA  ( 1 )%PsiF = 0._rp
      endif

      if (allocated (my_infiltration%SCS)) then
        allocate( infil%SCS ( infil%nland ) )
        do i = 1, size(my_infiltration%SCS)
          infil%SCS(i) = my_infiltration%SCS(i)
        enddo
      else
        allocate( infil%SCS ( 1 ) )
        infil%SCS  ( 1 )%lambda = 0._rp
        infil%SCS  ( 1 )%CN = 0._rp
      endif

      allocate( infil%x_min ( infil%nland) )
      allocate( infil%x_max ( infil%nland) )
      allocate( infil%y_min ( infil%nland) )
      allocate( infil%y_max ( infil%nland) )

      ! loop on all cells to define patch correspondance
      do i = 1, size(my_infiltration%land )
          infil%land( i )  =  my_infiltration%land( i )
      end do

      !define values for each paches
      do i = 1, nland
        infil%x_min(i) = my_infiltration%x_min(i)
        infil%x_max(i) = my_infiltration%x_max(i)
        infil%y_min(i) = my_infiltration%y_min(i)
        infil%y_max(i) = my_infiltration%y_max(i)
      end do
endif


!LEO: to check for infil, bathy seems to work
!#ifdef USE_MPI
!    if (bc_infil .ne. 0) then
!          call swap_vec_r  ( infil%land , swap_index( 1 : mesh%nc ) )
!          call reallocate_i( infil%land ,                 mesh%nc   )
!    endif
!#endif


END SUBROUTINE my_infiltration_2_fortran






! use variable param_model%bathy_cell (wrapped varible)
! to set up fortran variables (bathy_cell)

SUBROUTINE my_bathy_2_fortran() !(my_param_model)

implicit none

!    type( param_model ), intent(in   )  ::  my_param_model

   ! bathy cell allocated in build_mesh()
!      bathy_cell(:) = my_param_model%bathy_cell(:)



    !===== previously in subroutine



!#write(*,*) "mesh_type in Initial() = ", mesh_type
!#write(*,*) "Bathymetry Initialization"
!~    if ( mesh_type == 'basic' .or. mesh_type == 'gmsh' ) then
!~ 		write(*,*) "mesh_type basic or gmsh"
!~     if ( .not. allocated( bathy_cell ) ) then
!~        !===================================================================================================================!
!~ !#		write(*,*) ".not. allocated( bathy_cell ) "
!~ 		allocate( bathy_cell( mesh%nc + mesh%ncb ) )
!~     end if

!~ !		write(*,*) "z cells "
!~       do i = 1,mesh%nc
!~          bathy_cell(i)  =  bathy_user( mesh%cell(i)%grav%x , mesh%cell(i)%grav%y )

!~       end do



!~ !#		write(*,*) "cal mpi stuff "
!~       call com_var_r( bathy_cell , mesh )

!~ !#		write(*,*) "z ghost cells "
!~       do ie = 1,mesh%neb

!~          i = mesh%edge( mesh%edgeb(ie)%ind )%cell(2)
!~          j = mesh%edge( mesh%edgeb(ie)%ind )%cell(1)

!~          bathy_cell(i)  =  bathy_user( mesh%cellb(ie)%grav%x , mesh%cellb(ie)%grav%y )

!~       end do

!~    else if ( mesh_type == 'dassflow' ) then
!#		write(*,*) "mesh_type dassflow"

!---------------------- TOADD --------------------------

! my_bathy_to_fortran





!~    end if
!~ ! #  		write(*,*) "write"!
!~ !#write(*,*) size(bathy_cell), i, j, ie, mesh%edgeb(ie)%ind

!     do ie = 1,mesh%neb
!
!        i = mesh%edge( mesh%edgeb(ie)%ind )%cell(2)
!        j = mesh%edge( mesh%edgeb(ie)%ind )%cell(1)
!
! !       if ( mesh%edgeb(ie)%typlim == 'wall' ) print*, size(bathy_cell),i,j,ie,mesh%edgeb(ie)%ind
!       if ( mesh%edgeb(ie)%typlim == 'wall' ) bathy_cell(i)  =  bathy_cell(j)
!
!     end do

END SUBROUTINE my_bathy_2_fortran

SUBROUTINE my_phys_desc_2_fortran(my_phys_desc)

implicit none

   type( input_data ), intent(in   )  ::  my_phys_desc
!
     allocate( phys_desc%soil( size(my_phys_desc%soil) ) )
     allocate( phys_desc%surf( size(my_phys_desc%surf) ) )
     allocate( phys_desc%structures( size(my_phys_desc%structures) ) )

      do i = 1,size(my_phys_desc%soil)
          phys_desc%soil(i)%clay = my_phys_desc%soil(i)%clay
          phys_desc%soil(i)%silt = my_phys_desc%soil(i)%silt
          phys_desc%soil(i)%sand = my_phys_desc%soil(i)%sand
      enddo
      do i = 1,size(my_phys_desc%surf)
          phys_desc%surf(i)%imperm = my_phys_desc%surf(i)%imperm
          phys_desc%surf(i)%Dmax   = my_phys_desc%surf(i)%Dmax
      enddo
      do i = 1,size(my_phys_desc%structures)
          phys_desc%structures(i)%s_type = my_phys_desc%structures(i)%s_type
          phys_desc%structures(i)%true_x = my_phys_desc%structures(i)%true_x
          phys_desc%structures(i)%true_y = my_phys_desc%structures(i)%true_y
          phys_desc%structures(i)%name   = my_phys_desc%structures(i)%name
      end do

END SUBROUTINE my_phys_desc_2_fortran

END SUBROUTINE Initial
