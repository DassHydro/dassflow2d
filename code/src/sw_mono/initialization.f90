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
SUBROUTINE Initial( dof0, mesh, my_friction, my_infiltration, my_param_model, my_phys_desc, my_bc)

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
   type( bcs ), intent(in   )  ::  my_bc

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   character(len=lchar)  ::  filename

   ! array of filenames for error identification (20 file names of 30 characters max each)
   character(len=20), dimension(20) :: filenames
   ! error message for read errors / eof
   character(len=512) :: err_msg
   ! nb of opened file, line read
   integer(ip)  ::  file_nb, line_read

   integer(ip)  ::  icod, i_loc, l, spatial_index, num_bc

   integer(ip), dimension(:), allocatable  ::  index_bathy_min
   real(rp), dimension(:), allocatable  ::  bathy_min , bathy_min_glob, temp

!     call my_bathy_2_fortran() !(my_param_model)

     if (allocated(my_friction%manning))call my_friction_2_fortran(my_friction) ! propagate definition of friction from fortran to manning,

     if (bc_infil .ne. 0) call my_infiltration_2_fortran(my_infiltration)

     if (allocated(my_phys_desc%soil)) call my_phys_desc_2_fortran(my_phys_desc)

     if (allocated(my_bc%rain)) then
         call my_bc_2_fortran(my_bc)
     elseif (.not. allocated(bc%rain)) then 
        allocate( bc%rain(1) )
        allocate( bc%rain(1)%q(1) )
        allocate( bc%rain(1)%t(1) )

        bc%rain(1)%q = 0._rp
        bc%rain(1)%t = 0._rp
        bc%rain(1)%cumul = 0._rp
        bc%rain(1)%qin = 0._rp
     endif
     
#ifdef USE_MPI

   call swap_vec_i  ( land , swap_index( 1 : mesh%nc ) )
   call reallocate_i( land ,                 mesh%nc   )

   if (bc_rain == 1) then
         call swap_vec_i  ( bc%rain_land , swap_index( 1 : mesh%nc ) )
         call reallocate_i( bc%rain_land ,                 mesh%nc   )
   endif

   if (bc_infil .ne. 0) then
         call swap_vec_i  ( infil%land , swap_index( 1 : mesh%nc ) )
         call reallocate_i( infil%land ,                 mesh%nc   )
   endif
   
   if (allocated(phys_desc%soil_land)) then
         call swap_vec_i  ( phys_desc%soil_land , swap_index( 1 : mesh%nc ) )
         call reallocate_i( phys_desc%soil_land ,                 mesh%nc   )
         
!          if (use_ptf == 1) then
!             call swap_vec_i  ( phys_desc%ptf_land , swap_index( 1 : mesh%nc ) )
!             call reallocate_i( phys_desc%ptf_land ,                 mesh%nc   )
!          endif

   endif

   call swap_vec_r  ( bathy_cell , swap_index( 1 : mesh%nc + mesh%ncb ) )
   
#endif

   !===================================================================================================================!
   !  DOF Initialization --> dof pushed from python
   ! just keep parralel order com_dof
   !===================================================================================================================!

     filenames(1) = "ic.bin"
     inquire( file = 'ic.bin'      , exist = file_exist(1) )
     filenames(2) = "restart.bin"
     inquire( file = 'restart.bin' , exist = file_exist(2) )
     filenames(3) = "dof_init.txt"
     inquire( file = 'dof_init.txt', exist = file_exist(3) )
     filenames(4) = "zs_init.txt"
     inquire( file = 'zs_init.txt' , exist = file_exist(4) )
     
     if      ( file_exist(1) ) then

        file_nb = 1
        open(10,file='ic.bin',form='unformatted',status='old',access='direct',recl=3*length_real)
        line_read = 1

        read(10,rec=1, err=100) tc0

        tc0 = 0._rp

     else if (file_exist(3)) then
        file_nb = 3
        open(10,file='dof_init.txt',status='old',form='formatted')
        line_read = 1
        do i = 1,mesh%nc
  			read(10,*, err=100, end=100) dof0%h(i), dof0%u(i), dof0%v(i)
         line_read = line_read + 1
      end do


     else if (file_exist(4)) then

       file_nb = 4
       open(10,file='zs_init.txt',status='old',form='formatted')
  		 line_read = 1
       do i = 1,mesh%nc
  			read(10,*, err=100, end=100) dof0%h(i), dof0%u(i), dof0%v(i)
         line_read = line_read + 1
      end do
        
     else if ( file_exist(2) ) then

       file_nb = 2
       open(10,file='restart.bin',form='unformatted',status='old',access='direct',recl=3*length_real)
       line_read = 1
        read(10,rec=1, err=100) tc0
        line_read = line_read + 1
        if ( abs( ts - tc0 ) < zerom ) call Stopping_Program_Sub( 'End of simulation time reached' )

     end if

   call com_dof( dof0 , mesh )

   !===================================================================================================================!
   !  Boundary Condition Initialization
   !===================================================================================================================!

   allocate( bc%inflow ( 2 * mesh%neb ) )
   allocate( bc%outflow( 2 * mesh%neb ) )
   allocate( bc%sum_mass_flux( bc%nb ) )


#ifdef USE_HYDRO
   !===================================================================================================================!
   !  Reading input data for GR4 module
   !===================================================================================================================!

   filenames(5)="GR4params.txt"
   inquire( file = 'GR4params.txt' , exist = file_exist(1) )

   call mpi_wait_all

   if ( file_exist(1) ) then

      file_nb = 5
      open(10,file='GR4params.txt',status='old')
      line_read = 1

      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100) bc%nb_gr4in
      line_read = line_read + 1

      if (.not. allocated(bc%gr4)) allocate( bc%gr4( bc%nb_gr4in ))

      if (use_Qobs_gr4 == 1) then
         if (.not. allocated(stationQ)) allocate( stationQ(bc%nb_gr4in))
      else
        write(*,*) "Attempt to allocate already allocated variable stationQ when reading GR4params.txt"
      endif

      do i=1,bc%nb_gr4in

        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100) bc%gr4( i )%params(:)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100) bc%gr4( i )%state(:)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100) bc%gr4( i )%cell_id, bc%gr4( i )%surf
        line_read = line_read + 1

      enddo

      close(10)

      do i=1,bc%nb_gr4in

       write(filename,'("GR4forcings_",1I1.1,".txt")') i

      file_nb = 5
      open(10,file=,status='old',form='formatted')
      line_read = 0

            read(10,*, err=100, end=100)
            line_read = line_read + 1
            read(10,*, err=100, end=100) l
            line_read = line_read + 1

            allocate( bc%gr4( i )%t( l ))
            allocate( bc%gr4( i )%P( l ))
            allocate( bc%gr4( i )%E( l ))
            allocate( bc%gr4( i )%Q( l ))

        do j=1,l
            read(10,*, err=100, end=100) bc%gr4( i )%t( j ), bc%gr4( i )%P( j ), bc%gr4( i )%E( j )
            line_read = line_read + 1
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
   filenames(6)="GR4_PEwarmup.txt"
   inquire( file = 'GR4_PEwarmup.txt' , exist = file_exist(1) )

   call mpi_wait_all

   if ( file_exist(1) ) then

      file_nb = 6
      open(10,file='GR4_PEwarmup.txt',status='old')
      line_read = 0

      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100) l
      line_read = line_read + 1

      if (l .ne. bc%nb_gr4in) write(*,*) 'There are ',l,' GR4 warmup catchements and ',&
     bc%nb_gr4in, ' GR4 launch catchments'

      do i=1,bc%nb_gr4in

        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100) l
        line_read = line_read + 1

        allocate( bc%gr4( i )%P0( l ))
        allocate( bc%gr4( i )%E0( l ))

        do j=1,l
         read(10,*, err=100, end=100) bc%gr4( i )%P0( j ), bc%gr4( i )%E0( j )
         line_read = line_read + 1
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

   filenames(7)="hydrograph.txt"
   inquire( file = 'hydrograph.txt' , exist = file_exist(1) )

   call mpi_wait_all

   if ( file_exist(1) ) then

      file_nb = 7
      open(10,file='hydrograph.txt',status='old')
      line_read = 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)

      line_read = line_read + 1
      read(10,*, err=100, end=100) bc%nb_in

    allocate( bc%hyd( bc%nb_in ) )

      do i = 1,bc%nb_in
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)

        line_read = line_read + 1
        read(10,*, err=100, end=100) j

         allocate( bc%hyd( i )%t( j ) )
         allocate( bc%hyd( i )%q( j ) )


         do k = 1,j
            line_read = line_read + 1
            read(10,*, err=100, end=100) bc%hyd( i ) %t( k ) , bc%hyd( i )%q( k )
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
   filenames(8)="hpresc.txt"
   inquire( file = 'hpresc.txt' , exist = file_exist(1) )

   call mpi_wait_all

   if ( file_exist(1) ) then

      file_nb = 8
      open(10,file='hpresc.txt',status='old')
      line_read = 1

      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1

      read(10,*, err=100, end=100) bc%nb_out
      line_read = line_read + 1

    allocate( bc%hpresc( bc%nb_out ) )

      do i = 1,bc%nb_out
		! 3 comment lines of the file
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
		! size of hpresc
        read(10,*, err=100, end=100) j
        line_read = line_read + 1

         allocate( bc%hpresc( i )%t( j ) )
         allocate( bc%hpresc( i )%h( j ) )

		! read all available value
         do k = 1,j
            read(10,*, err=100, end=100) bc%hpresc( i ) %t( k ) , bc%hpresc( i )%h( k )
            line_read = line_read + 1
         end do
      end do
		close(10)

    end if


   !===================================================================================================================!
   !  Loading/Creating hpresc File
   !===================================================================================================================!

   filenames(9)="zspresc.txt"
   inquire( file = 'zspresc.txt' , exist = file_exist(2) )

   call mpi_wait_all

   if ( file_exist(2) ) then

      file_nb = 9
      open(10,file='zspresc.txt',status='old')
      line_read = 1

      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1

      read(10,*, err=100, end=100) bc%nb_out
      line_read = line_read + 1

    allocate( bc%zspresc( bc%nb_out ) )

      do i = 1,bc%nb_out
		! 3 comment lines of the file
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
		! size of zspresc
        read(10,*, err=100, end=100) j
        line_read = line_read + 1

         allocate( bc%zspresc( i )%t( j ) )
         allocate( bc%zspresc( i )%z( j ) )

		! read all available value
         do k = 1,j
            read(10,*, err=100, end=100) bc%zspresc( i ) %t( k ) , bc%zspresc( i )%z( k )
            line_read = line_read + 1
         end do
      end do
		close(10)

    end if

   
   !===================================================================================================================!
   !  Loading/Creating rain File
   !===================================================================================================================!

   if (bc_rain  == 1) then

   filenames(10)="rain.txt"
   inquire( file = 'rain.txt' , exist = file_exist(1) )

   call mpi_wait_all

   if ( file_exist(1) ) then

      file_nb = 10
      open(10,file='rain.txt',status='old')
      line_read = 1

      if (.not. allocated(bc%rain_land)) allocate(bc%rain_land(mesh%nc)) !in case rain.txt is read directly

      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100) bc%nb_rn
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1
      read(10,*, err=100, end=100)
      line_read = line_read + 1

      allocate(bc%rain(bc%nb_rn))

      do i = 1,bc%nb_rn
         read(10,*, err=100, end=100) bc%rain(i)%x_min, bc%rain(i)%x_max, bc%rain(i)%y_min, bc%rain(i)%y_max
         line_read = line_read + 1
      enddo

      do i = 1,bc%nb_rn

         read(10,*, err=100, end=100)
         line_read = line_read + 1
         read(10,*, err=100, end=100)
         line_read = line_read + 1
         read(10,*, err=100, end=100)
         line_read = line_read + 1
         read(10,*, err=100, end=100) j
         line_read = line_read + 1

         allocate( bc%rain(i)%t(j) )
         allocate( bc%rain(i)%q(j) )

         do k = 1,j
            read(10,*, err=100, end=100) bc%rain(i)%t(k) , bc%rain(i)%q(k)
            line_read = line_read + 1
            bc%rain(i)%q(k) = bc%rain(i)%q(k) / 1000._rp / 3600._rp
         enddo
      enddo
      close(10)
      
    !    ===================================================================================================================!
    !     Rain condition Initialization
    !    ===================================================================================================================!

      do i=1,mesh%nc

         do k=1,bc%nb_rn

              call spatial_index_fromxy(mesh, bc%rain(k)%x_min, bc%rain(k)%x_max,&
                                              bc%rain(k)%y_min, bc%rain(k)%y_max, spatial_index)

              if (spatial_index .ne. 0) exit

         enddo

      enddo

!         Output cell/rain attribution
        open(10,file='rain_post.dat',status='replace',form='formatted')

        write(10,*) '# Gnuplot DataFile Version'
        write(10,*) '# id x y rain_group'

        do i=1,mesh%nc !ADD ID_CELL AFTER MERGE

                    write(10,'(i5,2ES15.8,i3)') i, mesh%cell(i)%grav%x    , &
                                        mesh%cell(i)%grav%y    , &
                                         bc%rain_land(i)
        end do

        close(10)


    else
         mesh%cell(:)%rain = 0_ip
         write(*,*) "WARNING: Rain was not initialized from rain.txt. You need to initialize it through init_bc."
!          call Stopping_Program_Sub( 'File rain.txt not provided ...')
    endif

   end if

!    ===================================================================================================================!
!     Reading infiltration parameters and localization
!    ===================================================================================================================!

!    infil%nland = 0

   if (( mesh_type == 'dassflow' ) .and. ( .not. allocated(infil%land) )) then ! This is skipped if infiltration_initilise was called through python earlier

    allocate( infil%land( mesh%nc ) )

    filenames(11)="land_uses_GA.txt"
    inquire( file = 'land_uses_GA.txt' , exist = file_exist(1) )
    filenames(12)="land_uses_SCS.txt"
    inquire( file = 'land_uses_SCS.txt', exist = file_exist(2) )

    if ( file_exist(1) ) then

        file_nb = 11
        open(10,file='land_uses_GA.txt',status='old',form='formatted')
        line_read = 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100) infil%nland
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1

        allocate( infil%h_infil_max( infil%nland ) )
        infil%h_infil_max( : ) = 9999999._rp
        allocate( infil%GA( infil%nland ) )
        allocate( infil%SCS( 1 ) )
        infil%SCS( 1 )%lambdacn = 0._rp
        infil%SCS( 1 )%CN = 0._rp

        allocate (  infil%coord( 4, infil%nland )   )


        do i = 1,infil%nland
            read(10,*, err=100, end=100) k , infil%GA(k)%Ks , infil%GA(k)%PsiF ,&
                            infil%GA(k)%DeltaTheta
            line_read = line_read + 1
        end do

        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1

        do i = 1,infil%nland
!                     read(10,*, err=100, end=100) infil%x_min(i), infil%x_max(i), &
!                                infil%y_min(i), infil%y_max(i)
                      read(10,*, err=100, end=100) infil%coord(1,i), infil%coord(2,i), &
                                 infil%coord(3,i), infil%coord(4,i)
                      line_read = line_read + 1
        end do

        close(10)

    elseif ( file_exist(2) ) then

            file_nb = 12
            open(10,file='land_uses_SCS.txt',status='old',form='formatted')
            line_read = 1
            read(10,*, err=100, end=100)
            line_read = line_read + 1
            read(10,*, err=100, end=100)
            line_read = line_read + 1
            read(10,*, err=100, end=100)
            line_read = line_read + 1
            read(10,*, err=100, end=100) infil%nland
            line_read = line_read + 1
            read(10,*, err=100, end=100)
            line_read = line_read + 1
            read(10,*, err=100, end=100)
            line_read = line_read + 1
            read(10,*, err=100, end=100)
            line_read = line_read + 1

            allocate( infil%h_infil_max( infil%nland ) )
            infil%h_infil_max( : ) = 9999999._rp
            allocate( infil%SCS( infil%nland ) )
            allocate( infil%GA( 1 ) )

            infil%GA( 1 )%Ks = 0._rp
            infil%GA( 1 )%DeltaTheta = 0._rp
            infil%GA( 1 )%PsiF = 0._rp

            allocate (  infil%coord( 4, infil%nland )   )


            do i = 1,infil%nland
                read(10,*, err=100, end=100) k , infil%SCS(k)%lambdacn , infil%SCS(k)%CN
                line_read = line_read + 1
            end do

        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1
        read(10,*, err=100, end=100)
        line_read = line_read + 1

        do i = 1,infil%nland
                    read(10,*, err=100, end=100) infil%coord(1,i), infil%coord(2,i), &
                            infil%coord(3,i), infil%coord(4,i)
                    line_read = line_read + 1
        end do

        close(10)

    else

        write(*,*) "WARNING: Infiltration was not initialized from land_uses_XX.txt. You need to initialize it through init_infiltration."

    endif

    !===================================================================================================================!
    !  Handle land parameters attribution to cells from (x,y) tiles
    !===================================================================================================================!

    if ( file_exist(1) .or. file_exist(2)) then

       do i=1,mesh%nc
           do k=1,infil%nland

              call spatial_index_fromxy(mesh, infil%coord(1,k), infil%coord(2,k),&
                                              infil%coord(3,k), infil%coord(4,k), spatial_index)

              if (spatial_index .ne. 0) exit

           enddo
!            write(*,*) i,"spatial_index INFIL",spatial_index
           infil%land(i) = spatial_index

       enddo

        ! Output cell/infiltration attribution
        open(10,file='infil_post.dat',status='replace',form='formatted')

        write(10,*) '# Gnuplot DataFile Version'
        write(10,*) '# id x y infil_group'

        do i=1,mesh%nc !ADD ID_CELL AFTER MERGE
                    write(10,'(i5,2ES15.8,i3)') i, mesh%cell(i)%grav%x    , &
                                        mesh%cell(i)%grav%y    , &
                                        infil%land(i)
        end do

        close(10)
    endif

   endif

   
   !===================================================================================================================!
   !  Read geometry parameters if needed
   !===================================================================================================================!  

   if (use_xsshp == 1) then
     filenames(13)="geometry_params.txt"
     inquire( file = 'geometry_params.txt' , exist = file_exist(1) )
     if ( file_exist(1) ) then

       file_nb = 13
       open(20,file='geometry_params.txt',status='old',form='formatted')
       line_read = 1

       read(20,*, err=100, end=100)
       line_read = line_read + 1
       read(20,*, err=100, end=100) k
       line_read = line_read + 1
       read(20,*, err=100, end=100)
       line_read = line_read + 1

       allocate(XSshape(k))
       allocate(slope_x(k))
       allocate(slope_y(k))

       do i = 1,k
        read(20,*, err=100, end=100) XSshape(i)%xleft, XSshape(i)%xcenter, XSshape(i)%xright, XSshape(i)%s, XSshape(i)%hmax, &
        XSshape(i)%topz, slope_x(i), slope_y(i)
        line_read = line_read + 1
       enddo

       close(20)

     else

       allocate(XSshape(1))
       allocate(slope_x(1))
       allocate(slope_y(1))
        write(*,*) "WARNING: you used the option for parameterized bathymetry (use_xsshp == 1), but you did not provide the parameter file (geometry_params.txt)"

     endif

   else

    allocate(XSshape(1))
    allocate(slope_x(1))
    allocate(slope_y(1))
    slope_x(1) = 0_ip
    slope_y(1) = 0_ip

   endif
   
   
   !===================================================================================================================!
   !  Loading rating curve
   !===================================================================================================================!
   filenames(14)="rating_curve.txt"
   inquire( file = 'rating_curve.txt' , exist = file_exist(1) )

  call mpi_wait_all

   if ( file_exist(1) ) then

      file_nb = 14
		open(10,file='rating_curve.txt',status='old')
      line_read = 1

		read(10,*, err=100, end=100)
      line_read = line_read + 1
		read(10,*, err=100, end=100)
      line_read = line_read + 1
		read(10,*, err=100, end=100)
      line_read = line_read + 1

		read(10,*, err=100, end=100) bc%nb_out
      line_read = line_read + 1

		allocate( bc%rat( bc%nb_out ) )

		do i = 1,bc%nb_out

			read(10,*, err=100, end=100)
         line_read = line_read + 1
			read(10,*, err=100, end=100)
         line_read = line_read + 1
			read(10,*, err=100, end=100)
         line_read = line_read + 1

			read(10,*, err=100, end=100) j , bc%rat( i )%z_rat_ref
         line_read = line_read + 1

			allocate( bc%rat( i )%h( j ) )
			allocate( bc%rat( i )%q( j ) )

			do k = 1,j
				read(10,*, err=100, end=100) bc%rat( i )%h( k ) , bc%rat( i )%q( k )
            line_read = line_read + 1
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

!          call mpi_allgather_r( bathy_min(i) , bathy_min_glob )

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

   filenames(15)="obs.txt"
   inquire( file = 'obs.txt' , exist = file_exist(1) )

   if ( file_exist(1) ) then
      !================================================================================================================!
      !  Opening Data File concerning Stations and Sections Recording at prescribed frequency
      !================================================================================================================!

      file_nb = 15
      open(10,file='obs.txt',form='formatted',status='old')
      line_read = 1
      !================================================================================================================!
      !  Treating First Time Stations or Sections Data Record
      !================================================================================================================!

      buffer = '' ; icod = 1

      do while ( buffer(1:8) /= 'stations'          .and. &
                 buffer(1:17) /= 'stations_with_grp' .and. &
                 buffer(1:8) /= 'sections'          .and. &
                 buffer(1:10) /= 'stations_Q'          .and. icod >= 0 )



		read(10,'(A)',iostat=icod, err=100) buffer
      line_read = line_read + 1

      end do
      if ( icod >= 0 ) call read_obs_file(line_read)
		!================================================================================================================!
		!  Treating Second Time Stations or Sections Data Record
		!================================================================================================================!

		buffer = '' ; icod = 1

      do while ( buffer(1:8) /= 'stations'          .and. &
                 buffer(1:17) /= 'stations_with_grp' .and. &
                 buffer(1:8) /= 'sections'          .and.&
                 buffer(1:10) /= 'stations_Q'          .and. icod >= 0 )

		read(10,'(A)',iostat=icod, err=100) buffer
      line_read = line_read + 1

! 		write(*,*) buffer

		end do

		if ( icod >= 0 ) call read_obs_file(line_read)
		!================================================================================================================!
		!  Treating Thrid Time Stations or Sections Data Record
		!================================================================================================================!

		buffer = '' ; icod = 1

      do while ( buffer(1:8) /= 'stations'          .and. &
                buffer(1:17) /= 'stations_with_grp' .and. &
                 buffer(1:8) /= 'sections'          .and.&
                 buffer(1:10) /= 'stations_Q'          .and. icod >= 0 )

		 read(10,'(A)',iostat=icod, err=100) buffer
       line_read = line_read + 1

! 		write(*,*) buffer
		end do

		if ( icod >= 0 ) call read_obs_file(line_read)


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
   filenames(16)="param_obs.txt"
   inquire( file = 'param_obs.txt' , exist = file_exist(1) )

   if ( w_obs == 1 .and. file_exist(1) ) then

      file_nb = 16
      open(10,file='param_obs.txt',form='formatted',status='old')
      line_read = 1
      read(10,'(A)',iostat=icod, err=100, end=100) buffer
      line_read = line_read + 1

      do iobs = 1,size( station )
         read(10,*, err=100, end=100) station( iobs )%length , station(iobs)%dt_offset, temp,temp,temp,temp
         line_read = line_read + 1
      end do
      close(10)

   endif


   !===================================================================================================================!
   !  Array dt obs
   !===================================================================================================================!

   if (( w_obs == 1 ) .and. (( use_Zobs == 1 ) .or. ( use_UVobs == 1 )))then


     do iobs = 1,size(station)

        station( iobs )%nb_dt = int( ( ts - station(iobs)%dt_offset ) / station( iobs )%dt ) + 1_ip

        allocate( station( iobs )%dt_obs( station( iobs )%nb_dt ) )

        station( iobs )%dt_obs( 1 ) = station( iobs )%dt_offset

        do i = 2, station( iobs )%nb_dt

          station( iobs )%dt_obs( i ) = station( iobs )%dt_obs( i - 1 ) + station( iobs )%dt ! Create array of observation time ( offset + repetivity satellite)

        end do

        station( iobs )%ind_t = 1_ip

    enddo

  end if

   !===================================================================================================================!
   !  Initialize obs
   !===================================================================================================================!

   ! Allocate and initialize innovation vector for WSE observations
   ! innovW variable is unused

   if (( use_obs == 1 ) .and. ( use_Zobs == 1 )) then

      call read_stations()

      if (allocated( innovation )) deallocate( innovation )
      if (allocated( innovW ))     deallocate( innovW )

	  allocate( innovation ( size( station ) ) )
	  allocate( innovW (     size( station ) ) )

      do iobs = 1,size( station )

         innovation ( iobs )%nb_dt  =  size( station( iobs )%t(:) )
		 innovW (     iobs )%nb_dt  =  size( station( iobs )%t(:) )

         innovation ( iobs )%nb_dx  =  1
         innovW (     iobs )%nb_dx  =  1


		if ( allocated(innovation ( iobs )%diff) ) then
				write(*,*) "Already allocated innovation ( iobs )%diff"
		else
				allocate( innovation ( iobs )%diff( innovation( iobs )%nb_dt * &
                                            innovation( iobs )%nb_dx ) )
				allocate( innovW(      iobs )%diff( innovW( iobs )%nb_dt * &
                                            innovW(     iobs )%nb_dx ) )
		end if

         innovation ( iobs )%diff(:)  =  0._rp
         innovW(      iobs )%diff(:)  =  0._rp


      end do

  else

      if (allocated( innovation )) deallocate( innovation )
      if (allocated( innovW     )) deallocate( innovW )
      allocate( innovation( 1 ) )
      allocate( innovW(     1 ) )

 endif

  ! Allocate and initialize innovation vector for velocity observations
  if (( use_obs == 1 ) .and. ( use_UVobs == 1 )) then

      if ( use_Zobs == 0 ) call read_stations()

      if (allocated( innovUV )) deallocate( innovUV )

	  allocate( innovUV( size( station ) ) )

      do iobs = 1,size( station )

		 innovUV( iobs )%nb_dt  =  size( station( iobs )%t(:) )
         innovUV( iobs )%nb_dx  =  1

		if (allocated ( innovUV ( iobs )%diff ) ) then
				write(*,*) "Already allocated innovUV ( iobs )%diff"
		else
                allocate( innovUV( iobs )%diff( innovUV( iobs )%nb_dt * &
                                                innovUV( iobs )%nb_dx ) )
		end if

        innovUV( iobs )%diff(:)  =  0._rp

      end do

  else

      if ( allocated( innovUV ) ) deallocate( innovUV )
      allocate( innovUV( 1 ) )

 endif

   ! Allocate and initialize innovation vector for flow observations
   ! Either use_Qobs or usse_Qobs_gr4 can be used, not both. This should correspond to most inverse setups.
   if (( use_obs == 1 ) .and. ((use_Qobs == 1) .or. (use_Qobs_gr4 == 1))) then

      call read_stationsQ()

      if ( allocated( innovQ ) ) deallocate( innovQ )

      allocate( innovQ( size( stationQ ) * 2 ) ) ! The latter half of the innovQ vector is used only with use_NSE = 1.

      do iobs = 1,size( stationQ )

         innovQ( iobs )%nb_dt  =  stationQ( iobs )%nb_dt
         innovQ( iobs )%nb_dx  =  1

         if (allocated ( innovQ ( iobs )%diff ) ) then
				write(*,*) "Already allocated innovQ ( iobs )%diff"
		 else
                allocate( innovQ( iobs                    )%diff( innovQ( iobs )%nb_dt ))
                allocate( innovQ( iobs + size( stationQ ) )%diff( innovQ( iobs )%nb_dt ))
         endif

         innovQ( iobs )%diff(:)  =  0._rp

      end do

   else

      if ( allocated( innovQ ) ) deallocate( innovQ )
      allocate( innovQ( 1 ) )

   end if


   time(:) = 0._rp

   call display_mesh_cell(mesh)

   return

   100 write(err_msg, '(3A, I4)') "ERROR in file ", trim(filenames(file_nb))," at line ", line_read
       close(10)
       call f90wrap_abort(err_msg)
   101 write(err_msg, '(5A, I4)') "ERROR in file ", trim(filenames(1)), " or ", trim(filenames(2))," at line ", line_read
       close(10)
CONTAINS


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Reading Stations or Sections informations in obs.txt file and Creating associated structure
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE read_obs_file(line_read)

      implicit none

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      type( vec2d )  ::  pt_temp(2)

      real(rp)  ::  length

      integer(ip)  ::  nb_rec , nb_pt , pt , search_cell_inc_point , etq , cell

      ! read line number 
      integer(ip), intent(out)  ::  line_read
      integer(ip)  ::  line_station

      logical  ::  point_in_cell

      character(len=8)  ::  typ

      ! Name of obs file used
      character(len=20)  ::  obs_file  
      obs_file = "obs.txt"

      !================================================================================================================!
      !  Beginning Treatment of File 'obs.txt'
      !================================================================================================================!
      
      backspace 10
      line_read = line_read - 1

      read(10,*, err=100, end=100) buffer , nb_rec
      line_read = line_read + 1

      select case( trim(buffer) )

         !=============================================================================================================!
         !  Stations at One Point Case
         !=============================================================================================================!

		case( 'stations' )

            call alloc_or_realloc_station( station , nb_obs + nb_rec )

            do iobs = nb_obs+1 , nb_obs+nb_rec

               !=======================================================================================================!
               !  Allocate Point Coordinate (only one in this case)
               !=======================================================================================================!

               allocate( station( iobs )%pt(1) )

               !=======================================================================================================!
               !  Reading Stations Informations in 'obs.txt' File
               !=======================================================================================================!

               read(10,*, err=100, end=100) station( iobs )%pt(1)%coord%x , &
                          station( iobs )%pt(1)%coord%y , &
                          station( iobs )%dt            , &
                          station( iobs )%weight
               line_read = line_read + 1

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

         case( 'stations_Q' ) !Currently, stations are linked to a boundary condition index. Thus, (x,y) data defined below from obs.txt is not used.

          call alloc_or_realloc_stationQ( stationQ , nb_rec )

          do iobs = 1 , nb_rec

               !=======================================================================================================!
               !  Allocate Point Coordinate (only one in this case)
               !=======================================================================================================!

               allocate( stationQ( iobs )%pt(1) )

               !=======================================================================================================!
               !  Reading Statiostations_Qns Informations in 'obs.txt' File
               !=======================================================================================================!

               read(10,*, err=100, end=100) stationQ( iobs )%pt(1)%coord%x , &
                          stationQ( iobs )%pt(1)%coord%y , &
                          stationQ( iobs )%dt            , &
                          stationQ( iobs )%weight
               line_read = line_read + 1

               !=======================================================================================================!
               !  Searching Cells for Stations
               !=======================================================================================================!

               stationQ( iobs )%pt(1)%cell = search_cell_inc_point( mesh , stationQ( iobs )%pt(1)%coord )

               !=======================================================================================================!
               !  Writing Ouput Files to Visualize Stations
               !=======================================================================================================!

               write(buffer,'(A,I4.4)') 'res/obs/pos_stationQ_', iobs

               call write_station_pos_in_file( buffer , stationQ( iobs )%pt(1)%coord )

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

               read(10,*, err=100, end=100) etq , station( nb_obs + etq )%dt , station( nb_obs + etq )%weight
               line_read = line_read + 1

               write(buffer,'(A,I4.4,A)') 'station_grp_' , iobs - nb_obs , '.txt'

               nb_pt = count_lines( buffer ) - 1

               allocate( station( iobs )%pt( nb_pt ) )

               open(20,file=buffer,status='old',form='formatted')
               line_station = 1

               read(20,'(A)', err=101, end=101) typ
               line_station = line_station + 1

               do pt = 1,nb_pt

                  !====================================================================================================!
                  !  Reading Stations Informations in 'station_grp_xxxx.txt' File
                  !====================================================================================================!

                  if      ( trim(typ) == 'points' ) then
                     read(20,*, err=101, end=101) station( iobs )%pt( pt )%coord%x , &
                                station( iobs )%pt( pt )%coord%y
                     line_station = line_station + 1
                     !=================================================================================================!
                     !  Searching Cells for Stations
                     !=================================================================================================!

                     station( iobs )%pt( pt )%cell = search_cell_inc_point( mesh , station( iobs )%pt( pt )%coord )

                  else if ( trim(typ) == 'indexes' ) then

                     read(20,*, err=101, end=101) cell
                     line_station = line_station + 1

                     station( iobs )%pt( pt )%cell = cell
                     station( iobs )%pt( pt )%coord  =  mesh%cell( station( iobs )%pt( pt )%cell )%grav

                     if ( part( cell ) == proc ) then

!                         station( iobs )%pt( pt )%cell = mesh%inv_swap_index( cell )
write(*,*) " line commented temporarily ! station( iobs )%pt( pt )%cell = mesh%inv_swap_index( cell ) "

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

               read(10,*, err=100, end=100) pt_temp(1)%x , &
                          pt_temp(1)%y , &
                          pt_temp(2)%x , &
                          pt_temp(2)%y , &
                          nb_pt , &
                          section( iobs )%dt
               line_read = line_read + 1
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
      
      return
   
      
   100 write(err_msg, '(3A, I4)') "ERROR in file ", trim(obs_file)," at line ", line_read
       close(10)
       call f90wrap_abort(err_msg)
   101 write(err_msg, '(3A, I4)') "ERROR in file ", trim(buffer)," at line ", line_station
       close(20)
       call f90wrap_abort(err_msg)

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

!       allocate( land( size( my_friction%land ) ) )
      allocate( manning( my_friction%nland ) )
      allocate( manning_beta( my_friction%nland ) )

      ! loop on all cells to define patch correspondance
      do i = 1,size(my_friction%land)
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

   if( infil%nland  > 0) then

      allocate( infil%land( size( my_infiltration%land ) ) )
      allocate( infil%h_infil_max( size( my_infiltration%land ) ) )

      if (bc_infil == 1) then

        allocate( infil%GA ( size(my_infiltration%GA) ) )

        do i = 1, size(my_infiltration%GA)
          infil%GA(i)  = my_infiltration%GA(i)
        enddo

        allocate( infil%SCS ( 1 ) )
        infil%SCS  ( 1 )%lambdacn = 0._rp
        infil%SCS  ( 1 )%CN = 0._rp

      elseif (bc_infil == 2) then

        allocate( infil%SCS ( infil%nland ) )

        do i = 1, size(my_infiltration%SCS)
          infil%SCS(i) = my_infiltration%SCS(i)
        enddo

        allocate( infil%GA  ( 1 ) )
        infil%GA  ( 1 )%Ks = 0._rp
        infil%GA  ( 1 )%DeltaTheta = 0._rp
        infil%GA  ( 1 )%PsiF = 0._rp

      endif

!       allocate( infil%coord ( 4, infil%nland ))

      do i = 1, size(my_infiltration%land )
          infil%land( i )  =  my_infiltration%land( i )
          infil%h_infil_max( i )  =  my_infiltration%h_infil_max( i )
      end do

      do i = 1, size(my_infiltration%h_infil_max )
          infil%h_infil_max( i )  =  my_infiltration%h_infil_max( i )
      end do

      !define values for each patch
!       do i = 1, infil%nland
!         infil%coord(1,i) = my_infiltration%coord(1,i)
!         infil%coord(2,i) = my_infiltration%coord(2,i)
!         infil%coord(3,i) = my_infiltration%coord(3,i)
!         infil%coord(4,i) = my_infiltration%coord(4,i)
!       end do
   endif

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
     allocate( phys_desc%soil_land( size(my_phys_desc%soil_land) ) )

!      allocate( phys_desc%surf( size(my_phys_desc%surf) ) )
!      allocate( phys_desc%structures( size(my_phys_desc%structures) ) )


      do i = 1,size(my_phys_desc%soil)
          phys_desc%soil(i)%clay = my_phys_desc%soil(i)%clay
          phys_desc%soil(i)%silt = my_phys_desc%soil(i)%silt
          phys_desc%soil(i)%sand = my_phys_desc%soil(i)%sand
      enddo

      phys_desc%soil_land(:) = my_phys_desc%soil_land(:)
      
      if (use_ptf == 1) then
      
        allocate( PTF( size(my_phys_desc%ptf) ) )
        allocate( phys_desc%ptf_land( size(my_phys_desc%ptf_land) ) )
        
        do i = 1,size(my_phys_desc%ptf)
          PTF(i)%Kappa(:) = my_phys_desc%ptf(i)%Kappa(:)
        enddo
        
        phys_desc%ptf_nland = my_phys_desc%ptf_nland
        phys_desc%ptf_land(:) = my_phys_desc%ptf_land(:)
      
      endif
      
!       do i = 1,size(my_phys_desc%surf)
!       write(*,*) i
!           phys_desc%surf(i)%imperm = my_phys_desc%surf(i)%imperm
!           phys_desc%surf(i)%Dmax   = my_phys_desc%surf(i)%Dmax
!       enddo
!       do i = 1,size(my_phys_desc%structures)
!       write(*,*) i
!           phys_desc%structures(i)%s_type = my_phys_desc%structures(i)%s_type
!           phys_desc%structures(i)%true_x = my_phys_desc%structures(i)%true_x
!           phys_desc%structures(i)%true_y = my_phys_desc%structures(i)%true_y
!           phys_desc%structures(i)%name   = my_phys_desc%structures(i)%name
!       end do

END SUBROUTINE my_phys_desc_2_fortran

SUBROUTINE my_bc_2_fortran(my_bc)

    implicit none

    type( bcs ), intent(in   )  ::  my_bc

      allocate(bc%rain(my_bc%nb_rn))
      bc%nb_rn = my_bc%nb_rn

      do i = 1, my_bc%nb_rn
        allocate(bc%rain(i)%t(my_bc%nb_rn_t))
        allocate(bc%rain(i)%q(my_bc%nb_rn_t))
        bc%rain(i)%t(:) = my_bc%rain(i)%t(:)
        bc%rain(i)%q(:) = my_bc%rain(i)%q(:)
        bc%rain(i)%x_min = my_bc%rain(i)%x_min
        bc%rain(i)%y_min = my_bc%rain(i)%y_min
        bc%rain(i)%x_max = my_bc%rain(i)%x_max
        bc%rain(i)%y_max = my_bc%rain(i)%y_max
        bc%rain(i)%tile_index = my_bc%rain(i)%tile_index
        bc%rain(i)%cumul = my_bc%rain(i)%cumul
        bc%rain(i)%qin = 0._rp !my_bc%rain(i)%qin
      enddo

      allocate(bc%rain_land(size(my_bc%rain_land)))

      do i = 1,size(my_bc%rain_land)
        bc%rain_land(i) = my_bc%rain_land(i)
      enddo
    

END SUBROUTINE my_bc_2_fortran

END SUBROUTINE Initial
