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
!  Fill during simulation the Innovation Vector to compute at end the Cost Function
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE calc_innovation( dof,mesh )

   USE m_common
   USE m_mpi
   USE m_model
   USE m_obs

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type( msh ), intent(in)  ::  mesh

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  cell , searched_time , pt,N_average

   real(rp)  ::  h_mean, s_total, u_mean, v_mean

   !===================================================================================================================!
   !  Begin
   !===================================================================================================================!

    if ( use_obs == 1 ) then

        do iobs = 1, size( station )

            searched_time  =  innovation( iobs )%ind_t

            if ( searched_time > innovation( iobs )%nb_dt ) cycle

            if ( tc >= station( iobs )%t( searched_time ) ) then

                h_mean  = 0._rp
                s_total = 0._rp

                do pt = 1,size( station( iobs )%pt )

                    cell = station( iobs )%pt( pt )%cell

                    if ( cell < 0 ) cycle

                    !*****************************
                    ! H = 1/swet \int_swet Hdx
                    !*****************************
                    if ( dof%h( cell ) > 0 ) then  !test on water presence determining if cell is used for calculating h_average

                    h_mean  = h_mean  + ( dof%h( cell ) + bathy_cell( cell ) ) * mesh%cell(cell)%surf
                    s_total = s_total + mesh%cell(cell)%surf

                    endif


                    !*****************************
                    ! H = 1/sobs \int_sobs Hdx
                    !*****************************
        !             s_total= s_total + mesh%cell(cell)%surf
        !             h_mean = h_mean  + dof%h( cell )*mesh%cell(cell)%surf

                enddo

                call mpi_sum_r( h_mean )
                call mpi_sum_r( s_total )
                call mpi_sum_i( N_average )

                if (s_total>0) then

                    h_mean=h_mean / s_total

                end if

                innovation ( iobs )%diff( searched_time )    = (h_mean - station( iobs )%h( searched_time ))
                innovation ( iobs )%ind_t  =  innovation ( iobs )%ind_t + 1


            endif

        enddo

    endif

END SUBROUTINE calc_innovation


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Fill during simulation the Innovation Vector to compute at end the Cost Function
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

#ifdef USE_HYDRO
SUBROUTINE calc_innovQ_gr4( dof,mesh )

   USE m_common
   USE m_mpi
   USE m_model
   USE m_obs

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type( msh ), intent(in)  ::  mesh

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  cell , searched_time , pt,N_average

   real(rp)  ::  w_total!,w_mean

   !===================================================================================================================!
   !  Begin
   !===================================================================================================================!

   do iobs = 1,SIZE(stationQ)

      searched_time  =  innovQ( iobs )%ind_t

      if ( searched_time > innovQ( iobs )%nb_dt ) cycle

      if ( tc >= stationQ( iobs )%t( searched_time ) ) then

        !!!!!!!!!!!!!!!!!!!
        !Classical RMSE
        !!!!!!!!!!!!!!!!!!!        
! write(*,*) iobs, (bc%gr4( iobs )%Q( searched_time ) / bc%gr4(iobs)%surf * 3.6_rp) , stationQ( iobs )%Q( searched_time )
         innovQ ( iobs )%diff( searched_time )    = (bc%gr4( iobs )%Q( searched_time ) / bc%gr4(iobs)%surf * 3.6_rp) &
         - stationQ( iobs )%Q( searched_time ) !! mm/h, not m3/s !


        !!!!!!!!!!!!!!!!!!!
        !NSE
        !!!!!!!!!!!!!!!!!!

!         innovQ ( iobs )%diff( searched_time )    =  (bc%gr4( iobs )%Q( searched_time ) * 3.6_rp / bc%gr4( iobs )%surf &
!          - stationQ( iobs )%Q( searched_time )) !! mm/h, not m3/s !
!          
!         innovQ ( iobs + size(stationQ) )%diff( searched_time ) = stationQ( iobs )%Q( searched_time ) &
!         - sum(stationQ( iobs )%Q( : )) / size(stationQ( iobs )%Q( : ))
! ! 
!         innovQ ( iobs )%ind_t  =  innovQ ( iobs )%ind_t + 1

      end if

   end do
   
END SUBROUTINE calc_innovQ_gr4
#endif

SUBROUTINE calc_innovQ( dof,mesh )

   USE m_common
   USE m_mpi
   USE m_model
   USE m_obs

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type( msh ), intent(in)  ::  mesh

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  cell , searched_time , pt,N_average

   real(rp)  ::  w_total!,w_mean

   !===================================================================================================================!
   !  Begin
   !===================================================================================================================!

   do iobs = 1,SIZE(stationQ) !station number is also BCout number, tested for 1 BCout

      searched_time  =  innovQ( iobs )%ind_t

      if ( searched_time > innovQ( iobs )%nb_dt ) cycle

      if ( tc >= stationQ( iobs )%t( searched_time ) ) then

         innovQ ( iobs )%diff( searched_time )    = bc%sum_mass_flux( 3 ) - stationQ( iobs )%Q( searched_time )
         innovQ ( iobs )%ind_t  =  innovQ ( iobs )%ind_t + 1

      end if

   end do
   
END SUBROUTINE calc_innovQ




SUBROUTINE calc_innovUV( dof, mesh )

   USE m_common
   USE m_mpi
   USE m_model
   USE m_obs

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type( msh ), intent(in)  ::  mesh

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  cell , searched_time , pt,N_average

   real(rp)  ::  h_mean, s_total, u_mean, v_mean
character(50)  ::  filename
   !===================================================================================================================!
   !  Begin
   !===================================================================================================================!

   do iobs = 1, size( station )

      searched_time  =  innovUV( iobs )%ind_t

      if ( searched_time > innovUV( iobs )%nb_dt ) cycle

      if ( tc >= station( iobs )%t( searched_time ) ) then

         u_mean  = 0._rp
         v_mean  = 0._rp
         s_total = 0._rp

         do pt = 1,size( station( iobs )%pt )

            cell = station( iobs )%pt( pt )%cell

            if ( cell < 0 ) cycle

            if ( dof%h( cell ) > 0 ) then  !test on water presence determining if cell is used for calculating h_average

               u_mean  = u_mean  + dof%u( cell ) * mesh%cell(cell)%surf
               v_mean  = v_mean  + dof%v( cell ) * mesh%cell(cell)%surf
               s_total = s_total + mesh%cell(cell)%surf

            end if

         end do

         call mpi_sum_r( u_mean )
         call mpi_sum_r( v_mean )

         if (s_total>0) then
            u_mean = u_mean / s_total
            v_mean = v_mean / s_total
         end if

          innovUV ( iobs )%diff( searched_time ) =  (u_mean**2_rp + v_mean**2_rp)**(0.5_rp) - &
          (station( iobs )%u( searched_time )**2_rp + station( iobs )%v( searched_time )**2_rp)**(0.5_rp)

          innovUV ( iobs )%ind_t  =  innovUV ( iobs )%ind_t + 1

! Output x,y,obs_diff, cell_id
!
!         open(10,file='res/post/innovUV.dat',status='replace',form='formatted')
!
!         write(10,*) '# Gnuplot DataFile Version'
!         write(10,*) '# id_station x y id_cell innovUV'
!
!         do i=1,mesh%nc !ADD ID_CELL AFTER MERGE
!
!                     write(10,'(i5,2ES15.8,i3)') iobs, &
!                                                 mesh%cell(cell)%grav%x    , &
!                                                 mesh%cell(cell)%grav%y    , &
!                                                 cell, innovUV ( iobs )%diff( searched_time )
!         end do
!
!         close(10)

      endif

   enddo

END SUBROUTINE calc_innovUV


SUBROUTINE calc_innovW( dof,mesh )

   USE m_common
   USE m_mpi
   USE m_model
   USE m_obs

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type( msh ), intent(in)  ::  mesh

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  cell , searched_time , pt,N_average

   real(rp)  ::  w_total!,w_mean

   !===================================================================================================================!
   !  Begin
   !===================================================================================================================!

   do iobs = 1,size( station )

      searched_time  =  innovW( iobs )%ind_t

      if ( searched_time > innovW( iobs )%nb_dt ) cycle

      if ( tc >= station( iobs )%t( searched_time ) ) then


         w_total = 0._rp

         do pt = 1,size( station( iobs )%pt )

            cell = station( iobs )%pt( pt )%cell

            if ( cell < 0 ) cycle

            !if ( dof%h( cell ) > 0 ) then  !test on water presence determining if cell is used for calculating h_average
            !   w_total = w_total + mesh%cell(cell)%surf
            !end if
            w_total= w_total+ mesh%cell(cell)%surf*(max(dof%h(cell),0.0_rp)/max(dof%h(cell),0.0000000001_rp))

         end do

         call mpi_sum_r( w_total )



         if ( station( iobs )%length > 0) then
            w_total = w_total / station( iobs )%length
         else
            w_total= 0.0_rp
         endif

         innovW( iobs )%diff( searched_time )  = (w_total - station( iobs )%w( searched_time ))!
         innovW( iobs )%ind_t  =  innovW( iobs )%ind_t + 1

      end if

   end do

END SUBROUTINE calc_innovW


                                                                                                                 !<NOADJ

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Stations Output
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE write_stations( dof ,mesh )

   USE m_common
   USE m_mpi
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type( msh ), intent(in   )  ::  mesh

   !===================================================================================================================!
   !  Writing Stations Records in File
   !===================================================================================================================!


!    write(*,*) "INTO  write station"
   if ( .not. allocated( station ) .or. w_obs /= 1 ) then
!    write(*,*) allocated( station ), w_obs
!    write(*,*) ".not. allocated( station ) .or. w_obs /= 1 --> RETURN"
    return
   end if
   if ( w_tecplot == 1 ) then

      call sub_write( 'tecplot' ) !; return

   end if

   if ( w_gnuplot == 1 ) then

      call sub_write( 'gnuplot' ) !; return

   end if

!    write(*,*) "call  sub_write( 'bin' ) "
!   call sub_write( 'bin' )
!    write(*,*) "out  sub_write( 'gnuplot' ) "

!    WRITE(*,*) "out write_stations"

   CONTAINS


      SUBROUTINE sub_write( file_type )

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         character(len=*), intent(in)  ::  file_type

         !=============================================================================================================!
         !  Local Variables
         !=============================================================================================================!

         character(len=lchar)  ::  file_name

         integer(ip)  ::  iobs , cell , pt , N_average,searched_time

         real(rp)  ::  h_mean , u_mean, v_mean,s_total , w_meanb

         !=============================================================================================================!
         !  Begin Subroutine
         !=============================================================================================================!
        if ((use_Zobs == 1) .or. (use_UVobs == 1)) then

         do iobs = 1,size( station )

!         if (iobs .EQ. 1 ) then
! !            write(*,*) "iobs =", iobs
! !            write(*,*)"tc =",tc
!         end if
            !==========================================================================================================!
            !  Testing if Simulation Time match with Observation Time Step
            !==========================================================================================================!

            !if ( .not. test_dt_just_after( station( iobs )%dt ) ) cycle

            searched_time  =  station( iobs )%ind_t

            if ( searched_time > station( iobs )%nb_dt ) cycle
!write(*,*) "end_time_loop 3.3.2 -- into write_stations",end_time_loop



            if ( tc >= station( iobs )%dt_obs( searched_time ) ) then

               !==========================================================================================================!
               !  Creating File Name
               !==========================================================================================================!

               write(file_name,'(A,"_",I4.4)') 'res/obs/obs_station' , iobs

               file_name = file_name_ext( file_name , file_type )

               !==========================================================================================================!
               !
               !==========================================================================================================!

               ! creating wet cells file
               if ( iobs == 1 ) then
                    write( buffer,'(A,1ES23.15,A)') 'res/obs/wet_cells_',tc,'.txt'

                    inquire( file = buffer , exist = file_exist(1) )

                    if ( (.not. file_exist(1)) .and.(test_dt_just_after (dtw)) ) then

                       open(20,file=buffer,status='replace',form='formatted')

                    end if
               end if


               if ( station( iobs )%pt(1)%cell > 0 .and. all( is_file_open(:) /= file_name ) ) then
                  inquire( file = file_name , exist = file_exist(1) )
!write(*,*) "end_time_loop 3.3.3.2 ",end_time_loop

                  if ( .not. file_exist(1) ) then
                     if ( file_type /= 'bin' ) then
!write(*,*) "end_time_loop 3.3.3.4 ",end_time_loop

                        open(10,file=file_name,status='replace',form='formatted')
!write(*,*) "end_time_loop 3.3.3.5 ",end_time_loop

                        select case( file_type )

                           case( 'tecplot' )

                              write(10,'(A)') 'VARIABLES = "time" "h_mean" "u_mean" "v_mean" "w_mean"'

                           case default

                              write(10,'(A)') 'time h_mean u_mean v_mean w_mean'

                        end select
!write(*,*) "end_time_loop 3.3.3.6",end_time_loop

                     else
!write(*,*) "end_time_loop 3.3.3.7",end_time_loop

                        open(10,file=file_name,status='replace',form='unformatted')

                     end if
!write(*,*) "end_time_loop 3.3.3.8",end_time_loop

                  end if
!write(*,*) "end_time_loop 3.3.3.9",end_time_loop

                  close(10)

!write(*,*) "close file"
                  file_open_counter = file_open_counter + 1
!write(*,*) "end_time_loop 3.3.3.10",end_time_loop

                  is_file_open( file_open_counter ) = file_name
                  
!write(*,*) "file_open_counter", file_open_counter
!write(*,*) "file_name=", file_name
!write(*,*) "is_file_open(:)=", is_file_open(:)
!write(*,*) "end_time_loop 3.3.3.11",end_time_loop

               end if
!write(*,*) "end_time_loop 3.3.4 ",end_time_loop
               !==========================================================================================================!
               !
               !==========================================================================================================!

               h_mean = 0._rp
               u_mean = 0._rp
               v_mean = 0._rp
               w_meanb = 0._rp
               s_total= 0._rp
               N_average = 0

               do pt = 1,size( station( iobs )%pt )
					 
                  cell = station( iobs )%pt( pt )%cell

                  if ( cell < 0 ) cycle


                  !*****************************
                  ! H = 1/swet \int_swet Hdx
                  !*****************************
                  if ( dof%h( cell ) > 0 ) then  !test on water presence determining if cell is used for calculating h_average

                     u_mean = u_mean + dof%u( cell )
                     v_mean = v_mean + dof%v( cell )

                     h_mean = h_mean + (dof%h( cell )+bathy_cell( cell ) )*mesh%cell( cell)%surf
                     s_total = s_total+mesh%cell( cell)%surf
                  end if

                  !*****************************
                  ! H = 1/sobs \int_sobs Hdx
                  !*****************************
                  !s_total= s_total + mesh%cell(cell)%surf
                  !h_mean = h_mean  + dof%h( cell )*mesh%cell(cell)%surf


                  !*****************************
                  ! H = 1/sriver \int_sobs Hdx
                  !*****************************
                  !h_mean = h_mean  + dof%h( cell )*mesh%cell(cell)%surf


                  if (dof%h(cell)>0) then
                     w_meanb  =  w_meanb+mesh%cell(cell)%surf
                     N_average= N_average+1
                  end if

               end do
!write(*,*) "done loop on station( iobs )%pt( pt )"
               call mpi_sum_r( h_mean )
               call mpi_sum_r( u_mean )
               call mpi_sum_r( v_mean )
               call mpi_sum_r( w_meanb )
               call mpi_sum_r( s_total )
               call mpi_sum_i( N_average )

               !*****************************
               ! H = 1/sriver \int_sobs Hdx
               !*****************************

		if (N_average>0) then ! not divide by 0

                 h_mean = h_mean / s_total
! 		         u_mean = u_mean / real( N_average , 8) !real( size( station( iobs )%pt ) , 8 )
! 		         v_mean = v_mean / real( N_average , 8) !real( size( station( iobs )%pt ) , 8 )
		end if

               if (station(iobs)%length>0) then
                  w_meanb=w_meanb/(station(iobs)%length)
               else
                  w_meanb=0.0
               endif

               if ( station( iobs )%pt(1)%cell > 0 ) then

                  if ( file_type == 'bin' ) then

                     open(10,file=file_name,status='old',position='append',form='unformatted')

                     write(10) tc , h_mean , u_mean , v_mean ,w_meanb

                  else

                     open(10,file=file_name,status='old',position='append',form='formatted')

                     write(10,'(5ES23.15)') tc , h_mean , u_mean , v_mean ,w_meanb

                  end if

                  close(10)


                  !Begin

                  !write(file_name,'(A,"_",I4.4)') 'res/obs_station_wet' , iobs

                  !file_name = file_name_ext( file_name , file_type )

                  !open(10,file=file_name,status='replace',position='append',form='formatted')

                  !write(10,*) 'indexes'

                  !do pt = 1,size( station( iobs )%pt )
                  !   cell = station( iobs )%pt( pt )%cell
                  !   if (dof%h(cell)>0) then
                  !      write(10,*) cell
                  !   end if
                  !end do

                  !close(10)
                  !ENd


               end if

                station( iobs )%ind_t= station( iobs )%ind_t+1

            end if



          end do
        endif
       close(20)

      END SUBROUTINE sub_write

END SUBROUTINE write_stations


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Stations Input
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE read_stations

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_obs
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Reading Stations Records in File
   !===================================================================================================================!
! write(*,*) "in read_stations"
   if ( w_tecplot == 1 ) then

      call sub_read( 'tecplot' ) ; return

   end if

   if ( w_gnuplot == 1 ) then

      call sub_read( 'gnuplot' ) ; return

   end if

   # call sub_read( 'bin' )


   CONTAINS


      SUBROUTINE sub_read( file_type )

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         character(len=*), intent(in)  ::  file_type

         !=============================================================================================================!
         !  Local Variables
         !=============================================================================================================!

         character(len=lchar)  ::  file_name

         integer(ip)  ::  iobs , nbdt_obs , dt_obs , icod

         !=============================================================================================================!
         !  Begin Subroutine
         !=============================================================================================================!

         iobs = 0

         do

            !==========================================================================================================!
            !
            !==========================================================================================================!

            iobs = iobs + 1

            write(file_name,'(A,"_",I4.4)') 'obs/obs_station' , iobs

            file_name = file_name_ext( file_name , file_type )

            !==========================================================================================================!
            !
            !==========================================================================================================!

            inquire( file = file_name , exist = file_exist(1) )

            if ( .not. file_exist(1) ) then

                if (iobs <= size( station ) ) call Stopping_Program_Sub('Not enough number of observation files or not provided, check obs directory')

               exit

            else if ( iobs > size( station ) ) then

               call Stopping_Program_Sub( 'Mismatch between the number of files in obs directory and obs.txt file' )

            end if

          !  if ( file_type /= 'bin' ) then

               !=======================================================================================================!
               !
               !=======================================================================================================!

               nbdt_obs = count_lines( file_name ) - 1

               allocate( station( iobs )%t( nbdt_obs ) )
               allocate( station( iobs )%h( nbdt_obs ) )
               allocate( station( iobs )%u( nbdt_obs ) )
               allocate( station( iobs )%v( nbdt_obs ) )
               allocate( station( iobs )%w( nbdt_obs ) )

               !=======================================================================================================!
               !
               !=======================================================================================================!

               open(10,file=file_name,status='old',form='formatted')

               read(10,*)

               do dt_obs = 1,nbdt_obs

                  read(10,*) station( iobs )%t( dt_obs ) , &
                             station( iobs )%h( dt_obs ) , &
                             station( iobs )%u( dt_obs ) , &
                             station( iobs )%v( dt_obs ) , &
                             station( iobs )%w( dt_obs )

               end do

               close(10)

          !  end if

         end do

      END SUBROUTINE sub_read


END SUBROUTINE read_stations

SUBROUTINE read_stationsQ

   USE m_common
   USE m_mpi
   USE m_model

   implicit none


      !=============================================================================================================!
      !  Local Variables
      !=============================================================================================================!

         integer(ip)  ::  l
         character(len=lchar)  ::  filename

      !================================================================================================================!
      !  Opening Data File concerning Q Stations
      !================================================================================================================!
#ifdef USE_SW_MONO

      if ( use_Qobs == 1) then

            inquire( file = 'obs/flow_obs.txt'      , exist = file_exist(1) ) !NB : the flow_obs.txt file is obtained by manually converting data from the sum_mass_flux_out/inflow file from the res/post directory. It contains mass flow calculated at the Riemann interface. To be automated.

            if      ( file_exist(1) ) then
                open(10,file='obs/flow_obs.txt',status='old')

              do i=1,SIZE(stationQ)

                read(10,*)
                read(10,*) stationQ( i )%nb_dt, stationQ( i )%ind_bc

                allocate( stationQ( i )%t( stationQ( i )%nb_dt ))
                allocate( stationQ( i )%Q( stationQ( i )%nb_dt ))

                do j=1,stationQ( i )%nb_dt
                    read(10,*) stationQ( i )%t( j ), stationQ( i )%q( j )
                enddo

              enddo
              close(10)
            else

              call Stopping_Program_Sub("use_Qobs = 1 but no observation file can be found in the obs directory. Please provide obs/flow_obs.txt.")

            endif
      endif
#endif

#ifdef USE_HYDRO
      if (use_Qobs_gr4 == 1) then

        do i=1,size(stationQ)

       write(filename,'("obs/GR4obs_",1I1.1,".txt")') i
       inquire( file = filename     , exist = file_exist(1) )

            if ( file_exist(1) ) then

                open(10,file=filename,status='old')

                read(10,*)
                read(10,*)
                read(10,*)
                read(10,*) l, stationQ( i )%weight

                allocate( stationQ( i )%t( l ))
                allocate( stationQ( i )%Q( l ))

                do j=1,l
                    read(10,*) stationQ( i )%t( j ), stationQ( i )%q( j )
                enddo

                close(10)

            else

                call Stopping_Program_Sub("use_Qobs_gr4 = 1 but no observation file can be found in the obs directory. Please provide "//filename)

            endif

        enddo

    endif

#endif

END SUBROUTINE read_stationsQ
!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Sections Output
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE write_sections( dof,mesh )

   USE m_common
   USE m_mpi
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type( msh ), intent(in   )  ::  mesh

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  iobs , cell

   real(rp)  ::  h , q

   !===================================================================================================================!
   !  Writing Sections Records in File
   !===================================================================================================================!

   if ( .not. allocated( section ) ) return

   if ( w_tecplot == 1 ) then

      call sub_write( 'tecplot' ) ; return

   end if

   if ( w_gnuplot == 1 ) then

      call sub_write( 'gnuplot' ) ; return

   end if

   call sub_write( 'bin' )


   CONTAINS


      SUBROUTINE sub_write( file_type )

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         character(len=*), intent(in)  ::  file_type

         !=============================================================================================================!
         !  Local Variables
         !=============================================================================================================!

         character(len=lchar)  ::  file_name(2)

         integer(ip)  ::  pt

         !=============================================================================================================!
         !  Begin Subroutine
         !=============================================================================================================!
		!write(*,*) "in sub_write_file"

            do iobs = 1,size( section )
            if ( .not. test_dt_just_after( section( iobs )%dt ) ) cycle

            write(file_name(1),'(A,"_",I4.4)') 'res/obs/obs_section' , iobs
            write(file_name(2),'(A,"_",I4.4)') 'res/obs/obs_q_h_section' , iobs

            !==========================================================================================================!
            !
            !==========================================================================================================!

            select case( trim(file_type) )

               case( 'tecplot' )

                  file_name(1) = trim(file_name(1))//'.plt'
                  file_name(2) = trim(file_name(2))//'.plt'

               case( 'gnuplot' )

                  file_name(1) = trim(file_name(1))//'.dat'
                  file_name(2) = trim(file_name(2))//'.dat'

               case( 'bin' )

                  file_name(1) = trim(file_name(1))//'.bin'
                  file_name(2) = trim(file_name(2))//'.bin'

               case default

                  file_name(1) = trim(file_name(1))//'.txt'
                  file_name(2) = trim(file_name(2))//'.txt'

            end select

            !==========================================================================================================!
            !
            !==========================================================================================================!

            if ( section( iobs )%pt(1)%cell > 0 .and. all( is_file_open(:) /= file_name(1) ) ) then

               if ( file_type /= 'bin' ) then

                  inquire( file = file_name(1) , exist = file_exist(1) )

                  if ( .not. file_exist(1) ) then

                     open(10,file=file_name(1),status='replace',form='formatted')

                     select case( file_type )

                        case( 'tecplot' )

                           write(10,'(A)') 'VARIABLES = "x" "y" "bathy" "h" "zs" "u" "v"'

                        case default

                           write(10,'(A)') '# x y bathy h zs u v'

                     end select

                  end if

               else

                  inquire( file = file_name(1) , exist = file_exist(1) )

                  if ( .not. file_exist(1) ) open(10,file=file_name(1),status='replace',form='unformatted')

               end if

               close(10)

               file_open_counter = file_open_counter + 1

               is_file_open( file_open_counter ) = file_name(1)

            end if

            !==========================================================================================================!
            !
            !==========================================================================================================!

            if ( section( iobs )%pt(1)%cell > 0 .and. file_type /= 'bin' ) then

               open(10,file=file_name(1),status='old',position='append',form='formatted')

               select case( file_type )

                  case( 'tecplot' )
			
                     write(10,'(A,ES15.8,A,I6,A)') 'ZONE T = "' , tc , '" , I =' , &
                                                   size( section( iobs )%pt ) , ' , DATAPACKING=POINT'

                  case default

                     write(10,'(A,ES15.8,A,I6)') '# time = ', tc , ' nb_points = ' , size( section( iobs )%pt )

               end select

               close(10)

            end if

            call mpi_wait_all

            !==========================================================================================================!
            !
            !==========================================================================================================!

            do pt = 1,size( section( iobs )%pt )

               cell = section( iobs )%pt( pt )%cell

               if ( cell > 0 ) then

                  if ( file_type == 'bin' ) then

                     open(10,file=file_name(1),status='old',position='append',form='unformatted')

                     write(10) section( iobs )%pt( pt )%coord%x , &
                               section( iobs )%pt( pt )%coord%y , &
                               bathy_cell( cell ) , &
                               dof%h( cell ) , &
                               dof%u( cell ) , &
                               dof%v( cell )

                  else

                     open(10,file=file_name(1),status='old',position='append',form='formatted')

                     write(10,'(7ES16.8)') section( iobs )%pt( pt )%coord%x , &
                                           section( iobs )%pt( pt )%coord%y , &
                                           bathy_cell( cell ) , &
                                           dof%h( cell ) , &
                                           bathy_cell( cell ) + dof%h( cell ) , &
                                           dof%u( cell ) , &
                                           dof%v( cell )

                  end if

                  close(10)

               end if

               call mpi_wait_all

            end do

            !==========================================================================================================!
            !
            !==========================================================================================================!

            if ( proc == 0 .and. all( is_file_open(:) /= file_name(2) ) ) then

               if ( file_type /= 'bin' ) then

                  inquire( file = file_name(2) , exist = file_exist(1) )

                  if ( .not. file_exist(1) ) then

                     open(10,file=file_name(2),status='replace',form='formatted')

                     select case( file_type )

                        case( 'tecplot' )

                           write(10,'(A)') 'VARIABLES = "h_max" "Q"'

                        case default

                           write(10,'(A)') '# h_max Q'

                     end select

                  end if

               else

                  inquire( file = file_name(2) , exist = file_exist(1) )

                  if ( .not. file_exist(1) ) open(10,file=file_name(2),status='replace',form='unformatted')

               end if

               close(10)

               file_open_counter = file_open_counter + 1

               is_file_open( file_open_counter ) = file_name(2)

            end if

            !==========================================================================================================!
            !
            !==========================================================================================================!

            call calc_h_q_in_section

            !==========================================================================================================!
            !
            !==========================================================================================================!

            if ( proc == 0 ) then

               if ( file_type == 'bin' ) then

                  open(10,file=file_name(2),status='old',position='append',form='unformatted')

                  write(10) h , q

               else

                  open(10,file=file_name(2),status='old',position='append',form='formatted')

                  write(10,'(7ES16.8)') h , q

               end if

               close(10)

            end if

         end do

      END SUBROUTINE sub_write


      SUBROUTINE calc_h_q_in_section

         !=============================================================================================================!
         !  Local Variables
         !=============================================================================================================!

         integer(ip)  ::  pt , nb_step

         !=============================================================================================================!
         !  Begin Subroutine
         !=============================================================================================================!

         nb_step = size( section( iobs )%pt ) - 1

         h = zero
         q = zero

         do pt = 1,nb_step + 1

            cell = section( iobs )%pt( pt )%cell

            if ( cell > 0 ) then

               h = max( h , dof%h( cell ) )

               if ( pt == 1 .or. &
                    pt == nb_step + 1 ) then

                  q = q + demi * dof%h( cell ) * ( dof%u( cell ) * section( iobs )%normal%x + &
                                                   dof%v( cell ) * section( iobs )%normal%y )

               else

                  q = q +        dof%h( cell ) * ( dof%u( cell ) * section( iobs )%normal%x + &
                                                   dof%v( cell ) * section( iobs )%normal%y )

               end if

            end if

         end do

         q = q * section( iobs )%dx
! write(*,*) "call mpi_max_r( h )", proc
         call mpi_max_r( h )
!          write(*,*) "call mpi_sum_r( q )", proc
         call mpi_sum_r( q )

      END SUBROUTINE

END SUBROUTINE write_sections                                                                                    !>NOADJ
