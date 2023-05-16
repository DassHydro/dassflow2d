SUBROUTINE calc_innovation( dof,mesh )
   USE m_common
   USE m_model
   USE m_obs
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( unk ), intent(in) :: dof
   type( msh ), intent(in) :: mesh
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   integer(ip) :: cell , searched_time , pt, N_average
   real(rp) :: h_mean,s_total
   !===================================================================================================================!
   ! Begin
   !===================================================================================================================!
    !write(*,*) "size( station )",size( station )
   do iobs = 1,size( station )
      searched_time = innovation( iobs )%ind_t
      if ( searched_time > innovation( iobs )%nb_dt ) cycle
      if ( tc >= station( iobs )%t( searched_time ) ) then
    !write(*,*) " tc >= station( iobs )%t( searched_time ) "
         h_mean = 0._rp
         s_total= 0._rp
         do pt = 1,size( station( iobs )%pt )
            cell = station( iobs )%pt( pt )%cell
            if ( cell < 0 ) cycle
            !*****************************
            ! H = 1/swet \int_swet Hdx
            !*****************************
            if ( dof%h( cell ) > 0 ) then !test on water presence determining if cell is used for calculating h_average
               h_mean = h_mean + (dof%h( cell )+bathy_cell( cell ) )*mesh%cell(cell)%surf
               s_total= s_total + mesh%cell(cell)%surf
            end if
            !*****************************
            ! H = 1/sobs \int_sobs Hdx
            !*****************************
            !*****************************
            ! H = 1/sriver \int_sobs Hdx
            !*****************************
         end do
         !*****************************
         ! H = 1/sriver \int_sobs Hdx
         !*****************************
         !s_total=100._rp
         if (s_total>0) then
            h_mean=h_mean/ s_total
         end if
         innovation ( iobs )%diff( searched_time ) = (h_mean - station( iobs )%h( searched_time ))
         innovation ( iobs )%ind_t = innovation ( iobs )%ind_t + 1
      end if
   end do
END SUBROUTINE calc_innovation
                                                                                                                 !<NOADJ
SUBROUTINE write_stations( dof ,mesh )
   USE m_common
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( unk ), intent(in) :: dof
   type( msh ), intent(in ) :: mesh
   !===================================================================================================================!
   ! Writing Stations Records in File
   !===================================================================================================================!
   if ( .not. allocated( station ) .or. w_obs /= 1 ) then
    return
   end if
   if ( w_tecplot == 1 ) then
      call sub_write( 'tecplot' ) !; return
   end if
   if ( w_gnuplot == 1 ) then
      call sub_write( 'gnuplot' )
   end if
   CONTAINS
      SUBROUTINE sub_write( file_type )
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         character(len=*), intent(in) :: file_type
         !=============================================================================================================!
         ! Local Variables
         !=============================================================================================================!
         character(len=lchar) :: file_name
         integer(ip) :: iobs , cell , pt , N_average,searched_time
         real(rp) :: h_mean , u_mean, v_mean,s_total , w_meanb
         !=============================================================================================================!
         ! Begin Subroutine
         !=============================================================================================================!
         do iobs = 1,size( station )
            !==========================================================================================================!
            ! Testing if Simulation Time match with Observation Time Step
            !==========================================================================================================!
            !if ( .not. test_dt_just_after( station( iobs )%dt ) ) cycle
            searched_time = station( iobs )%ind_t
            if ( searched_time > station( iobs )%nb_dt ) cycle
            if ( tc >= station( iobs )%dt_obs( searched_time ) ) then
               !==========================================================================================================!
               ! Creating File Name
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
                  if ( .not. file_exist(1) ) then
                     if ( file_type /= 'bin' ) then
                        open(10,file=file_name,status='replace',form='formatted')
                        select case( file_type )
                           case( 'tecplot' )
                              write(10,'(A)') 'VARIABLES = "time" "h_mean" "u_mean" "v_mean" "w_mean"'
                           case default
                              write(10,'(A)') '# time h_mean u_mean v_mean w_mean'
                        end select
                     else
                        open(10,file=file_name,status='replace',form='unformatted')
                     end if
                  end if
                  close(10)
                  file_open_counter = file_open_counter + 1
                  is_file_open( file_open_counter ) = file_name
               end if
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
                  if ( dof%h( cell ) > 0 ) then !test on water presence determining if cell is used for calculating h_average
                     u_mean = u_mean + dof%u( cell )
                     v_mean = v_mean + dof%v( cell )
                     h_mean = h_mean + (dof%h( cell )+bathy_cell( cell ) )*mesh%cell( cell)%surf
                     s_total = s_total+mesh%cell( cell)%surf
                  end if
                  !*****************************
                  ! H = 1/sobs \int_sobs Hdx
                  !*****************************
                  !s_total= s_total + mesh%cell(cell)%surf
                  !h_mean = h_mean + dof%h( cell )*mesh%cell(cell)%surf
                  !*****************************
                  ! H = 1/sriver \int_sobs Hdx
                  !*****************************
                  !h_mean = h_mean + dof%h( cell )*mesh%cell(cell)%surf
                  if (dof%h(cell)>0) then
                     w_meanb = w_meanb+mesh%cell(cell)%surf
                     N_average= N_average+1
                  end if
               end do
               !*****************************
               ! H = 1/sriver \int_sobs Hdx
               !*****************************
               !s_total=100
  if (N_average>0) then ! not divide by 0
                 h_mean = h_mean / s_total
  end if
               if (station(iobs)%length>0) then
                  w_meanb=w_meanb/(station(iobs)%length)
               else
                  w_meanb=0.0
               endif
               !write(*,*) 'h_mean', h_mean, 'cell' ,station( iobs )%pt(1)%cell
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
                  ! cell = station( iobs )%pt( pt )%cell
                  ! if (dof%h(cell)>0) then
                  ! write(10,*) cell
                  ! end if
                  !end do
                  !close(10)
                  !ENd
               end if
                station( iobs )%ind_t= station( iobs )%ind_t+1
            end if
       end do
       close(20)
      END SUBROUTINE sub_write
END SUBROUTINE write_stations
SUBROUTINE read_stations
   USE m_common
   USE m_mesh
   USE m_obs
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Reading Stations Records in File
   !===================================================================================================================!
   if ( w_tecplot == 1 ) then
      call sub_read( 'tecplot' ) ; return
   end if
   if ( w_gnuplot == 1 ) then
      call sub_read( 'gnuplot' ) ; return
   end if
   CONTAINS
      SUBROUTINE sub_read( file_type )
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         character(len=*), intent(in) :: file_type
         !=============================================================================================================!
         ! Local Variables
         !=============================================================================================================!
         character(len=lchar) :: file_name
         integer(ip) :: iobs , nbdt_obs , dt_obs , icod
         !=============================================================================================================!
         ! Begin Subroutine
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
          ! if (iobs <= size( station ) ) call Stopping_Program_Sub('Not enough number of observation files or not provided, check obs directory')
               exit
            else if ( iobs > size( station ) ) then
          ! call Stopping_Program_Sub( 'Mismatch between the number of files in obs directory and obs.txt file' )
            end if
          ! if ( file_type /= 'bin' ) then
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
          ! end if
         end do
      END SUBROUTINE sub_read
END SUBROUTINE read_stations
SUBROUTINE read_stationsQ
   USE m_common
   USE m_model
   implicit none
      !=============================================================================================================!
      ! Local Variables
      !=============================================================================================================!
         integer(ip) :: l
         character(len=lchar) :: filename
      !================================================================================================================!
      ! Opening Data File concerning Q Stations
      !================================================================================================================!
      write(*,*) "was only implemented for GR4, which is removed because higly bugged"
END SUBROUTINE read_stationsQ
SUBROUTINE write_sections( dof,mesh )
   USE m_common
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( unk ), intent(in) :: dof
   type( msh ), intent(in ) :: mesh
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   integer(ip) :: iobs , cell
   real(rp) :: h , q
   !===================================================================================================================!
   ! Writing Sections Records in File
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
         ! Interface Variables
         !=============================================================================================================!
         character(len=*), intent(in) :: file_type
         !=============================================================================================================!
         ! Local Variables
         !=============================================================================================================!
         character(len=lchar) :: file_name(2)
         integer(ip) :: pt
         !=============================================================================================================!
         ! Begin Subroutine
         !=============================================================================================================!
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
         ! Local Variables
         !=============================================================================================================!
         integer(ip) :: pt , nb_step
         !=============================================================================================================!
         ! Begin Subroutine
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
                  q = q + dof%h( cell ) * ( dof%u( cell ) * section( iobs )%normal%x + &
                                                   dof%v( cell ) * section( iobs )%normal%y )
               end if
            end if
         end do
         q = q * section( iobs )%dx
      END SUBROUTINE
END SUBROUTINE write_sections !>NOADJ
