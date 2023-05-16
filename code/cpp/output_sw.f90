SUBROUTINE write_results( dof , mesh )
   USE m_common
   USE m_time_screen
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   TYPE( unk ), intent(in) :: dof
   TYPE( msh ), intent(in) :: mesh
   !===================================================================================================================!
   ! Testing if it is time to write a result file
   !===================================================================================================================!
   if ( .not. test_dt_just_after( dtw ) .and. &
              tc > zerom .and. &
              abs( tc - ts ) > zerom .and. &
              nt /= max_nt_for_direct ) return
   !===================================================================================================================!
   ! Output Screen
   !===================================================================================================================!
   call Print_Screen( 'result' )
   !===================================================================================================================!
   ! Writing Result File
   !===================================================================================================================!
   call system('mkdir -p res')
   call write_result_file( dof , mesh , 'res/result' )
   !call write_restart_with_tc0( dof , mesh )
   !===================================================================================================================!
   ! Writing Exact Solution
   !===================================================================================================================!
END SUBROUTINE write_results
SUBROUTINE write_result_file( dof , mesh , namefile )
   USE m_common
   USE m_time_screen
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   TYPE( unk ), intent(in) :: dof
   TYPE( msh ), intent(in) :: mesh
   character(len=*), intent(in) :: namefile
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   character(50) :: filename
   !===================================================================================================================!
   ! Creating result file name
   !===================================================================================================================!
   if ( tc < zerom ) then
      write(filename,'(A,"_initial")') namefile
   else if ( abs( tc - ts ) < zerom ) then
      write(filename,'(A,"_final")') namefile
   else
      write(filename,'(A,"_",ES12.6)') namefile , tc
   end if
   !===================================================================================================================!
   ! Gnuplot result file output
   !===================================================================================================================!
   if ( w_gnuplot == 1 ) then
      call v_gnuplot( dof , mesh , trim(filename)//'.dat' )
    ! call v_gnuplot_with_ghostcells( dof , mesh , trim(filename)//'_with_ghostcells.dat' )
   end if
   !===================================================================================================================!
   ! Infiltration write , TODO make its own routine
   !===================================================================================================================!
   if (bc_infil/=0) then
    open(10,file='infil_temp.dat',status='replace',form='formatted')
    write(10,*) '# id_cell id_infil infil_qty'
    do i=1,mesh%nc
                write(10,'(2i4, ES15.8)') i, infil%land(i), dof%infil(i)
    end do
        close(10)
   endif
END SUBROUTINE write_result_file
SUBROUTINE write_restart( dof , mesh )
   USE m_common
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( unk ), intent(in) :: dof
   type( msh ), intent(in) :: mesh
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   character(len=lchar) :: file_name
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   call system('mkdir -p restart')
   write(file_name,'(A,I4.4,A)') 'restart/restart' , proc , '.bin'
   open(10,file=trim(file_name),status='replace',form='unformatted')
   write(10) tc
   write(10) nt
   write(10) dof%h
   write(10) dof%u
   write(10) dof%v
   close(10)
END SUBROUTINE write_restart
SUBROUTINE write_restart_direct( dof , mesh )
   USE m_common
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( unk ), intent(in) :: dof
   type( msh ), intent(in) :: mesh
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   if ( proc == 0 ) then
      open(10,file='restart.bin',status='replace',form='unformatted',access='direct',recl=3*length_real)
      write(10,*) tc
      close(10)
   end if
if ( proc == 0) then
   open(10,file='restart.bin',status='old',form='unformatted',access='direct',recl=3*length_real)
   do i = 1,mesh%nc
      write(10,*) dof%h(i) , dof%u(i) , dof%v(i)
   end do
   close(10)
end if
END SUBROUTINE write_restart_direct
SUBROUTINE v_gnuplot( dof , mesh , filename )
   USE m_common
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   TYPE(unk), intent(inout) :: dof
   TYPE(msh), intent(in) :: mesh
   character(len=*), intent(in) :: filename
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   integer(ip) :: rec_index, mesh_total_size
   real(rp) :: h , u_temp , v_temp , c , hx , hy , dtx , dty , dmin
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   mesh_total_size = mesh%nc
   !===================================================================================================================!
   ! Opening/Creating Result File and Writing Header
   !===================================================================================================================!
   if ( proc == 0 ) then
      if ( all( is_file_open(:) /= filename ) ) then
         open(10,file=filename,status='replace',form='formatted')
         write(10,*) '# Gnuplot DataFile Version'
         write(10,*) '# i x y bathy h zs Manning u v'
         close(10)
         file_open_counter = file_open_counter + 1
         is_file_open( file_open_counter ) = filename
      end if
   end if
      if ( proc == 0 ) then
         open(10,file=filename,status='old',position='append',form='formatted')
         do i=1,mesh%nc
if ( abs( dof%u( i ) ) < 1E-16 ) then
         u_temp = 0._rp ! dirty force to zero values below 10**-99
         write(*,*) "u_temp= dof%u( i )-modif", u_temp
else
         u_temp = dof%u(i)
end if
if (abs(dof%v(i)) < 1E-16 ) then
         v_temp =0._rp ! dirty force to zero values below 10**-100
         write(*,*) "v= dof%v( i )", v_temp
else
        v_temp = dof%v(i)
end if
            write(10,'(I5,8(" ",ES15.8))') i , &
        mesh%cell(i)%grav%x , &
                                 mesh%cell(i)%grav%y , &
                                 bathy_cell(i) , &
                                 dof%h(i) , &
                                 bathy_cell(i)+dof%h(i) , &
                                 manning( land(i) ) , &
                                 u_temp , &
                                 v_temp
         end do
         close(10)
      end if
END SUBROUTINE v_gnuplot
