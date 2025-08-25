SUBROUTINE write_scalar_in_time( var , file_name )
   USE m_common
   USE m_mpi
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   real(rp), intent(in) :: var
   character(len=*), intent(in) :: file_name
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   character(len=lchar) :: var_name
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   var_name = trim(file_name)
   if ( w_tecplot == 1 ) then
      call sub_write( 'tecplot' ) ; return
   end if
   if ( w_gnuplot == 1 ) then
      call sub_write( 'gnuplot' ) ; return
   end if
   call sub_write( 'default' )
   CONTAINS
      !> Generic routine to write a real in time with dtp time step (input.txt)
      !! \details NOT CHECKED
      !! \details This subroutine writes a real var in the file named 'file_name' from the type of the file wanted.
      !! \param[in] file_type Type of the file wanted ('gnuplot' and 'tecplot' oly avaible).
      SUBROUTINE sub_write( file_type )
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         character(len=*), intent(in) :: file_type
         !=============================================================================================================!
         ! Local Variables
         !=============================================================================================================!
         character(len=lchar) :: complete_file_name
         !=============================================================================================================!
         ! Begin Subroutine
         !=============================================================================================================!
         if ( proc == 0 .and. test_dt_nearest( dtp ) ) then
            !==========================================================================================================!
            !
            !==========================================================================================================!
            select case( file_type )
               case( 'tecplot' )
                  complete_file_name = 'res/post/'//trim(file_name)//'.plt'
               case( 'gnuplot' )
                  complete_file_name = 'res/post/'//trim(file_name)//'.dat'
               case default
                  complete_file_name = 'res/post/'//trim(file_name)//'.txt'
            end select
            !==========================================================================================================!
            !
            !==========================================================================================================!
            ! create the target file if not existed
            if ( all( is_file_open(:) /= complete_file_name ) ) then
               inquire( file = complete_file_name , exist = file_exist(1) )
               if ( .not. file_exist(1) ) then
               call system('mkdir -p res/post')
                  open(10,file=complete_file_name,status='replace',form='formatted')
                  select case( file_type )
                     case( 'tecplot' )
                        write(10,'(A,A,A)') 'VARIABLES = "time" "' , trim(var_name) , '"'
                     case default
                        write(10,'(A,A)') '# time ' , trim(var_name)
                  end select
                  close(10)
               end if
               file_open_counter = file_open_counter + 1
               is_file_open( file_open_counter ) = complete_file_name
            end if
            !==========================================================================================================!
            !
            !==========================================================================================================!
            ! write the value in target file
            open(10,file=complete_file_name,status='old',position='append',form='formatted')
            write(10,'(2ES23.15)') tc , var
            close(10)
         end if
      END SUBROUTINE sub_write
END SUBROUTINE write_scalar_in_time
SUBROUTINE write_integer_for_all_nt( var , file_name )
   USE m_common
   USE m_mpi
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   integer(ip), intent(in) :: var
   character(len=*), intent(in) :: file_name
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   character(len=lchar) :: var_name
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   var_name = trim(file_name)
   if ( w_tecplot == 1 ) then
      call sub_write( 'tecplot' ) ; return
   end if
   if ( w_gnuplot == 1 ) then
      call sub_write( 'gnuplot' ) ; return
   end if
   call sub_write( 'default' )
   CONTAINS
      SUBROUTINE sub_write( file_type )
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         character(len=*), intent(in) :: file_type
         !=============================================================================================================!
         ! Local Variables
         !=============================================================================================================!
         character(len=lchar) :: complete_file_name
         !=============================================================================================================!
         ! Begin Subroutine
         !=============================================================================================================!
         select case( file_type )
            case( 'tecplot' )
               complete_file_name = 'res/post/'//trim(file_name)//'.plt'
            case( 'gnuplot' )
               complete_file_name = 'res/post/'//trim(file_name)//'.dat'
            case default
               complete_file_name = 'res/post/'//trim(file_name)//'.txt'
         end select
         !=============================================================================================================!
         !
         !=============================================================================================================!
   ! create the file if doesnot exist
         if ( all( is_file_open(:) /= complete_file_name ) ) then
            inquire( file = complete_file_name , exist = file_exist(1) )
            if ( .not. file_exist(1) ) then
    call system('mkdir -p res/post')
               open(10,file=complete_file_name,status='replace',form='formatted')
               select case( file_type )
                  case( 'tecplot' )
                     write(10,'(A,A,A)') 'VARIABLES = "time" "' , trim(var_name) , '"'
                  case default
                     write(10,'(A,A)') '# time ' , trim(var_name)
               end select
               close(10)
            end if
            file_open_counter = file_open_counter + 1
            is_file_open( file_open_counter ) = complete_file_name
         end if
         !=============================================================================================================!
         !
         !=============================================================================================================!
         open(10,file=complete_file_name,status='old',position='append',form='formatted')
         write(10,'(2I12)') nt , var
         close(10)
         !=============================================================================================================!
         !
         !=============================================================================================================!
         call mpi_wait_all
      END SUBROUTINE sub_write
END SUBROUTINE write_integer_for_all_nt
SUBROUTINE write_pscalar_in_time( var , var_names , file_name , nbvars )
   USE m_common
   USE m_mpi
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   integer(ip), intent(in) :: nbvars
   real(rp) , dimension(nbvars), intent(in) :: var
   character(len=*), dimension(nbvars), intent(in) :: var_names
   character(len=*), intent(in) :: file_name
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   if ( w_tecplot == 1 ) then
      call sub_write( 'tecplot' ) ; return
   end if
   if ( w_gnuplot == 1 ) then
      call sub_write( 'gnuplot' ) ; return
   end if
   call sub_write( 'default' )
   CONTAINS
      SUBROUTINE sub_write( file_type )
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         character(len=*), intent(in) :: file_type
         !=============================================================================================================!
         ! Local Variables
         !=============================================================================================================!
         character(len=lchar) :: complete_file_name
         integer(ip) :: ivar
         !=============================================================================================================!
         ! Begin Subroutine
         !=============================================================================================================!
         if ( proc == 0 .and. test_dt_nearest( dtp ) ) then
            !==========================================================================================================!
            !
            !==========================================================================================================!
            select case( file_type )
               case( 'tecplot' )
                  complete_file_name = 'res/post/'//trim(file_name)//'.plt'
               case( 'gnuplot' )
                  complete_file_name = 'res/post/'//trim(file_name)//'.dat'
               case default
                  complete_file_name = 'res/post/'//trim(file_name)//'.txt'
            end select
            !==========================================================================================================!
            !
            !==========================================================================================================!
            if ( all( is_file_open(:) /= complete_file_name ) ) then
               inquire( file = complete_file_name , exist = file_exist(1) )
               if ( .not. file_exist(1) ) then
               call system('mkdir -p res/post')
                  open(10,file=complete_file_name,status='replace',form='formatted')
                  select case( file_type )
                     case( 'tecplot' )
                        write(10,'(A)',advance='no') 'VARIABLES = "time" "'
                        do ivar = 1,nbvars
                           if ( ivar == nbvars ) then
                              write(10,'(A)' ) trim( var_names( ivar ) ) // '"'
                           else
                              write(10,'(A)',advance='no') trim( var_names( ivar ) ) // '" "'
                           end if
                        end do
                     case default
                        write(10,'(A)',advance='no') '# time '
                        do ivar = 1,nbvars
                           if ( ivar == nbvars ) then
                              write(10,'(A)' ) trim( var_names( ivar ) )
                           else
                              write(10,'(A)',advance='no') trim( var_names( ivar ) ) // ' '
                           end if
                        end do
                  end select
                  close(10)
               end if
               file_open_counter = file_open_counter + 1
               is_file_open( file_open_counter ) = complete_file_name
            end if
            !==========================================================================================================!
            !
            !==========================================================================================================!
            open(10,file=complete_file_name,status='old',position='append',form='formatted')
            write(10,'(ES23.15)',advance='no') tc
            do ivar = 1,nbvars
               if ( ivar == nbvars ) then
                  write(10,'(ES23.15)' ) var( ivar )
               else
                  write(10,'(ES23.15)',advance='no') var( ivar )
               end if
            end do
            close(10)
         end if
         call mpi_wait_all
      END SUBROUTINE sub_write
END SUBROUTINE write_pscalar_in_time
SUBROUTINE write_pscalar( var , var_names , file_name , nbvars )
   USE m_common
   USE m_mpi
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   integer(ip), intent(in) :: nbvars
   real(rp) , dimension(nbvars), intent(in) :: var
   character(len=*), dimension(nbvars), intent(in) :: var_names
   character(len=*), intent(inout) :: file_name
   !===================================================================================================================!
   ! Calling Generic sub_write with Output Options defined in 'input.txt'
   !===================================================================================================================!
   if ( w_tecplot == 1 ) then
      call sub_write( 'tecplot' ) ; return
   end if
   if ( w_gnuplot == 1 ) then
      call sub_write( 'gnuplot' ) ; return
   end if
   call sub_write( 'default' )
   CONTAINS
      SUBROUTINE sub_write( file_type )
         !=============================================================================================================!
         ! Interface Variables
         !=============================================================================================================!
         character(len=*), intent(in) :: file_type
         !=============================================================================================================!
         ! Local Variables
         !=============================================================================================================!
         character(len=lchar) :: complete_file_name
         integer(ip) :: ivar
         !=============================================================================================================!
         ! Begin Subroutine
         !=============================================================================================================!
         if ( proc == 0 .and. test_dt_nearest( dtp ) ) then
            !==========================================================================================================!
            !
            !==========================================================================================================!
            select case( file_type )
               case( 'tecplot' )
                  complete_file_name = 'res/post/'//trim(file_name)//'.plt'
               case( 'gnuplot' )
                  complete_file_name = 'res/post/'//trim(file_name)//'.dat'
               case default
                  complete_file_name = 'res/post/'//trim(file_name)//'.txt'
            end select
            !==========================================================================================================!
            !
            !==========================================================================================================!
            ! create the file if doesnot exist
            if ( all( is_file_open(:) /= complete_file_name ) ) then
               inquire( file = complete_file_name , exist = file_exist(1) )
               if ( .not. file_exist(1) ) then
    call system('mkdir -p res/post')
                  open(10,file=complete_file_name,status='replace',form='formatted')
                  select case( file_type )
                     case( 'tecplot' )
                        write(10,'(A)',advance='no') 'VARIABLES = "'
                        do ivar = 1,nbvars
                           if ( ivar == nbvars ) then
                              write(10,'(A)' ) trim( var_names( ivar ) ) // '"'
                           else
                              write(10,'(A)',advance='no') trim( var_names( ivar ) ) // '" "'
                           end if
                        end do
                     case default
                        write(10,'(A)',advance='no') '# '
                        do ivar = 1,nbvars
                           if ( ivar == nbvars ) then
                              write(10,'(A)' ) trim( var_names( ivar ) )
                           else
                              write(10,'(A)',advance='no') trim( var_names( ivar ) ) // ' '
                           end if
                        end do
                  end select
                  close(10)
               end if
               file_open_counter = file_open_counter + 1
               is_file_open( file_open_counter ) = complete_file_name
            end if
            !==========================================================================================================!
            !
            !==========================================================================================================!
            open(10,file=complete_file_name,status='old',position='append',form='formatted')
            do ivar = 1,nbvars
               if ( ivar == nbvars ) then
                  write(10,'(ES23.15)' ) var( ivar )
               else
                  write(10,'(ES23.15)',advance='no') var( ivar )
               end if
            end do
            close(10)
         end if
         call mpi_wait_all
      END SUBROUTINE sub_write
END SUBROUTINE write_pscalar
SUBROUTINE write_scalar_field( var , mesh , file_name )
   USE m_common
   USE m_mpi
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type(msh), intent(in) :: mesh
   real(rp), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
   character(len=*), intent(in) :: file_name
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   character(len=lchar) :: complete_file_name
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   if ( w_tecplot == 1 ) then
      write(complete_file_name,'(A,".plt")') trim(file_name)
      call sub_write_tecplot
   end if
   CONTAINS
      SUBROUTINE sub_write_tecplot
         !=============================================================================================================!
         ! Local Variables
         !=============================================================================================================!
            integer(ip) :: mesh_total_size
         !=============================================================================================================!
         ! Begin Subroutine
         !=============================================================================================================!
               mesh_total_size = mesh%nc
         !=============================================================================================================!
         ! Opening/Creating result File and Writing Header
         !=============================================================================================================!
         if ( proc == 0 ) then
            if ( all( is_file_open(:) /= complete_file_name ) ) then
               open(10,file=complete_file_name,status='replace',form='formatted')
               write(10,'(A)') 'TITLE = "DassFlow result File"'
               write(10,'(A)') 'VARIABLES = "x","y","var"'
               close(10)
               file_open_counter = file_open_counter + 1
               is_file_open( file_open_counter ) = complete_file_name
            end if
            open(10,file=complete_file_name,status='old',position='append',form='formatted')
            write(10,'(A,ES15.8,A)') 'ZONE T = "' , tc , '"'
            write(10,'(A        )') 'DATAPACKING = BLOCK'
            write(10,'(A,I10    )') 'N = ' , mesh%nn
            write(10,'(A,I10    )') 'E = ' , mesh_total_size
            write(10,'(A        )') 'ZONETYPE = FEQUADRILATERAL'
            write(10,'(A        )') 'VARLOCATION = ([3-8]=CELLCENTERED)'
            !==========================================================================================================!
            ! Writing Tecplot file x and y nodes coordinates
            !==========================================================================================================!
            write(10,'(ES15.8)') mesh%node(:)%coord%x
            write(10,'(ES15.8)') mesh%node(:)%coord%y
            close(10)
         end if
         !=============================================================================================================!
         ! Writing Tecplot file var cell data
         !=============================================================================================================!
         do k = 0,np-1
            if ( proc == k ) then
               open(10,file=complete_file_name,status='old',position='append',form='formatted')
               write(10,'(ES15.8)') var(1:mesh%nc)
               close(10)
            end if
            call mpi_wait_all
         end do
         !=============================================================================================================!
         ! Writing Tecplot file nodes connectivity
         !=============================================================================================================!
         do k = 0,np-1
            if ( proc == k ) then
               open(10,file=complete_file_name,status='old',position='append',form='formatted')
               do i = 1,mesh%nc
                  write(10,'(I10)') mesh%cell(i)%node(:)
               end do
               close(10)
            end if
            call mpi_wait_all
         end do
      END SUBROUTINE sub_write_tecplot
END SUBROUTINE write_scalar_field
SUBROUTINE write_vector_field( var , mesh , file_name )
   USE m_common
   USE m_mpi
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type(msh), intent(in) :: mesh
   type(vec2d), dimension( mesh%nc + mesh%ncb ), intent(in) :: var
   character(len=*), intent(in) :: file_name
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   character(len=lchar) :: complete_file_name
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   if ( w_tecplot == 1 ) then
      write(complete_file_name,'(A,".plt")') trim(file_name)
      call sub_write_tecplot
   end if
   CONTAINS
      SUBROUTINE sub_write_tecplot
         !=============================================================================================================!
         ! Local Variables
         !=============================================================================================================!
            integer(ip) :: mesh_total_size
         !=============================================================================================================!
         ! Begin Subroutine
         !=============================================================================================================!
               mesh_total_size = mesh%nc
         !=============================================================================================================!
         ! Opening/Creating result File and Writing Header
         !=============================================================================================================!
         if ( proc == 0 ) then
            if ( all( is_file_open(:) /= complete_file_name ) ) then
               open(10,file=complete_file_name,status='replace',form='formatted')
               write(10,'(A)') 'TITLE = "DassFlow result File"'
               write(10,'(A)') 'VARIABLES = "x","y","vx","vy"'
               close(10)
               file_open_counter = file_open_counter + 1
               is_file_open( file_open_counter ) = complete_file_name
            end if
            open(10,file=complete_file_name,status='old',position='append',form='formatted')
            write(10,'(A,ES15.8,A)') 'ZONE T = "' , tc , '"'
            write(10,'(A        )') 'DATAPACKING = BLOCK'
            write(10,'(A,I10    )') 'N = ' , mesh%nn
            write(10,'(A,I10    )') 'E = ' , mesh_total_size
            write(10,'(A        )') 'ZONETYPE = FEQUADRILATERAL'
            write(10,'(A        )') 'VARLOCATION = ([3-8]=CELLCENTERED)'
            !==========================================================================================================!
            ! Writing Tecplot file x and y nodes coordinates
            !==========================================================================================================!
            write(10,'(ES15.8)') mesh%node(:)%coord%x
            write(10,'(ES15.8)') mesh%node(:)%coord%y
            close(10)
         end if
         !=============================================================================================================!
         ! Writing Tecplot file var cell data
         !=============================================================================================================!
         do k = 0,np-1
            if ( proc == k ) then
               open(10,file=complete_file_name,status='old',position='append',form='formatted')
               write(10,'(ES15.8)') var(1:mesh%nc)%x
               close(10)
            end if
            call mpi_wait_all
         end do
         do k = 0,np-1
            if ( proc == k ) then
               open(10,file=complete_file_name,status='old',position='append',form='formatted')
               write(10,'(ES15.8)') var(1:mesh%nc)%y
               close(10)
            end if
            call mpi_wait_all
         end do
         !=============================================================================================================!
         ! Writing Tecplot file nodes connectivity
         !=============================================================================================================!
         do k = 0,np-1
            if ( proc == k ) then
               open(10,file=complete_file_name,status='old',position='append',form='formatted')
               do i = 1,mesh%nc
                  write(10,'(I10)') mesh%cell(i)%node(:)
               end do
               close(10)
            end if
            call mpi_wait_all
         end do
      END SUBROUTINE sub_write_tecplot
END SUBROUTINE write_vector_field
