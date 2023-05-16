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
!> \file output.f90
!! \brief This file includes routines of creation of ouput file.

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Generic routine to write a real in time with dtp time step (input.txt)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


!>  Generic routine to write a real in time with dtp time step (input.txt)
!!  \details NOT CHECKED
!! \details This subroutine write a real var in the file named 'file_name'.
!! this routine should only be used during the main time loop
!! \param[in]    var Real to print in the output file.
!! \param[in]    file_name Name of the file.
SUBROUTINE write_scalar_in_time( var , file_name )

   USE m_common
   USE m_mpi

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   real(rp), intent(in)  ::  var

   character(len=*), intent(in)  ::  file_name

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   character(len=lchar)  ::  var_name

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   var_name  =  trim(file_name)

   if ( w_tecplot == 1 ) then

      call sub_write( 'tecplot' ) ; return

   end if

   if ( w_gnuplot == 1 ) then

      call sub_write( 'gnuplot' ) ; return

   end if

   call sub_write( 'default' )


   CONTAINS

      !>  Generic routine to write a real in time with dtp time step (input.txt)
      !! \details NOT CHECKED
      !! \details This subroutine writes a real var in the file named 'file_name' from the type of the file wanted.
      !! \param[in]    file_type Type of the file wanted ('gnuplot' and 'tecplot' oly avaible).
      SUBROUTINE sub_write( file_type )

         !=============================================================================================================!
         !  Interface Variables
         !=============================================================================================================!

         character(len=*), intent(in)  ::  file_type

         !=============================================================================================================!
         !  Local Variables
         !=============================================================================================================!

         character(len=lchar)  ::  complete_file_name

         !=============================================================================================================!
         !  Begin Subroutine
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

               file_open_counter  =  file_open_counter  +  1

               is_file_open( file_open_counter )  =  complete_file_name

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


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Generic routine to write a real in time with dtp time step (input.txt)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!>  Generic routine to write a real in time with dtp time step (input.txt)
!!  \details NOT CHECKED
SUBROUTINE write_integer_for_all_nt( var , file_name )

   USE m_common
   USE m_mpi

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   integer(ip), intent(in)  ::  var

   character(len=*), intent(in)  ::  file_name

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   character(len=lchar)  ::  var_name

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   var_name  =  trim(file_name)

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
         !  Interface Variables
         !=============================================================================================================!

         character(len=*), intent(in)  ::  file_type

         !=============================================================================================================!
         !  Local Variables
         !=============================================================================================================!

         character(len=lchar)  ::  complete_file_name

         !=============================================================================================================!
         !  Begin Subroutine
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

            file_open_counter  =  file_open_counter  +  1

            is_file_open( file_open_counter )  =  complete_file_name

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


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Generic routine to write n scalars in time with dtp time step (input.txt)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!
! write nbvars values in time
! this routine should also be used inside the main time loop

!>   Generic routine to write n scalars in time with dtp time step (input.txt)
!!
!! \details  NOT CHECKED
!! This subroutine write a real var in the file named 'file_name'.
!! write nbvars values in time
!! this routine should also be used inside the main time loop
!! \param[in]    var Real to print in the output file.
!! \param[in]    var_names Array of name variable to print in the file.
!! \param[in]    file_name Name of the file.
!! \param[in]    nbvars Number of the variable to print.
SUBROUTINE write_pscalar_in_time( var , var_names , file_name , nbvars )

   USE m_common
   USE m_mpi

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   integer(ip), intent(in)  ::  nbvars

   real(rp)        , dimension(nbvars), intent(in)  ::  var
   character(len=*), dimension(nbvars), intent(in)  ::  var_names

   character(len=*), intent(in)  ::  file_name

   !===================================================================================================================!
   !  Begin Subroutine
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
         !  Interface Variables
         !=============================================================================================================!

         character(len=*), intent(in)  ::  file_type

         !=============================================================================================================!
         !  Local Variables
         !=============================================================================================================!

         character(len=lchar)  ::  complete_file_name

         integer(ip)  ::  ivar

         !=============================================================================================================!
         !  Begin Subroutine
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

                              write(10,'(A)'             ) trim( var_names( ivar ) ) // '"'

                           else

                              write(10,'(A)',advance='no') trim( var_names( ivar ) ) // '" "'

                           end if

                        end do

                     case default

                        write(10,'(A)',advance='no') '# time '

                        do ivar = 1,nbvars

                           if ( ivar == nbvars ) then

                              write(10,'(A)'             ) trim( var_names( ivar ) )

                           else

                              write(10,'(A)',advance='no') trim( var_names( ivar ) ) // ' '

                           end if

                        end do

                  end select

                  close(10)

               end if

               file_open_counter  =  file_open_counter  +  1

               is_file_open( file_open_counter )  =  complete_file_name

            end if

            !==========================================================================================================!
            !
            !==========================================================================================================!

            open(10,file=complete_file_name,status='old',position='append',form='formatted')

            write(10,'(ES23.15)',advance='no') tc

            do ivar = 1,nbvars

               if ( ivar == nbvars ) then

                  write(10,'(ES23.15)'             ) var( ivar )

               else

                  write(10,'(ES23.15)',advance='no') var( ivar )

               end if

            end do

            close(10)

         end if

         call mpi_wait_all

      END SUBROUTINE sub_write


END SUBROUTINE write_pscalar_in_time


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Generic routine to write n scalars with dtp time step (input.txt)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!>   Generic routine to write n scalars with dtp time step (input.txt)
!!
!! \details NOT CHECKED 
!! This subroutine write a real var in the file named 'file_name'.
!! write nbvars values in time
!! this routine should also be used inside the main time loop
!! \param[in]    var Real to print in the output file.
!! \param[in]    var_names Array of name variable to print in the file.
!! \param[in]    file_name Name of the file.
!! \param[in]    nbvars Number of the variable to print.
SUBROUTINE write_pscalar( var , var_names , file_name , nbvars )

   USE m_common
   USE m_mpi

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   integer(ip), intent(in)  ::  nbvars

   real(rp)        , dimension(nbvars), intent(in)  ::  var
   character(len=*), dimension(nbvars), intent(in)  ::  var_names

   character(len=*), intent(inout)  ::  file_name

   !===================================================================================================================!
   !  Calling Generic sub_write with Output Options defined in 'input.txt'
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
         !  Interface Variables
         !=============================================================================================================!

         character(len=*), intent(in)  ::  file_type

         !=============================================================================================================!
         !  Local Variables
         !=============================================================================================================!

         character(len=lchar)  ::  complete_file_name

         integer(ip)  ::  ivar

         !=============================================================================================================!
         !  Begin Subroutine
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

                              write(10,'(A)'             ) trim( var_names( ivar ) ) // '"'

                           else

                              write(10,'(A)',advance='no') trim( var_names( ivar ) ) // '" "'

                           end if

                        end do

                     case default

                        write(10,'(A)',advance='no') '# '

                        do ivar = 1,nbvars

                           if ( ivar == nbvars ) then

                              write(10,'(A)'             ) trim( var_names( ivar ) )

                           else

                              write(10,'(A)',advance='no') trim( var_names( ivar ) ) // ' '

                           end if

                        end do

                  end select

                  close(10)

               end if

               file_open_counter  =  file_open_counter  +  1

               is_file_open( file_open_counter )  =  complete_file_name

            end if

            !==========================================================================================================!
            !
            !==========================================================================================================!

            open(10,file=complete_file_name,status='old',position='append',form='formatted')

            do ivar = 1,nbvars

               if ( ivar == nbvars ) then

                  write(10,'(ES23.15)'             ) var( ivar )

               else

                  write(10,'(ES23.15)',advance='no') var( ivar )

               end if

            end do

            close(10)

         end if

         call mpi_wait_all

      END SUBROUTINE sub_write


END SUBROUTINE write_pscalar


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Generic routine to write a scalar field
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!



!>  Generic routine to write a scalar field
!!
!! \details NOT CHECKED 
!! This subroutine write a real var in the file named 'file_name'.
!! write nbvars values in time
!! this routine should also be used inside the main time loop
!! \param[in]    var Real to print in the output file.
!! \param[in]    mesh msh mesh data
!! \param[in]    file_name Name of the file.
SUBROUTINE write_scalar_field( var , mesh , file_name )

   USE m_common
   USE m_mpi

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type(msh), intent(in)  ::  mesh

   real(rp), dimension( mesh%nc + mesh%ncb ), intent(in)  ::  var

   character(len=*), intent(in)  ::  file_name

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   character(len=lchar)  ::  complete_file_name

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   if ( w_tecplot == 1 ) then

      write(complete_file_name,'(A,".plt")') trim(file_name)

      call sub_write_tecplot

   end if

!   if ( w_vtk == 1 ) then

!      write(complete_file_name,'(A,".vtk")') trim(file_name)

!      call sub_write_vtk

!   end if


   CONTAINS


      SUBROUTINE sub_write_tecplot

         !=============================================================================================================!
         !  Local Variables
         !=============================================================================================================!

            integer(ip)  ::  mesh_total_size

         !=============================================================================================================!
         !  Begin Subroutine
         !=============================================================================================================!

!             #ifdef USE_MPI
!             write(*,*) "MPI_ALLREDUCE in sub_write_tecplot + write_scalar_field"
!                call MPI_ALLREDUCE( mesh%nc , mesh_total_size , 1 , inttype , MPI_SUM , MPI_COMM_WORLD , code )
!             #else
!                mesh_total_size = mesh%nc
!             #endif

         !=============================================================================================================!
         !  Opening/Creating result File and Writing Header
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
            !  Writing Tecplot file x and y nodes coordinates
            !==========================================================================================================!

            write(10,'(ES15.8)') mesh%node(:)%coord%x
            write(10,'(ES15.8)') mesh%node(:)%coord%y

            close(10)

         end if

         !=============================================================================================================!
         !  Writing Tecplot file var cell data
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
         !   Writing Tecplot file nodes connectivity
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

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Generic routine to write a scalar field
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!>  Generic routine to write a scalar field
!!
!! \details NOT CHECKED 
!! This subroutine write a real var in the file named 'file_name'.
!! write nbvars values in time
!! this routine should also be used inside the main time loop
!! \param[in]    var Real to print in the output file.
!! \param[in]    mesh msh mesh data
!! \param[in]    file_name Name of the file.
SUBROUTINE write_vector_field( var , mesh , file_name )

   USE m_common
   USE m_mpi

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type(msh), intent(in)  ::  mesh

   type(vec2d), dimension( mesh%nc + mesh%ncb ), intent(in)  ::  var

   character(len=*), intent(in)  ::  file_name

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   character(len=lchar)  ::  complete_file_name

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   if ( w_tecplot == 1 ) then

      write(complete_file_name,'(A,".plt")') trim(file_name)

      call sub_write_tecplot

   end if

!   if ( w_vtk == 1 ) then

!      write(complete_file_name,'(A,".vtk")') trim(file_name)

!      call sub_write_vtk

!   end if


   CONTAINS


      SUBROUTINE sub_write_tecplot

         !=============================================================================================================!
         !  Local Variables
         !=============================================================================================================!

            integer(ip)  ::  mesh_total_size

         !=============================================================================================================!
         !  Begin Subroutine
         !=============================================================================================================!

            #ifdef USE_MPI
            write(*,*) "MPI_ALLREDUCE in sub_write_tecplot + write_vector_field"
               call MPI_ALLREDUCE( mesh%nc , mesh_total_size , 1 , inttype , MPI_SUM , MPI_COMM_WORLD , code )
            #else
               mesh_total_size = mesh%nc
            #endif

         !=============================================================================================================!
         !  Opening/Creating result File and Writing Header
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
            !  Writing Tecplot file x and y nodes coordinates
            !==========================================================================================================!

            write(10,'(ES15.8)') mesh%node(:)%coord%x
            write(10,'(ES15.8)') mesh%node(:)%coord%y

            close(10)

         end if

         !=============================================================================================================!
         !  Writing Tecplot file var cell data
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
         !   Writing Tecplot file nodes connectivity
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
