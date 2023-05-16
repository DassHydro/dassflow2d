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
!  Main Subroutine Managing Output Files
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE write_results( dof , mesh )

   USE m_common
   USE m_mpi
   USE m_time_screen
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   TYPE( unk ), intent(in)  ::  dof
   TYPE( msh ), intent(in)  ::  mesh

   !===================================================================================================================!
   !  Testing if it is time to write a result file
   !===================================================================================================================!
!write(*,*) "Within write results ---"
!#	write(*,*) 'tc', tc
!#	write(*,*)  'dt', dt
   if ( .not. test_dt_just_after( dtw ) .and. &
              tc > zerom                .and. &
              abs( tc - ts ) > zerom    .and. &
              nt /= max_nt_for_direct ) return

   !===================================================================================================================!
   !  Output Screen
   !===================================================================================================================!

   call Print_Screen( 'result' )

   !===================================================================================================================!
   !  Writing Result File
   !===================================================================================================================!

   call system('mkdir -p res')
!write(*,*) "call write results files"
    call write_result_file( dof , mesh , 'res/result' )
!write(*,*) "out results files"

    ! call write_restart_direct( dof , mesh ) >>> replaced with tc0 : TC=0
    call write_restart_with_tc0( dof , mesh )
   !===================================================================================================================!
   !  Writing Exact Solution
   !===================================================================================================================!

   #ifdef USE_VALID

      if      ( w_exact == 1 ) then

         call write_exact_solution( mesh )

      else if ( w_exact == 2 .and. nt == 0 ) then

         call write_exact_solution( mesh )

      end if

   #endif

!write(*,*) "out write_results routine"
END SUBROUTINE write_results


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Main Subroutine to Write a Result File
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE write_result_file( dof , mesh , namefile )

   USE m_common
   USE m_mpi
   USE m_time_screen
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   TYPE( unk ), intent(in)  ::  dof
   TYPE( msh ), intent(in)  ::  mesh

   character(len=*), intent(in)  ::  namefile

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   character(50)  ::  filename

   !===================================================================================================================!
   !  Creating result file name
   !===================================================================================================================!
!write(*,*) "whithin write result file"
   if      ( tc < zerom ) then

      write(filename,'(A,"_initial")') namefile

   else if ( abs( tc - ts ) < zerom ) then

      write(filename,'(A,"_final")') namefile

   else

      write(filename,'(A,"_",ES12.6)') namefile , tc

   end if

   !===================================================================================================================!
   !  VTK result file output
   !===================================================================================================================!

   if      ( w_vtk == 1 ) then

      call v_vtk    ( dof , mesh , trim(filename)//'.vtk' )

   else if ( w_vtk == 2 ) then

      call v_vtk_bin( dof , mesh , trim(filename)//'.vtk' )

   end if
! write(*,*) "--> vtk done"
   !===================================================================================================================!
   !  Tecplot result file output
   !===================================================================================================================!

   if ( w_tecplot == 1 ) then

      call v_tecplot( dof , mesh , trim(filename)//'.plt' )

   end if
! write(*,*) "--> Tecplot done"
   !===================================================================================================================!
   !  Gnuplot result file output
   !===================================================================================================================!

   if ( w_gnuplot == 1 ) then

      call v_gnuplot( dof , mesh , trim(filename)//'.dat' )

    !  call v_gnuplot_with_ghostcells( dof , mesh , trim(filename)//'_with_ghostcells.dat' )
   end if

! write(*,*) "--> Gnuplot done"

   !===================================================================================================================!
   !  Infiltration write , TODO make its own routine
   !===================================================================================================================!
   if (bc_infil/=0) then
    open(10,file='infil_temp.dat',status='replace',form='formatted')

    write(10,*) '# id_cell id_infil infil_qty'

    do i=1,mesh%nc
                write(10,'(2i4, ES15.8)') i, infil%land(i), dof%infil(i)
    end do

        close(10)
   endif
! write(*,*) "--> infil done"


! write(*,*) "OUT Write_result_file"
END SUBROUTINE write_result_file


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Write a restart file
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE write_restart( dof , mesh )

   USE m_common
   USE m_mpi
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type( msh ), intent(in)  ::  mesh

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   character(len=lchar)  ::  file_name

   !===================================================================================================================!
   !  Begin Subroutine
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
   USE m_mpi
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type( msh ), intent(in)  ::  mesh

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   if ( proc == 0 ) then

      open(10,file='restart.bin',status='replace',form='unformatted',access='direct',recl=3*length_real)

      write(10,rec=1) tc

      close(10)

   end if

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file='restart.bin',status='old',form='unformatted',access='direct',recl=3*length_real)

         do i = 1,mesh%nc

            write(10,rec=1+swap_index(i)) dof%h(i) , dof%u(i) , dof%v(i)

         end do

         close(10)

      end if

      call mpi_wait_all

   end do

END SUBROUTINE write_restart_direct



! Same as write_restart_direct but we set tc = 0 instead of tc = ts, it allow
! to generate initial conditions into a restart.bin file.
! this way new simulation can be ran with appropriate initial conditions
! from tc = 0.0
SUBROUTINE write_restart_with_tc0( dof , mesh )

   USE m_common
   USE m_mpi
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type( msh ), intent(in)  ::  mesh

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   if ( proc == 0 ) then

      open(10,file='restart.bin',status='replace',form='unformatted',access='direct',recl=3*length_real)

      write(10,rec=1) 0.0

      close(10)

   end if

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file='restart.bin',status='old',form='unformatted',access='direct',recl=3*length_real)

         do i = 1,mesh%nc

            write(10,rec=1+swap_index(i)) dof%h(i) , dof%u(i) , dof%v(i)

         end do

         close(10)

      end if

      call mpi_wait_all

   end do

END SUBROUTINE write_restart_with_tc0


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Write a Tecplot Output Result File in standard mode
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE v_tecplot( dof , mesh , filename )

   USE m_common
   USE m_mpi
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type( msh ), intent(in)  ::  mesh

   character(len=*), intent(in)  ::  filename

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  mesh_total_size

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   mesh_total_size = mesh%nc

   call mpi_sum_i( mesh_total_size )

   !===================================================================================================================!
   !  Opening/Creating Result File and Writing Header
   !===================================================================================================================!

   if ( proc == 0 ) then

      if ( all( is_file_open(:) /= filename ) ) then

         open(10,file=filename,status='replace',form='formatted')

         write(10,'(A)') 'TITLE = "DassFlow Result File in Time"'
         write(10,'(A)') 'VARIABLES = "x","y","bathy","h","zs","Manning","u","v"'

         close(10)

         file_open_counter = file_open_counter + 1

         is_file_open( file_open_counter ) = filename

      end if

      open(10,file=filename,status='old',position='append',form='formatted')

      write(10,'(A,ES15.8,A)') 'ZONE T = "' , dof%t_display , '"'
      write(10,'(A        )') 'DATAPACKING = BLOCK'
      write(10,'(A,I10    )') 'N = ' , mesh%nn
      write(10,'(A,I10    )') 'E = ' , mesh_total_size
      write(10,'(A        )') 'ZONETYPE = FEQUADRILATERAL'
      write(10,'(A        )') 'VARLOCATION = ([3-8]=CELLCENTERED)'

      !================================================================================================================!
      !  Writing Tecplot file x and y nodes coordinates
      !================================================================================================================!

      write(10,'(ES15.8)') mesh%node(:)%coord%x
      write(10,'(ES15.8)') mesh%node(:)%coord%y

      close(10)

   end if

   !===================================================================================================================!
   !  Writing Tecplot file bathy cell data
   !===================================================================================================================!

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file=filename,status='old',position='append',form='formatted')

         write(10,'(ES15.8)') bathy_cell(1:mesh%nc)

         close(10)

      end if

      call mpi_wait_all

   end do


   !===================================================================================================================!
   !  Writing Tecplot file h cell data
   !===================================================================================================================!

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file=filename,status='old',position='append',form='formatted')

         write(10,'(ES15.8)') dof%h(1:mesh%nc)

         close(10)

      end if

      call mpi_wait_all

   end do

   !===================================================================================================================!
   !  Writing Tecplot file zs cell data
   !===================================================================================================================!

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file=filename,status='old',position='append',form='formatted')

         write(10,'(ES15.8)') bathy_cell(1:mesh%nc) + dof%h(1:mesh%nc)

         close(10)

      end if

      call mpi_wait_all

   end do

   !===================================================================================================================!
   !  Writing Tecplot file Manning cell data
   !===================================================================================================================!

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file=filename,status='old',position='append',form='formatted')

         write(10,'(ES15.8)') manning( land(1:mesh%nc) )

         close(10)

      end if

      call mpi_wait_all

   end do

   !===================================================================================================================!
   !  Writing Tecplot file u cell data
   !===================================================================================================================!

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file=filename,status='old',position='append',form='formatted')

         write(10,'(ES15.8)') dof%u(1:mesh%nc)

         close(10)

      end if

      call mpi_wait_all

   end do

   !===================================================================================================================!
   !  Writing Tecplot file v cell data
   !===================================================================================================================!

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file=filename,status='old',position='append',form='formatted')

         write(10,'(ES15.8)') dof%v(1:mesh%nc)

         close(10)

      end if

      call mpi_wait_all

   end do

   !===================================================================================================================!
   !   Writing Tecplot file nodes connectivity
   !===================================================================================================================!

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file=filename,status='old',position='append',form='formatted')

         do i = 1,mesh%nc

            write(10,'(I10)') mesh%cell(i)%node(:)

         end do

         close(10)

      end if

      call mpi_wait_all

   end do

END SUBROUTINE v_tecplot


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Write a Tecplot Output Result File in direct mode
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE v_tecplot_direct( dof , mesh , filename )

   USE m_common
   USE m_mpi
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   TYPE( unk ), intent(in)  ::  dof
   TYPE( msh ), intent(in)  ::  mesh

   character(len=*), intent(in)  ::  filename

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  rec_index , mesh_total_size

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   mesh_total_size = mesh%nc

   call mpi_sum_i( mesh_total_size )

   open(10,file=filename,status='replace',form='formatted',access='direct',recl=16)

   !===================================================================================================================!
   !  Creating Tecplot file and header
   !===================================================================================================================!

   write(10,rec= 1,fmt='(A16   )') 'TITLE  =  "time '
   write(10,rec= 2,fmt='(ES14.7,A1,A1)') dof%t_display , '"' , char(10)

   write(10,rec= 3,fmt='(A16   )') 'VARIABLES = "x",'
   write(10,rec= 4,fmt='(A16   )') '"y","bathy","h",'
   write(10,rec= 5,fmt='(A16   )') '"zs","Manning","'
   write(10,rec= 6,fmt='(A15,A1)') 'u","v"         ' , char(10)

   write(10,rec= 7,fmt='(A16   )') 'ZONE NODES =    '
   write(10,rec= 8,fmt='(I15,A1)')  mesh%nn , char(10)
   write(10,rec= 9,fmt='(A16   )') 'ELEMENTS =      '
   write(10,rec=10,fmt='(I15,A1)')  mesh_total_size , char(10)
   write(10,rec=11,fmt='(A16   )') 'DATAPACKING = BL'
   write(10,rec=12,fmt='(A15,A1)') 'OCK            ' , char(10)
   write(10,rec=13,fmt='(A16   )') 'ZONETYPE = FEQUA'
   write(10,rec=14,fmt='(A15,A1)') 'DRILATERAL     ' , char(10)
   write(10,rec=15,fmt='(A16   )') 'VARLOCATION=([3-'
   write(10,rec=16,fmt='(A16   )') '8]=CELLCENTERED)'
   write(10,rec=17,fmt='(A15,A1)') '               ' , char(10)

   rec_index = 17

   do i = 1,mesh%nn

      if ( proc == 0 ) write(10,rec=rec_index+i,fmt='(ES15.8,A1)') mesh%node(i)%coord%x , char(10)

   end do

   rec_index = rec_index + mesh%nn

   do i = 1,mesh%nn

      if ( proc == 0 ) write(10,rec=rec_index+i,fmt='(ES15.8,A1)') mesh%node(i)%coord%y , char(10)

   end do

   rec_index = rec_index + mesh%nn

   !===================================================================================================================!
   !  Writing Tecplot file bathy cell data
   !===================================================================================================================!

   do i = 1,mesh%nc

      write(10,rec=rec_index+swap_index(i),fmt='(ES15.8,A1)') bathy_cell(i) , char(10)

   end do

   rec_index = rec_index + mesh_total_size

   !===================================================================================================================!
   !  Writing Tecplot file h cell data
   !===================================================================================================================!

   do i = 1,mesh%nc

      write(10,rec=rec_index+swap_index(i),fmt='(ES15.8,A1)') dof%h(i) , char(10)

   end do

   rec_index = rec_index + mesh_total_size

   !===================================================================================================================!
   !  Writing Tecplot file zs cell data
   !===================================================================================================================!

   do i = 1,mesh%nc

      write(10,rec=rec_index+swap_index(i),fmt='(ES15.8,A1)') dof%h(i) + bathy_cell(i) , char(10)

   end do

   rec_index = rec_index + mesh_total_size

   !===================================================================================================================!
   !  Writing Tecplot file Manning cell data
   !===================================================================================================================!

   do i = 1,mesh%nc

      write(10,rec=rec_index+swap_index(i),fmt='(ES15.8,A1)') manning( land(i) ) , char(10)

   end do

   rec_index = rec_index + mesh_total_size

   !===================================================================================================================!
   !  Writing Tecplot file u cell data
   !===================================================================================================================!

   do i = 1,mesh%nc

      write(10,rec=rec_index+swap_index(i),fmt='(ES15.8,A1)') dof%u(i) , char(10)

   end do

   rec_index = rec_index + mesh_total_size

   !===================================================================================================================!
   !  Writing Tecplot file v cell data
   !===================================================================================================================!

   do i = 1,mesh%nc

      write(10,rec=rec_index+swap_index(i),fmt='(ES15.8,A1)') dof%v(i) , char(10)

   end do

   rec_index = rec_index + mesh_total_size

   !===================================================================================================================!
   !   Writing Tecplot file nodes connectivity
   !===================================================================================================================!

   do i = 1,mesh%nc
      write(10,rec=rec_index+4*swap_index(i)-3,fmt='(I16   )') mesh%cell(i)%node(1)
      write(10,rec=rec_index+4*swap_index(i)-2,fmt='(I16   )') mesh%cell(i)%node(2)
      write(10,rec=rec_index+4*swap_index(i)-1,fmt='(I16   )') mesh%cell(i)%node(3)
      write(10,rec=rec_index+4*swap_index(i)  ,fmt='(I15,A1)') mesh%cell(i)%node(4) , char(10)
   end do

   close(10)

   call mpi_wait_all

END SUBROUTINE v_tecplot_direct


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Write Exact Solution File
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


#ifdef USE_VALID

   SUBROUTINE write_exact_solution( mesh )

      USE m_common
      USE m_mpi
      USE m_time_screen
      USE m_model

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      TYPE( msh ), intent(in)  ::  mesh

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      TYPE( unk )  ::  dof_exact

      !================================================================================================================!
      !  Calculating exact solution
      !================================================================================================================!

      call alloc_dof( dof_exact , mesh )

      do i = 1,mesh%nc

         dof_exact%h(i) = max( 0._rp , zs_exact( mesh%cell(i)%grav%x , mesh%cell(i)%grav%y , tc ) - bathy_cell(i) )

         if ( dof_exact%h(i) <= heps ) then

            dof_exact%u(i) = 0._rp
            dof_exact%v(i) = 0._rp

         else

            dof_exact%u(i) = u_exact( mesh%cell(i)%grav%x , mesh%cell(i)%grav%y , tc )
            dof_exact%v(i) = v_exact( mesh%cell(i)%grav%x , mesh%cell(i)%grav%y , tc )

         end if

      end do

      !================================================================================================================!
      !  Writing exact solution
      !================================================================================================================!

      call write_result_file( dof_exact , mesh , 'result_exact' )

      call dealloc_dof( dof_exact )

   END SUBROUTINE write_exact_solution

#endif


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Write an Gnuplot Output Result File
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE v_gnuplot( dof , mesh , filename )

   USE m_common
   USE m_mpi
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   TYPE(unk), intent(inout)  ::  dof
   TYPE(msh), intent(in)  ::  mesh

   character(len=*), intent(in)  ::  filename

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  rec_index, mesh_total_size

   real(rp)  ::  h , u_temp , v_temp , c , hx , hy , dtx , dty , dmin

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   mesh_total_size = mesh%nc

   call mpi_sum_i( mesh_total_size )

   !===================================================================================================================!
   !  Opening/Creating Result File and Writing Header
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


   do k = 0,np-1

      if ( proc == k ) then

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

            write(10,'(I5,8(" ",ES15.8))') swap_index(i)					, &
								mesh%cell(i)%grav%x    , &
                                 mesh%cell(i)%grav%y    , &
                                 bathy_cell(i)          , &
                                 dof%h(i)               , &
                                 bathy_cell(i)+dof%h(i) , &
                                 manning( land(i) )     , &
                                 u_temp             , &
                                 v_temp

         end do


         close(10)

      end if

      call mpi_wait_all

   end do


END SUBROUTINE v_gnuplot




!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Write an ASCII VTK Output Result File
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE v_vtk( dof , mesh , filename )

   USE m_common
   USE m_mpi
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   TYPE( unk ), intent(in)  ::  dof
   TYPE( msh ), intent(in)  ::  mesh

   character(len=*), intent(in)  ::  filename

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  rec_index

   real(rp)  ::  h , u , v , c , dtx , dty , dmin

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   open(10,file=filename,status='replace',form='formatted',access='direct',recl=16)

   !===================================================================================================================!
   !  Creating VTK file and header
   !===================================================================================================================!

   write(10,rec= 1,fmt='(A16   )') '# vtk DataFile V'
   write(10,rec= 2,fmt='(A15,A1)') 'ersion 3.0     ', char(10)
   write(10,rec= 3,fmt='(A16   )') 'DassFlow Output '
   write(10,rec= 4,fmt='(A15,A1)') 'File           ', char(10)
   write(10,rec= 5,fmt='(A15,A1)') 'ASCII          ', char(10)
   write(10,rec= 6,fmt='(A16   )') 'DATASET UNSTRUCT'
   write(10,rec= 7,fmt='(A15,A1)') 'URED_GRID      ', char(10)
   write(10,rec= 8,fmt='(A16   )') 'POINTS          '
   write(10,rec= 9,fmt='(I16   )')  mesh%nn
   write(10,rec=10,fmt='(A15,A1)') ' double        ', char(10)

   rec_index = 10

   do i = 1,mesh%nn
      write(10,rec=rec_index+3*i-2,fmt='(ES15.8   )') mesh%node(i)%coord%x
      write(10,rec=rec_index+3*i-1,fmt='(ES15.8   )') mesh%node(i)%coord%y
      write(10,rec=rec_index+3*i  ,fmt='(ES15.8,A1)') 0._rp , char(10)
   end do

   rec_index = rec_index + 3 * mesh%nn

   write(10,rec=rec_index+1,fmt='(A16   )') 'CELLS          '
   write(10,rec=rec_index+2,fmt='(I15   )')  mesh%nc
   write(10,rec=rec_index+3,fmt='(I15,A1)')  mesh%nc*5 , char(10)

   rec_index = rec_index + 3

   do i = 1,mesh%nc
      write(10,rec=rec_index+5*swap_index(i)-4,fmt='(I15   )') 4
      write(10,rec=rec_index+5*swap_index(i)-3,fmt='(I15   )') mesh%cell(i)%node(1)-1
      write(10,rec=rec_index+5*swap_index(i)-2,fmt='(I15   )') mesh%cell(i)%node(2)-1
      write(10,rec=rec_index+5*swap_index(i)-1,fmt='(I15   )') mesh%cell(i)%node(3)-1
      write(10,rec=rec_index+5*swap_index(i)  ,fmt='(I15,A1)') mesh%cell(i)%node(4)-1, char(10)
   end do

   rec_index = rec_index + 5 * mesh%nc

   write(10,rec=rec_index+1,fmt='(A16   )') 'CELL_TYPES      '
   write(10,rec=rec_index+2,fmt='(I15,A1)')  mesh%nc , char(10)

   rec_index = rec_index + 2

   do i = 1,mesh%nc
      write(10,rec=rec_index+i,fmt='(I15,A1)') 9 , char(10)
   end do

   rec_index = rec_index + mesh%nc

   write(10,rec=rec_index+1,fmt='(A16   )')  'CELL_DATA       '
   write(10,rec=rec_index+2,fmt='(I15,A1)')  mesh%nc , char(10)

   rec_index = rec_index + 2

   !===================================================================================================================!
   !  Writing VTK file bathy cell data
   !===================================================================================================================!

   write(10,rec=rec_index+1,fmt='(A16   )') 'SCALARS         '
   write(10,rec=rec_index+2,fmt='(A16   )') 'bathy           '
   write(10,rec=rec_index+3,fmt='(A15,A1)') 'double 1       ' , char(10)
   write(10,rec=rec_index+4,fmt='(A16   )') 'LOOKUP_TABLE def'
   write(10,rec=rec_index+5,fmt='(A15,A1)') 'ault           ' , char(10)

   rec_index = rec_index + 5

   do i = 1,mesh%nc
      write(10,rec=rec_index+swap_index(i),fmt='(ES15.8,A1)') bathy_cell(i) , char(10)
   end do

   rec_index = rec_index + mesh%nc

   !===================================================================================================================!
   !  Writing VTK file h cell data
   !===================================================================================================================!

   write(10,rec=rec_index+1,fmt='(A16   )') 'SCALARS         '
   write(10,rec=rec_index+2,fmt='(A16   )') 'h               '
   write(10,rec=rec_index+3,fmt='(A15,A1)') 'double 1       ' , char(10)
   write(10,rec=rec_index+4,fmt='(A16   )') 'LOOKUP_TABLE def'
   write(10,rec=rec_index+5,fmt='(A15,A1)') 'ault           ' , char(10)

   rec_index = rec_index + 5

   do i = 1,mesh%nc
      write(10,rec=rec_index+swap_index(i),fmt='(ES15.8,A1)') dof%h(i) , char(10)
   end do

   rec_index = rec_index + mesh%nc

   !===================================================================================================================!
   !  Writing VTK file zs cell data
   !===================================================================================================================!

   write(10,rec=rec_index+1,fmt='(A16   )') 'SCALARS         '
   write(10,rec=rec_index+2,fmt='(A16   )') 'zs              '
   write(10,rec=rec_index+3,fmt='(A15,A1)') 'double 1       ' , char(10)
   write(10,rec=rec_index+4,fmt='(A16   )') 'LOOKUP_TABLE def'
   write(10,rec=rec_index+5,fmt='(A15,A1)') 'ault           ' , char(10)

   rec_index = rec_index + 5

   do i = 1,mesh%nc
      write(10,rec=rec_index+swap_index(i),fmt='(ES15.8,A1)') bathy_cell(i) + dof%h(i) , char(10)
   end do

   rec_index = rec_index + mesh%nc

   !===================================================================================================================!
   !   Writing VTK file Manning cell data
   !===================================================================================================================!

   write(10,rec=rec_index+1,fmt='(A16   )') 'SCALARS         '
   write(10,rec=rec_index+2,fmt='(A16   )') 'Manning         '
   write(10,rec=rec_index+3,fmt='(A15,A1)') 'double 1       ' , char(10)
   write(10,rec=rec_index+4,fmt='(A16   )') 'LOOKUP_TABLE def'
   write(10,rec=rec_index+5,fmt='(A15,A1)') 'ault           ' , char(10)

   rec_index = rec_index + 5

   do i = 1,mesh%nc
      write(10,rec=rec_index+swap_index(i),fmt='(ES15.8,A1)') manning( land(i) ) , char(10)
   end do

   rec_index = rec_index + mesh%nc

   !===================================================================================================================!
   !   Writing VTK file u cell data
   !===================================================================================================================!

   write(10,rec=rec_index+1,fmt='(A16   )') 'SCALARS         '
   write(10,rec=rec_index+2,fmt='(A16   )') 'u               '
   write(10,rec=rec_index+3,fmt='(A15,A1)') 'double 1       ' , char(10)
   write(10,rec=rec_index+4,fmt='(A16   )') 'LOOKUP_TABLE def'
   write(10,rec=rec_index+5,fmt='(A15,A1)') 'ault           ' , char(10)

   rec_index = rec_index + 5

   do i = 1,mesh%nc
      write(10,rec=rec_index+swap_index(i),fmt='(ES15.8,A1)') dof%u(i) , char(10)
   end do

   rec_index = rec_index + mesh%nc

   !===================================================================================================================!
   !   Writing VTK file v cell data
   !===================================================================================================================!

   write(10,rec=rec_index+1,fmt='(A16   )') 'SCALARS         '
   write(10,rec=rec_index+2,fmt='(A16   )') 'v               '
   write(10,rec=rec_index+3,fmt='(A15,A1)') 'double 1       ' , char(10)
   write(10,rec=rec_index+4,fmt='(A16   )') 'LOOKUP_TABLE def'
   write(10,rec=rec_index+5,fmt='(A15,A1)') 'ault           ' , char(10)

   rec_index = rec_index + 5

   do i = 1,mesh%nc
      write(10,rec=rec_index+swap_index(i),fmt='(ES15.8,A1)') dof%v(i) , char(10)
   end do

   rec_index = rec_index + mesh%nc

   call mpi_wait_all

END SUBROUTINE v_vtk


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Write a Binary VTK Output Result File
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE v_vtk_bin( dof , mesh , filename )

   USE m_common
   USE m_mpi
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   TYPE( unk ), intent(in)  ::  dof
   TYPE( msh ), intent(in)  ::  mesh

   character(len=*), intent(in)  ::  filename

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  rec_index, mesh_total_size

   !===================================================================================================================!
   !   Creating VTK file and header
   !===================================================================================================================!

   mesh_total_size = mesh%nc

   call mpi_sum_i( mesh_total_size )

   if ( proc == 0 ) then

      open(10,file=filename,status='replace',form= 'unformatted',access='stream',convert='big_endian')

      write(10) "# vtk DataFile Version 3.0"//char(10)
      write(10) "Results file generated by DassFlow"//char(10)
      write(10) "BINARY"//char(10)
      write(10) "DATASET UNSTRUCTURED_GRID"//char(10)

      write(buffer,fmt='(A,I12,A)') 'POINTS ', mesh%nn ,' double'

      write(10) trim(buffer)//char(10)

      write(10) ( mesh%node(i)%coord%x , mesh%node(i)%coord%y , 0.d0 , i = 1,mesh%nn )

      write(10) char(10)

      write(buffer,fmt='(A,2I12)') 'CELLS ' , mesh_total_size , 5*mesh_total_size

      write(10) trim(buffer)//char(10)

      close(10)

   end if

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file=filename,status='old',form= 'unformatted',access='stream',position='append',convert='big_endian')

         write(10) ( 4 , mesh%cell(i)%node(1)-1 , &
                         mesh%cell(i)%node(2)-1 , &
                         mesh%cell(i)%node(3)-1 , &
                         mesh%cell(i)%node(4)-1 , i = 1,mesh%nc )

         if ( proc == np-1 ) write(10) char(10)

         close(10)

      end if

      call mpi_wait_all

   end do

   if ( proc == 0 ) then

      open(10,file=filename,status='old',form= 'unformatted',access='stream',position='append',convert='big_endian')

      write(buffer,fmt='(A,I12)') 'CELL_TYPES ' , mesh_total_size

      write(10) trim(buffer)//char(10)

      write(10) ( 9 , i = 1,mesh_total_size )

      write(10) char(10)

      write(buffer,fmt='(A,I12)') 'CELL_DATA ' , mesh_total_size

      write(10) trim(buffer)//char(10)

      close(10)

   end if

   !===================================================================================================================!
   !   Writing VTK file bathy cell data
   !===================================================================================================================!

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file=filename,status='old',form= 'unformatted',access='stream',position='append',convert='big_endian')

   	   if ( proc == 0    ) write(10) 'SCALARS '//'bathy'//' double 1'//char(10)
         if ( proc == 0    ) write(10) 'LOOKUP_TABLE default'//char(10)
                             write(10) bathy_cell(1:mesh%nc)
         if ( proc == np-1 ) write(10) char(10)

         close(10)

      end if

      call mpi_wait_all

   end do

   !===================================================================================================================!
   !   Writing VTK file h cell data
   !===================================================================================================================!

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file=filename,status='old',form= 'unformatted',access='stream',position='append',convert='big_endian')

   	   if ( proc == 0    ) write(10) 'SCALARS '//'h'//' double 1'//char(10)
         if ( proc == 0    ) write(10) 'LOOKUP_TABLE default'//char(10)
                             write(10) dof%h(1:mesh%nc)
         if ( proc == np-1 ) write(10) char(10)

         close(10)

      end if

      call mpi_wait_all

   end do

   !===================================================================================================================!
   !   Writing VTK file u cell data
   !===================================================================================================================!

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file=filename,status='old',form= 'unformatted',access='stream',position='append',convert='big_endian')

   	   if ( proc == 0    ) write(10) 'SCALARS '//'u'//' double 1'//char(10)
         if ( proc == 0    ) write(10) 'LOOKUP_TABLE default'//char(10)
                             write(10) dof%u(1:mesh%nc)
         if ( proc == np-1 ) write(10) char(10)

         close(10)

      end if

      call mpi_wait_all

   end do

   !===================================================================================================================!
   !   Writing VTK file v cell data
   !===================================================================================================================!

   do k = 0,np-1

      if ( proc == k ) then

         open(10,file=filename,status='old',form= 'unformatted',access='stream',position='append',convert='big_endian')

   	   if ( proc == 0    ) write(10) 'SCALARS '//'v'//' double 1'//char(10)
         if ( proc == 0    ) write(10) 'LOOKUP_TABLE default'//char(10)
                             write(10) dof%v(1:mesh%nc)
         if ( proc == np-1 ) write(10) char(10)

         close(10)

      end if

      call mpi_wait_all

   end do

END SUBROUTINE v_vtk_bin



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! AJOUT LILIAN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE write_norms_end_timeloop( var , var_names , file_name , nbvars )
!
!    USE m_common
!    USE m_mpi
!
!    implicit none
!
!    !===================================================================================================================!
!    !  Interface Variables
!    !===================================================================================================================!
!
!    integer(ip), intent(in)  ::  nbvars
!
!    real(rp)        , dimension(nbvars), intent(in)  ::  var
!    character(len=*), dimension(nbvars), intent(in)  ::  var_names
!
!    character(len=*), intent(inout)  ::  file_name
!
!    !===================================================================================================================!
!    !  Calling Generic sub_write with Output Options defined in 'input.txt'
!    !===================================================================================================================!
! !#    write(*,*) "w_gnuplot==",  w_gnuplot
!    if ( w_tecplot == 1 ) then
!
!       call sub_write( 'tecplot' ) ; return
!
!    end if
!
!    if ( w_gnuplot == 1 ) then
!
!       call sub_write( 'gnuplot' ) ; return
!
!    end if
!
!    call sub_write( 'default' )
!
!
!    CONTAINS
!
!
!       SUBROUTINE sub_write( file_type )
!
!          !=============================================================================================================!
!          !  Interface Variables
!          !=============================================================================================================!
!
!          character(len=*), intent(in)  ::  file_type
!
!          !=============================================================================================================!
!          !  Local Variables
!          !=============================================================================================================!
!
!          character(len=lchar)  ::  complete_file_name
!
!          integer(ip)  ::  ivar
!
!          !=============================================================================================================!
!          !  Begin Subroutine
!          !=============================================================================================================!
!
!         if(end_time_loop) then
!
!             ! /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!             ! not checked in parallel version, it could give issues ? (i guess not because we are in "end_time_loop" case but beware...
!             !/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
!
!                 select case( file_type )
!
!                      case( 'gnuplot' )
!  ! INITIALIZE NORM RESULT FILE
!         complete_file_name = 'res/'//trim(file_name)//'.dat'
!       if ( all( is_file_open(:) /= complete_file_name ) ) then
!
!          open(10,file=complete_file_name,status='replace',form='formatted')
!
!          write(10,*) '# Gnuplot DataFile Version'
!          write(10,*) '# file generated within subroutine write_norms_end_timeloop in src/common/output'
!          close(10)
!
!          file_open_counter = file_open_counter + 1
!
!          is_file_open( file_open_counter ) = complete_file_name
!
!       end if
!          open(10,file=complete_file_name,status='old',position='append',form='formatted') !
!              write(10,'(10(Ax))') '#', &
!                                  var_names(1) , &
!                                  var_names(2) , &
!                                  var_names(3) , &
!                                  var_names(4) , &
!                                  var_names(5) , &
!                                  var_names(6)  , &
!                                  var_names(7)  , &
!                                  var_names(8)  , &
!                                  var_names(9)
!               close(10)
!
!          open(10,file=complete_file_name,status='old',position='append',form='formatted') !
!
!             write(10,'(9ES15.8)') var(1) , & ! LILIAN : '(8ES15.8)' modifié vers '(8ES16.8)' pour régler des pb de moins "-" !'(9(F16.8x))'
!                                  var(2) , &
!                                  var(3) , &
!                                  var(4) , &
!                                  var(5) , &
!                                  var(6)  , &
!                                  var(7)  , &
!                                  var(8)  , &
!                                  var(9)
!
!
!          close(10)
!
!                      case default
!       write(*,*) 'only case gnuplot is implemented for writting norms'
!                   end select
!
!
!          end if ! end if time_ime_loop
!
!          call mpi_wait_all
!
!       END SUBROUTINE sub_write
!
!
! END SUBROUTINE write_norms_end_timeloop
