PROGRAM perturb_obs

   implicit none

   integer(4)  ::  i , j , k

   integer(4), parameter  ::  nb_dx = 4
   integer(4), parameter  ::  nb_dt = 240

   real(8), parameter  ::  max_perturb = 0.1

   real(8)  ::  r

   real(8), dimension(:,:), allocatable  ::  t , h , u , v

   character(len=128)  ::  file_out , file_in , buffer

   allocate( t( nb_dt , nb_dx ) )
   allocate( h( nb_dt , nb_dx ) )
   allocate( u( nb_dt , nb_dx ) )
   allocate( v( nb_dt , nb_dx ) )

   do j = 1,nb_dx

      write(file_in,'(A,I4.4,A)') 'obs_station_' , j , '.plt'

      open(10,file=file_in,status='old',form='formatted')

      read(10,'(A)') buffer

      do i = 1,nb_dt

         read(10,*) t(i,j) , h(i,j) , u(i,j) , v(i,j)

      end do

      close(10)

   end do

   do j = 1,nb_dx

      write(file_out,'(A,I4.4,A)') 'obs_station_' , j , '.plt'

      open(10,file=file_out,status='replace',form='formatted')

      write(10,'(A)') buffer

      do i = 1,nb_dt

         call random_number( r )

         r = - 1._8 + 2._8 * r

         if ( h(i,j) < 1.d-8 ) then

            write(10,'(4ES23.15)') t(i,j) , 0._8                                   , u(i,j) , v(i,j)

         else

            write(10,'(4ES23.15)') t(i,j) , max( 0._8 , h(i,j) + max_perturb * r ) , u(i,j) , v(i,j)

         end if

      end do

      close(10)

   end do

END PROGRAM perturb_obs
