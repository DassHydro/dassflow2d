subroutine gen_basic_channel(nx, ny, lx, ly) !, bc_E, bc_W, nb_timestep_in, nb_timestep_out

   implicit none

   !===================================================================================================================!
   !  SW model parameters
   !===================================================================================================================!

   real(8), parameter  	  ::  slope         =  0.0025
   real(8), parameter     ::  h_n  =  1.
   real(8), parameter     ::  pi  =  3.14159265358979_8
   real(8), parameter     ::  g  =  10.
   real(8), parameter     ::  q  =  2.

   !===================================================================================================================!
   !  Parameters defining problem
   !===================================================================================================================!
   integer(4), intent(in) :: nx,ny
   real(8), intent(in)    :: lx, ly 


   integer(4)  ::  i , j , k
   real(8)  ::  dx , dy


   !===================================================================================================================!
   ! GEN CHANNEL  f(nx,ny,lx,ly)
   !===================================================================================================================!

   open(10,file='files/channel.geo',status='replace',form='formatted')

   write(10,'(A)') '# Simple channel mesh'
   write(10,'(2I8,E15.7)')  nx * ny , (nx-1)*(ny-1) , 0._8

   write(10,'(A)') '# Nodes'

   dx = lx / (nx - 1) 
   dy = ly / (ny - 1 )
   
   k = 0

   do j = 1,ny
      do i = 1,nx

         k = k + 1
         
         write(10,'(I8,3E15.7)') k , &
                                 ( i - 1 ) * dx , &
                                 ( j - 1 ) * dy , &
                                 0._8

      end do
   end do

   write(10,'(A)') '# Cells'

   k = 0

   do j = 1,ny-1
      do i = 1,nx-1

         k = k + 1

         write(10,'(6I8,E15.7)') k , &
                                 i + ( j - 1 ) * nx , &
                                 i + ( j - 1 ) * nx + 1 , &
                                 i + ( j - 1 ) * nx + nx + 1 , &
                                 i + ( j - 1 ) * nx + nx , &
                                 1 , &
                                 bathy_user( 0.5_8 * dx + ( i - 1 ) * dx , &
                                             0.5_8 * dy + ( j - 1 ) * dy )

      end do
   end do

   write(10,'(A)') '# Boundary conditions'

   write(10,'(A,2I8)') 'INLET' , ny - 1 , 0

   do j = 1,ny-1

      write(10,'(3I8,E15.7)') 1 + ( j - 1 ) * ( nx - 1 ) , 4 , 1 , &
                              bathy_user( - 0.5_8 * dx , 0.5_8 * dy + ( j - 1 ) * dy )

   end do

   write(10,'(A,2I8)') 'OUTLET' , ny - 1 , 0

   do j = 1,ny-1

      write(10,'(3I8,E15.7)') j * ( nx - 1 ) , 2 , 2 , &
                              bathy_user( lx + 0.5_8 * dx , 0.5_8 * dy + ( j - 1 ) * dy )

   end do

   close(10)

CONTAINS

	real(8) FUNCTION bathy_user( x , y )

      implicit none

      real(8), intent(in)  ::  x , y

      integer(4)  ::  ite

      real(8)  ::  xx , sum_b

      bathy_user  =  0.5 * q**2 / g * ( 1 / h_ex(lx,y)**2 - 1 / h_ex(x,y)**2 ) + h_ex(lx,y) - h_ex(x,y)

      bathy_user  =  bathy_user - 0.03333333333**2 * q**2 * FUNSimpson( h_ex_p , lx , x , int( 200 * abs( lx - x ) / dx ) )

	END FUNCTION bathy_user
!======================================================================================================================!
!~ !  Exact water depth / derivate
!~ !======================================================================================================================!

	real(8) FUNCTION h_ex( x , y )

      implicit none

		real(8), intent(in)  ::  x , y

      h_ex  =  (4._8/g)**0.33333333333 * ( 1 + 0.5 * exp( - 16._8 * ( x / lx - 0.5 )**2 ) )

	END FUNCTION h_ex

	real(8) FUNCTION h_ex_p( x )

      implicit none

		real(8), intent(in)  ::  x

      h_ex_p  =  1 / h_ex(x,0._8)**3.33333333333

	END FUNCTION h_ex_p

	real(8) FUNCTION dh_ex( x , y )

      implicit none

		real(8), intent(in)  ::  x , y

      dh_ex  =  - 16._8 / lx * (4._8/g)**0.33333333333 * ( x / lx - 0.5 ) * exp( - 16._8 * ( x / lx - 0.5 )**2 )

	END FUNCTION dh_ex
	
	
	
   
!> \brief  Numerical integration using Simpson rule
!! \details  USED for macdonnald's benchmark with perturbated topography
real(8) FUNCTION FUNSimpson( f , xmin , xmax , nbp )

  implicit none
  
  real(8)  ::  f ; external f

  real(8), intent(in)  ::  xmin , xmax

  integer(4), intent(in)  ::  nbp

  real(8)  ::  step , sum_s , x_stage

  integer(4)  ::  i_stage , nb_stage

  if ( modulo( nbp , 2 ) == 0 ) then

	 nb_stage = nbp

  else

	 nb_stage = nbp + 1

  end if

  step  =  ( xmax - xmin ) / real( nb_stage , 8 )

  sum_s  =  f( xmin )

  do i_stage = 1,nb_stage,2

	 x_stage  =  xmin  +  real( i_stage , 8 ) * step

	 sum_s  =  sum_s  +  4._8 * f( x_stage )

  end do

  do i_stage = 2,nb_stage-1,2

	 x_stage  =  xmin  +  real( i_stage , 8 ) * step

	 sum_s  =  sum_s  +  2._8 * f( x_stage )

  end do

  sum_s  =  sum_s  +  f( xmax )

  FUNSimpson = step * sum_s / 3._8

END FUNCTION FUNSimpson

END subroutine gen_basic_channel



subroutine gen_bc(in_type, out_type)

   implicit none
   
   character(10), intent(in) :: in_type, out_type

   open(10,file='files/bc.txt',status='replace',form='formatted')

   write(10,'(A)') '!=================================================================================================!'
   write(10,'(A)') '!  Number of boundary conditions'
   write(10,'(A)') '!=================================================================================================!'

   write(10,'(I1)') 2

   write(10,'(A)') '!=================================================================================================!'
   write(10,'(A)') '!   List of boundary conditions'
   write(10,'(A)') '!=================================================================================================!'

   write(10,'(I1,A,A,A)') 1 ,' ', in_type, ' file'
   write(10,'(I1,A,A,A)') 2 ,' ', out_type  ,' file'

   close(10)
end subroutine gen_bc







subroutine gen_land_use(manning_alpha, manning_beta)

   implicit none
   
   real(8), intent(in)  ::  manning_alpha, manning_beta

   open(10,file='files/land_uses.txt',status='replace',form='formatted')

   write(10,'(A)') '!=================================================================================================!'
   write(10,'(A)') '!  Number of Land Uses'
   write(10,'(A)') '!=================================================================================================!'
   write(10,'(I1)') 1
   write(10,'(A)') '!=================================================================================================!'
   write(10,'(A)') '!  List of Land Uses'
   write(10,'(A)') '!=================================================================================================!'
   write(10,'(I4,E22.15,E22.15)') 1 , manning_alpha, manning_beta

   close(10)


end subroutine gen_land_use



!~ subroutine gen_ratcurve(nrow)

!~    implicit none
   
!~    integer(4), intent(in) :: nrow  
!~    integer(4)  ::  i
   
!~    open(10,file='rating_curve.txt',status='replace',form='formatted')

!~    write(10,'(A)') '!=================================================================================================!'
!~    write(10,'(A)') '!  Number of ratcurve'
!~    write(10,'(A)') '!=================================================================================================!'
!~    write(10,'(I1)') 1
!~    write(10,'(A)') '!=================================================================================================!'
!~    write(10,'(A)') '!  Number of row of ratcurve'
!~    write(10,'(A)') '!=================================================================================================!'
!~    write(10,'(I4,E22.15)') nrow+1 , 0.

!~    do i = 0,nrow
!~       write(10,'(3E15.7)') i*2./nrow , 100._8 * (i*2./nrow)**(5./3.)
!~    end do
   
!~    close(10)
!~ end subroutine gen_ratcurve





!~ subroutine gen_hydrograph(nrow, q, stop_time)

!~    implicit none
   
!~    integer(4), intent(in) :: nrow
!~    real(8), intent(in)  ::  q, stop_time
!~    real(8) ::  time, dt
!~    integer(4)  ::  i
   
!~    dt = stop_time /nrow
   
!~    open(10,file='hydrograph.txt',status='replace',form='formatted')

!~    write(10,'(A)') '!=================================================================================================!'
!~    write(10,'(A)') '!  Number of hydrograph'
!~    write(10,'(A)') '!=================================================================================================!'
!~    write(10,'(I1)') 1
!~    write(10,'(A)') '!=================================================================================================!'
!~    write(10,'(A)') '!  Number of row of hydrograph'
!~    write(10,'(A)') '!=================================================================================================!'
!~    write(10,'(I4,E22.15)') nrow+1 , 0.

!~    do i = 0,nrow
!~       time = i*dt
!~       write(10,'(3E15.7)') time, q
!~    end do
   
!~    close(10)
!~ end subroutine gen_hydrograph




subroutine gen_bc_data(bc_typ, nrow, var1, var2)

   implicit none
   
   character(10), intent(in) :: bc_typ
   integer(4), intent(in) :: nrow
   real, dimension(nrow), intent(in) :: var1, var2
   integer(4)  ::  i
   character(30) :: file_name
   
   if (bc_typ(1:8) .eq. 'discharg') then
        file_name = "files/hydrograph.txt"
   else if (bc_typ(1:8) .eq. 'ratcurve') then
        file_name = "files/rating_curve.txt"
   else if (bc_typ(1:6) .eq. 'hpresc') then
        file_name = "files/hpresc.txt"
   else if (bc_typ(1:6) .eq. 'zpresc') then
        file_name = "files/zpresc.txt"
   end if
   
   open(10,file=file_name,status='replace',form='formatted')

   write(10,'(A)') '!=================================================================================================!'
   write(10,'(2A)') '!  Number of ', file_name
   write(10,'(A)') '!=================================================================================================!'
   write(10,'(I1)') 1
   write(10,'(A)') '!=================================================================================================!'
   write(10,'(2A)') '!  Number of row'
   write(10,'(A)') '!=================================================================================================!'
   write(10,'(I4,E22.15)') nrow, 0.

   do i = 1,nrow
      write(10,'(2E15.7)') var1(i), var2(i)
   end do   
   close(10)
   
end subroutine gen_bc_data



subroutine gen_obs(nx_obs, ny_obs, xmax_obs, ymax_obs, xmin_obs, ymin_obs, dt_obs)

   implicit none
   
   integer(4),  intent(in)   ::  nx_obs, ny_obs
   real(8),  intent(in)  ::  dt_obs

   real(8),  intent(in)  ::  xmin_obs, xmax_obs
   real(8),  intent(in)  ::  ymin_obs, ymax_obs

   real(8), dimension(nx_obs)  ::  x_obs
   real(8), dimension(ny_obs)  ::  y_obs
   real(8) :: dx, dy  
   integer(4)  ::  i , j
   
   dx = ( xmax_obs - xmin_obs ) /  (nx_obs - 1 )
   dy = ( ymax_obs - ymin_obs )
   
   do i = 1,nx_obs

      x_obs(i)  =  xmin_obs + dx * ( i - 1 )

   end do

   do j = 1,ny_obs

      y_obs(j)  =  ymin_obs + dy * ( j - 1 )

   end do

   open(10,file='files/obs.txt',status='replace',form='formatted')

   write(10,'(A)') '!=================================================================================================!'
   write(10,'(A)') '! File Defining Points and/or Lines to output result in files at a prescribed temporal frequency'
   write(10,'(A)') '!=================================================================================================!'
   write(10,'(A)')
   write(10,'(A,I6)')  'stations' , nx_obs * ny_obs
   write(10,'(A)')

   do j = 1,ny_obs
      do i = 1,nx_obs
               write(10,'(4E15.7)') x_obs(i) , y_obs(j) , dt_obs, 0.
      end do
   end do

   write(10,'(A)') 
   write(10,'(A)') 'stations_with_grp		0'
   write(10,'(A)') 
   write(10,'(A)') 'sections		0'
   close(10)
   
end subroutine gen_obs




subroutine  h_true_macdo( x , lx, g, h_true)

      implicit none

		real(8), intent(in)  ::  x , lx, g
		real(8), intent(out) :: h_true

      h_true  =  (4._8/g)**0.33333333333 * ( 1 + 0.5 * exp( - 16._8 * ( x / lx - 0.5 )**2 ) )

END subroutine h_true_macdo
