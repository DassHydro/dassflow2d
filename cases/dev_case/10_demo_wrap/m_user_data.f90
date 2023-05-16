MODULE m_user_data

   USE m_mesh
   USE m_model

   implicit none

!======================================================================================================================!
!  Parameters defining problem
!======================================================================================================================!

   real(rp), parameter  ::  manning_coef  =  0.033_rp
   real(rp), parameter  ::  q  =  2._rp


   CONTAINS

!======================================================================================================================!
!  Exact water depth / derivate
!======================================================================================================================!

	real(rp) FUNCTION h_ex( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      h_ex  =  (4._rp/g)**d1p3 * ( one + demi * exp( - 16._rp * ( x / lx - demi )**2 ) )

	END FUNCTION h_ex

	real(rp) FUNCTION h_ex_p( x )

      implicit none

		real(rp), intent(in)  ::  x

      h_ex_p  =  one / h_ex(x,0._rp)**d10p3

	END FUNCTION h_ex_p

	real(rp) FUNCTION dh_ex( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      dh_ex  =  - 16._rp / lx * (4._rp/g)**d1p3 * ( x / lx - demi ) * exp( - 16._rp * ( x / lx - demi )**2 )

	END FUNCTION dh_ex

!======================================================================================================================!
!  Bed elevation
!======================================================================================================================!

	real(rp) FUNCTION bathy_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      integer(ip)  ::  ite

      real(rp)  ::  xx , sum_b

      bathy_user  =  0!demi * q**2 / g * ( one / h_ex(lx,y)**2 - one / h_ex(x,y)**2 ) + h_ex(lx,y) - h_ex(x,y)

	END FUNCTION bathy_user

!======================================================================================================================!
!  Free surface elevation ( zs = h + bathy )
!======================================================================================================================!

	real(rp) FUNCTION zs0_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      zs0_user  =  bathy_user(x,y)  +  h_ex(x,y)

	END FUNCTION zs0_user

!======================================================================================================================!
!  Manning coefficient
!======================================================================================================================!

	real(rp) FUNCTION manning_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      manning_user  =  manning_coef

	END FUNCTION manning_user

!======================================================================================================================!
!  Initial x velocity
!======================================================================================================================!

	real(rp) FUNCTION u0_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      u0_user  =  q / h_ex(x,y)

	END FUNCTION u0_user

!======================================================================================================================!
!  Initial y velocity
!======================================================================================================================!

	real(rp) FUNCTION v0_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

		v0_user  =  0._rp

	END FUNCTION v0_user

!======================================================================================================================!
!  Inflow boundary condition
!======================================================================================================================!

	real(rp) FUNCTION inflow_user( t , x , y )

      implicit none

		real(rp), intent(in)  ::  t , x , y

		inflow_user  =  ly * q

	END FUNCTION inflow_user

!======================================================================================================================!
!  Outflow boundary condition
!======================================================================================================================!

	real(rp) FUNCTION outflow_user( t , x , y )

      implicit none

		real(rp), intent(in)  ::  t , x , y

		outflow_user  =  h_ex( lx , zero )

	END FUNCTION outflow_user

!======================================================================================================================!
!  Exact water elevation
!======================================================================================================================!

	real(rp) FUNCTION zs_exact( x , y , t )

      implicit none

		real(rp), intent(in)  ::  x , y , t

      zs_exact  =  bathy_user(x,y)  +  h_ex(x,y)

	END FUNCTION zs_exact

!======================================================================================================================!
!  Exact x velocity
!======================================================================================================================!

	real(rp) FUNCTION u_exact( x , y , t )

      implicit none

		real(rp), intent(in)  ::  x , y , t

      u_exact  =  q / h_ex(x,y)

	END FUNCTION u_exact

!======================================================================================================================!
!  Exact y velocity
!======================================================================================================================!

	real(rp) FUNCTION v_exact( x , y , t )

      implicit none

		real(rp), intent(in)  ::  x , y , t

      v_exact  =  0._rp

	END FUNCTION v_exact

END MODULE m_user_data
