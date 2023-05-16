MODULE m_user_data

   USE m_common
   USE m_mesh
   USE m_sw

   implicit none

   real(rp), parameter  ::  h0  =  10._rp
   real(rp), parameter  ::  a   =  3000._rp
   real(rp), parameter  ::  b   =  5._rp

   CONTAINS

!===========================================================================================================!
!  Bed elevation
!===========================================================================================================!

	real(rp) FUNCTION bathy_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      real(rp)  ::  xx

      xx  =  x  -  demi * lx

      bathy_user  =  h0 * ( xx / a )**2

	END FUNCTION bathy_user

!===========================================================================================================!
!  Free surface elevation ( zs = h + bathy )
!===========================================================================================================!

	real(rp) FUNCTION zs0_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      zs0_user  =  zs_exact( x , y , zero )

	END FUNCTION zs0_user

!===========================================================================================================!
!  Manning coefficient
!===========================================================================================================!

	real(rp) FUNCTION manning_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      manning_user  =  0.001_rp

	END FUNCTION manning_user

!===========================================================================================================!
!  Initial x velocity
!===========================================================================================================!

	real(rp) FUNCTION u0_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      u0_user  =  0._rp

	END FUNCTION u0_user

!===========================================================================================================!
!  Initial y velocity
!===========================================================================================================!

	real(rp) FUNCTION v0_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

		v0_user  =  0._rp

	END FUNCTION v0_user

!===========================================================================================================!
!  Inflow boundary condition
!===========================================================================================================!

	real(rp) FUNCTION inflow_user( t , x , y )

      implicit none

		real(rp), intent(in)  ::  t , x , y

		inflow_user  =  0.0_rp

	END FUNCTION inflow_user

!===========================================================================================================!
!  Outflow boundary condition
!===========================================================================================================!

	real(rp) FUNCTION outflow_user( t , x , y )

      implicit none

		real(rp), intent(in)  ::  t , x , y

		outflow_user  =  0.0_rp

	END FUNCTION outflow_user

!===========================================================================================================!
!  Test Constants
!===========================================================================================================!

	real(rp) FUNCTION om( a , h0 )

      implicit none

		real(rp), intent(in)  ::  a , h0

      om  =  sqrt( 8._rp * g * h0 / a**2 )

	END FUNCTION om

	real(rp) FUNCTION s( om , tau )

      implicit none

		real(rp), intent(in)  ::  om , tau

      s  =  demi * sqrt( om**2 - tau**2 )

	END FUNCTION s

!===========================================================================================================!
!  Exact water elevation
!===========================================================================================================!

	real(rp) FUNCTION zs_exact( x , y , t )

      implicit none

		real(rp), intent(in)  ::  x , y , t

      real(rp)  ::  tau , ss , xx

      xx   =  x - demi * lx

      tau  =  manning_user( xx , y )

      ss   =  s( om( a , h0 ) , tau )

      zs_exact  =  ( a**2 * b**2 * exp( - tau * t ) ) / ( 8._rp * g**2 * h0 )

      zs_exact  =  zs_exact * ( - ss * tau * sin( two * ss * t ) + &
                                ( d1p4 * tau**2 - ss**2 ) * cos( two * ss * t ) )

      zs_exact  =  zs_exact - ( b**2 * exp( - tau * t ) ) / ( 4._rp * g )

      zs_exact  =  zs_exact - exp( - demi * tau * t ) * &
                              b * ( ss * cos( ss * t ) + tau * demi * sin( ss * t ) ) * xx / g

      zs_exact  =  zs_exact + h0

	END FUNCTION zs_exact

!===========================================================================================================!
!  Exact x velocity
!===========================================================================================================!

	real(rp) FUNCTION u_exact( x , y , t )

      implicit none

		real(rp), intent(in)  ::  x , y , t

      real(rp)  ::  tau , ss

      tau  =  manning_user( x , y )

      ss   =  s( om( a , h0 ) , tau )

      u_exact  =  b * exp( - demi * tau * t ) * sin( ss * t )

	END FUNCTION u_exact

!===========================================================================================================!
!  Exact y velocity
!===========================================================================================================!

	real(rp) FUNCTION v_exact( x , y , t )

      implicit none

		real(rp), intent(in)  ::  x , y , t

      v_exact  =  0._rp

	END FUNCTION v_exact

END MODULE m_user_data
