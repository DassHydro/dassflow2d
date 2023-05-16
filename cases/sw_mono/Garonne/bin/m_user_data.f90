MODULE m_user_data

   USE m_mesh
   USE m_common

   implicit none

   CONTAINS

!======================================================================================================================!
!  Bed elevation
!======================================================================================================================!

	real(rp) FUNCTION bathy_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      bathy_user = 0.

	END FUNCTION bathy_user

!======================================================================================================================!
!  Free surface elevation ( zs = h + bathy )
!======================================================================================================================!

	real(rp) FUNCTION zs0_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      zs0_user = 0.

	END FUNCTION zs0_user

!======================================================================================================================!
!  Manning coefficient
!======================================================================================================================!

	real(rp) FUNCTION manning_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      manning_user  =  0.05_rp

	END FUNCTION manning_user
	
	real(rp) FUNCTION manning_beta_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      manning_beta_user  =  0._rp

	END FUNCTION manning_beta_user

!======================================================================================================================!
!  Initial x velocity
!======================================================================================================================!

	real(rp) FUNCTION u0_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      u0_user  =  0._rp

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

		inflow_user  =  10._rp * ( 1._rp + 7._rp * exp( - ( t - 14400._rp )**2 / ( 3600._rp )**2 ) )

	END FUNCTION inflow_user

!======================================================================================================================!
!  Outflow boundary condition
!======================================================================================================================!

	real(rp) FUNCTION outflow_user( t , x , y )

      implicit none

		real(rp), intent(in)  ::  t , x , y

		outflow_user  =  0._rp

	END FUNCTION outflow_user

END MODULE m_user_data
