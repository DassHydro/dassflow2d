MODULE m_user_data

   USE m_mesh
   USE m_common
   USE m_model

   implicit none

   !===================================================================================================================!
   !  SW model parameters
   !===================================================================================================================!

   real(rp), parameter  ::  slope         =  0.0025_rp
   real(rp), parameter  ::  manning_coef  =  0.065_rp

   !===================================================================================================================!
   !  Parameters defining problem
   !===================================================================================================================!

   real(rp), parameter  ::  h_n  =  1._rp


CONTAINS


   !===================================================================================================================!
   !  Bed elevation
   !===================================================================================================================!

	real(rp) FUNCTION bathy_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      bathy_user  =  slope * ( lx - x )
      
      if ( x >= 0.25 * lx .and. x <= 0.75_rp * lx ) then

         bathy_user  =  bathy_user  +  0.5_rp * h_n * sin( 4.0_rp * pi * x / lx ) * &
                                                      sin( 1.0_rp * pi * y / ly )

      end if

	END FUNCTION bathy_user

   !===================================================================================================================!
   !  Free surface elevation ( zs = h + bathy )
   !===================================================================================================================!

	real(rp) FUNCTION zs0_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      zs0_user  =  bathy_user(x,y)  +  h_n

	END FUNCTION zs0_user

   !===================================================================================================================!
   !  Manning coefficient
   !===================================================================================================================!

	real(rp) FUNCTION manning_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      manning_user  =  manning_coef

	END FUNCTION manning_user
	
	real(rp) FUNCTION manning_beta_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      manning_beta_user  =  0._rp

	END FUNCTION manning_beta_user

   !===================================================================================================================!
   !  Initial x velocity
   !===================================================================================================================!

	real(rp) FUNCTION u0_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

      u0_user  =  sqrt( slope ) * h_n**d2p3 / manning_coef

	END FUNCTION u0_user

   !===================================================================================================================!
   !  Initial y velocity
   !===================================================================================================================!

	real(rp) FUNCTION v0_user( x , y )

      implicit none

		real(rp), intent(in)  ::  x , y

		v0_user  =  0._rp

	END FUNCTION v0_user

   !===================================================================================================================!
   !  Inflow boundary condition
   !===================================================================================================================!

	real(rp) FUNCTION inflow_user( t , x , y )

      implicit none

		real(rp), intent(in)  ::  t , x , y

      real(rp)  ::  Qn

      Qn = 100._rp * sqrt( slope ) * h_n**d5p3 / manning_coef

      inflow_user  =  Qn * ( 1.0_rp + 1.0_rp * exp( - ( t - 7200._rp )**2 / ( 1800._rp )**2 ) )

	END FUNCTION inflow_user

   !===================================================================================================================!
   !  Outflow boundary condition
   !===================================================================================================================!

	real(rp) FUNCTION outflow_user( t , x , y )

      implicit none

		real(rp), intent(in)  ::  t , x , y

      outflow_user  =  h_n

	END FUNCTION outflow_user


END MODULE m_user_data
