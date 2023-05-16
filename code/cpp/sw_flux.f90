SUBROUTINE sw_lf( hL , uL , vL , hR , uR , vR , flux, surfL, surfR, periL, periR,bathyL,bathyR )
   USE m_common
   USE m_model
   implicit none
   real(rp), intent(in) :: hL , uL , vL ! dof for left cell of edge ie
   real(rp), intent(in) :: hR , uR , vR ! dof for right cell of edge ie
   real(rp), intent(in) :: surfL , surfR ! surfaces of legt and girht cells of the edge ie
   real(rp), intent(in) :: periL , periR ! surfaces of legt and girht cells of the edge ie
   real(rp), intent(in) :: bathyL , bathyR ! bathymethry of legt and girht cells of the edge ie
   real(rp), dimension(3), intent(out) :: flux
   real(rp) :: hL_rusanov ! h value for left cell of edge ie
   real(rp) :: hR_rusanov ! h value for right cell of edge ie
   real(rp) :: c1 , c2 ! wave speeds
   real(rp) :: lambda_rusanov ! true_wave speed
   real(rp) :: gamma_lf ! gamma value for low froude scheme
   real(rp) :: alpha_lf ! alpha value for low froude scheme
   real(rp) :: phiR_lf, phiL_lf ! entropy values for low froude scheme
   real(rp) :: phiR_rusanov, phiL_rusanov ! entropy values for rusanov scheme
   real(rp), dimension(3) :: fL_rusanov , fR_rusanov ! intermediate values of flux
   real(rp), dimension(3) :: flux_rusanov
   real(rp), dimension(3) :: fL_lf , fR_lf ! intermediate values of flux
   real(rp), dimension(3) :: flux_lf
   real(rp) :: tmp ! temporary variable for calculus
   gamma_lf = 1._rp !
   alpha_lf = 1._rp
   phiR_lf = 0._rp
   phiL_lf = 0._rp
   ! "famous" hydrostatic reconstruction of andusse et al
   tmp = max(bathyL, bathyR)
   hL_rusanov= max(0._rp, hL + bathyL - tmp )
   tmp = max(bathyL, bathyR)
   hR_rusanov= max(0._rp, hR + bathyR - tmp)
   ! reconstruction of pressure potential to find balance when water is below bathymetry of a dry neighbouring cell
   if(hL + bathyL < bathyR .AND. hR<hL) then
        phiR_rusanov = phiL_lf
   else
        phiR_rusanov = phiR_lf
   end if
   if(hR + bathyR < bathyL .AND. hL<hR) then
        phiL_rusanov = phiR_lf
   else
        phiL_rusanov = phiL_lf
   end if
   !===================================================================================================================!
   ! rusanov flux computation
   !===================================================================================================================!
   ! ==================== Wave speed computation ================================================================= !
   c1 = ABS(UL) + sqrt( g * hL )
   c2 = ABS(UR) + sqrt( g * hR )
   lambda_rusanov = max( c1, c2 )
    ! ==================== Left and Right flux computation ========================================== !
   fL_rusanov(1) = hL_rusanov * uL
   fL_rusanov(2) = hL_rusanov * uL * uL + 0.5_rp * g * hL_rusanov * hL_rusanov
   fR_rusanov(1) = hR_rusanov * uR
   fR_rusanov(2) = hR_rusanov * uR * uR + 0.5_rp * g * hR_rusanov * hR_rusanov
      ! ==================== Flux computation ================================================================= !
   flux_rusanov(1) = 0.5_rp * ( fL_rusanov(1) + fR_rusanov(1) ) - lambda_rusanov * ( hR_rusanov - hL_rusanov ) !FIX lambda calculation
   flux_rusanov(2) = 0.5_rp * ( fL_rusanov(2) + fR_rusanov(2) ) - lambda_rusanov * ( fR_rusanov(1) - fL_rusanov(1)) !FIX lambda calculation
   ! fix phi
   flux_rusanov(3) = 0.5_rp * ( phiR_rusanov + phiL_rusanov)
   !===================================================================================================================!
   ! LOW FROUDE flux computation
   !===================================================================================================================!
    ! ==================== Left and Right flux computation ========================================== !
   fL_lf(1) = hL * uL
   fR_lf(1) = hR * uR
      ! ==================== Flux computation ================================================================== !
 ! 1/8 = 0.125
   flux_lf(1) = 0.5_rp * ( fL_lf(1) + fR_lf(1) ) - 0.125_rp * gamma_lf * ( periL/surfL + periR/surfR ) * (phiL_lf- phiR_lf)
   ! ?????? flux_lf(1) + flux_lf(1)- ?????????
   flux_lf(2) = uL * flux_lf(1) +uR * flux_lf(1)
   ! 1/4 = 0.25 ! g in module ???
   flux_lf(3) = 0.5_rp * ( phiR_lf + phiL_lf ) - 0.25_rp * alpha_lf * g * ( periL/surfL + periR/surfR ) * (hL * vL - hR * vR ) !FIX lambda calculation
END SUBROUTINE sw_lf
SUBROUTINE sw_hll( hL , uL , vL , hR , uR , vR , flux )
   USE m_common
   USE m_model
   implicit none
   real(rp), intent(in) :: hL , uL , vL
   real(rp), intent(in) :: hR , uR , vR
   real(rp), dimension(3), intent(out) :: flux
   real(rp) :: sL , sR , cL , cR
   real(rp), dimension(3) :: fL , fR
   !===================================================================================================================!
   ! Wave speed computation
   !===================================================================================================================!
   cL = sqrt( g * hL )
   cR = sqrt( g * hR )
   sL = min( 0._rp , uL - cL , uR - 2._rp * cR + cL )
   sR = max( 0._rp , uR + cR , uL + 2._rp * cL - cR )
   if ( ( sL > - zerom .and. sR < zerom ) .or. &
        ( hL < zerom .and. hR < zerom ) ) then
      flux(1:3) = 0._rp ; return
   end if
   !===================================================================================================================!
   ! Left and Right flux computation
   !===================================================================================================================!
   fL(1) = hL * uL
   fL(2) = hL * uL * uL + 0.5_rp * g * hL * hL
   fL(3) = hL * uL * vL
   fR(1) = hR * uR
   fR(2) = hR * uR * uR + 0.5_rp * g * hR * hR
   fR(3) = hR * uR * vR
   !===================================================================================================================!
   ! LF flux computation
   !===================================================================================================================!
   flux(1) = sR * fL(1) - sL * fR(1) + sL * sR * ( hR - hL )
   flux(2) = sR * fL(2) - sL * fR(2) + sL * sR * ( hR * uR - hL * uL )
   flux(3) = sR * fL(3) - sL * fR(3) + sL * sR * ( hR * vR - hL * vL )
   flux(1:3) = flux(1:3) / ( sR - sL )
END SUBROUTINE sw_hll
SUBROUTINE sw_hllc( hL , uL , vL , hR , uR , vR , flux )
   USE m_common
   USE m_model
   implicit none
   real(rp), intent(in) :: hL , uL , vL
   real(rp), intent(in) :: hR , uR , vR
   real(rp), dimension(3), intent(out) :: flux
   real(rp) :: sL , sR , sM , cL , cR
   real(rp), dimension(3) :: fL , fR
   !===================================================================================================================!
   ! Wave speed computation
   !===================================================================================================================!
   cL = sqrt( g * hL )
   cR = sqrt( g * hR )
   sL = min( 0._rp , uL - cL , uR - 2._rp * cR + cL )
   sR = max( 0._rp , uR + cR , uL + 2._rp * cL - cR )
   if ( ( sL > - zerom .and. sR < zerom ) .or. &
        ( hL < zerom .and. hR < zerom ) ) then
      flux(1:3) = 0._rp ; return
   end if
   !===================================================================================================================!
   ! Left and Right flux computation
   !===================================================================================================================!
   fL(1) = hL * uL
   fL(2) = hL * uL * uL + 0.5_rp * g * hL * hL
   fR(1) = hR * uR
   fR(2) = hR * uR * uR + 0.5_rp * g * hR * hR
   !===================================================================================================================!
   ! hll flux computation
   !===================================================================================================================!
   flux(1) = sR * fL(1) - sL * fR(1) + sL * sR * ( hR - hL )
   flux(2) = sR * fL(2) - sL * fR(2) + sL * sR * ( fR(1) - fL(1) )
   flux(1:2) = flux(1:2) / ( sR - sL )
   !===================================================================================================================!
   ! hllC flux computation
   !===================================================================================================================!
   sM = ( sL * hR * uR - sR * hL * uL - sL * sR * ( hR - hL ) ) / ( hR * ( uR - sR ) - hL * ( uL - sL ) )
   if ( sM < 0._rp ) then
      flux(3) = flux(1) * vR
   else if ( sM > 0._rp ) then
      flux(3) = flux(1) * vL
   else
      flux(3) = flux(1) * demi * ( vL + vR )
   end if
END SUBROUTINE sw_hllc
SUBROUTINE sw_god( hL , uL , vL , hR , uR , vR , flux )
   USE m_common
   USE m_model
   implicit none
   real(rp), intent(in) :: hL , uL , vL
   real(rp), intent(in) :: hR , uR , vR
   real(rp), dimension(3), intent(out) :: flux
   real(rp) :: cL , cR , cS , cC , hS , uS1 , uS2 , uu , hC , uC , vC
   logical :: L0 , L1 , L2 , L3 , L4 , L5
   !===================================================================================================================!
   ! Wave speed computation
   !===================================================================================================================!
   cL = sqrt( g * hL )
   cR = sqrt( g * hR )
   !===================================================================================================================!
   !
   !===================================================================================================================!
   cS = ( uL - uR + 2._rp * ( cL + cR ) ) * d1p4
   if ( cS < 0._rp ) then
      cS = 0._rp
      hS = 0._rp
      uS1 = uL
      uS2 = uR
   else
      hS = cS**2 / g
      uS1 = ( uL + uR ) * demi + cL - cR
      uS2 = uS1
   end if
   L0 = cS > cL
   L1 = cS > cR
   !===================================================================================================================!
   ! Case uC = uL
   !===================================================================================================================!
   uu = uL - cL
   L2 = uu < 0._rp
   L3 = uS1 - cS + uu < 0._rp
   if ( ( (.not.L0) .and. (.not.L2) ) .or. &
        ( L0 .and. (.not.L3) ) ) then
      flux(1) = hL * uL
      flux(2) = hL * uL * uL + 0.5_rp * g * hL * hL
      flux(3) = hL * uL * vL
      return
   end if
   !===================================================================================================================!
   ! Case uC = uR
   !===================================================================================================================!
   uu = uR + cR
   L4 = uu > 0._rp
   L5 = uS2 + cS + uu > 0._rp
   if ( ( (.not.L1) .and. (.not.L4) ) .or. &
        ( L1 .and. (.not.L5) ) ) then
      flux(1) = hR * uR
      flux(2) = hR * uR * uR + 0.5_rp * g * hR * hR
      flux(3) = hR * uR * vR
      return
   end if
   !===================================================================================================================!
   ! Case uC = OMEGA L inter LAMBDA L
   !===================================================================================================================!
   uu = uS1 - cS
   if ( uu >= 0._rp .and. L2 ) then
      cC = ( uL + 2.0 * cL ) * d1p3
      hC = cC * cC / g
      uC = cC
      vC = vL
      flux(1) = hC * uC
      flux(2) = hC * uC * uC + 0.5_rp * g * hC * hC
      flux(3) = hC * uC * vC
      return
   end if
   !===================================================================================================================!
   ! Case uC = OMEGA R inter LAMBDA R
   !===================================================================================================================!
   uu = uS2 + cS
   if ( uu <= 0._rp .and. L4 ) then
      cC = ( uR - 2._rp * cR ) * d1p3
      hC = cC * cC / g
      uC = cC
      vC = vR
      flux(1) = hC * uC
      flux(2) = hC * uC * uC + 0.5_rp * g * hC * hC
      flux(3) = hC * uC * vC
      return
   end if
   !===================================================================================================================!
   ! Case uC = uS
   !===================================================================================================================!
   hC = hS
   uC = uS1
   if ( UC < 0._rp ) then
      vC = vR
   else
      vC = vL
   end if
   flux(1) = hC * uC
   flux(2) = hC * uC * uC + 0.5_rp * g * hC * hC
   flux(3) = hC * uC * vC
END SUBROUTINE sw_god
SUBROUTINE sw_balanced_hllc( hL , uL , vL , zbL , hR , uR , vR , zbR , flux , source )
   USE m_common
   USE m_model
   implicit none
   real(rp), intent(in) :: hL , uL , vL , zbL
   real(rp), intent(in) :: hR , uR , vR , zbR
   real(rp), dimension(3), intent(out) :: flux
   real(rp), intent(out) :: source
   real(rp) :: sL , sR , sM , cL , cR
   real(rp), dimension(3) :: fL , fR
   !===================================================================================================================!
   ! Wave speed computation
   !===================================================================================================================!
   cL = sqrt( g * hL )
   cR = sqrt( g * hR )
   sL = min( 0._rp , uL - cL , uR - 2._rp * cR + cL )
   sR = max( 0._rp , uR + cR , uL + 2._rp * cL - cR )
   !===================================================================================================================!
   ! Left and Right flux computation
   !===================================================================================================================!
   fL(1) = hL * uL
   fL(2) = hL * uL * uL + demi * g * hL * hL
   fR(1) = hR * uR
   fR(2) = hR * uR * uR + demi * g * hR * hR
   !===================================================================================================================!
   ! hll flux computation
   !===================================================================================================================!
   if ( ( sL > - zerom .and. sR < zerom ) .or. &
        ( hL < zerom .and. hR < zerom ) ) then
      flux(1:3) = 0._rp
      source = 0._rp
      return
   else
      flux(1) = sR * fL(1) - sL * fR(1) + sL * sR * ( ( hR + zbR ) - ( hL + zbL ) )
      flux(2) = sR * fL(2) - sL * fR(2) + sL * sR * ( fR(1) - fL(1) ) - &
                  0.25_rp * ( sL + sR ) * g * ( hL + hR ) * ( zbR - zbL )
      flux(1:2) = flux(1:2) / ( sR - sL )
      source = 0.25_rp * g * ( hL + hR ) * ( zbR - zbL )
   end if
   !===================================================================================================================!
   ! hllC flux computation
   !===================================================================================================================!
   sM = ( sL * hR * uR - sR * hL * uL - sL * sR * ( ( hR + zbR ) - ( hL + zbL ) ) ) / &
          ( hR * ( uR - sR ) - hL * ( uL - sL ) )
   if ( sM < 0._rp ) then
      flux(3) = flux(1) * vR
   else if ( sM > 0._rp ) then
      flux(3) = flux(1) * vL
   else
      flux(3) = flux(1) * demi * ( vL + vR )
   end if
END SUBROUTINE sw_balanced_hllc
SUBROUTINE sw_balanced_hllc_src_out( hL , uL , vL , zbL , hR , uR , vR , zbR , flux , source )
   USE m_common
   USE m_model
   implicit none
   real(rp), intent(in) :: hL , uL , vL , zbL
   real(rp), intent(in) :: hR , uR , vR , zbR
   real(rp), dimension(3), intent(out) :: flux
   real(rp), dimension(2), intent(out) :: source
   real(rp) :: sL , sR , sM , cL , cR
   real(rp), dimension(3) :: fL , fR
   !===================================================================================================================!
   ! Wave speed computation
   !===================================================================================================================!
   cL = sqrt( g * hL )
   cR = sqrt( g * hR )
   sL = min( 0._rp , uL - cL , uR - 2._rp * cR + cL )
   sR = max( 0._rp , uR + cR , uL + 2._rp * cL - cR )
   !===================================================================================================================!
   ! Left and Right flux computation
   !===================================================================================================================!
   fL(1) = hL * uL
   fL(2) = hL * uL * uL + demi * g * hL * hL
   fR(1) = hR * uR
   fR(2) = hR * uR * uR + demi * g * hR * hR
   !===================================================================================================================!
   ! hll flux computation
   !===================================================================================================================!
   if ( ( sL > - zerom .and. sR < zerom ) .or. &
        ( hL < zerom .and. hR < zerom ) ) then
      flux(1:3) = 0._rp
      source(:) = 0._rp
      return
   else
      flux(1) = sR * fL(1) - sL * fR(1) + sL * sR * ( ( hR + zbR ) - ( hL + zbL ) )
      flux(2) = sR * fL(2) - sL * fR(2) + sL * sR * ( fR(1) - fL(1) )
      flux(1:2) = flux(1:2) / ( sR - sL )
      source(1) = 0.25_rp * g * ( hL + hR ) * ( zbR - zbL )
      source(2) = - 0.25_rp * g * ( hL + hR ) * ( zbR - zbL ) * ( sL + sR ) / ( sR - sL )
   end if
   !===================================================================================================================!
   ! hllC flux computation
   !===================================================================================================================!
   sM = ( sL * hR * uR - sR * hL * uL - sL * sR * ( ( hR + zbR ) - ( hL + zbL ) ) ) / &
          ( hR * ( uR - sR ) - hL * ( uL - sL ) )
   if ( sM < 0._rp ) then
      flux(3) = flux(1) * vR
   else if ( sM > 0._rp ) then
      flux(3) = flux(1) * vL
   else
      flux(3) = flux(1) * demi * ( vL + vR )
   end if
END SUBROUTINE sw_balanced_hllc_src_out
SUBROUTINE sw_hllc_balanced_muscl_src_out( hL , uL , vL , zbL , hR , uR , vR , zbR , flux , source )
   USE m_common
   USE m_model
   implicit none
   real(rp), intent(in) :: hL(2) , uL , vL , zbL
   real(rp), intent(in) :: hR(2) , uR , vR , zbR
   real(rp), intent(out) :: flux(3) , source(2)
   real(rp) :: sL , sR , sM , cL , cR
   real(rp) :: fL(3) , fR(3)
   real(rp) :: alpha(2) , h_eq(2)
   !===================================================================================================================!
   ! Wave speed computation
   !===================================================================================================================!
   cL = sqrt( g * hL(2) )
   cR = sqrt( g * hR(2) )
   sL = min( 0._rp , uL - cL , uR - 2._rp * cR + cL )
   sR = max( 0._rp , uR + cR , uL + 2._rp * cL - cR )
   !===================================================================================================================!
   ! Left and Right flux computation
   !===================================================================================================================!
   fL(1) = hL(2) * uL
   fL(2) = hL(2) * uL * uL + demi * g * hL(2) * hL(2)
   fR(1) = hR(2) * uR
   fR(2) = hR(2) * uR * uR + demi * g * hR(2) * hR(2)
   !===================================================================================================================!
   ! MUSCL alpha coefficient computation
   !===================================================================================================================!
   if ( abs( hL(1) - hR(1) ) < zerom ) then
      alpha(:) = 0._rp
   else
      alpha(1) = ( hL(2) - hL(1) ) / ( hR(1) - hL(1) )
      alpha(2) = ( hR(2) - hR(1) ) / ( hL(1) - hR(1) )
   end if
   h_eq(1) = hL(2) + ( 1._rp - alpha(1) - alpha(2) ) * zbL
   h_eq(2) = hR(2) + ( 1._rp - alpha(1) - alpha(2) ) * zbR
   !===================================================================================================================!
   ! hll flux computation
   !===================================================================================================================!
   if ( ( sL > - zerom .and. sR < zerom ) ) then
      flux(1:3) = 0._rp
      source(:) = 0._rp
      return
   else
      flux(1) = sR * fL(1) - sL * fR(1) + sL * sR * ( h_eq(2) - h_eq(1) )
      flux(2) = sR * fL(2) - sL * fR(2) + sL * sR * ( fR(1) - fL(1) )
      flux(1:2) = flux(1:2) / ( sR - sL )
      source(1) = 0.25_rp * g * ( hL(1) + hR(1) ) * ( zbR - zbL )
      source(2) = - 0.25_rp * g * ( hL(1) + hR(1) ) * ( zbR - zbL ) * ( sL + sR ) / ( sR - sL )
      source(2) = source(2) + 0.5_rp * g * ( sR * alpha(1) * ( hL(1) + hL(2) ) + &
                                                  sL * alpha(2) * ( hR(1) + hR(2) ) ) * ( zbR - zbL ) / ( sR - sL )
   end if
   !===================================================================================================================!
   ! hllC flux computation
   !===================================================================================================================!
   sM = ( sL * hR(2) * uR - sR * hL(2) * uL - sL * sR * ( hR(2) - hL(2) ) ) / &
          ( hR(2) * ( uR - sR ) - hL(2) * ( uL - sL ) )
   if ( sM < 0._rp ) then
      flux(3) = flux(1) * vR
   else if ( sM > 0._rp ) then
      flux(3) = flux(1) * vL
   else
      flux(3) = flux(1) * demi * ( vL + vR )
   end if
END SUBROUTINE sw_hllc_balanced_muscl_src_out
