SUBROUTINE friction_imex( dof_imp , dof_exp , rk1 , rk2 , mesh )
   USE m_common
   USE m_mesh
   USE m_time_screen
   USE m_model
   implicit none
   TYPE( msh ), intent(in ) :: mesh
   TYPE( unk ), intent(inout) :: dof_imp
   TYPE( unk ), intent(in ) :: dof_exp
   real(rp), intent(in) :: rk1 , rk2
   real(rp) :: vel ! Velocity norm
   real(rp) :: Ks
   real(rp) :: imp_part , exp_part_u , exp_part_v
    do i = 1,mesh%nc
      !================================================================================================================!
      ! Positivity cut-off
      !================================================================================================================!
      if ( dof_imp%h(i) <= heps ) then
         dof_imp%u(i) = 0._rp
         dof_imp%v(i) = 0._rp
      else
         !=============================================================================================================!
         ! Semi-Implicit Treatment of Friction Source Term (Manning/Strickler Formula)
         !=============================================================================================================!
         if ( friction == 1 ) then
            Ks = g * manning( land(i) )**2
            vel = sqrt( dof_exp%u(i)**2 + dof_exp%v(i)**2 )
            if ( dof_exp%h(i) <= heps ) then
               exp_part_u = 0._rp
               exp_part_v = 0._rp
            else
               exp_part_u = rk2 * dt * Ks * dof_exp%u(i) * vel / dof_exp%h(i)**d1p3
               exp_part_v = rk2 * dt * Ks * dof_exp%v(i) * vel / dof_exp%h(i)**d1p3
            end if
            exp_part_u = dof_imp%u(i) - exp_part_u / dof_imp%h(i)
            exp_part_v = dof_imp%v(i) - exp_part_v / dof_imp%h(i)
            vel = sqrt( exp_part_u**2 + exp_part_v**2 )
            imp_part = dof_imp%h(i)**d2p3 + sqrt( dof_imp%h(i)**d4p3 + 4._rp * rk1 * dt * Ks * vel )
            imp_part = 2._rp * dof_imp%h(i)**d2p3 / imp_part
         else if ( friction == 2 ) then
            Ks = manning( land(i) )
            exp_part_u = rk2 * dt * Ks * dof_exp%h(i) * dof_exp%u(i)
            exp_part_v = rk2 * dt * Ks * dof_exp%h(i) * dof_exp%v(i)
            exp_part_u = dof_imp%u(i) - exp_part_u / dof_imp%h(i)
            exp_part_v = dof_imp%v(i) - exp_part_v / dof_imp%h(i)
            imp_part = one / ( one + rk1 * dt * Ks )
         end if
         dof_imp%u(i) = exp_part_u * imp_part
         dof_imp%v(i) = exp_part_v * imp_part
      end if
   end do
END SUBROUTINE friction_imex
