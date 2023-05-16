SUBROUTINE imex_time_step( dof , mesh )
   USE m_common
   USE m_mesh
   USE m_time_screen
   USE m_model
   implicit none
   TYPE( msh ), intent(in ) :: mesh
   TYPE( unk ), intent(inout) :: dof
   TYPE( unk ) :: dof_1 , dof_2 , dof_3 , dof_4
   integer(ip) :: icopy
   real(rp) :: fric2 , fric4
   call alloc_dof( dof_1 , mesh )
   call alloc_dof( dof_2 , mesh )
   call alloc_dof( dof_3 , mesh )
   call alloc_dof( dof_4 , mesh )
   icopy = mesh%nc
   dof_1%h( 1 : icopy ) = dof%h( 1 : icopy )
   dof_1%u( 1 : icopy ) = dof%u( 1 : icopy )
   dof_1%v( 1 : icopy ) = dof%v( 1 : icopy )
   dof_2%h( 1 : icopy ) = dof%h( 1 : icopy )
   dof_3%h( 1 : icopy ) = dof%h( 1 : icopy )
   call friction_imex( dof_1 , dof , demi , zero , mesh )
   dof_2%u( 1 : icopy ) = 2._rp * dof%u( 1 : icopy ) - dof_1%u( 1 : icopy )
   dof_2%v( 1 : icopy ) = 2._rp * dof%v( 1 : icopy ) - dof_1%v( 1 : icopy )
   call friction_imex( dof_2 , dof , demi , zero , mesh )
   dof_3%u( 1 : icopy ) = dof_2%u( 1 : icopy )
   dof_3%v( 1 : icopy ) = dof_2%v( 1 : icopy )
    call muscl_aud_flux_n ( dof_3 , dof , mesh )
   dof_4%h( 1 : icopy ) = dof_3%h( 1 : icopy )
   dof_4%u( 1 : icopy ) = dof_1%u( 1 : icopy ) + dof_2%u( 1 : icopy ) + dof_3%u( 1 : icopy ) - 2._rp * dof%u( 1 : icopy )
   dof_4%v( 1 : icopy ) = dof_1%v( 1 : icopy ) + dof_2%v( 1 : icopy ) + dof_3%v( 1 : icopy ) - 2._rp * dof%v( 1 : icopy )
   call friction_imex( dof_4 , dof , demi , zero , mesh )
   dof_1%h( 1 : icopy ) = dof_4%h( 1 : icopy )
   dof_1%u( 1 : icopy ) = dof_4%u( 1 : icopy )
   dof_1%v( 1 : icopy ) = dof_4%v( 1 : icopy )
 ! remove the condition on the spatial scheme choice as only muscl_b1_b is used
   call muscl_aud_flux_n ( dof_4 , dof , mesh )
   dof%u( 1 : icopy ) = 0.5_rp * ( dof_4%h( 1 : icopy ) * dof_4%u( 1 : icopy ) - &
                                     dof_3%h( 1 : icopy ) * dof_3%u( 1 : icopy ) ) + &
                                     dof_1%h( 1 : icopy ) * dof_1%u( 1 : icopy )
   dof%v( 1 : icopy ) = 0.5_rp * ( dof_4%h( 1 : icopy ) * dof_4%v( 1 : icopy ) - &
                                     dof_3%h( 1 : icopy ) * dof_3%v( 1 : icopy ) ) + &
                                     dof_1%h( 1 : icopy ) * dof_1%v( 1 : icopy )
   dof%h( 1 : icopy ) = 0.5_rp * ( dof_3%h( 1 : icopy ) + dof_4%h( 1 : icopy ) )
   where( dof%h( 1 : icopy ) <= heps )
      dof%u( 1 : icopy ) = 0._rp
      dof%v( 1 : icopy ) = 0._rp
   elsewhere
      dof%u( 1 : icopy ) = dof%u( 1 : icopy ) / dof%h( 1 : icopy )
      dof%v( 1 : icopy ) = dof%v( 1 : icopy ) / dof%h( 1 : icopy )
   end where
   call dealloc_dof( dof_1 )
   call dealloc_dof( dof_2 )
   call dealloc_dof( dof_3 )
   call dealloc_dof( dof_4 )
END SUBROUTINE imex_time_step
