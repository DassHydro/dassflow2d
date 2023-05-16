SUBROUTINE boundary_post( mass_flux , index_ghost , mesh )
   USE m_common
   USE m_mesh
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(in) :: mesh
   real(rp), intent(in) :: mass_flux
   integer(ip), intent(in) :: index_ghost
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   integer(ip) :: group
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   ib = mesh%edge(ie)%lim
   group = mesh%edgeb(ib)%group
   if ( mesh%edgeb(ib)%typlim(1:8) == 'discharg' ) then
      bc%sum_mass_flux( group ) = bc%sum_mass_flux( group ) - mass_flux * mesh%edge(ie)%length
      if ( feedback_inflow == 1 ) then
         bathy_cell( index_ghost ) = bathy_cell( index_ghost ) + &
         coef_feedback * ( mass_flux - bc%inflow( mesh%neb + ib ) )
      end if
   end if
   if ( mesh%edgeb(ib)%typlim(1:6) == 'transm' .or. &
        mesh%edgeb(ib)%typlim(1:7) == 'neumann' .or. &
        mesh%edgeb(ib)%typlim(1:8) == 'ratcurve' .or. &
        mesh%edgeb(ib)%typlim(1:7) == 'zspresc' .or. &
        mesh%edgeb(ib)%typlim(1:6) == 'hpresc' ) then
      bc%sum_mass_flux( group ) = bc%sum_mass_flux( group ) + mass_flux * mesh%edge(ie)%length
   end if
   if ( mesh%edgeb(ib)%typlim(1:11) == 'internal_2D') then
        bc%sum_mass_flux( group ) = bc%sum_mass_flux( group ) + mass_flux * mesh%edge(ie)%length
   endif
END SUBROUTINE boundary_post
                                                                                                                 !<NOADJ
SUBROUTINE spread_mass_added( dof , mesh )
   USE m_common
   USE m_mesh
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( msh ), intent(in) :: mesh
   type( unk ), intent(inout) :: dof
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
   real(rp) :: coef_mass_balance
   integer(ip), dimension( maxed ) :: icut
   !===================================================================================================================!
   ! Begin Subroutine
   !===================================================================================================================!
   k = count( mesh%cell(i)%cell(:) >= 1 .and. mesh%cell(i)%cell(:) <= mesh%nc )
   icut(1:k) = pack ( mesh%cell(i)%cell(:) , mesh%cell(i)%cell(:) >= 1 .and. mesh%cell(i)%cell(:) <= mesh%nc )
   coef_mass_balance = sum( dof%h( icut(1:k) ) * mesh%cell( icut(1:k) )%surf )
   if ( coef_mass_balance > - dof%h(i) * mesh%cell(i)%surf ) then
      coef_mass_balance = 1._rp + dof%h(i) * mesh%cell(i)%surf / coef_mass_balance
      dof%h( icut(1:k) ) = coef_mass_balance * dof%h( icut(1:k) )
   else
      mass_cut = mass_cut - dof%h(i) * mesh%cell(i)%surf - coef_mass_balance
      dof%h( icut(1:k) ) = 0._rp
   end if
   dof%h(i) = 0._rp
END SUBROUTINE spread_mass_added !>NOADJ
