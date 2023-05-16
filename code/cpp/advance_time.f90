SUBROUTINE advance_time( dof , mesh )
   USE m_common
   USE m_mesh
   USE m_model
   USE m_time_screen !NOADJ
   implicit none
   TYPE( unk ), intent(in) :: dof
   TYPE( msh ), intent(in) :: mesh
   real(rp) :: c !> wave celerity ? (sqrt(g*h) )
   real(rp) :: vel !> flow velocity (sqrt(u2+v2) )
   real(rp) :: dist !> distance calculated as if it was a squared mesh dist = 2*surface/perimeter
   real(rp) :: s !> la réponse D
   real(rp) :: dt_min !> smalest timestep , temporary allocation for loop on mesh cells
   real(rp) :: test ! store velocity value in this variable in some cases
   integer(ip) :: iL , iR
   real(rp) :: hL , hR , hM , uM , dL , dR
   integer(ip) :: imin
   nt = nt + 1
   if ( nt == max_nt_for_direct ) end_time_loop = .true.
   if ( adapt_dt == 1 ) then
      dt = dt
      if ( fix_time_step_serie == 2 ) then
         read(80) tc , dt , end_time_loop
         return
      end if
      dt = hugem ! hugem is Machine overflow - precision limits numbers
      do i = 1,mesh%nc
         if ( dof%h(i) > heps ) then
            dist = 2._rp * mesh%cell(i)%surf / mesh%cell(i)%peri
            c = sqrt( g * dof%h(i) )
            vel = sqrt( dof%u(i)**2 + dof%v(i)**2 )
            dt_min = min( dt , dist / ( vel + c ) )
            if (dt_min <= dt) then
    dt = dt_min
    imin = i
    test = vel
            end if
         end if
      end do
      dt = cfl * dt
 if(dt >100000) dt = 0.1
      if ( tc + dt < ts ) then
         tc = tc + dt
      else
         dt = ts - tc ; tc = ts
         end_time_loop = .true.
      end if
      call write_scalar_in_time( dt , 'time_step' )
      if ( fix_time_step_serie == 1 ) write(80) tc , dt , end_time_loop
   else if ( adapt_dt == 2 ) then
      dt = dt
      dt = hugem
      do i = 1,mesh%nc
         if ( dof%h(i) > heps ) then
            dist = 2._rp * mesh%cell(i)%surf / mesh%cell(i)%peri
            c = sqrt( g * dof%h(i) )
            vel = sqrt( dof%u(i)**2 + dof%v(i)**2 )
            dt = min( dt , dist / ( vel + c ) )
         end if
      end do
      dt = cfl * dt
     call write_scalar_in_time( dt , 'time_step_exp' )
      dt = hugem
      do ie = 1,mesh%ne
         if ( .not. mesh%edge(ie)%boundary ) then
            iL = mesh%edge(ie)%cell(1)
            iR = mesh%edge(ie)%cell(2)
            hL = dof%h( iL )
            hR = dof%h( iR )
            dL = mesh%cell( iL )%peri * mesh%cell( iL )%invsurf
            dR = mesh%cell( iR )%peri * mesh%cell( iR )%invsurf
            uM = 0.5_rp * ( mesh%edge(ie)%normal%x * ( dof%u( iL ) + dof%u( iR ) ) + &
                              mesh%edge(ie)%normal%y * ( dof%v( iL ) + dof%v( iR ) ) )
            uM = abs( uM )
            hM = max( dL * hR / hL , dR * hL / hR )
            dt = min ( dt , 1._rp / max( zerom , uM * hM ) )
         end if
      end do
      dt = cfl * dt
      if ( tc + dt < ts ) then
         tc = tc + dt
      else
         dt = ts - tc ; tc = ts
         end_time_loop = .true.
      end if
      call write_scalar_in_time( dt , 'time_step_imp' )
   else
      tc = real(nt,rp) * dt
      if ( tc + dt + zerom > ts ) end_time_loop = .true.
   end if
   call Print_Screen( 'dt' ) !NOADJ
END SUBROUTINE advance_time
