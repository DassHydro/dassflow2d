SUBROUTINE sw_post_treatment( dof , mesh )
   USE m_common
   USE m_mesh
   USE m_time_screen !NOADJ
   USE m_model
   implicit none
   !===================================================================================================================!
   ! Interface Variables
   !===================================================================================================================!
   type( unk ), intent(in) :: dof
   type( msh ), intent(in) :: mesh
   !===================================================================================================================!
   ! Local Variables
   !===================================================================================================================!
                                                                                                                 !<NOADJ
   real(rp) :: mass
   real(rp) :: norm_inf_e(3) , norm_L1_e(3) , norm_L2_e(3)
   real(rp) :: norm_inf_s(3) , norm_L1_s(3) , norm_L2_s(3)
   real(rp) :: h_e , u_e , v_e , vel , diff , surf_total
   real(rp) :: discharg , mass_cut_tot !>NOADJ
   integer(ip) :: num_bc
   !===================================================================================================================!
   ! Total flow rates at boundaries inputs and outputs
   !===================================================================================================================!
   do num_bc = 1,bc%nb
      if ( temp_scheme(1:2) == 'rk' .or. &
           temp_scheme(1:4) == 'imex' ) then
         bc%sum_mass_flux( num_bc ) = 0.5_rp * bc%sum_mass_flux( num_bc )
      end if
   end do
                                                                                                           !<NOADJ
   !===================================================================================================================!
   ! Total of water volume
   !===================================================================================================================!
   mass = 0._rp
   do i = 1,mesh%nc
      mass = mass + dof%h(i) * mesh%cell(i)%surf
   end do
   call write_scalar_in_time( mass , 'water_vol' )
   !===================================================================================================================!
   ! Eventual water volume added cutting to zero negative water depths
   !===================================================================================================================!
   mass_cut_tot = mass_cut
   if ( temp_scheme(1:2) == 'rk' .or. &
        temp_scheme(1:4) == 'imex' ) then
      mass_cut_tot = demi * mass_cut_tot
   end if
   call write_scalar_in_time( mass_cut_tot , 'water_vol_num_add' )
   !===================================================================================================================!
   ! Output in File prescribed Stations/Sections at given dtp frequency
   !===================================================================================================================!
   call write_stations( dof, mesh )
   !===================================================================================================================!
   ! Output Total flow rates at boundaries inputs and outputs
   !===================================================================================================================!
   do num_bc = 1,bc%nb
      buffer = ''
      if ( bc%typ(num_bc,1)(1:8) == 'discharg') then
            write(buffer,'(A,I3.3)') 'sum_mass_flux_inflow_' , num_bc!bc%grpf( num_bc )
      elseif ( bc%typ(num_bc,1)(1:6) == 'transm' .or. &
                bc%typ(num_bc,1)(1:8) == 'ratcurve' .or. &
                bc%typ(num_bc,1)(1:7) == 'zspresc' .or. &
                bc%typ(num_bc,1)(1:6) == 'hpresc' ) then
    write(buffer,'(A,I3.3)') 'sum_mass_flux_outflow_' , num_bc!bc%grpf( num_bc )
       elseif ( bc%typ(num_bc,1)(1:11) == 'internal_2D') then
           write(buffer,'(A,I3.3)') 'sum_mass_flux_internalflow_' , num_bc!bc%grpf( num_bc )
      endif
      if ( buffer /= '' ) call write_scalar_in_time( bc%sum_mass_flux( num_bc ) , buffer )
   end do
   !===================================================================================================================!
   ! Test qin/qout
   !===================================================================================================================!
   do num_bc = 1,bc%nb
      discharg = zero
      do ib = 1,mesh%neb
         if ( mesh%edgeb(ib)%group == num_bc ) then
            ie = mesh%edgeb(ib)%ind
            i = mesh%edge(ie)%cell(1)
            if ( mesh%edgeb(ib)%typlim(1:8) == 'discharg' ) then
               discharg = discharg - dof%h(i) * ( dof%u(i) * mesh%edge(ie)%normal%x + &
                                                      dof%v(i) * mesh%edge(ie)%normal%y ) * mesh%edge(ie)%length
            else if ( mesh%edgeb(ib)%typlim(1:6) == 'transm' .or. &
                      mesh%edgeb(ib)%typlim(1:8) == 'ratcurve' .or. &
                      mesh%edgeb(ib)%typlim(1:7) == 'zspresc' .or. &
                      mesh%edgeb(ib)%typlim(1:6) == 'hpresc' ) then
               discharg = discharg + dof%h(i) * ( dof%u(i) * mesh%edge(ie)%normal%x + &
                                                      dof%v(i) * mesh%edge(ie)%normal%y ) * mesh%edge(ie)%length
             end if
         end if
      end do
      buffer = ''
      if ( bc%typ(num_bc,1)(1:8) == 'discharg') then
         if ( mesh_type == 'basic' ) then
            write(buffer,'(A,I3.3)') 'sum_q_inflow_' , num_bc
         else
            write(buffer,'(A,I3.3)') 'sum_q_inflow_' , bc%grpf( num_bc )
         end if
      else if ( bc%typ(num_bc,1)(1:6) == 'transm' .or. &
                bc%typ(num_bc,1)(1:8) == 'ratcurve' .or. &
                bc%typ(num_bc,1)(1:7) == 'zspresc' .or. &
                bc%typ(num_bc,1)(1:6) == 'hpresc' ) then
         if ( mesh_type == 'basic' ) then
            write(buffer,'(A,I3.3)') 'sum_q_outflow_' , num_bc
         else
            write(buffer,'(A,I3.3)') 'sum_q_outflow_' , bc%grpf( num_bc )
         end if
      end if
      if ( buffer /= '' ) call write_scalar_in_time( discharg , buffer )
   end do ! end loop on bc%nb
   !===================================================================================================================!
   ! Rain write
   !===================================================================================================================!
   if (bc_rain /= 0) then
      do k=1,bc%nb_rn
         buffer = ''
         write(buffer,'(A,I3.3)') 'Rain_' , k
         if (buffer /= '') then
          call write_scalar_in_time(bc%rain(k)%qin , buffer )
         end if
      enddo
   endif
END SUBROUTINE sw_post_treatment
