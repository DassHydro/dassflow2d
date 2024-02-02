!======================================================================================================================!
!
!                    DassFlow Version 3.0
!
!======================================================================================================================!
!
!  Copyright University of Toulouse-INSA, Univ. of Strasbourg, INRAE & CNRS (France)
!
!  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
!  DassFlow is a computational software whose purpose is to simulate geophysical free surface flows,
!  designed for variational sensitivities and data assimilation (4D-var). Inverse capabilities are
!  based on the adjoint code generation by a source-to-source algorithmic differentiation (Tapenade software used).
!
!  DassFlow software includes few mostly independent "modules" with common architectures and structures.
!  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
!
!  Many people have contributed to the DassFlow development from the initial version to the latest ones.
!  Current contributions:
!               L. Pujol (PhD Unistra)
!               L. Villenave (PhD student)
!               P.-A. Garambois (INRAE Aix-en-Provence)
!               J. Monnier (INSA & Mathematics Institute of Toulouse IMT).
!               K. Larnier (CS group - IMT-INSA).
!  Former scientific or programming contributions of:
!               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
!               J.-P. Vila (INSA & Mathematics Institute of Toulouse IMT).
!               R. Madec   (Mathematics Institute of Toulouse IMT).
!  plus less recent other developers (M. Honnorat and J. Marin).
!
!  Contact : see the DassFlow webpage
!
!  This software is governed by the CeCILL license under French law and abiding by the rules of distribution
!  of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license
!  as circulated by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
!
!  As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the
!  license, users are provided only with a limited warranty and the software's author, the holder of the economic
!  rights, and the successive licensors have only limited liability.
!
!  In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or
!  developing or reproducing the software by the user in light of its specific status of free software, that may
!  mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
!  experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
!  software's suitability as regards their requirements in conditions enabling the security of their systems and/or
!  data to be ensured and, more generally, to use and operate it in the same conditions as regards security.
!
!  The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you
!  accept its terms.
!
!======================================================================================================================!


SUBROUTINE sw_post_treatment( dof , mesh )

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_time_screen                                                                                              !NOADJ
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( unk ), intent(in)  ::  dof
   type( msh ), intent(in)  ::  mesh

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!
                                                                                                                 !<NOADJ
   real(rp)  ::  mass

   real(rp)  ::  norm_inf_e(3) , norm_L1_e(3) , norm_L2_e(3)
   real(rp)  ::  norm_inf_s(3) , norm_L1_s(3) , norm_L2_s(3)

   real(rp)  ::  h_e , u_e , v_e , vel , diff , surf_total

   real(rp)  ::  discharg , mass_cut_tot                                                                         !>NOADJ

   integer(ip)  ::  num_bc

   !===================================================================================================================!
   !  Total flow rates at boundaries inputs and outputs
   !===================================================================================================================!

   do num_bc = 1,bc%nb

      call mpi_sum_r( bc%sum_mass_flux( num_bc ) )

      if ( temp_scheme(1:2) == 'rk' .or. &
           temp_scheme(1:4) == 'imex' ) then

         bc%sum_mass_flux( num_bc )  =  0.5_rp * bc%sum_mass_flux( num_bc )

      end if

   end do
                                                                                                           !<NOADJ
   !===================================================================================================================!
   !  Total of water volume
   !===================================================================================================================!

   mass  =  0._rp

   do i = 1,mesh%nc

      mass  =  mass  +  dof%h(i) * mesh%cell(i)%surf

   end do

   call mpi_sum_r( mass )

   call write_scalar_in_time( mass , 'water_vol' )

   !===================================================================================================================!
   !  Eventual water volume added cutting to zero negative water depths
   !===================================================================================================================!

   mass_cut_tot  =  mass_cut

   call mpi_sum_r( mass_cut_tot )

   if ( temp_scheme(1:2) == 'rk' .or. &
        temp_scheme(1:4) == 'imex' ) then

      mass_cut_tot  =  demi * mass_cut_tot

   end if

   call write_scalar_in_time( mass_cut_tot , 'water_vol_num_add' )

   !===================================================================================================================!
   !  Output in File prescribed Stations/Sections at given dtp frequency
   !===================================================================================================================!

   call write_stations( dof, mesh )

! call write_sections( dof )

   !===================================================================================================================!
   !  Output Total flow rates at boundaries inputs and outputs
   !===================================================================================================================!

   do num_bc = 1,bc%nb

      buffer = ''

      if ( bc%typ(num_bc,1)(1:8) == 'discharg'  .or. &
      bc%typ(num_bc,1)(1:3) == 'gr4' ) then ! .or. &
!    mesh%edgeb(ib)%typlim(1:17) == 'internal_discharg') DEPRECATED

!          if ( mesh_type == 'basic' ) then

!             write(buffer,'(A,I3.3)') 'sum_mass_flux_inflow_' , num_bc

!          else

            write(buffer,'(A,I3.3)') 'sum_mass_flux_inflow_' , num_bc!bc%grpf( num_bc )

!          end if

      elseif ( bc%typ(num_bc,1)(1:6) == 'transm'   .or. &
                bc%typ(num_bc,1)(1:8) == 'ratcurve' .or. &
                bc%typ(num_bc,1)(1:7) == 'zspresc'  .or. &
                bc%typ(num_bc,1)(1:6) == 'hpresc' ) then
!               bc%typ(num_bc,1)(1:17) == 'internal_ratcurve' .or. & !DEPRECATED

!          if ( mesh_type == 'basic' ) then
!
!             write(buffer,'(A,I3.3)') 'sum_mass_flux_outflow_' , num_bc
!
!          else

            write(buffer,'(A,I3.3)') 'sum_mass_flux_outflow_' , num_bc!bc%grpf( num_bc )

!          endif

       elseif ( bc%typ(num_bc,1)(1:11) == 'internal_2D') then

           write(buffer,'(A,I3.3)') 'sum_mass_flux_internalflow_' , num_bc!bc%grpf( num_bc )

      endif

      if ( buffer /= '' ) call write_scalar_in_time( bc%sum_mass_flux( num_bc )  , buffer )

   end do

   !===================================================================================================================!
   !  Test qin/qout
   !===================================================================================================================!

   do num_bc = 1,bc%nb

      discharg  = zero

      do ib = 1,mesh%neb

         if ( mesh%edgeb(ib)%group  == num_bc ) then

            ie  =  mesh%edgeb(ib)%ind

            i  =  mesh%edge(ie)%cell(1)

            if ( mesh%edgeb(ib)%typlim(1:8) == 'discharg' &
                 .or. bc%typ(num_bc,1)(1:4) == 'gr4') then

               discharg  =  discharg  -  dof%h(i) * ( dof%u(i) * mesh%edge(ie)%normal%x + &
                                                      dof%v(i) * mesh%edge(ie)%normal%y ) * mesh%edge(ie)%length

            else if ( mesh%edgeb(ib)%typlim(1:6) == 'transm'   .or. &
                      mesh%edgeb(ib)%typlim(1:8) == 'ratcurve' .or. &
                      mesh%edgeb(ib)%typlim(1:7) == 'zspresc'  .or. &
                      mesh%edgeb(ib)%typlim(1:6) == 'hpresc' ) then

               discharg  =  discharg  +  dof%h(i) * ( dof%u(i) * mesh%edge(ie)%normal%x + &
                                                      dof%v(i) * mesh%edge(ie)%normal%y ) * mesh%edge(ie)%length

             end if

         end if

      end do

      call mpi_sum_r( discharg )

      buffer = ''

      if ( bc%typ(num_bc,1)(1:8) == 'discharg'  &
      .or. bc%typ(num_bc,1)(1:4) == 'gr4' ) then

         if ( mesh_type == 'basic' ) then

            write(buffer,'(A,I3.3)') 'sum_q_inflow_' , num_bc

         else

            write(buffer,'(A,I3.3)') 'sum_q_inflow_' , bc%grpf( num_bc )

         end if

      else if ( bc%typ(num_bc,1)(1:6) == 'transm'   .or. &
                bc%typ(num_bc,1)(1:8) == 'ratcurve' .or. &
                bc%typ(num_bc,1)(1:7) == 'zspresc'  .or. &
                bc%typ(num_bc,1)(1:6) == 'hpresc' ) then

         if ( mesh_type == 'basic' ) then

            write(buffer,'(A,I3.3)') 'sum_q_outflow_' , num_bc

         else

            write(buffer,'(A,I3.3)') 'sum_q_outflow_' , bc%grpf( num_bc )

         end if

      end if

      if ( buffer /= '' ) call write_scalar_in_time( discharg , buffer )

   end do

   !===================================================================================================================!
   !  Rain write
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

      !===================================================================================================================!
   !  Infil write
   !===================================================================================================================!

   if (bc_infil .ne. 0) then

         buffer = ''

         write(buffer,'(A,I3.3)') 'Mean_infil'

         if (buffer /= '') then

          call write_scalar_in_time(sum(dof%infil)/size(dof%infil) , buffer )

         end if

   endif

!======================================================================================================================!
!  Norms infinity, L1 and L2 in case of validation Makefile flag and a prescribed exact solution
!======================================================================================================================!

#ifdef USE_VALID

      !================================================================================================================!
      !  Output Norms in time or integrated in time or at the end of the simulation (t = ts)
      !================================================================================================================!

      if ( ( w_norm == 1                     ) .or. & ! relative norm in time
           ( w_norm == 2                     ) .or. & ! relative norms integrated in time
           ( w_norm == 3 .and. end_time_loop ) .or. & ! relative norms at the end of the simulation
           ( w_norm == 4                     ) .or. & ! absolute norm in time
           ( w_norm == 5                     ) .or. & ! absolute norm integrated in time
           ( w_norm == 6 .and. end_time_loop ) ) then ! absolute norm  at the end of the simulation

         if ( w_norm == 2 .or. & ! integrated in time
              w_norm == 5 ) then

            norm_inf_s(:)  =  norm_inf(:)
            norm_L1_s (:)  =  norm_L1 (:)
            norm_L2_s (:)  =  norm_L2 (:)

         end if
        ! ALWAYS initialise for this timestep
         norm_inf(:)  =  0._rp ; norm_inf_e(:)  =  0._rp
         norm_L1 (:)  =  0._rp ; norm_L1_e (:)  =  0._rp
         norm_L2 (:)  =  0._rp ; norm_L2_e (:)  =  0._rp

         surf_total  =  0._rp
        !NORM ON EACH CELL
         do i = 1,mesh%nc

            !==========================================================================================================!
            !                           GET EXACT SOLUTIONS                                                          !
            h_e  =  max( 0._rp , zs_exact( mesh%cell(i)%grav%x , mesh%cell(i)%grav%y , tc ) - bathy_cell(i) )

            u_e  =  u_exact( mesh%cell(i)%grav%x , mesh%cell(i)%grav%y , tc )
            v_e  =  v_exact( mesh%cell(i)%grav%x , mesh%cell(i)%grav%y , tc )

            vel  =  sqrt( u_e**2 + v_e**2 )

            !==========================================================================================================!
            !                           update h norm                                                        !
            diff  =  abs( dof%h(i) - h_e )

            norm_inf  (1)  =  max( norm_inf  (1) , diff )
            norm_inf_e(1)  =  max( norm_inf_e(1) , h_e  )

            norm_L1   (1)  =  norm_L1  (1)  +  mesh%cell(i)%surf * diff
            norm_L1_e (1)  =  norm_L1_e(1)  +  mesh%cell(i)%surf * h_e

            norm_L2   (1)  =  norm_L2  (1)  +  mesh%cell(i)%surf * diff**2
            norm_L2_e (1)  =  norm_L2_e(1)  +  mesh%cell(i)%surf * h_e**2

            !==========================================================================================================!
            !                           update u norm                                                        !
            diff  =  abs( sqrt( dof%u(i)**2 + dof%v(i)**2 ) - vel )

            norm_inf  (2)  =  max( norm_inf  (2) , diff )
            norm_inf_e(2)  =  max( norm_inf_e(2) , vel  )

            norm_L1   (2)  =  norm_L1  (2)  +  mesh%cell(i)%surf * diff
            norm_L1_e (2)  =  norm_L1_e(2)  +  mesh%cell(i)%surf * vel

            norm_L2   (2)  =  norm_L2  (2)  +  mesh%cell(i)%surf * diff**2
            norm_L2_e (2)  =  norm_L2_e(2)  +  mesh%cell(i)%surf * vel**2

            !==========================================================================================================!
            !                           update q norm                                                        !
            diff  =  abs( dof%h(i) * sqrt( dof%u(i)**2 + dof%v(i)**2 ) - h_e * vel )

            norm_inf  (3)  =  max( norm_inf  (3) , diff      )
            norm_inf_e(3)  =  max( norm_inf_e(3) , h_e * vel )

            norm_L1   (3)  =  norm_L1  (3)  +  mesh%cell(i)%surf * diff
            norm_L1_e (3)  =  norm_L1_e(3)  +  mesh%cell(i)%surf * h_e * vel

            norm_L2   (3)  =  norm_L2  (3)  +  mesh%cell(i)%surf * diff**2
            norm_L2_e (3)  =  norm_L2_e(3)  +  mesh%cell(i)%surf * ( h_e * vel )**2
            !==========================================================================================================!
            ! Update total surface
            surf_total  =  surf_total  +  mesh%cell(i)%surf

         end do
        ! update norm on all cells for this timestep
         call mpi_sum_r( surf_total )

         call mpi_max_r( norm_inf(1) ) ; call mpi_max_r( norm_inf_e(1) )
         call mpi_sum_r( norm_L1 (1) ) ; call mpi_sum_r( norm_L1_e (1) )
         call mpi_sum_r( norm_L2 (1) ) ; call mpi_sum_r( norm_L2_e (1) )

         call mpi_max_r( norm_inf(2) ) ; call mpi_max_r( norm_inf_e(2) )
         call mpi_sum_r( norm_L1 (2) ) ; call mpi_sum_r( norm_L1_e (2) )
         call mpi_sum_r( norm_L2 (2) ) ; call mpi_sum_r( norm_L2_e (2) )

         call mpi_max_r( norm_inf(3) ) ; call mpi_max_r( norm_inf_e(3) )
         call mpi_sum_r( norm_L1 (3) ) ; call mpi_sum_r( norm_L1_e (3) )
         call mpi_sum_r( norm_L2 (3) ) ; call mpi_sum_r( norm_L2_e (3) )

         if ( w_norm <= 3 ) then ! for relative norms

            norm_inf(1)  =        div_by_except_0( norm_inf(1) , norm_inf_e(1) )
            norm_L1 (1)  =        div_by_except_0( norm_L1 (1) , norm_L1_e (1) )
            norm_L2 (1)  =  sqrt( div_by_except_0( norm_L2 (1) , norm_L2_e (1) ) )

            norm_inf(2)  =        div_by_except_0( norm_inf(2) , norm_inf_e(2) )
            norm_L1 (2)  =        div_by_except_0( norm_L1 (2) , norm_L1_e (2) )
            norm_L2 (2)  =  sqrt( div_by_except_0( norm_L2 (2) , norm_L2_e (2) ) )

            norm_inf(3)  =        div_by_except_0( norm_inf(3) , norm_inf_e(3) )
            norm_L1 (3)  =        div_by_except_0( norm_L1 (3) , norm_L1_e (3) )
            norm_L2 (3)  =  sqrt( div_by_except_0( norm_L2 (3) , norm_L2_e (3) ) )

         else ! for absolute norms

            norm_L1 (:)  =        norm_L1 (:) / surf_total
            norm_L2 (:)  =  sqrt( norm_L2 (:) / surf_total )

         end if

         if ( w_norm == 2 .or. &
              w_norm == 5 ) then !integrated norm in time at this time step

            norm_inf(:)  =  norm_inf(:)  +  norm_inf_s(:)
            norm_L1 (:)  =  norm_L1 (:)  +  norm_L1_s (:)
            norm_L2 (:)  =  norm_L2 (:)  +  norm_L2_s (:)

            if ( end_time_loop ) then  ! mean on the number of timestep

               norm_inf(:)  =  norm_inf(:) / real( nt , 8 )
               norm_L1 (:)  =  norm_L1 (:) / real( nt , 8 )
               norm_L2 (:)  =  norm_L2 (:) / real( nt , 8 )

            end if

         end if
!REMOVED 21/08 LEO TESTS
!          if ( w_norm == 1 ) then ! norm at this timestep (relative)
!
!             call write_scalar_in_time( norm_inf(1) , 'norm_rel_h_inf' )
!             call write_scalar_in_time( norm_L1 (1) , 'norm_rel_h_L1'  )
!             call write_scalar_in_time( norm_L2 (1) , 'norm_rel_h_L2'  )
!
!             call write_scalar_in_time( norm_inf(2) , 'norm_rel_u_inf' )
!             call write_scalar_in_time( norm_L1 (2) , 'norm_rel_u_L1'  )
!             call write_scalar_in_time( norm_L2 (2) , 'norm_rel_u_L2'  )
!
!             call write_scalar_in_time( norm_inf(3) , 'norm_rel_q_inf' )
!             call write_scalar_in_time( norm_L1 (3) , 'norm_rel_q_L1'  )
!             call write_scalar_in_time( norm_L2 (3) , 'norm_rel_q_L2'  )
!
!          else if ( w_norm == 4 ) then ! norm at this timestep (absolute)
!
!             call write_scalar_in_time( norm_inf(1) , 'norm_abs_h_inf' )
!             call write_scalar_in_time( norm_L1 (1) , 'norm_abs_h_L1'  )
!             call write_scalar_in_time( norm_L2 (1) , 'norm_abs_h_L2'  )
!
!             call write_scalar_in_time( norm_inf(2) , 'norm_abs_u_inf' )
!             call write_scalar_in_time( norm_L1 (2) , 'norm_abs_u_L1'  )
!             call write_scalar_in_time( norm_L2 (2) , 'norm_abs_u_L2'  )
!
!             call write_scalar_in_time( norm_inf(3) , 'norm_abs_q_inf' )
!             call write_scalar_in_time( norm_L1 (3) , 'norm_abs_q_L1'  )
!             call write_scalar_in_time( norm_L2 (3) , 'norm_abs_q_L2'  )
!
!          else if ( w_norm == 2 .and. end_time_loop ) then ! norm integrated in time untill this timestep (relative)
! ! write_pscalar(var ,
! !               var_names ,
! !               file_name ,
! !               nbvars)
!
!             call write_pscalar( (/ norm_L1(1) , norm_L2(1) , norm_inf(1) , &                        ! var
!                                    norm_L1(3) , norm_L2(3) , norm_inf(3) /) , &
!                                 (/ 'norm_rel_L1__h' , 'norm_rel_L2__h' , 'norm_rel_inf_h' , &       ! var_names
!                                    'norm_rel_L1__q' , 'norm_rel_L2__q' , 'norm_rel_inf_q' /) ,&
!                                    'norms_rel_int_time' ,&                                          ! file_name
!                                    6 )                                                              ! nbvars
!
!
!             ! AJOUT LILIAN
!             call write_norms_end_timeloop( (/ norm_L1(1) , norm_L2(1) , norm_inf(1) , &
!                                    norm_L1(2) , norm_L2(2) , norm_inf(2), &
!                                    norm_L1(3) , norm_L2(3) , norm_inf(3) /) , &
!                                 (/ 'norm_rel_L1__h' , 'norm_rel_L2__h' , 'norm_rel_inf_h' , &
!                                    'norm_rel_L1__u' , 'norm_rel_L2__u' , 'norm_rel_inf_u' , &
!                                    'norm_rel_L1__q' , 'norm_rel_L2__q' , 'norm_rel_inf_q' /) ,&
!                                    'norms_rel_int_time' ,&
!                                    9 )
!             ! end ajout lilian
!
!
!          else if ( w_norm == 5 .and. end_time_loop ) then ! norm integrated in time untill this timestep (absolute)
!
!             call write_pscalar( (/ norm_L1(1) , norm_L2(1) , norm_inf(1) , &
!                                    norm_L1(3) , norm_L2(3) , norm_inf(3) /) , &
!                                 (/ 'norm_abs_L1__h' , 'norm_abs_L2__h' , 'norm_abs_inf_h' , &
!                                    'norm_abs_L1__q' , 'norm_abs_L2__q' , 'norm_abs_inf_q' /) ,&
!                                    'norms_abs_int_time' ,&
!                                    6 )
!
!             ! AJOUT LILIAN
!             call write_norms_end_timeloop( (/ norm_L1(1) , norm_L2(1) , norm_inf(1) , &
!                                    norm_L1(2) , norm_L2(2) , norm_inf(2), &
!                                    norm_L1(3) , norm_L2(3) , norm_inf(3) /) , &
!                                 (/ 'norm_abs_L1__h' , 'norm_abs_L2__h' , 'norm_abs_inf_h' , &
!                                    'norm_abs_L1__u' , 'norm_abs_L2__u' , 'norm_abs_inf_u' , &
!                                    'norm_abs_L1__q' , 'norm_abs_L2__q' , 'norm_abs_inf_q' /) ,&
!                                    'norms_abs_int_time' ,&
!                                    9 )
!             ! end ajout lilian
!
!          else if ( w_norm == 3 ) then ! norm at the end of the simulation (relative)
!
!             call write_pscalar( (/ norm_L1(1) , norm_L2(1) , norm_inf(1) , &
!                                    norm_L1(3) , norm_L2(3) , norm_inf(3) /) , &
!                                 (/ 'norm_rel_L1__h' , 'norm_rel_L2__h' , 'norm_rel_inf_h' , &
!                                    'norm_rel_L1__q' , 'norm_rel_L2__q' , 'norm_rel_inf_q' /) ,&
!                                    'norms_rel_end_time' ,&
!                                    6 )
!
!             ! AJOUT LILIAN
!             call write_norms_end_timeloop( (/ norm_L1(1) , norm_L2(1) , norm_inf(1) , &
!                                    norm_L1(2) , norm_L2(2) , norm_inf(2), &
!                                    norm_L1(3) , norm_L2(3) , norm_inf(3) /) , &
!                                 (/ 'norm_rel_L1__h' , 'norm_rel_L2__h' , 'norm_rel_inf_h' , &
!                                    'norm_rel_L1__u' , 'norm_rel_L2__u' , 'norm_rel_inf_u' , &
!                                    'norm_rel_L1__q' , 'norm_rel_L2__q' , 'norm_rel_inf_q' /) ,&
!                                    'norms_rel_end_time' ,&
!                                    9 )
!             ! end ajout lilian
!
!          else if ( w_norm == 6 ) then ! norm at the end of the simulation (absolute)
!
!             call write_pscalar( (/ norm_L1(1) , norm_L2(1) , norm_inf(1) , &
!                                    norm_L1(3) , norm_L2(3) , norm_inf(3) /) , &
!                                 (/ 'norm_abs_L1__h' , 'norm_abs_L2__h' , 'norm_abs_inf_h' , &
!                                    'norm_abs_L1__q' , 'norm_abs_L2__q' , 'norm_abs_inf_q' /) ,&
!                                    'norms_abs_end_time' ,&
!                                    6 )
!
!             ! AJOUT LILIAN
!             call write_norms_end_timeloop( (/ norm_L1(1) , norm_L2(1) , norm_inf(1) , &
!                                    norm_L1(2) , norm_L2(2) , norm_inf(2), &
!                                    norm_L1(3) , norm_L2(3) , norm_inf(3) /) , &
!                                 (/ 'norm_abs_L1__h' , 'norm_abs_L2__h' , 'norm_abs_inf_h' , &
!                                    'norm_abs_L1__u' , 'norm_abs_L2__u' , 'norm_abs_inf_u' , &
!                                    'norm_abs_L1__q' , 'norm_abs_L2__q' , 'norm_abs_inf_q' /) ,&
!                                    'norms_abs_end_time' ,&
!                                    9 )
!             ! end ajout lilian
!
!          end if

         if ( end_time_loop ) then

            call Print_Screen( 'norms_h' )
            call Print_Screen( 'norms_q' )
            call Print_Screen( 'norms_u' )

         end if

      end if

#endif   																													!>NOADJ

END SUBROUTINE sw_post_treatment
