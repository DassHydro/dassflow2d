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


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Module managing minimisation process
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


MODULE m_minimization

   USE m_adjoint

#ifdef USE_SW_MONO
      USE m_model
#endif

#ifdef USE_HYDRO
!       USE m_gr4
#endif

   implicit none

   integer(ip)  ::  ite_line_search


CONTAINS


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Main subroutine of minmization
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE minimize_cost( mesh , dof0 , dof )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( msh ), intent(in   )  ::  mesh
      type( unk ), intent(inout)  ::  dof0
      type( unk ), intent(inout)  ::  dof

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      external  simul_rc , euclid , ctonbe , ctcabe

      integer(ip)  ::  n            ! Control vector size
      real(rp)     ::  dxmin        ! Resolution for x in Linf norm
      real(rp)     ::  dJ_expect    ! Estimation of the expected decrease in cost during the first iteration
      real(rp)     ::  epsg         ! Stopping criterion
      character*3  ::  normtype     ! specifies the norm used to test optimality
      integer(ip)  ::  impres       ! Print option for m1qn3
      integer(ip)  ::  io           ! File output label
      integer(ip)  ::  imode(3)     ! Input mode
      integer(ip)  ::  omode        ! Output mode
      integer(ip)  ::  niter        ! Maximum number of minimization iterations
      integer(ip)  ::  nsim         ! Maximum number of simulator calls
!       integer(ip)  ::  iz(5)        ! Adress of working array for m1qn3
      integer(ip)  ::  ndz          ! Dimension of dz
      integer(ip)  ::  reverse      ! Specifies direct or reverse mode
      integer(ip)  ::  indic        ! m1qn3 indicates it needs cost/gradient

      integer(ip)  ::  izs(1)       !
      real(rp)     ::  rzs(1)       ! Working areas not needed in reverse mode
      real(rp)     ::  dzs(1)       !

      integer(ip)  :: max_ite_line_search
      real(rp), dimension(:), allocatable  ::  dz     ! Adress of working array for m1qn3
      integer(ip), dimension(:), allocatable  ::  iz        ! Adress of working array for m1qn3 (size 5) or lbfgsb3 (varied size)

      integer(ip)                             ::  m         ! m for LBFGSB-3.0
      integer(ip), dimension(:), allocatable  ::  nbd       ! nbd for LBFGSB-3.0
      real(rp)                                ::  factr     ! factr for LBFGSB-3.0
      character(len=60)                       ::  task      ! task for LBFGSB-3.0
      character(len=60)                       ::  csave     ! csave for LBFGSB-3.0
      logical, dimension(29)                  ::  lsave     ! lsave for LBFGSB-3.0
      real(rp), dimension(29)                 ::  dsave     ! dsave for LBFGSB-3.0
      integer(ip), dimension(44)              ::  isave     ! isave for LBFGSB-3.0

      logical   ::  normalize_gradJ_and_J
      real(rp)  ::  norm_gradJ

      real(rp)  ::  cost_ini , norm_grad_cost_ini(100) , norm_grad_cost(100)

!~       real(rp), dimension(:), allocatable :: control_temp

      !================================================================================================================!
      !  Initialization
      !================================================================================================================!

      call system('mkdir -p min')
      call dealloc_back_vars()
      call alloc_back_vars(dof0_back, dof_back, mesh ) ! allocate _back variables (both dof0_back and other control vector components, to be call only on time

      call write_control( dof0 , mesh )

      if (size(control) < 20_ip) then
        write(*,*) "Current control vector is ", control
      else
        write(*,*) "Current control vector is too large to show"
      endif

#ifndef USE_M1QN3
#ifndef USE_LBFGSB3
      write (*,'(A)') 'Wrong MINMETHOD in Makefile (must be 1 or 2)'
      stop
#endif
#endif

#ifdef USE_M1QN3
    Write(*,*) 'Calling M1QN3'
      !================================================================================================================!
      !  m1qn3 arguments
      !================================================================================================================!

      n          =   size( control )
      dxmin      =   1.d-5
      epsg       =   eps_min
      normtype   =   'two'
      impres     =   5
      io         =   40
      imode(1)   =   0
      imode(2)   =   0
      imode(3)   =   1
      ndz        =   4 * n + 100 * ( 2 * n + 1 )
      reverse    =   1
      indic      =   4
      max_ite_line_search = 25

      allocate( dz( ndz ) ) ; dz(:) = 0._rp
      allocate( iz(  5  ) )

      nsim  = 100 !This should be accessible via Python
      write(*,*) "---------------------------------------------------------------------------------"
      write(*,*) "WARNING : Maximum number of inverse model interations (including in-line iterations) is hardcoded to 100 in m_minimization!"
      write(*,*) "---------------------------------------------------------------------------------"

      if ( restart_min == 0 ) then

         niter  =  nsim

      else

         niter  =  restart_min

      end if

      !================================================================================================================!
      !  Initialization / Restarting
      !================================================================================================================!

      inquire( file = 'min/min_cost.txt' , exist = file_exist(1) )

      ! call mpi_wait_all

!~       if ( .not. file_exist(1) ) then

         ite_min         = 0
         ite_line_search = 0

         open(30,file='min/min_cost.txt'    ,status='replace',form = 'formatted')
         open(40,file='min/m1qn3_output.txt',status='replace',form = 'formatted')

         write(30,'(A)') '# ite cost norm_grad_cost'

!~       else

!~          open(30,file='min/min_cost.txt'    ,status='old',position = 'append',form='formatted')
!~          open(40,file='min/m1qn3_output.txt',status='old',position = 'append',form='formatted')

!~          allocate( control_back( size(control) ) )

!~          call read_restart_m1qn3

!~          imode(2)   =   1
!~          epsg       =   eps_min / epsg
!~          reverse    =   1

!~       end if

      !================================================================================================================!
      !  Minimization ... M1QN3 loop in reverse mode
      !================================================================================================================!


      do while( reverse >= 0 .and. indic /= 0 )
!             write(*,*) "-----------------------------------------------------------------------------"
!             write(*,*) "New loop ", "indic=",indic, "reverse =", reverse, "ite_min", ite_min, "ite_line_search", ite_line_search
!             write(*,*) "-----------------------------------------------------------------------------"
         !=============================================================================================================!
         !  Case indic = 4 -> M1QN3 needs new values of the cost and its gradient to compute
         !  either the step descent (line-search) or the new iterate
         !=============================================================================================================!

         if ( indic == 4 ) then

            cost_back = one

            verbose = -1
            call adjoint_model( mesh , dof0 , dof )
            ite_line_search = ite_line_search + 1

         end if

!~          !=============================================================================================================!
!~          !  As M1QN3 is sequential, only master thread perform the minimization
!~          !  The control vector is managed in the module m_adjoint
!~          !=============================================================================================================!

          if ( proc == 0 ) then

            !==========================================================================================================!
            !  Store the initial cost and its gradient norm
            !==========================================================================================================!

            if ( ite_min == 0 .and. indic == 4 ) then

               cost_ini  =  cost
write(*,*) "cost_ini = ", cost_ini
! write(*,*) "control_back = ", control_back
               i = 1

               do k = 1,nb_vars_in_control
! write(*,*) k, i, i - 1 + dim_vars_in_control(k)
                  norm_grad_cost_ini(k) = sqrt( sum( control_back( i : i - 1 + dim_vars_in_control(k) ) * &
                                                     control_back( i : i - 1 + dim_vars_in_control(k) ) ) )
                  i = i + dim_vars_in_control(k)

               end do

            end if

            !==========================================================================================================!
            !  Normalization of the control vector
            !==========================================================================================================!

            i = 1

            do k = 1,nb_vars_in_control
! write(*,*) "control( i : i - 1 + dim_vars_in_control(k) )", control( i : i - 1 + dim_vars_in_control(k) )
               control( i : i - 1 + dim_vars_in_control(k) ) = &
               control( i : i - 1 + dim_vars_in_control(k) ) * norm_grad_cost_ini(k) / cost_ini
! write(*,*) "control( i : i - 1 + dim_vars_in_control(k) )", control( i : i - 1 + dim_vars_in_control(k) )
               i = i + dim_vars_in_control(k)

            end do

            !==========================================================================================================!
            !  Normalization of the cost, its gradient and the control vector
            !==========================================================================================================!

            if ( indic == 4 ) then

               cost = cost / cost_ini

               i = 1

               do k = 1,nb_vars_in_control

                  control_back( i : i - 1 + dim_vars_in_control(k) ) = &
                  control_back( i : i - 1 + dim_vars_in_control(k) ) / control_back(k)

                  norm_grad_cost(k) = sqrt( sum( control_back( i : i - 1 + dim_vars_in_control(k) ) * &
                                                 control_back( i : i - 1 + dim_vars_in_control(k) ) ) )
! write(*,*) "norm_grad_cost(k)", norm_grad_cost(k)
! write(*,*) "control_back( i : i - 1 + dim_vars_in_control(k) )", control_back( i : i - 1 + dim_vars_in_control(k) )
                  i = i + dim_vars_in_control(k)

               end do

               norm_gradJ  =  sqrt( sum( control_back(:) * control_back(:) ) )
               write(6,'("ite ",I3," , J =",ES13.6," , |grad J| =",ES13.6, " , time =",ES13.6)') &
               ite_min , cost , norm_gradJ , time(1)


            end if

            !==========================================================================================================!
            !  dJ_expect is df1 in the M1QN3 documentation
            !==========================================================================================================!

            dJ_expect = 0.5_rp * cost

            !==========================================================================================================!
            !  Call of M1QN3
!           !==========================================================================================================!
!~             write(*,*) "-----------------------------------------------------------------------------"
!~             write(*,*) "control vector size:            "
!~             write(*,*) n
!~             write(*,*) "   control         "
!~             write(*,*) control(:)
!~             write(*,*) "control back"
!~             write(*,*) control_back(:)
!~             write(*,*) " input mode            "
!~             write(*,*) imode
!~             write(*,*) " output mode            "
!~             write(*,*) omode
!~             write(*,*) " cost            "
!~             write(*,*) cost
!~             write(*,*) "indic"
!~             write(*,*) indic
!~             write(*,*) "reverse"
!~             write(*,*) reverse
!~             write(*,*) "ndz"
!~             write(*,*) ndz
!~             write(*,*) "iz(:)"
!~             write(*,*) iz(:)
!~             write(*,*) "nsim(:)"
!~             write(*,*) nsim
!~             write(*,*) "io(:)"
!~             write(*,*) io
!~             write(*,*) "normtype(:)"
!~             write(*,*) normtype
!~             write(*,*) "epsg(:)"
!~             write(*,*) epsg
!~             write(*,*) "dJ_expect(:)"
!~             write(*,*) dJ_expect
!~             write(*,*) "-----------------------------------------------------------------------------"
!~             control_temp(:) = control(:)
            call m1qn3( simul_rc     , &
                        euclid       , &
                        ctonbe       , &
                        ctcabe       , &
                        n            , &
                        control      , &
                        cost         , &
                        control_back , &
                        dxmin        , &
                        dJ_expect    , &
                        epsg         , &
                        normtype     , &
                        impres       , &
                        io           , &
                        imode        , &
                        omode        , &
                        niter        , &
                        nsim         , &
                        iz           , &
                        dz           , &
                        ndz          , &
                        reverse      , &
                        indic        , &
                        izs          , &
                        rzs          , &
                        dzs )
!~             write(*,*) ">>>> DONE CALL M1QN3  "
!~             write(*,*) "control vector size:            "
!~             write(*,*) n
!~             write(*,*) "   control         "
!~             write(*,*) control(:)
!~             write(*,*) "control back"
!~             write(*,*) control_back(:)
!~             write(*,*) " input mode            "
!~             write(*,*) imode
!~             write(*,*) " output mode            "
!~             write(*,*) omode
!~             write(*,*) " cost            "
!~             write(*,*) cost
!~             write(*,*) "indic"
!~             write(*,*) indic
!~             write(*,*) "reverse"
!~             write(*,*) reverse
!~             write(*,*) "ndz"
!~             write(*,*) ndz
!~             write(*,*) "iz(:)"
!~             write(*,*) iz(:)
!~             write(*,*) "nsim(:)"
!~             write(*,*) nsim
!~             write(*,*) "io(:)"
!~             write(*,*) io
!~             write(*,*) "normtype(:)"
!~             write(*,*) normtype
!~             write(*,*) "epsg(:)"
!~             write(*,*) epsg
!~             write(*,*) "dJ_expect(:)"
!~             write(*,*) dJ_expect
!~             if (all(control(:) == control_temp(:) )) print *, "control is not changed here l344"
            !==========================================================================================================!
            !  Back normalization of the control vector in order to output it
            !==========================================================================================================!

            i = 1

            do k = 1,nb_vars_in_control

               control( i : i - 1 + dim_vars_in_control(k) ) = &
               control( i : i - 1 + dim_vars_in_control(k) ) / norm_grad_cost_ini(k) * cost_ini

               i = i + dim_vars_in_control(k)

            end do

         end if

         call mpi_wait_all

         call mpi_bcast_i( reverse , 0 )
         call mpi_bcast_i( indic   , 0 )
         call mpi_bcast_i( omode   , 0 )

!          call write_restart_m1qn3

         !=============================================================================================================!
         !  Case indic = 1 -> M1QN3 has finished an iterate ( imode(3) = 1 )
         !=============================================================================================================!

         if      ( indic == 1 ) then

!             write(*,*) "Case indic = 1 -> M1QN3 has finished an iterate ( imode(3) = 1 )"
            call output_control( dof0 , mesh )

            if ( proc == 0 ) write(30,'(I4,102ES15.7)') ite_min , cost , norm_grad_cost(1:nb_vars_in_control)

            ite_min = ite_min + 1

            ite_line_search = 0

         else if ( indic == 4 .and. ite_line_search > max_ite_line_search) then

            write(*,*) " indic == 4 .and. ite_line_search > max_ite_line_search"
            call Stopping_Program_Sub( 'Stopping Minimization, too much iterations in line-search' )

            write(*,*) "Warning : dirty stop calibration lilian in m_minimization"
            reverse = -1 ! dirty stop calibration lilian

         end if

      end do

      close(30)
      close(40)
#endif

#ifdef USE_LBFGSB3
      !================================================================================================================!
      !  Parameters for LBFGSB 3.0
      !================================================================================================================!
      n          =   size( control )
      m          =   5
      factr      =   eps_min / epsilon(one)
      impres     =   -1
      epsg       =   eps_min
      allocate( iz( 3 * n ) )
      allocate( dz( 2 * m * n + 5 * n + 11 * m * m + 8 * m ) )
      allocate( nbd( n ) )
!       nbd(1) = 2
!       nbd(2) = 2
      ite_line_search = 0
!       allocate( dx( n ) )
!       allocate( gradJ( n ) )

      !================================================================================================================!
      !  Compute bounds
      !================================================================================================================!
      call write_control_bounds( dof0 , mesh )
      do k = 1, size(control)
         if (control_ubound(k) - control_lbound(k) < 1e-9) then
            nbd(k) = 0
         else
            nbd(k) = 2
         end if
      end do

      !================================================================================================================!
      !  Rerun
      !================================================================================================================!
      inquire( file = 'min/rerun.txt' , exist = file_exist(1) )

      if ( ite_min == 0 .and. file_exist(1) ) then

         call Print_Screen('rerun_min')

         open(60, file='min/rerun.txt',status='old')
         read(60, *) control
         close(60)

      end if

      !================================================================================================================!
      !  Initialize LBFGSB
      !================================================================================================================!
      task       = 'START'
      normalize_gradJ_and_J = .True. !Modif Leo
      call setulb ( n, &
                    m, &
                    control, &            ! x
                    control_lbound, &     ! l
                    control_ubound,&      ! u
                    nbd, &
                    cost, &               ! f
                    control_back, &       ! g
                    factr, & !eps_min, &            ! factr
                    epsg, &               ! pgtol
                    dz, &                 ! wa
                    iz, &                 ! iwa
                    task, &
                    impres, &             !iprint
                    csave, &
                    lsave, &
                    isave, &
                    dsave )




      !================================================================================================================!
      !  Initialization (TODO : Restarting ?)
      !================================================================================================================!
      open(30,file='min/min_cost.txt'    ,status='replace',form = 'formatted')
!       open(40,file='min/m1qn2_output.txt',status='replace',form = 'formatted')
      write(30,'(A)') '# ite cost norm_grad_cost'



      do while(task(1:2) == "FG" .or. task == "NEW_X" .or. &
               task == "START")

         !=============================================================================================================!
         !  task = FG -> LBFGSB-3.0 needs new values of the cost and its gradient
         !=============================================================================================================!
         if (task(1:2) == "FG") then

            if (isave(36) > 25) then
              write(*,*) 'Stopping Minimization, too much iterations in line-search'
              task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
            end if

            cost_back = one

            verbose = -1

            call adjoint_model( mesh , dof0 , dof )


            ! make plots !
            if ( c_gr4params == 1 ) then
            open(1234,file="min/gr4_params_states_current",form='formatted')
            do i=1, size(bc%gr4)
               write(1234,*) bc%gr4(i)%params, bc%gr4(i)%state
            enddo
            close(1234)
            endif

            open(1234,file="min/control_back_current",form='formatted')
            do i=1, size(control)
               write(1234,*) control_back(i)
            enddo
            close(1234)

         end if


         !=============================================================================================================!
         !  Store the initial cost and its gradient norm
         !=============================================================================================================!

         if ( isave(30) == 0 .and. isave(36) == 0 .and. task(1:2) == "FG" ) then

            cost_ini  =  cost

            i = 1

            do k = 1, nb_vars_in_control

              norm_grad_cost_ini(k) = sqrt( sum( control_back( i : i - 1 + dim_vars_in_control(k) ) * &
                                                  control_back( i : i - 1 + dim_vars_in_control(k) ) ) )

              i = i + dim_vars_in_control(k)

            end do

         end if

         call write_restart_lbfgsb3


         !=============================================================================================================!
         !  Normalization of the control vector and bounds
         !=============================================================================================================!

         if (task == "FG_START") then
            call output_control( dof0 , mesh )
         end if


         if ( normalize_gradJ_and_J ) then
            i = 1
            do k = 1,nb_vars_in_control
                control( i : i - 1 + dim_vars_in_control(k) ) = &
                control( i : i - 1 + dim_vars_in_control(k) ) * norm_grad_cost_ini(k) / cost_ini
                control_lbound( i : i - 1 + dim_vars_in_control(k) ) = &
                control_lbound( i : i - 1 + dim_vars_in_control(k) ) * norm_grad_cost_ini(k) / cost_ini
                control_ubound( i : i - 1 + dim_vars_in_control(k) ) = &
                control_ubound( i : i - 1 + dim_vars_in_control(k) ) * norm_grad_cost_ini(k) / cost_ini
                i = i + dim_vars_in_control(k)
!write(*,*) 'lbfgsb3 - normalized control and bounds', control, control_lbound, control_ubound
            end do
         end if

         !=============================================================================================================!
         !  Normalization of the cost, its gradient and the control vector
         !=============================================================================================================!

         if ( task(1:2) == "FG" ) then

            cost = cost / cost_ini

            i = 1

            do k = 1,nb_vars_in_control

              if ( normalize_gradJ_and_J ) then
                 control_back( i : i - 1 + dim_vars_in_control(k) ) = &
                 control_back( i : i - 1 + dim_vars_in_control(k) ) / norm_grad_cost_ini(k)
              end if

              norm_grad_cost(k) = sqrt( sum( control_back( i : i - 1 + dim_vars_in_control(k) ) * &
                                              control_back( i : i - 1 + dim_vars_in_control(k) ) ) )

              i = i + dim_vars_in_control(k)

            end do

            norm_gradJ  =  sqrt( sum( control_back(:) * control_back(:) ) )

            write(6,'("ite ",I3," , J =",ES13.6," , |grad J| =",ES13.6, " , time =",ES13.6)') &
            isave(30) , cost , norm_gradJ , time(1)

            if (task == "FG_START") then
            write(30,'(I4,102ES15.7)') 0, cost , norm_grad_cost(1:nb_vars_in_control)
            end if

         end if

         !=============================================================================================================!
         !  Call of LBFGSB
         !=============================================================================================================!
         call setulb ( n, &
                       m, &
                       control, &            ! x
                       control_lbound, &     ! l
                       control_ubound,&      ! u
                       nbd, &
                       cost, &               ! f
                       control_back, &       ! g
                       factr, & !eps_min, &            ! factr
                       epsg, &               ! pgtol
                       dz, &                 ! wa
                       iz, &                 ! iwa
                       task, &
                       impres, &             !iprint
                       csave, &
                       lsave, &
                       isave, &
                       dsave )

         print *, "TASK:", trim(task), cost

         !=============================================================================================================!
         !  Back normalization of the control vector in order to output it
         !=============================================================================================================!

         if ( normalize_gradJ_and_J ) then
            i = 1

            do k = 1,nb_vars_in_control

                control( i : i - 1 + dim_vars_in_control(k) ) = &
                control( i : i - 1 + dim_vars_in_control(k) ) / norm_grad_cost_ini(k) * cost_ini
                control_lbound( i : i - 1 + dim_vars_in_control(k) ) = &
                control_lbound( i : i - 1 + dim_vars_in_control(k) ) / norm_grad_cost_ini(k) * cost_ini
                control_ubound( i : i - 1 + dim_vars_in_control(k) ) = &
                control_ubound( i : i - 1 + dim_vars_in_control(k) ) / norm_grad_cost_ini(k) * cost_ini

                i = i + dim_vars_in_control(k)

            end do
         end if

         !=============================================================================================================!
         !  Case task = NEW_X -> LBFGSB-3.0 has finished an iterate
         !=============================================================================================================!
         if (task == "NEW_X") then

            ite_min = ite_min + 1
            ite_line_search = 0

            call output_control( dof0 , mesh )
            write(30,'(I4,102ES15.7)') ite_min , cost , norm_grad_cost(1:nb_vars_in_control)

         end if

      end do
      close(30)

      if (task(1:5) == "ERROR") then
        omode = 2
         write(*,*) 'LBFGSB-3: '//trim(task)
      else if (task(1:4) == "ABNO") then
        omode = 6
      else
        omode = 1
      end if
      print *, 'FINAL_TASK=', trim(task)
#endif

      call mpi_wait_all

      deallocate( dz )

      if ( omode == 1 .or. &
           omode == 3 .or. &
           omode == 4 .or. &
           omode == 5 .or. &
           omode == 6 ) then

         verbose = 1

			write(*,*) "Direct simulation within minimization script was using last control vector and not "
			write(*,*) "optimal control vector, until this is fixed, the direct simulation associated "
			write(*,*) "to the optimal control vector must be done manualy "
!~          call model_direct( mesh , dof0 , dof )
!~          call Print_Screen( 'cost' , cost )

      end if

      if      ( omode == 1 ) then
         call Stopping_Program_Sub( 'Successful Minimization, the test on the gradient norm (eps_min) is satisfied' )

      else if ( omode == 2 ) then
         call Stopping_Program_Sub( 'Stopping Minimization, input parmeter not well initialized (omode=2 in the M1QN3 documentation)' )

      else if ( omode == 3 ) then
         call Stopping_Program_Sub( 'Stopping Minimization, the line-search is blocked (omode=3 in the M1QN3 documentation)' )

      else if ( omode == 4 ) then
         call Stopping_Program_Sub( 'Stopping Minimization, maximal number of iterations reached (omode=4 in the M1QN3 documentation)' )

      else if ( omode == 5 ) then
         call Stopping_Program_Sub( 'Stopping Minimization, maximal number of simulations reached (omode=5 in the M1QN3 documentation)' )

      else if ( omode == 6 ) then
         call Stopping_Program_Sub( 'Stopping Minimization, stop on dxmin during the line-search (omode=6 in the M1QN3 documentation)' )

      else

         call Stopping_Program_Sub( 'Stopping Minimization, unknown reason !?)' )

      end if


   CONTAINS


      SUBROUTINE write_restart_m1qn3

         implicit none

         if ( proc == 0 ) then

            open(50,file='min/restart_min.bin',status='replace',form='unformatted')

            write(50) ite_min , ite_line_search , reverse , indic , epsg , iz , dz, &
                      control , cost , control_back , cost_ini , norm_grad_cost_ini

            close(50)

         end if

      END SUBROUTINE write_restart_m1qn3

      SUBROUTINE write_restart_lbfgsb3

         implicit none

         if ( proc == 0 ) then

            open(50,file='min/restart_min.bin',status='replace',form='unformatted')

            write(50) ite_min , ite_line_search , reverse , indic , epsg , iz , dz, &
                      control , cost , control_back , cost_ini , norm_grad_cost_ini

            close(50)

         end if

      END SUBROUTINE write_restart_lbfgsb3


      SUBROUTINE read_restart_m1qn3

         implicit none

         open(50,file='min/restart_min.bin',status='old',form='unformatted')

         read(50) ite_min , ite_line_search , reverse , indic , epsg , iz , dz, &
                  control , cost , control_back , cost_ini , norm_grad_cost_ini

         close(50)

      END SUBROUTINE read_restart_m1qn3


   END SUBROUTINE minimize_cost


END MODULE m_minimization
