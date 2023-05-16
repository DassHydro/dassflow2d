MODULE m_minimization
   USE m_adjoint
   USE m_model
   implicit none
   integer(ip) :: ite_line_search
CONTAINS
   SUBROUTINE minimize_cost( mesh , dof0 , dof )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in ) :: mesh
      type( unk ), intent(inout) :: dof0
      type( unk ), intent(inout) :: dof
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      external simul_rc , euclid , ctonbe , ctcabe
      integer(ip) :: n ! Control vector size
      real(rp) :: dxmin ! Resolution for x in Linf norm
      real(rp) :: dJ_expect ! Estimation of the expected decrease in cost during the first iteration
      real(rp) :: epsg ! Stopping criterion
      character*3 :: normtype ! specifies the norm used to test optimality
      integer(ip) :: impres ! Print option for m1qn3
      integer(ip) :: io ! File output label
      integer(ip) :: imode(3) ! Input mode
      integer(ip) :: omode ! Output mode
      integer(ip) :: niter ! Maximum number of minimization iterations
      integer(ip) :: nsim ! Maximum number of simulator calls
      integer(ip) :: ndz ! Dimension of dz
      integer(ip) :: reverse ! Specifies direct or reverse mode
      integer(ip) :: indic ! m1qn3 indicates it needs cost/gradient
      integer(ip) :: izs(1) !
      real(rp) :: rzs(1) ! Working areas not needed in reverse mode
      real(rp) :: dzs(1) !
      integer(ip) :: max_ite_line_search
      real(rp), dimension(:), allocatable :: dz ! Adress of working array for m1qn3
      integer(ip), dimension(:), allocatable :: iz ! Adress of working array for m1qn3 (size 5) or lbfgsb3 (varied size)
      integer(ip) :: m ! m for LBFGSB-3.0
      integer(ip), dimension(:), allocatable :: nbd ! nbd for LBFGSB-3.0
      real(rp) :: factr ! factr for LBFGSB-3.0
      character(len=60) :: task ! task for LBFGSB-3.0
      character(len=60) :: csave ! csave for LBFGSB-3.0
      logical, dimension(29) :: lsave ! lsave for LBFGSB-3.0
      real(rp), dimension(29) :: dsave ! dsave for LBFGSB-3.0
      integer(ip), dimension(44) :: isave ! isave for LBFGSB-3.0
      logical :: normalize_gradJ_and_J
      real(rp) :: norm_gradJ
      real(rp) :: cost_ini , norm_grad_cost_ini(100) , norm_grad_cost(100)
      !================================================================================================================!
      ! Initialization
      !================================================================================================================!
      call system('mkdir -p min')
      call dealloc_back_vars()
      call alloc_back_vars(dof0_back, dof_back, mesh ) ! allocate _back variables (both dof0_back and other control vector components, to be call only on time
      call write_control( dof0 , mesh )
write(*,*) "control in", control
    Write(*,*) 'Calling M1QN3'
      !================================================================================================================!
      ! m1qn3 arguments
      !================================================================================================================!
      n = size( control )
      dxmin = 1.d-5
      epsg = eps_min
      normtype = 'two'
      impres = 5
      io = 40
      imode(1) = 0
      imode(2) = 0
      imode(3) = 1
      ndz = 4 * n + 100 * ( 2 * n + 1 )
      reverse = 1
      indic = 4
      max_ite_line_search = 50
      allocate( dz( ndz ) ) ; dz(:) = 0._rp
      allocate( iz( 5 ) )
      nsim = 100
      if ( restart_min == 0 ) then
         niter = nsim
      else
         niter = restart_min
      end if
      !================================================================================================================!
      ! Initialization / Restarting
      !================================================================================================================!
      inquire( file = 'min/min_cost.txt' , exist = file_exist(1) )
         ite_min = 0
         ite_line_search = 0
         open(30,file='min/min_cost.txt' ,status='replace',form = 'formatted')
         open(40,file='min/m1qn3_output.txt',status='replace',form = 'formatted')
         write(30,'(A)') '# ite cost norm_grad_cost'
      !================================================================================================================!
      ! Minimization ... M1QN3 loop in reverse mode
      !================================================================================================================!
      do while( reverse >= 0 .and. indic /= 0 )
            write(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
            write(*,*) "New loop", "indic=",indic, "reverse=", reverse, "ite_min", ite_min, "ite_line_search", ite_line_search
            write(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
         !=============================================================================================================!
         ! Case indic = 4 -> M1QN3 needs new values of the cost and its gradient to compute
         ! either the step descent (line-search) or the new iterate
         !=============================================================================================================!
         if ( indic == 4 ) then
            cost_back = one
            verbose = -1
            call adjoint_model( mesh , dof0 , dof )
            ite_line_search = ite_line_search + 1
         end if
          if ( proc == 0 ) then
            !==========================================================================================================!
            ! Store the initial cost and its gradient norm
            !==========================================================================================================!
            if ( ite_min == 0 .and. indic == 4 ) then
               cost_ini = cost
               i = 1
               do k = 1,nb_vars_in_control
                  norm_grad_cost_ini(k) = sqrt( sum( control_back( i : i - 1 + dim_vars_in_control(k) ) * &
                                                     control_back( i : i - 1 + dim_vars_in_control(k) ) ) )
                  i = i + dim_vars_in_control(k)
               end do
            end if
            !==========================================================================================================!
            ! Normalization of the control vector
            !==========================================================================================================!
            i = 1
            do k = 1,nb_vars_in_control
               control( i : i - 1 + dim_vars_in_control(k) ) = &
               control( i : i - 1 + dim_vars_in_control(k) ) * norm_grad_cost_ini(k) / cost_ini
               i = i + dim_vars_in_control(k)
            end do
            !==========================================================================================================!
            ! Normalization of the cost, its gradient and the control vector
            !==========================================================================================================!
            if ( indic == 4 ) then
               cost = cost / cost_ini
               i = 1
               do k = 1,nb_vars_in_control
                  control_back( i : i - 1 + dim_vars_in_control(k) ) = &
                  control_back( i : i - 1 + dim_vars_in_control(k) ) / norm_grad_cost_ini(k)
                  norm_grad_cost(k) = sqrt( sum( control_back( i : i - 1 + dim_vars_in_control(k) ) * &
                                                 control_back( i : i - 1 + dim_vars_in_control(k) ) ) )
                  i = i + dim_vars_in_control(k)
               end do
               norm_gradJ = sqrt( sum( control_back(:) * control_back(:) ) )
               write(6,'("ite ",I3," , J =",ES13.6," , |grad J| =",ES13.6, " , time =",ES13.6)') &
               ite_min , cost , norm_gradJ , time(1)
            end if
            !==========================================================================================================!
            ! dJ_expect is df1 in the M1QN3 documentation
            !==========================================================================================================!
            dJ_expect = 0.5_rp * cost
            !==========================================================================================================!
            ! Call of M1QN3
Write(*,*) "control before   call m1qn3 =", control
            call m1qn3( simul_rc , &
                        euclid , &
                        ctonbe , &
                        ctcabe , &
                        n , &
                        control , &
                        cost , &
                        control_back , &
                        dxmin , &
                        dJ_expect , &
                        epsg , &
                        normtype , &
                        impres , &
                        io , &
                        imode , &
                        omode , &
                        niter , &
                        nsim , &
                        iz , &
                        dz , &
                        ndz , &
                        reverse , &
                        indic , &
                        izs , &
                        rzs , &
                        dzs )
Write(*,*) "control  after lbfgsb =", control
            !==========================================================================================================!
            ! Back normalization of the control vector in order to output it
            !==========================================================================================================!
            i = 1
            do k = 1,nb_vars_in_control
               control( i : i - 1 + dim_vars_in_control(k) ) = &
               control( i : i - 1 + dim_vars_in_control(k) ) / norm_grad_cost_ini(k) * cost_ini
               i = i + dim_vars_in_control(k)
            end do
         end if
         call write_restart_m1qn3
         !=============================================================================================================!
         ! Case indic = 1 -> M1QN3 has finished an iterate ( imode(3) = 1 )
         !=============================================================================================================!
         if ( indic == 1 ) then
            write(*,*) "Case indic = 1 -> M1QN3 has finished an iterate ( imode(3) = 1 )"
            write(*,*) "Call output_control"
            call output_control( dof0 , mesh )
            if ( proc == 0 ) write(30,'(I4,102ES15.7)') ite_min , cost , norm_grad_cost(1:nb_vars_in_control)
            ite_min = ite_min + 1
            ite_line_search = 0
         else if ( indic == 4 .and. ite_line_search > max_ite_line_search) then
            write(*,*) " indic == 4 .and. ite_line_search > max_ite_line_search"
        ! call Stopping_Program_Sub( 'Stopping Minimization, too much iterations in line-search' )
            reverse = -1 ! dirty stop calibration lilian
         end if
      end do
      close(30)
      close(40)
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
      end if
      if ( omode == 1 ) then
      ! call Stopping_Program_Sub( 'Successful Minimization, the test on the gradient norm (eps_min) is satisfied' )
      else if ( omode == 2 ) then
    ! call Stopping_Program_Sub( 'Stopping Minimization, input parmeter not well initialized (omode=2 in the M1QN3 documentation)' )
      else if ( omode == 3 ) then
    ! call Stopping_Program_Sub( 'Stopping Minimization, the line-search is blocked (omode=3 in the M1QN3 documentation)' )
      else if ( omode == 4 ) then
  ! call Stopping_Program_Sub( 'Stopping Minimization, maximal number of iterations reached (omode=4 in the M1QN3 documentation)' )
      else if ( omode == 5 ) then
  ! call Stopping_Program_Sub( 'Stopping Minimization, maximal number of simulations reached (omode=5 in the M1QN3 documentation)' )
      else if ( omode == 6 ) then
    ! call Stopping_Program_Sub( 'Stopping Minimization, stop on dxmin during the line-search (omode=6 in the M1QN3 documentation)' )
      else
  ! call Stopping_Program_Sub( 'Stopping Minimization, unknown reason !?)' )
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
