MODULE m_time_screen
   USE m_common
   USE m_mesh
   implicit none
   integer, parameter :: max_time_index = 100
   real(rp), dimension( max_time_index ) :: t_ini
   real(rp), dimension( max_time_index ) :: t_end
   real(rp), dimension( max_time_index ) :: time
CONTAINS
   SUBROUTINE Time_Init( i )
      implicit none
      intrinsic cpu_time
      integer(ip), intent(in) :: i
      time( i ) = 0._rp
         call CPU_TIME( t_ini( i ) )
   END SUBROUTINE Time_Init
   SUBROUTINE Time_End( i )
      implicit none
      intrinsic cpu_time
      integer(ip), intent(in) :: i
         call CPU_TIME( t_end( i ) )
      time( i ) = t_end( i ) - t_ini( i )
   END SUBROUTINE Time_End
   SUBROUTINE Time_Init_Part( i )
      implicit none
      intrinsic cpu_time
      integer(ip), intent(in) :: i
         call CPU_TIME( t_ini( i ) )
   END SUBROUTINE Time_Init_Part
   SUBROUTINE Time_End_Part( i )
      implicit none
      intrinsic cpu_time
      integer(ip), intent(in) :: i
         call CPU_TIME( t_end( i ) )
      time( i ) = time( i ) + t_end( i ) - t_ini( i )
   END SUBROUTINE Time_End_Part
   SUBROUTINE Time_Screen( mesh )
      implicit none
      TYPE( msh ), intent(in) :: mesh
      integer(ip) :: mesh_total_size
      mesh_total_size = mesh%nc
         write(6,'(A)' )
         write(6,'(A)' ) '********************************************************************************'
         write(6,'(A,F15.2)' ) ' Time of simulation   =   ' , time(1)
         write(6,'(A)' ) '********************************************************************************'
         write(6,'(A,F15.2,A)') ' Performance          =   ' , time(1) * np / mesh_total_size / nt * 1.d6 , &
                                    ' microsecond / dx / dt / proc'
         write(6,'(A)' ) '********************************************************************************'
         if ( spatial_scheme(1:5) == 'muscl' ) then
            write(6,'(A,F5.2)' ) ' MUSCL reconstruction in percent  =  ' , 100._rp * time(2) / time(1)
            write(6,'(A,F5.2)' ) ' MUSCL limitation     in percent  =  ' , 100._rp * time(3) / time(1)
            write(6,'(A,F5.2)' ) ' Fluxes computation   in percent  =  ' , 100._rp * time(4) / time(1)
            write(6,'(A,F5.2)' ) ' Variables update     in percent  =  ' , 100._rp * time(5) / time(1)
            write(6,'(A)' ) '********************************************************************************'
         end if
   END SUBROUTINE Time_Screen
 !**********************************************************************************************************************!
 !**********************************************************************************************************************!
 !
 ! Stopping Program in case of bad parameter or numeric convergence
 !
 !**********************************************************************************************************************!
 !**********************************************************************************************************************!
 !> \brief PRINT WHY Stopping Program in case of bad parameter or numeric convergence
 !> \details print error message for different error cases
 !! - unknow_mesh : say unknown mesh
 !! - default : say stop dassflow run
 !! f90wrap_abort is called here to message wrapping error
    SUBROUTINE Stopping_Program( screen_case )
       implicit none
       character(len=*), intent(in) :: screen_case
          select case( screen_case )
             case( 'unknow_mesh' )
                write(6,*)
                write(6,'(A80)') '================================================================================'
                write(6,'(A80)') '*  Unknow mesh type in input.txt'
                write(6,'(A80)') '================================================================================'
             case default
          end select
          write(6,*)
          write(6,'(A)') '================================================================================'
          write(6,'(A)') '*'
          write(6,'(A)') '*  STOPPING DASSFLOW RUN'
          write(6,'(A)') '*'
          write(6,'(A)') '================================================================================'
   call f90wrap_abort(screen_case)
    END SUBROUTINE Stopping_Program
    !**********************************************************************************************************************!
    !**********************************************************************************************************************!
    !
    ! Subroutine dedicated to Screen Output Control
    !
    !**********************************************************************************************************************!
    !**********************************************************************************************************************!
    !> \brief Subroutine dedicated to Screen Output Control)
    !> \details following tags are treated
    !! following options are treated
    !! - proc == 0
    !! - verbose >= 0
    !! - test_dt_just_after( 0.01_rp * ts )
    !! following tags are treated
    !! following case can be proposed
    !! - start_mesh
    !! - end_mesh
    !! - start_direct
    !! - end_direct
    !! - start_testadj
    !! - end_testadj
    !! - start_grad_cost
    !! - end_grad_cost
    !! - start_minimize
    !! - end_minimize
    !! - result (if verbose >=0)
    !! - dt ( if verbose >=0 && test_dt_just_after( 0.01_rp * ts ))
    !! - number_limits
    !! - norms_h
    !! - norms_u
    !! - norms_q
    !! - tangent_gradient_test
    !! - backward_gradient_test
    !! - backward_scalar_product_test
    !! - cost
    !! - cost_diff
    !! - cost_back
       SUBROUTINE Print_Screen( screen_case , var )
          implicit none
          intrinsic range , precision
          character(len=*), intent(in) :: screen_case
          real(rp), optional, intent(in) :: var
             select case( screen_case )
                case( 'start_mesh' )
                   write(6,'(A)')
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A80)') '*  Mesh Loading                                                                *'
                   write(6,'(A80)') '================================================================================'
                case( 'end_mesh' )
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A80)') '*  Mesh Loading is OK                                                          *'
                   write(6,'(A80)') '================================================================================'
                case( 'start_direct' )
                   write(6,'(A)')
                      write(6,'(A80)') '================================================================================'
                      write(6,'(A80)') '*  Running Shallow-Water Model in direct mode                                  *'
                      write(6,'(A80)') '================================================================================'
                case( 'end_direct' )
                      write(6,'(A80)') '================================================================================'
                      write(6,'(A80)') '*  End of run of the Shallow-Water Model in direct mode                        *'
                      write(6,'(A80)') '================================================================================'
                case( 'start_testadj' )
                   write(6,'(A)')
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A80)') '*  Testing the generated Adjoint to Shallow-Water Model                        *'
                   write(6,'(A80)') '================================================================================'
                case( 'end_testadj' )
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A80)') '*  End of test of the generated Adjoint to Shallow-Water Model                 *'
                   write(6,'(A80)') '================================================================================'
                case( 'start_grad_cost' )
                   write(6,'(A)')
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A80)') '*  Running the Shallow-Water Model Adjoint to calculate a cost-gradient        *'
                   write(6,'(A80)') '================================================================================'
                case( 'end_grad_cost' )
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A80)') '*  End of run of the Shallow-Water Model Adjoint                               *'
                   write(6,'(A80)') '================================================================================'
                case( 'start_minimize' )
                   write(6,'(A)')
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A80)') '*  Running the Shallow-Water loop to minimize the cost function                *'
                   write(6,'(A80)') '================================================================================'
                case( 'end_minimize' )
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A80)') '*  End of run of the minimization Shallow-Water Model loop                     *'
                   write(6,'(A80)') '================================================================================'
                case( 'result' )
                   if ( verbose >= 0 ) then
                      write(6,'(A80)') '================================================================================'
                      write(6,'(A80)') '*  Writing Result File                                                         *'
                      write(6,'(A80)') '================================================================================'
                   end if
                case( 'dt' )
                   if ( test_dt_just_after( 0.01_rp * ts ) .and. verbose >= 0 ) then
                      write(6,'("nt = ",I8," t = ",ES12.5," / ",ES12.5," ( ",F5.1," % ) , dt = ",ES13.6)') &
                      nt , tc , ts , 100._rp * tc / ts , dt
                   end if
                case( 'number_limits' )
                   write(6,'(A)')
                   write(6,'(A80)') '********************************************************************************'
                   write(6,'(A,I15  )') ' range     = ' , range(1._rp)
                   write(6,'(A,I15  )') ' precision = ' , precision(1._rp)
                   write(6,'(A,ES15.7)') ' zerom     = ' , zerom
                   write(6,'(A,ES15.7)') ' pinfm     = ' , pinfm
                   write(6,'(A,ES15.7)') ' minfm     = ' , minfm
                   write(6,'(A,ES15.7)') ' hugem     = ' , hugem
                   write(6,'(A,ES15.7)') ' tinym     = ' , tinym
                   write(6,'(A80)') '********************************************************************************'
                case( 'norms_h' )
                   write(6,'(A80)') '********************************************************************************'
                   write(6,'(A,ES15.7)') ' norm_h_inf  =  ' , norm_inf(1)
                   write(6,'(A,ES15.7)') ' norm_h_L1   =  ' , norm_L1 (1)
                   write(6,'(A,ES15.7)') ' norm_h_L2   =  ' , norm_L2 (1)
                   write(6,'(A80)') '********************************************************************************'
                case( 'norms_u' )
                   write(6,'(A80)') '********************************************************************************'
                   write(6,'(A,ES15.7)') ' norm_u_inf  =  ' , norm_inf(2)
                   write(6,'(A,ES15.7)') ' norm_u_L1   =  ' , norm_L1 (2)
                   write(6,'(A,ES15.7)') ' norm_u_L2   =  ' , norm_L2 (2)
                   write(6,'(A80)') '********************************************************************************'
                case( 'norms_q' )
                   write(6,'(A80)') '********************************************************************************'
                   write(6,'(A,ES15.7)') ' norm_q_inf  =  ' , norm_inf(3)
                   write(6,'(A,ES15.7)') ' norm_q_L1   =  ' , norm_L1 (3)
                   write(6,'(A,ES15.7)') ' norm_q_L2   =  ' , norm_L2 (3)
                   write(6,'(A80)') '********************************************************************************'
                case( 'tangent_gradient_test' )
                   write(6,'(A)')
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A80)') '*  Tangent Gradient Test                                                       *'
                   write(6,'(A80)') '================================================================================'
                case( 'backward_gradient_test' )
                   write(6,'(A)')
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A80)') '*  Backward Gradient Test                                                      *'
                   write(6,'(A80)') '================================================================================'
                case( 'backward_scalar_product_test' )
                   write(6,'(A)')
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A80)') '*  Backward Scalar Product Test                                                *'
                   write(6,'(A80)') '================================================================================'
                case( 'cost' )
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A,ES24.15)') ' cost function       =  ' , var
                   write(6,'(A80)') '================================================================================'
                case( 'cost_diff' )
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A,ES24.15)') ' diff cost function  =  ' , var
                   write(6,'(A80)') '================================================================================'
                case( 'cost_back' )
                   write(6,'(A80)') '================================================================================'
                   write(6,'(A,ES24.15)') ' back cost function  =  ' , var
                   write(6,'(A80)') '================================================================================'
             end select
       END SUBROUTINE Print_Screen
END MODULE m_time_screen
