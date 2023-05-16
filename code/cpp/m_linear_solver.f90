MODULE m_linear_solver
   USE m_common
   USE m_mesh
   USE m_linear_algebra
   USE m_model
   implicit none
CONTAINS
   SUBROUTINE Init_Linear_Solver(mesh)
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type(msh), intent(in) :: mesh
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      integer(ip) :: iL , iR
      !================================================================================================================!
      ! Begin Subroutine
      !================================================================================================================!
   END SUBROUTINE Init_Linear_Solver
   SUBROUTINE matrix_COO_to_CSR( mat )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( sys_lin_sparse ), intent(inout) :: mat
      !================================================================================================================!
      ! Local Variables
      !================================================================================================================!
      type( sys_lin_sparse ) :: copy
      integer(ip) :: elt , row , ct( mat%nz )
      !================================================================================================================!
      ! Begin Subroutine
      !================================================================================================================!
      allocate( copy%ia ( mat%nz ) ) ; copy%ia = mat%ia
      allocate( copy%ja ( mat%nz ) ) ; copy%ja = mat%ja
      call reallocate_i( mat%ia , mat%n + 1 )
      mat%ia(:) = 1
      do elt = 1,mat%nz
         row = copy%ia( elt )
         mat%ia( row + 1 ) = mat%ia( row + 1 ) + 1
      end do
      do elt = 1,mat%n
         mat%ia( elt + 1 ) = mat%ia( elt + 1 ) + mat%ia( elt ) - 1
      end do
      ct(:) = 0
      do elt = 1,mat%nz
         row = mat%ia( copy%ia( elt ) )
         mat%ja( row + ct( row ) ) = copy%ja( elt )
         mat%swap( row + ct( row ) ) = elt
         ct( row ) = ct( row ) + 1
      end do
      deallocate( copy%ia )
      deallocate( copy%ja )
   END SUBROUTINE matrix_COO_to_CSR
END MODULE m_linear_solver
