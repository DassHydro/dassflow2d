MODULE m_mpi_diff
   USE m_mpi
      USE m_model
   implicit none
CONTAINS
   SUBROUTINE com_var_i_diff( var , var_diff , mesh )
      USE m_mesh
      implicit none
      type(msh), intent(in) :: mesh
      integer(ip), dimension( mesh%nc + mesh%ncb ), intent(inout) :: var
      integer(ip), dimension( mesh%nc + mesh%ncb ), intent(inout) :: var_diff
   END SUBROUTINE com_var_i_diff
   SUBROUTINE com_var_r_diff( var , var_diff , mesh )
      USE m_mesh
      implicit none
      type(msh), intent(in) :: mesh
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(inout) :: var
      real(rp), dimension( mesh%nc + mesh%ncb ), intent(inout) :: var_diff
   END SUBROUTINE com_var_r_diff
   SUBROUTINE mpi_sum_r_diff( val , val_diff )
      implicit none
      real(rp), intent(inout) :: val
      real(rp), intent(inout) :: val_diff
   END SUBROUTINE mpi_sum_r_diff
   SUBROUTINE mpi_sum_i_diff( val , val_diff )
      implicit none
      integer(ip), intent(inout) :: val
      integer(ip), intent(inout) :: val_diff
   END SUBROUTINE mpi_sum_i_diff
   SUBROUTINE mpi_max_r_diff( val , val_diff )
      implicit none
      real(rp), intent(inout) :: val
      real(rp), intent(inout) :: val_diff
      real(rp), dimension(np) :: temp , temp_diff
      integer(ip) :: index_m
   END SUBROUTINE mpi_max_r_diff
   SUBROUTINE mpi_max_i_diff( val , val_diff )
      implicit none
      integer(ip), intent(inout) :: val
      integer(ip), intent(inout) :: val_diff
      integer(ip), dimension(np) :: temp , temp_diff
      integer(ip) :: index_m
   END SUBROUTINE mpi_max_i_diff
   SUBROUTINE mpi_min_r_diff( val , val_diff )
      implicit none
      real(rp), intent(inout) :: val
      real(rp), intent(inout) :: val_diff
      real(rp), dimension(np) :: temp , temp_diff
      integer(ip) :: index_m
   END SUBROUTINE mpi_min_r_diff
   SUBROUTINE mpi_min_i_diff( val , val_diff )
      implicit none
      integer(ip), intent(inout) :: val
      integer(ip), intent(inout) :: val_diff
      integer(ip), dimension(np) :: temp , temp_diff
      integer(ip) :: index_m
   END SUBROUTINE mpi_min_i_diff
   SUBROUTINE com_dof_diff( dof , dof_diff , mesh )
      implicit none
      !================================================================================================================!
      ! Interface Variables
      !================================================================================================================!
      type( msh ), intent(in ) :: mesh
      type( unk ), intent(inout) :: dof
      type( unk ), intent(inout) :: dof_diff
      !================================================================================================================!
      !
      !================================================================================================================!
         call com_var_r_diff( dof%h(:) , dof_diff%h(:) , mesh )
         call com_var_r_diff( dof%u(:) , dof_diff%u(:) , mesh )
         call com_var_r_diff( dof%v(:) , dof_diff%v(:) , mesh )
   END SUBROUTINE com_dof_diff
END MODULE m_mpi_diff
