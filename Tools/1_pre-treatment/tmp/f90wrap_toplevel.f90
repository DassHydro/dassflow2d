subroutine f90wrap_gen_basic_channel(nx, ny, lx, ly)
    implicit none
    external gen_basic_channel
    
    integer(4), intent(in) :: nx
    integer(4), intent(in) :: ny
    real(4), intent(in) :: lx
    real(4), intent(in) :: ly
    call gen_basic_channel(nx, ny, lx, ly)
end subroutine f90wrap_gen_basic_channel

subroutine f90wrap_gen_bc(in_type, out_type)
    implicit none
    external gen_bc
    
    character(10), intent(in) :: in_type
    character(10), intent(in) :: out_type
    call gen_bc(in_type, out_type)
end subroutine f90wrap_gen_bc

subroutine f90wrap_gen_land_use(manning_alpha, manning_beta)
    implicit none
    external gen_land_use
    
    real(4), intent(in) :: manning_alpha
    real(4), intent(in) :: manning_beta
    call gen_land_use(manning_alpha, manning_beta)
end subroutine f90wrap_gen_land_use

subroutine f90wrap_gen_bc_data(bc_typ, nrow, var1, var2, n0, n1)
    implicit none
    external gen_bc_data
    
    character(10), intent(in) :: bc_typ
    integer(4), intent(in) :: nrow
    real, intent(in), dimension(n0) :: var1
    real, intent(in), dimension(n1) :: var2
    integer :: n0
    !f2py intent(hide), depend(var1) :: n0 = shape(var1,0)
    integer :: n1
    !f2py intent(hide), depend(var2) :: n1 = shape(var2,0)
    call gen_bc_data(bc_typ, nrow, var1, var2)
end subroutine f90wrap_gen_bc_data

subroutine f90wrap_gen_obs(nx_obs, ny_obs, xmax_obs, ymax_obs, xmin_obs, ymin_obs, dt_obs)
    implicit none
    external gen_obs
    
    integer(4), intent(in) :: nx_obs
    integer(4), intent(in) :: ny_obs
    real(4), intent(in) :: xmax_obs
    real(4), intent(in) :: ymax_obs
    real(4), intent(in) :: xmin_obs
    real(4), intent(in) :: ymin_obs
    real(4), intent(in) :: dt_obs
    call gen_obs(nx_obs, ny_obs, xmax_obs, ymax_obs, xmin_obs, ymin_obs, dt_obs)
end subroutine f90wrap_gen_obs

subroutine f90wrap_h_true_macdo(x, lx, g, h_true)
    implicit none
    external h_true_macdo
    
    real(4), intent(in) :: x
    real(4), intent(in) :: lx
    real(4), intent(in) :: g
    real(4), intent(out) :: h_true
    call h_true_macdo(x, lx, g, h_true)
end subroutine f90wrap_h_true_macdo

