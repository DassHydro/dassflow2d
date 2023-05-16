subroutine f90wrap_read_input(filename)
    implicit none
    external read_input
    
    character*(*), intent(in) :: filename
    call read_input(filename)
end subroutine f90wrap_read_input

subroutine f90wrap_print_all
    implicit none
    external print_all
    
    call print_all()
end subroutine f90wrap_print_all

subroutine f90wrap_read_bc_file
    implicit none
    external read_bc_file
    
    call read_bc_file()
end subroutine f90wrap_read_bc_file

