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

subroutine f90wrap_read_bc_file(line_read)
    implicit none
    external read_bc_file
    
    integer(4), intent(out) :: line_read
    call read_bc_file(line_read)
end subroutine f90wrap_read_bc_file

