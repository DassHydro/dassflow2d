!===============================================================================================================================
!
! Driver for using the M1QN3 solver (version 3.3) on the CUTEr collection
! =======================================================================
!
! M1qn3 solves an unconstrained minimization problem:
!
!    minimize f(x)
!
! using a limited memory quasi-Newton method (l-BFGS). More information on m1qn3 can be obtained at the address:
!
!    http://www-rocq.inria.fr/estime/modulopt/optimization-routines/m1qn3/m1qn3.html
!
! This file contains the main program (m1qn3ma) and a simulator (cutersim, which calls the simulator 'uofg' of the
! CUTEr environment). This driver is written in Fortran 90/95.
!
! Memory allocation
! ^^^^^^^^^^^^^^^^^
! In order to adapt properly the size of the variables of this main program to the dimension of the selected CUTEr
! problem, this dimension (= number of variables) is first read in the data file 'OUTSDIF.d' (it is assumed here
! that this dimension is the first data of the file), then the variables are allocated, next the subroutine 'usetup'
! is called. There can still be a memory-size problem due to the size of the variables defined in the CUTEr
! environment, which cannot be fixed at the level of this program.
!
! Author: J.Ch. Gilbert, November, 2007
!
!===============================================================================================================================

  program m1qn3ma

  implicit none

! Parameters

  double precision, parameter :: zero = 0.d0

! General variables

  integer :: i, j

! Variables related to CUTEr

  double precision, parameter :: minf = -1.0d20
  double precision, parameter :: pinf =  1.0d20

  character(len=10) :: pname
  integer :: ic, n, naux

  character(len=10), allocatable :: vnames(:)

! parameters and variables related to M1QN3

  character(len=3) :: normtype
  integer :: reverse, indic, mupdts, impres, io, imode(3), omode, niter, nsim, iz(5), ndz, izs
  real :: rzs
  double precision :: f, df1, df1_factor, prec, epsg, gn1, dzs, dxmin

  double precision, allocatable :: x(:), lx(:), ux(:), g(:), dz(:)

  external cutersim, euclid, ctonbe, ctcabe

! parameters and variables related to LIBOPT

  integer :: info

! variables and functions for CPU time

  real :: cput(2), actual_time, cput_init, total_time
  real :: etime
! external etime

!-----------------------------------------------------------------------
! initialization
!-----------------------------------------------------------------------

  io     =  6   ! input-output channel
  ic     = 15   ! cuter data file

!-----------------------------------------------------------------------
! get the dimensions from CUTEr
!-----------------------------------------------------------------------

  open (ic, file='OUTSDIF.d', form='formatted')
  call udimen (ic, n)

  if (n <= 0) then
    write (io,'(a)') "(qpalma) >>> error: udimen returns a nonpositive 'n'"
    stop
  end if

!-----------------------------------------------------------------------
! allocate variables depending on n
!-----------------------------------------------------------------------

  allocate (x(n), stat=i)
  if (i /= 0) then
    write (6,'(a)') "(m1qn3ma) >>> error in the allocation of 'x'"
    stop
  end if
  allocate (g(n), stat=i)
  if (i /= 0) then
    write (6,'(a)') "(m1qn3ma) >>> error in the allocation of 'g'"
    stop
  end if
  allocate (lx(n), stat=i)
  if (i /= 0) then
    write (6,'(a)') "(m1qn3ma) >>> error in the allocation of 'lx'"
    stop
  end if
  allocate (ux(n), stat=i)
  if (i /= 0) then
    write (6,'(a)') "(m1qn3ma) >>> error in the allocation of 'ux'"
    stop
  end if
  allocate (vnames(n), stat=i)
  if (i /= 0) then
    write (6,'(a)') "(m1qn3ma) >>> error in the allocation of 'vnames'"
    stop
  end if

!-----------------------------------------------------------------------
! get other data from CUTEr
!-----------------------------------------------------------------------

  rewind ic
  call usetup (ic, io, naux, x, lx, ux, n)
  close(ic)

!-----------------------------------------------------------------------
! check that the problem is suitable for M1QN3 (it must be unconstrained)
!-----------------------------------------------------------------------

  do i = 1,n
    if (lx(i) == ux(i)) then
      write (6,'(a,i0,a)') "(m1qn3ma) >>> variable ", i, " is fixed"
      stop
    else
      if ((lx(i) > minf) .or. (ux(i) < pinf)) then
        write (6,'(a,i0,a)') "(m1qn3ma) >>> variable ", i, " is bounded"
        stop
      end if
    end if
  end do

!-----------------------------------------------------------------------
! get the problem name
!-----------------------------------------------------------------------

  call unames (n, pname, vnames)

!-----------------------------------------------------------------------
! call the simulator (this is required by m1qn3)
!-----------------------------------------------------------------------

  indic = 4
  call cutersim (indic, n, x, f, g, izs, rzs, dzs)

!-----------------------------------------------------------------------
! call the optimization routine
!-----------------------------------------------------------------------

  ! set the default values of the solver parameters

  impres   =     3      ! print level (=0: no output, =1: 1st and last iteration, =3: 1 line per iteration, =5: talkative)
  imode(1) =     0      ! running mode (=0: diagonal, =1: scalar)
  imode(2) =     0      ! starting mode (=0: cold start, =1: warm restart)
  imode(3) =     0      ! call with indic=1 (=0: never, >0: every imode(3) iteration)
  mupdts   =    20      ! nbr of l-BFGS updates to form the inverse Hessian approximation
  niter    = 10000      ! max number of iterations
  nsim     = 20000      ! max number of simultations (computations of function and gradient)
  df1 = 1.d-1*abs(f)    ! looks better than df1 = 1.d-8; if (f > zero) df1 = 1.d-1*f
  df1_factor =   0.d+0  ! (<=0.d0: ignored; >0: expected decrease in f at the first iteration is df1_factor*abs(f(x_init))
  dxmin    =     1.d-20 ! resolution in x for the sup-norm (> 0.d0)
  epsg     =     1.d-7  ! required relative gradient accuracy in Euclidean norm (>0)
  reverse  =     0

  ! read the specification file, possibly modifying the default values

  call m1qn3_spc_read (io,mupdts,df1_factor,prec,dxmin,impres,imode,niter,nsim,info)
  if (info < 0) stop
  if (df1_factor > zero) df1 = df1_factor*abs(f)

  ! allocate the workspace for M1QN3

  ndz = 4*n+mupdts*(2*n+1)
  allocate (dz(ndz), stat=i)
  if (i /= 0) then
    write (6,'(a)') "(m1qn3ma) >>> error in the allocation of 'dz'"
    stop
  end if

  ! determine relative precision

  normtype = 'sup'      ! 'sup' for sup-norm, 'two' for 2-norm, 'dfn' for the norm defined by prosca
  gn1 = zero
  do i=1,n
    gn1 = max(gn1,abs(g(i)))
  end do
  epsg = prec/gn1

  ! call M1QN3, measuring the CPU time

  actual_time = etime (cput)
  cput_init = cput(1)

  call m1qn3 (cutersim, euclid, ctonbe, ctcabe, &
              n, x, f, g, dxmin, df1, epsg, normtype, impres, io, imode, omode, niter, nsim, iz, dz, ndz, &
              reverse, indic, izs, rzs, dzs)

  actual_time = etime (cput)
  total_time = max(0.e0, cput(1)-cput_init)
  write (6,'(/a,es12.5)') "total_time = ", total_time

!-----------------------------------------------------------------------
! Recover the results
!-----------------------------------------------------------------------

  ! print the solution found in case this is useful

  write (io,'(a,5(es21.14,1x))') "x = ", (x(i), i=1,n)

  ! the libopt line

  if (omode <= 1) omode = omode-1       ! in order to have 0 when the solver succeeds

  write (6,'(/a,a)',advance='no') "libopt%m1qn3%cuter%", trim(pname)
  write (6,'("%n=",i0,"%prec=",es8.2)',advance='no') n, prec
  write (6,'("%nfc=",i0,"%nga=",i0)',advance='no') nsim, nsim
  write (6,'("%time=",es8.2)',advance='no') total_time
  write (6,'("%objf=",es14.7,"%glag=",es11.5)',advance='no') f, epsg*gn1
  write (6,'("%info=",i0/)') omode

  stop

  end program m1qn3ma

!===============================================================================================================================

  subroutine cutersim (indic, n, x, f, g, izs, rzs, dzs)

  implicit none
  integer :: indic, n, izs
  real :: rzs
  double precision :: x(n), f, g(n), dzs

!-----------------------------------------------------------------------
!
! Simulator for m1qn3, which calls the CUTEr simulator 'uofg'.
!
!-----------------------------------------------------------------------

! parameters

  integer, parameter :: io = 6

! local variable

  logical :: grad_required

! start computation ------------------------------------------------

  if (indic.eq.1) then

    ! do nothing

  elseif (indic.eq.4) then

    ! evaluate the objective function and its gradient

    grad_required = .true.
    call uofg (n, x, f, g, grad_required)

  else

    ! unexpected value of indic

    write(io,'(/a,i5/)') '>>> CUTERSIM: unexpected value of indic = ',indic
    indic = 0   ! stop

  end if

! there is no reason to reset indic to a positive value here

  return
  end

!===============================================================================================================================
! M1QN3 specification file reader
!
! on return
! - info = 0: fine reading
!        < 0: problem
!        = 1: no specification file found
!===============================================================================================================================

  subroutine m1qn3_spc_read (io,mupdts,df1_factor,epsg,dxmin,impres,imode,niter,nsim,info)

  implicit none

  integer :: io, mupdts, impres, imode(3), niter, nsim, info
  double precision :: df1_factor, dxmin, epsg

  intent(in)  :: io
  intent(out) :: mupdts, impres, imode, niter, nsim, info, df1_factor, dxmin, epsg

! local variables

  character(len=128) :: spcfile_name    ! specification file name
  character(len=128) :: line            ! one line of the specification file
  logical :: spcfile_exists
  integer :: i, l, spcfile_id, stat, length

! default is a fine reading

  info = 0

! open the specification file

  spcfile_name = 'm1qn3.spc'
  inquire (file=spcfile_name, exist=spcfile_exists)
  if (.not. spcfile_exists) then
    write (io,'(/a)') "(m1qn3ma) >>> specification file 'm1qn3.spc' not found"
    info = 1
    return
  end if
  write (io,'(a,a,a)') "(m1qn3ma) reading specification file '", trim(spcfile_name), "'" ! should only be printed in verbose

  spcfile_id = 11
  open (unit=spcfile_id, file=spcfile_name, action='read', iostat=stat)
  if (stat /= 0) then
    write (io,'(a,a,a)') "(m1qn3ma) >>> cannot open specification file '", trim(spcfile_name), "'"
    info = -1
    return
  end if

! read the specification file

  l = 0 ! line number

  do 

    l = l+1

    ! read the line

    read (spcfile_id, '(a)', end=1, iostat=stat) line
    if (stat /= 0) then
      write (io,'(a,i0,a)') "(m1qn3ma) >>> specification file, line ", l, ": not readable"
      info = -1
      return
    end if

    ! remove the comment and skip an empty line

    line = adjustl(line)        ! remove leading blanks
    i = index(line,'!')         ! a comment starts after an exclamation mark
    if (i == 0) then            ! no exclamation mark
      length = len_trim(line)
    else                        ! there is an exclamation mark at position i
      line = line(1:i-1)
      length = i-1
    end if
    if (length <= 0) cycle

    ! decode the line

    if (index(line,"mupdts") == 1) then
      line = adjustl(line(7:length))
      read (line,*) mupdts
    elseif (index(line,"df1_factor") == 1) then
      line = adjustl(line(11:length))
      read (line,*) df1_factor
    elseif (index(line,"impres") == 1) then
      line = adjustl(line(7:length))
      read (line,*) impres
    elseif (index(line,"imode1") == 1) then
      line = adjustl(line(7:length))
      read (line,*) imode(1)
    elseif (index(line,"imode2") == 1) then
      line = adjustl(line(7:length))
      read (line,*) imode(2)
    elseif (index(line,"imode3") == 1) then
      line = adjustl(line(7:length))
      read (line,*) imode(3)
    elseif (index(line,"niter") == 1) then
      line = adjustl(line(6:length))
      read (line,*) niter
    elseif (index(line,"nsim") == 1) then
      line = adjustl(line(5:length))
      read (line,*) nsim
    elseif (index(line,"dxmin") == 1) then
      line = adjustl(line(6:length))
      read (line,*) dxmin
    elseif (index(line,"epsg") == 1) then
      line = adjustl(line(5:length))
      read (line,*) epsg
    else
      i = 1
      do while ((i <= 128) .and. (line(i:i) /= ' '))
        i = i+1
      end do
      i = i-1
      write (io,'(a,i0,a,a,a)') "(m1qn3ma) >>> specification file, line ", l, ": unrecognized keyword '", &
        trim(line(1:i)), "', line ignored"
    end if

  end do

! close the specification file

1 close (spcfile_id)
  return

  end subroutine m1qn3_spc_read

