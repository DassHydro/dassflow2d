!=======================================================================
! Specification file for the solver M1QN3
!=======================================================================
! To use this specification file, place a *copy* in the working
! directory and adapt it to your needs (the present file, however,
! should not be changed). Not all the specifications of this file are
! used by all the programs running M1QN3, since for some collections of
! problems (in particular Modulopt and Modulopttoys), the specifications
! are obtained from other sources.
!=======================================================================
! - Each specification is written on a single line, starts whith a
!   'keyword', is followed by a 'value', and ends with an optional
!   comment.
! - An absent specification is given a default value.
! - The position of the keywords and values are free (multiple blanks
!   are considered as simple blank); a comment goes from an exclamation
!   mark to the end of the line.
! - The only known keywords are those used below.
! - The type of the value must be the same as the one given below.
!=======================================================================

imode1     0      ! running mode (=0: diagonal, =1: scalar)
imode2     0      ! starting mode (=0: cold start, =1: warm restart)
imode3     0      ! call with indic=1 (=0: never, >0: every imode(3) iteration)
mupdts    20      ! nbr of l-BFGS updates to form the inverse Hessian approximation
impres     3      ! printing level (0: no output, 1: at the 1st and last iteration, 3: 1 line per iteration, 5: talkative)
niter  10000      ! max number of iterations
nsim   20000      ! max number of simultations (computations of function and gradient)
df1_factor 0.d+0  ! (<=0.d0: ignored; >0: expected decrease in f at the first iteration is df1_factor*abs(f(x_init))
dxmin      1.d-20 ! resolution in x for the sup-norm (> 0.d0)
epsg       1.d-5  ! required absolute gradient accuracy in sup-norm (>0)
