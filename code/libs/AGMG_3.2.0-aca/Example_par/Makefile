# MPIopt    = -I/... (where to find mpif.h)
# MUMPSPopt = -I/... (where to find files to be included by 
#                     applications using MUMPS)
# MUMPSlib  = -l...  (link reference for MUMPS)
# SCALAP    = -l...  (link reference for SCALAPACK, needed by MUMPS)
# BLASLAPACK= -l...  (link reference for LAPACK & BLAS)
# MPIlib    = -l...  (link reference for MPI)

# your FORTRAN compiler
F90=ifort

# compilation options
opt = -O4

# directory where AGMG sources are located
agmgdir  = ../SRC


incl=$(MPIopt) $(MUMPSopt)

lib = $(MUMPSlib) $(SCALAP) $(MPIlib) $(BLASLAPACK)

export F90
export opt
export incl

list=Example_par.o $(agmgdir)/dagmg_par.o

Example_par: Example_par.o $(agmgdir)/dagmg_par.f90
	cd $(agmgdir);make dpar
	$(F90) $(opt) -o Example_par $(list) $(lib)

Example_par.o: Example_par.f90 
	$(F90) $(opt) $(incl) -c Example_par.f90
