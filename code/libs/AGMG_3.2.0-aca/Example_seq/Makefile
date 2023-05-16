# your FORTRAN compiler
F90=ifort

# set here link reference needed for BLAS & LAPACK
BLASLAPACK=-L/usr/lib -llapack -lblas  

# compilation options
opt = -O2 #-fdefault-integer-8

# directory where AGMG sources are located
agmgdir=../SRC

lib=$(BLASLAPACK)

export F90
export opt

list=Example_seq.o $(agmgdir)/dagmg.o $(agmgdir)/dagmg_mumps.o

Example_seq: Example_seq.o $(agmgdir)/dagmg.f90 $(agmgdir)/dagmg_mumps.f90
	cd $(agmgdir);make dseq
	$(F90) $(opt) -o Example_seq $(list) $(lib)

Example_seq.o: Example_seq.f90
	$(F90) $(opt) -c Example_seq.f90
