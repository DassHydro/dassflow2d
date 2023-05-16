dseq:dagmg.o dagmg_mumps.o
sseq:sagmg.o sagmg_mumps.o
cseq:cagmg.o cagmg_mumps.o
zseq:zagmg.o zagmg_mumps.o
dpar:dagmg_par.o
spar:sagmg_par.o
cpar:cagmg_par.o
zpar:zagmg_par.o
dparNM:dagmg_parNM.o dagmg_mumps.o
sparNM:sagmg_parNM.o sagmg_mumps.o
cparNM:cagmg_parNM.o cagmg_mumps.o
zparNM:zagmg_parNM.o zagmg_mumps.o
lapack:lapack_agmg.o
blas:blas_agmg.o

dagmg.o: dagmg.f90
	$(F90) $(opt) -c dagmg.f90
dagmg_mumps.o: dagmg_mumps.f90
	$(F90) $(opt) -c dagmg_mumps.f90
dagmg_par.o: dagmg_par.f90
	$(F90) $(opt) $(incl) -c dagmg_par.f90
dagmg_parNM.o: dagmg_parNM.f90
	$(F90) $(opt) $(incl) -c dagmg_parNM.f90

sagmg.o: sagmg.f90
	$(F90) $(opt) -c sagmg.f90
sagmg_mumps.o: sagmg_mumps.f90
	$(F90) $(opt) -c sagmg_mumps.f90
sagmg_par.o: sagmg_par.f90
	$(F90) $(opt) $(incl) -c sagmg_par.f90
sagmg_parNM.o: sagmg_parNM.f90
	$(F90) $(opt) $(incl) -c sagmg_parNM.f90

cagmg.o: cagmg.f90
	$(F90) $(opt) -c cagmg.f90
cagmg_mumps.o: cagmg_mumps.f90
	$(F90) $(opt) -c cagmg_mumps.f90
cagmg_par.o: cagmg_par.f90
	$(F90) $(opt) $(incl) -c cagmg_par.f90
cagmg_parNM.o: cagmg_parNM.f90
	$(F90) $(opt) $(incl) -c cagmg_parNM.f90

zagmg.o: zagmg.f90
	$(F90) $(opt) -c zagmg.f90
zagmg_mumps.o: zagmg_mumps.f90
	$(F90) $(opt) -c zagmg_mumps.f90
zagmg_par.o: zagmg_par.f90
	$(F90) $(opt) $(incl) -c zagmg_par.f90
zagmg_parNM.o: zagmg_parNM.f90
	$(F90) $(opt) $(incl) -c zagmg_parNM.f90

lapack_agmg.o: lapack_agmg.f
	$(F90) $(opt) -c lapack_agmg.f
blas_agmg.o: blas_agmg.f
	$(F90) $(opt) -c blas_agmg.f
