#-----------------------------------------------------------------------
# This Makefile generates the binaries to run 'm1qn3' on the
# test-problems 'essai', 'essai_rc', 'u1mt1', and 'u1mt1_rc'. See the
# README file for more details.
#-----------------------------------------------------------------------

# The full pathname of the M1QN3 directory is set by removing the
# string '/examples' from the end of the present directory pathname

m1qn3_dir = $(PWD:/examples=)

# Some useful directories

m1qn3_bin_dir = $(m1qn3_dir)/bin
m1qn3_src_dir = $(m1qn3_dir)/src
m1qn3_bin     = $(m1qn3_bin_dir)/m1qn3.$(PLAT).o
m1qn3_src     = $(m1qn3_src_dir)/m1qn3.f

# Include file

include $(m1qn3_dir)/make.$(PLAT)

# Determine which BLAS library to link

m1qn3_bins    = $(shell ls $(m1qn3_bin_dir))

blas_lib      = $(findstring libblas.$(PLAT).a, $(m1qn3_bins))
ifeq ($(blas_lib),)
blas_lib      = -lblas
else
blas_lib      = $(m1qn3_bin_dir)/libblas.$(PLAT).a
endif

#-----------------------------------------------------------------------

# Make nothing in case 'make' is typed uncautiously

nothing:

#-----------------------------------------------------------------------

essai: $(m1qn3_bin) essai.o
	$(F77) -o essai essai.o \
	    $(m1qn3_bin) \
	    $(blas_lib)

essai_wr: $(m1qn3_bin) essai_wr.o
	$(F77) -o essai_wr essai_wr.o \
	    $(m1qn3_bin) \
	    $(blas_lib)

essai_rc: $(m1qn3_bin) essai_rc.o
	$(F77) -o essai_rc essai_rc.o \
	    $(m1qn3_bin) \
	    $(blas_lib)

essai_rc_wr: $(m1qn3_bin) essai_rc_wr.o
	$(F77) -o essai_rc_wr essai_rc_wr.o \
	    $(m1qn3_bin) \
	    $(blas_lib)

u1mt1: $(m1qn3_bin) u1mt1.o
	$(F77) -o u1mt1 u1mt1.o \
	    $(m1qn3_bin) \
	    $(blas_lib)

u1mt1_rc: $(m1qn3_bin) u1mt1_rc.o
	$(F77) -o u1mt1_rc u1mt1_rc.o \
	    $(m1qn3_bin) \
	    $(blas_lib)

$(m1qn3_bin): $(m1qn3_src)
	$(F77) $(F77FLAGS) $(F77NOLD) \
	    $(m1qn3_src) \
	    -o $(m1qn3_bin)

clean:
	rm -f *.o essai essai_wr essai_rc essai_rc_wr u1mt1 u1mt1_rc res

purge:
	rm -f *.o essai essai_wr essai_rc essai_rc_wr u1mt1 u1mt1_rc res \
	    $(m1qn3_bin)

#-----------------------------------------------------------------------

%.o: %.f
	$(F77) -c $<
