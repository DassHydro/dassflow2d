#-----------------------------------------------------------------------
# This makefile compiles all the BLAS objects in this directory and
# store them in the archive "../bin/libblas.$(PLAT).a", where $(PLAT)
# identifies the current platform. This one must be described in the
# file "../make.$(PLAT)".
#-----------------------------------------------------------------------

# The full pathname of the M1QN3 directory is set by removing the
# string '/examples' from the end of the present directory pathname

m1qn3_dir = $(PWD:/blas=)

# Include file

include $(m1qn3_dir)/make.$(PLAT)

# Some useful directories and files

m1qn3_bin_dir = $(m1qn3_dir)/bin
blas_lib_plat = $(m1qn3_bin_dir)/libblas.$(PLAT).a

#-----------------------------------------------------------------------

$(blas_lib_plat): \
            $(blas_lib_plat)(ddot.o)
	ranlib $(blas_lib_plat)

clean:
	/bin/rm -f *.o $(blas_lib_plat)

#-----------------------------------------------------------------------
# $@ is the file name of the target of the rule
# $< (if the target got its commands from an implicit rule) this will be
#    the first prerequisite added by the implicit rule

%.o: %.f
	$(F77) $(F77FLAGS) -c $<
