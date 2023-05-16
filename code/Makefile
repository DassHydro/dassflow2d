#======================================================================================================================#
#
#                    DassFlow Version 3.0
#
#======================================================================================================================#
#
#  Copyright University of Toulouse-INSA, Univ. of Strasbourg, INRAE & CNRS (France)
#
#  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
#  DassFlow is a computational software whose purpose is to simulate geophysical free surface flows,
#  designed for variational sensitivities and data assimilation (4D-var). Inverse capabilities are
#  based on the adjoint code generation by a source-to-source algorithmic differentiation (Tapenade software used).
#
#  DassFlow software includes few mostly independent "modules" with common architectures and structures.
#  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
#
#  Many people have contributed to the DassFlow development from the initial version to the latest ones.
#  Current contributions:
#               L. Pujol (PhD Unistra)
#               L. Villenave (PhD student)
#               P.-A. Garambois (INRAE Aix-en-Provence)
#               J. Monnier (INSA & Mathematics Institute of Toulouse IMT).
#               K. Larnier (CS group - IMT-INSA).
#  Former scientific or programming contributions of:
#               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
#               J.-P. Vila (INSA & Mathematics Institute of Toulouse IMT).
#               R. Madec   (Mathematics Institute of Toulouse IMT).
#  plus less recent other developers (M. Honnorat and J. Marin).
#
#  Contact : see the DassFlow webpage
#
#  This software is governed by the CeCILL license under French law and abiding by the rules of distribution
#  of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license
#  as circulated by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
#
#  As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the
#  license, users are provided only with a limited warranty and the software's author, the holder of the economic
#  rights, and the successive licensors have only limited liability.
#
#  In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or
#  developing or reproducing the software by the user in light of its specific status of free software, that may
#  mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
#  experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
#  software's suitability as regards their requirements in conditions enabling the security of their systems and/or
#  data to be ensured and, more generally, to use and operate it in the same conditions as regards security.
#
#  The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you
#  accept its terms.
#
#======================================================================================================================#

#======================================================================================================================#
#   User Options
#======================================================================================================================#

include Makefile.inc


#======================================================================================================================#
#   Executable name created in the /bin directory
#======================================================================================================================#

executable  =  exe

#======================================================================================================================#
#   Module name and package option for wrapper part
#======================================================================================================================#

wrapping_name =  wrapping
DOC_PLUGIN    = doc_plugin.py
move         = -M

string       = -S 128

pkg = -P

debug_wrapping =

ifeq ($(OPTIM),0)
	ifeq ($(DEBUG_WRAPPING), 1)
		debug_wrapping += --debug-capi
	endif
endif

#======================================================================================================================#
#   Compilator Flags
#======================================================================================================================#

ifeq ($(COMPILO),1)

  ifeq ($(MPI),1)
    FC = mpif90
    CC = mpicc
  else
    FC = ifort
    CC = icc
  endif

  ifeq ($(OPTIM),0)
    export CFLAGS = fpp -module ./obj -g -traceback -check all -check noarg_temp_created -fpe0 -nofor_main -diag-disable 8291 -fPIC
    export FFLAGS = -nofor_main -fPIC

  else

    export CFLAGS = fpp -module ./obj -fast -nofor_main -diag-disable 8291
    export FFLAGS = -nofor_main -fPIC

  endif

else


  ifeq ($(MPI),1)
    FC = mpif90
    CC = mpicc
  else
    FC = gfortran
    CC = gcc
  endif

  ifeq ($(OPTIM),0)

    export CFLAGS = -g -fbounds-check -ffree-line-length-0 -I./obj -J./obj -fPIC -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
    export FFLAGS = -fPIC

  else

    export CFLAGS = -O3 -fno-range-check -ffree-line-length-0 -I./obj -J./obj -fPIC
    export FFLAGS = -fPIC

  endif

endif

export CLIBS =
export FLIBS =

#======================================================================================================================#
#   Scotch
#======================================================================================================================#

SCOTCH = $(PWD)/libs/scotch_5.1.12_esmumps

ifeq ($(MPI),1)

  export CLIBS  +=  -I$(SCOTCH)/include

  export FLIBS  +=    $(SCOTCH)/lib/libscotch.a
  export FLIBS  +=    $(SCOTCH)/lib/libscotcherr.a

endif

#======================================================================================================================#
#   MUMPS
#======================================================================================================================#

MUMPS = $(PWD)/libs/MUMPS_4.10.0

ifeq ($(SOLVER),1)

  export CLIBS  +=  -I$(SCOTCH)/include
  export FLIBS  +=    $(SCOTCH)/lib/libscotch.a
  export FLIBS  +=    $(SCOTCH)/lib/libscotcherr.a

  export CLIBS  +=  -I$(MUMPS)/include
  export FLIBS  +=    $(MUMPS)/lib/libdmumps.a
  export FLIBS  +=    $(MUMPS)/lib/libmumps_common.a
  export FLIBS  +=  -L$(MUMPS)/PORD/lib
  export FLIBS  +=  -L$(MUMPS)/libseq
	
  export FLIBS  +=  -lpord -lmpiseq -lblas -lpthread -L$(PWD)/libs/scotch_5.1.12_esmumps/lib -lesmumps -lscotch -lscotcherr

endif

#======================================================================================================================#
#   AGMG
#======================================================================================================================#

AGMG = $(PWD)/libs/AGMG_3.2.0-aca/SRC

ifeq ($(SOLVER),2)

  export FLIBS  +=  $(AGMG)/dagmg.o
  export FLIBS  +=  $(AGMG)/dagmg_mumps.o
  export FLIBS  +=  -L/usr/lib
  export FLIBS  +=  -llapack -lblas

endif

#======================================================================================================================#
#   Minimization method: M1QN3, N2QN1 (not implemented), LBFGSB-3
#======================================================================================================================#

ifeq ($(MINMETHOD),1)
  CPP_MIN = -DUSE_M1QN3
  OBJS_MINMETHOD=obj/m1qn3.o\
                 obj/ddot.o
endif

ifeq ($(MINMETHOD),2)
  CPP_MIN = -DUSE_LBFGSB3
  OBJS_MINMETHOD=obj/lbfgsb.o\
                 obj/lbfgsb_timer.o\
                 obj/lbfgsb_blas.o\
                 obj/lbfgsb_linpack.o
endif

#======================================================================================================================#
#	CPP options according the options choices above
#======================================================================================================================#

CPP_MODEL     =
CPP_MPI       =
CPP_MPI_ADJ   =
CPP_ADJ       =
CPP_SOL       =
CPP_VAL       =
CPP_COMP      =
CPP_PASS      =
CPP_HYD       =

ifeq ($(COMPILO),1)
  CPP_COMP = -DUSE_INTEL
endif

ifeq ($(MODEL),1)
  CPP_MODEL = -DUSE_SW_MONO
endif

ifeq ($(MODEL),2)
  CPP_MODEL = -DUSE_SW_MULTI
endif

ifeq ($(MODEL),3)
  CPP_MODEL = -DUSE_NS_MULTIFLUID
endif

ifeq ($(HYDRO),1)
  CPP_HYD = -DUSE_HYDRO
endif

ifeq ($(MPI),1)

  CPP_MPI = -DUSE_MPI
  ifeq ($(ADJOINT),1)
  CPP_MPI_ADJ = -DUSE_MPI_ADJ
  endif
endif

ifeq ($(ADJOINT),1)
  CPP_ADJ = -DUSE_ADJ
endif

ifeq ($(SOLVER),1)
  CPP_SOL = -DUSE_MUMPS
endif

ifeq ($(SOLVER),2)
  CPP_SOL = -DUSE_AGMG
endif

ifeq ($(VALID),1)
  CPP_VAL = -DUSE_VALID
endif

CPP_FLAGS		=  $(CPP_COMP) $(CPP_MODEL)  $(CPP_HYD) $(CPP_ADJ) $(CPP_MPI) $(CPP_SOL) $(CPP_VAL) $(CPP_MIN)
CPP_FLAGS_WRAP		=  $(CPP_COMP) $(CPP_MODEL) -DCPP_PASS  $(CPP_HYD) $(CPP_ADJ) $(CPP_MPI) $(CPP_SOL) $(CPP_VAL)
CPP_FLAGS_ADJ  		=  $(CPP_COMP) $(CPP_MODEL) $(CPP_HYD) $(CPP_MPI_ADJ)

#======================================================================================================================#
#
#======================================================================================================================#

MODELD =

ifeq ($(MODEL),1)
  MODELD = src/sw_mono
endif

ifeq ($(MODEL),2)
  MODELD = src/sw_multi
endif

ifeq ($(MODEL),3)
  MODELD = src/ns_multifluid
endif

ADJD =

ifeq ($(ADJOINT),1)
  ADJD = src/adjoint
endif

HYDROD   =
HYDROTAP =
HYDROVAR =

ifeq ($(HYDRO),1)
  HYDROD   = src/hydro
  #HYDROTAP =  tap/gr4_main.f90\
#	      tap/mod_launch.f90\
#	      tap/gr4_stsp.f90\
#	      tap/crit_calc.f90
  HYDROVAR =  bc%gr4%params bc%gr4%Q bc%gr4%t
endif


DIRS      =  src src/common src/base $(MODELD) $(HYDROD) $(ADJD) $(CASEDIR)
DIRS_INC  =  src/base/include $(MODELD)/include

POST_SRC  =  $(foreach DIR, $(DIRS)    , $(patsubst $(DIR)/%.f90,cpp/%.f90,$(wildcard $(DIR)/*.f90)))
POST_INC  =  $(foreach DIR, $(DIRS_INC), $(patsubst $(DIR)/%.f90,cpp/%.inc,$(wildcard $(DIR)/*.f90)))

VPATH  =  src:src/common:src/base:$(MODELD):$(HYDROD):$(ADJD):$(CASEDIR):src/base/include:$(MODELD)/include


ifeq ($(ADJOINT),1)

 ifeq ($(HYDRO),1)
  OBJS_TAP  =  obj/m_tap_vars.o\
    obj/m_mpi_back.o\
    obj/m_mpi_diff.o\
    obj/m_gr4.o\
    obj/m_numeric_back.o\
    obj/m_adjoint.o\
    obj/m_numeric_diff.o\
    $(patsubst src/adjoint/%.f90,obj/%.o,$(wildcard src/adjoint/*.f90))\
    $(patsubst src/hydro/%.f90,obj/%.o,$(wildcard src/hydro/*.f90))\
    obj/adBuffer.o\
    obj/adStack.o\
    obj/lbfgsb_blas.o\
    obj/lbfgsb_linpack.o\
    obj/lbfgsb_timer.o\
    obj/lbfgsb.o\
    obj/m1qn3.o\
    obj/ddot.o\
    $(patsubst tap/%.f90,obj/%.o,$(wildcard tap/m_*.f90))\
    $(patsubst tap/%.f90,obj/%.o,$(wildcard tap/*.f90))
 else
  OBJS_TAP  =  obj/m_tap_vars.o\
			   obj/m_mpi_back.o\
			   obj/m_mpi_diff.o\
			   obj/m_gr4.o\
			   obj/m_numeric_back.o\
			   obj/m_adjoint.o\
			   obj/m_numeric_diff.o\
			   $(patsubst src/adjoint/%.f90,obj/%.o,$(wildcard src/adjoint/*.f90))\
			   obj/adBuffer.o\
			   obj/adStack.o\
			   obj/lbfgsb_blas.o\
		   obj/lbfgsb_linpack.o\
		   obj/lbfgsb_timer.o\
		   obj/lbfgsb.o\
			   obj/m1qn3.o\
			   obj/ddot.o\
			   $(patsubst tap/%.f90,obj/%.o,$(wildcard tap/m_*.f90))\
			   $(patsubst tap/%.f90,obj/%.o,$(wildcard tap/*.f90))
endif


else
  OBJS_TAP  =
endif


OBJS_MOD  =  obj/m_common.o\
             obj/m_linear_algebra.o\
             obj/m_mesh.o\
	     obj/m_mpi.o \
	     obj/m_time_screen.o\
	     obj/m_sw_mono.o \
	     obj/m_numeric.o \
	     obj/m_obs.o \
             $(patsubst $(MODELD)/%.f90,obj/%.o,$(wildcard $(MODELD)/m_*.f90))\
             $(patsubst $(HYDROD)/%.f90,obj/%.o,$(wildcard $(HYDROD)/m_*.f90))\
             obj/m_linear_solver.o

OBJS_SRC  =  obj/m_common.o\
             obj/m_linear_algebra.o\
             obj/m_mesh.o\
             obj/m_mpi.o\
	     obj/m_time_screen.o\
	     obj/m_sw_mono.o \
	     obj/m_numeric.o\
	     obj/m_obs.o \
             $(patsubst $(MODELD)/%.f90,obj/%.o,$(wildcard $(MODELD)/m_*.f90))\
             $(patsubst $(HYDROD)/%.f90,obj/%.o,$(wildcard $(HYDROD)/m_*.f90))\
             $(patsubst src/common/%.f90,obj/%.o,$(wildcard src/common/*.f90))\
             $(patsubst src/base/%.f90,obj/%.o,$(wildcard src/base/*.f90))\
             $(OBJS_USER)\
             $(patsubst $(MODELD)/%.f90,obj/%.o,$(wildcard $(MODELD)/*.f90))\
             $(patsubst $(HYDROD)/%.f90,obj/%.o,$(wildcard $(HYDROD)/*.f90))\
             $(OBJS_TAP)\
             obj/main.o

#======================================================================================================================#
#
#======================================================================================================================#

ifeq ($(MODEL),1)

   FILES_TAP  =  tap/m_common.f90\
                 tap/m_linear_algebra.f90\
                 tap/m_mesh.f90\
                 tap/m_mpi.f90\
                 tap/m_numeric.f90\
                 tap/m_sw_mono.f90\
                 tap/m_obs.f90\
                 tap/run_model_sw_mono.f90\
                 tap/advance_time.f90\
                 tap/boundary.f90\
		  tap/euler_time_step_first_b1.f90\
                 tap/sw_flux.f90\
                 tap/read_write_obs.f90\
                 tap/numeric_sw.f90\
                 tap/sw_post.f90

endif

ifeq ($(HYDRO),1)
	   FILES_TAP  =  tap/m_common.f90\
		         tap/m_linear_algebra.f90\
		         tap/m_mesh.f90\
		         tap/m_mpi.f90\
		         tap/m_numeric.f90\
		         tap/m_sw_mono.f90\
		         tap/m_obs.f90\
		         tap/m_gr4.f90\
		         tap/run_model_sw_mono.f90\
		         tap/advance_time.f90\
		         tap/boundary.f90\
			 tap/euler_time_step_first_b1.f90\
		         tap/sw_flux.f90\
		         tap/read_write_obs.f90\
		         tap/numeric_sw.f90\
		         tap/sw_post.f90\
			 tap/gr4_main.f90\
			 tap/mod_launch.f90\
			 tap/gr4_stsp.f90\
			 tap/crit_calc.f90

endif

   TAP_INVARS=dof0%h dof0%u dof0%v manning manning_beta bathy_cell bc%hyd%t bc%hyd%q bc%rat%h bc%rat%q bc%rain%t bc%rain%q infil%GA%Ks infil%GA%PsiF infil%GA%DeltaTheta infil%SCS%lambda infil%SCS%CN $(HYDROVAR)
   TAP_OUTVARS=dof%h dof%u dof%v innovation%diff innovW%diff innovQ%diff cost

ifeq ($(MODEL),0)
   FILES_TAP  =
   TAP_INVARS =
   TAP_OUTVARS=
endif


#======================================================================================================================#
#
#======================================================================================================================#

all: preproc $(executable) mes1
wrap: wrappers/_${wrapping_name}.so post_wrap mes2

doc:
	sphinx-build  ../doc/DOC_SPHINX/source  ../doc/DOC_SPHINX/build
#======================================================================================================================#
#  Preprocessor
#======================================================================================================================#

preproc: $(POST_SRC) $(POST_INC)

mes1:
	@echo "============================================================================================================"
	@echo "  DassFlow Compilation and Link Done"
	@echo "============================================================================================================"

#======================================================================================================================#
#
#======================================================================================================================#

$(executable): $(OBJS_MOD) $(OBJS_SRC)
	@echo "============================================================================================================"
	$(FC) $(FFLAGS)   -o $(CASEDIR)/$(executable) ./obj/*.o $(FLIBS)

#======================================================================================================================#
#  Tapenade
#======================================================================================================================#

tap_files: copy tgt adj post

tap_files_no_post: copy tgt adj

copy: $(FILES_TAP)

tgt: $(FILES_TAP)
	@echo "============================================================================================================"

	tapenade -d -fixinterface -splitoutputfiles \
	-tgtfuncname _diff -tgtvarname _diff -tgtmodulename _diff -copyname _c_diff \
	-head "run_model($(TAP_INVARS))\($(TAP_OUTVARS))" -O ./tap  $(FILES_TAP)

adj: $(FILES_TAP)
	@echo "============================================================================================================"

	tapenade -b -fixinterface -splitoutputfiles \
	-adjfuncname _back -adjvarname _back -adjmodulename _back -copyname _c_back \
	-head "run_model($(TAP_INVARS))\($(TAP_OUTVARS))" -O ./tap  $(FILES_TAP)


post:
	@echo "============================================================================================================"
	cp ./src/adjoint/finish_to_gen_adjoint.pl ./tap/
	cd ./tap ; perl finish_to_gen_adjoint.pl
	@echo "============================================================================================================"

#======================================================================================================================#
#  Dependencies
#======================================================================================================================#

./cpp/%.f90: %.f90
	@echo "============================================================================================================"
	@echo $(<)
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

./cpp/%.f90: $(CASEDIR)/%.f90
	@echo "============================================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

./cpp/%.inc: %.f90
	@echo "============================================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

./tap/%.f90: %.f90
	@echo "============================================================================================================"
	sed -e '/<NOADJ/,/>NOADJ/d;/NOADJ/d;/^ *!.*/d' $< | cpp -C -P $(CPP_FLAGS_ADJ) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

./tap/%.f90: $(CASEDIR)/%.f90
	@echo "============================================================================================================"
	sed -e '/<NOADJ/,/>NOADJ/d;/NOADJ/d;/^ *!.*/d' $< | cpp -C -P $(CPP_FLAGS_ADJ) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

./obj/%.o: ./cpp/%.f90 $(POST_INC)
	@echo "============================================================================================================"
	$(FC) $(CFLAGS) $(CLIBS) -c $< -o $@

./obj/%.o: ./tap/%.f90
	@echo "============================================================================================================"
	$(FC) $(CFLAGS) -c $< -o $@

./obj/%.o: %.f
	@echo "============================================================================================================"
	$(FC) $(CFLAGS) -c $< -o $@

./obj/%.o: %.c
	@echo "============================================================================================================"
	$(CC) -fPIC -O2 -c $< -o $@

wrappers:
	mkdir -p wrappers

wrappers/${wrapping_name}: wrappers
	mkdir -p wrappers/${wrapping_name}

wrappers/${wrapping_name}/utils:
	mkdir -p wrappers/${wrapping_name}/utils

wrappers/${wrapping_name}/post:
	mkdir -p wrappers/${wrapping_name}/post

wrappers/${wrapping_name}/post/%.py: src/wrappers/post/%.py wrappers/${wrapping_name}/post
	@echo "================================================================================"
	cp $< $@

wrappers/${wrapping_name}/utils/%.py: src/wrappers/utils/%.py wrappers/${wrapping_name}/utils
	@echo "================================================================================"
	cp $< $@



#======================================================================================================================#
#     Libraries
#======================================================================================================================#

lib:
ifeq ($(MPI),1)
	rm -rf $(SCOTCH)/src/Makefile.inc
ifeq ($(COMPILO),0)
	cd $(SCOTCH)/src ; ln -s Make.inc/Makefile.inc.x86-64_pc_linux2 Makefile.inc
else
	cd $(SCOTCH)/src ; ln -s Make.inc/Makefile.inc.x86-64_pc_linux2.icc Makefile.inc
endif
	cd $(PWD)/libs/scotch_5.1.12_esmumps/src ;\
	make
endif
ifeq ($(SOLVER),1)
	rm -rf $(SCOTCH)/src/Makefile.inc
ifeq ($(COMPILO),0)
	cd $(SCOTCH)/src ; ln -s Make.inc/Makefile.inc.x86-64_pc_linux2 Makefile.inc
else
	cd $(SCOTCH)/src ; ln -s Make.inc/Makefile.inc.x86-64_pc_linux2.icc Makefile.inc
endif
	cd $(PWD)/libs/scotch_5.1.12_esmumps/src ;\
	make
	cd $(PWD)/libs/MUMPS_4.10.0 ;\
	make
endif
ifeq ($(SOLVER),2)
	cd $(AGMG) ; make dseq
endif

cleanlib:
	cd $(PWD)/libs/scotch_5.1.12_esmumps/src ;\
	make clean ;\
	cd .. ;\
	rm -rf $(CASEDIR) lib include ;\
	cd $(PWD)/libs/MUMPS_4.10.0 ;\
	make clean



#======================================================================================================================#
#  Wrappers Python
#======================================================================================================================#

WRAPPERS_CPP 	=  wrappers/m_common.f90 \
		   wrappers/m_linear_algebra.f90 \
		   wrappers/m_mpi.f90 \
		   wrappers/m_mesh.f90 \
		   wrappers/m_sw_mono.f90 \
		   wrappers/input.f90

 ifeq ($(HYDRO),1)
OBJ_TMP = 	./obj/m_common.o \
			./obj/m_random.o \
			./obj/m_linear_algebra.o \
			./obj/m_mesh.o \
			./obj/m_mpi.o  \
		./obj/m_time_screen.o \
		./obj/m_gr4.o \
	  ./obj/m_sw_mono.o \
	  ./obj/sw_pre.o \
	  ./obj/m_linear_solver.o \
	  ./obj/m_numeric.o \
	  ./obj/m_obs.o \
	  ./obj/boundary.o \
	  ./obj/advance_time.o \
	  ./obj/numeric_sw.o \
	  ./obj/sw_flux.o \
	  ./obj/euler_time_step_first_b1.o \
	  ./obj/low_froude_time_step.o \
	  ./obj/friction_imex.o \
  	  ./obj/muscl_aud_flux_n.o \
	  ./obj/imex_time_step.o \
	  ./obj/sw_post.o \
	  ./obj/run_model_sw_mono.o \
	  ./obj/geometry.o \
	  ./obj/initialization.o \
	  ./obj/input.o \
	  ./obj/output_sw.o \
	  ./obj/output.o \
	  ./obj/read_write_obs.o
else
OBJ_TMP = 	./obj/m_common.o \
			./obj/m_random.o \
			./obj/m_linear_algebra.o \
			./obj/m_mesh.o \
			./obj/m_mpi.o  \
		./obj/m_time_screen.o \
	  ./obj/m_sw_mono.o \
	  ./obj/sw_pre.o \
	  ./obj/m_linear_solver.o \
	  ./obj/m_numeric.o \
	  ./obj/m_obs.o \
	  ./obj/boundary.o \
	  ./obj/advance_time.o \
	  ./obj/numeric_sw.o \
	  ./obj/sw_flux.o \
	  ./obj/euler_time_step_first_b1.o \
	  ./obj/low_froude_time_step.o \
	  ./obj/friction_imex.o \
  	  ./obj/muscl_aud_flux_n.o \
	  ./obj/imex_time_step.o \
	  ./obj/sw_post.o \
	  ./obj/run_model_sw_mono.o \
	  ./obj/geometry.o \
	  ./obj/initialization.o \
	  ./obj/input.o \
	  ./obj/output_sw.o \
	  ./obj/output.o \
	  ./obj/read_write_obs.o


endif

ifeq ($(ADJOINT),1)
	WRAPPERS_CPP += wrappers/m_adjoint.f90 \
					wrappers/m_tap_vars.f90
	OBJ_TMP += $(OBJS_TAP)
endif

WRAPPERS_CPP	+= 	 wrappers/call_run_model.f90
OBJ_TMP			+=	./obj/call_run_model.o

POST_PY  	= ./wrappers/${wrapping_name}/post/__init__.py \
           	  ./wrappers/${wrapping_name}/post/metrics.py \
           	  ./wrappers/${wrapping_name}/post/results.py \
           	  ./wrappers/${wrapping_name}/post/minimization.py \
           	  ./wrappers/${wrapping_name}/post/read_results.py

UTILS_PY 	= ./wrappers/${wrapping_name}/utils/__init__.py \
           	  ./wrappers/${wrapping_name}/utils/SwotObs.py \
           	  ./wrappers/${wrapping_name}/utils/PepsiNetCDF.py \
           	  ./wrappers/${wrapping_name}/utils/mesh_from_swot_observations.py \
           	  ./wrappers/${wrapping_name}/utils/pepsi_nodes_observations.py \

#ASSIM_PY	= ./wrappers/${wrapping_name}/assim/__init__.py \
#				./wrappers/$(wrapping_name)/assim/lbfgs.py \
#				./wrappers/$(wrapping_name)/assim/utils.py \
#				./wrappers/$(wrapping_name)/assim/lbfgs_mpi.py
#


WRAPPER_LIST    = $(patsubst wrappers/%.f90, %.f90, ${WRAPPERS_CPP})
WRAPPED_OBJ 	= $(patsubst ./obj/%.o, ../obj/%.o, $(OBJ_TMP))


wrappers/%.f90: src/common/%.f90 obj/%.o wrappers/kind_map
	@echo "============================================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS_WRAP) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

wrappers/%.f90: src/adjoint/%.f90
	@echo "============================================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS_WRAP) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

wrappers/%.f90: src/base/%.f90
	@echo "============================================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS_WRAP) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

wrappers/%.f90: src/sw_mono/%.f90
	@echo "============================================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS_WRAP) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

wrappers/%.f90: $(CASEDIR)/%.f90
	@echo "============================================================================================================"
	sed -e 's/^!.*//g;s/^ *#\(.*\)/#\1/g' $< | cpp -C -P $(CPP_FLAGS_WRAP) | sed -e '/^\/\*.*\*\/$$/d;/^\/\*/,/\*\/$$/d' > $@

wrappers/_${wrapping_name}.so: $(OBJ_TMP) $(WRAPPERS_CPP)
	cp ./src/wrappers/kind_map ./wrappers/
	cp ./src/wrappers/include.py ./wrappers/
	cp  ./src/wrappers/doc_plugin.py ./wrappers/
	cd wrappers ; f90wrap -m ${wrapping_name} ${WRAPPER_LIST} -k kind_map --documentation-plugin doc_plugin.py  --package

ifeq ($(MPI),1)
	cd wrappers ; f2py-f90wrap --f90exec=mpif90 --f77exec=gfortran ${debug_wrapping} --opt="-O3" --build-dir . -c -m _${wrapping_name} -I../obj -L. f90wrap*.f90 ../obj/*.o ${FLIBS}
else
	cd wrappers ; f2py-f90wrap --f90exec=gfortran --f77exec=gfortran ${debug_wrapping} --opt="-O3" --build-dir . -c -m _${wrapping_name} -I../obj -L. f90wrap*.f90 ../obj/*.o ${FLIBS}
endif

post_wrap: $(POST_PY) $(UTILS_PY)
# $(ASSIM_PY)

mes2:
	@echo "============================================================================================================"
	@echo "  DassFlow Wrappers Done"
	@echo "============================================================================================================"

#======================================================================================================================#
#     Copy files to /bin directory
#======================================================================================================================#

bin_files:
	cp -r ./simu/sw_mono/channel_U/* $(CASEDIR)

#======================================================================================================================#
#     Clean
#======================================================================================================================#

clean: cleanwrap
	rm -rf ./obj/* ./cpp/* $(CASEDIR)/exe $(CASEDIR)/input.post $(CASEDIR)/screen_output

cleanres:
	rm -rf $(CASEDIR)/res/*
	rm -rf $(CASEDIR)/graph/*

cleanmsh:
	rm -rf $(CASEDIR)/msh/*

cleanobs:
	rm -rf $(CASEDIR)/obs/*

cleanmin:
	rm -rf $(CASEDIR)/min/*

cleancase:
	rm -rf $(CASEDIR)/*

cleantap:
	rm -rf ./tap/*

cleanwrap:
	rm -rf ./wrappers/*
	rm -rf ./wrappers_true/*

cleanall: clean cleanlib cleantap
	rm -rf $(CASEDIR)/*


#======================================================================================================================#
#     Run
#======================================================================================================================#

rundirect: cleanres
ifeq ($(MPI),1)
	cp ./src/wrappers/main_direct_parallel.py $(CASEDIR)
	#cd $(CASEDIR) ; python3 main.py direct; rm main.py
	#cd $(CASEDIR) ; mpirun -np $(NB_PROC) python3 main.py direct; rm main.py
	cd $(CASEDIR) ; mpirun -np $(NB_PROC) python3 main_direct_parallel.py
else
	cp ./src/wrappers/main.py $(CASEDIR)
	cd $(CASEDIR) ; python3 main.py direct; rm main.py
endif

runtestadj:
ifeq ($(MPI),1)
	cp ./src/wrappers/main.py $(CASEDIR)
	cd $(CASEDIR) ; mpirun -np $(NB_PROC) python3 main.py testadj; rm main.py
else
	cp ./src/wrappers/main.py $(CASEDIR)
	cd $(CASEDIR) ; python3 main.py testadj; rm main.py
endif


rungrad:
ifeq ($(MPI),1)
	cp ./src/wrappers/main.py $(CASEDIR)
	cd $(CASEDIR) ; mpirun -np $(NB_PROC) python3 main.py grad; rm main.py
else
	cp ./src/wrappers/main.py $(CASEDIR)
	cd $(CASEDIR) ; python3 main.py grad; rm main.py
endif

runmin:
ifeq ($(MPI),1)
	cp ./src/wrappers/main_min_parallel.py $(CASEDIR)
	#cd $(CASEDIR) ; python3 main.py min; rm main.py
	# cd $(CASEDIR) ; mpirun -np $(NB_PROC) python3 main.py min; rm main.py
	cd $(CASEDIR) ; mpirun -np $(NB_PROC) python3 main_min_parallel.py
else
	cp ./src/wrappers/main.py $(CASEDIR)
	cd $(CASEDIR) ; python3 main.py min; rm main.py
endif


runtwin:
ifeq ($(MPI),1)
	@echo " RUNTWIN NOT IMPLEMENTED IN PARRALLEL   $(ARG_HYDRO) $(ARG_bathy) $(ARG_MANNING_ALPHA)"
else
	@echo "RUN TWIN EXPERIMENT  ARG_HYDRO=$(ARG_HYDRO) ARG_BATHY=$(ARG_BATHY) ARG_MANNING_ALPHA=$(ARG_MANNING_ALPHA)"
	cp ./src/wrappers/main.py $(CASEDIR)
	cd $(CASEDIR) ; python3 runtwin.py $(ARG_HYDRO) $(ARG_BATHY) $(ARG_MANNING_ALPHA); rm main.py
endif

runlbfgs:
ifeq ($(MPI),1)
	cp ./src/wrappers/main.py $(CASEDIR)
	cd $(CASEDIR) ; mpirun -np $(NB_PROC) python3 main.py lbfgs; rm main.py
else
	cp ./src/wrappers/main.py $(CASEDIR)
	cd $(CASEDIR) ; python3 main.py lbfgs; rm main.py
endif

testwrap:
ifeq ($(MPI),1)
	cp ./src/wrappers/main.py $(CASEDIR)
	cd $(CASEDIR) ; mpirun -np $(NB_PROC) python3 main.py; rm main.py
else
	cp ./src/wrappers/main.py $(CASEDIR)
	cd $(CASEDIR) ; python3 main.py; rm main.py
endif

testlbfgs:
	cp test/test_lbfgs.py $(CASEDIR)
ifeq ($(MPI),1)
	cd $(CASEDIR) ; mpirun -np $(NB_PROC) python3 test_lbfgs.py 1; rm test_lbfgs.py
else
	cd $(CASEDIR) ; python3 test_lbfgs.py 0; rm test_lbfgs.py
endif

testdirect:
	cp test/test_direct.py $(CASEDIR)
	cd $(CASEDIR) ; python3 test_direct.py ; rm test_direct.py

testadjoint:
	cp test/test_adjoint.py $(CASEDIR)
	cd $(CASEDIR) ; python3 test_adjoint.py ; rm test_adjoint.py

testread:
	cp test/test_read_package.py $(CASEDIR)
	cd $(CASEDIR) ; python3 test_read_package.py ; rm test_read_package.py

testcontrolmpi:
	cp test/test_control_mpi.py $(CASEDIR)
	cd $(CASEDIR) ; mpirun -np $(NB_PROC) python3 test_control_mpi.py; rm test_control_mpi.py

#======================================================================================================================#
#     Gathered make command
#					- regenerates ADJOINT no matter makefile options (still can make direct simulation without problem)
#								-> the problem can be solved usin manually make command below, and skip make tap_files step.
#					- activate MPI if required in adjoint
#				 - generate fotran/python I/O (f90wrap-f2py), the pluggin "doc_pluggin.py" allows for python parsing for documentation
#         			--> in wrappers
#				 - finaly generates dassflow2d package gathering necessary information.
#							--> in wrappers_true
#======================================================================================================================#
lilian_dir = src/wrappers/package
module_name= dassflow2d
true_wrapper_dir = wrappers_true
install:
	make cleanlib
	make cleantap
	make cleanwrap
	make clean
	sleep 2
	make lib
	make tap_files
	make wrap
	make package


package:
		# create new directory, where we can write everything as necessary without modifying existing Compilation
		mkdir -p $(true_wrapper_dir)                               		# base working dir (necessary to story the setup.py)
		mkdir -p $(true_wrapper_dir)/$(module_name)								 		# within there is all python stuff
		mkdir -p $(true_wrapper_dir)/$(module_name)/$(wrapping_name)  # here we store generated files by f90wrap
		mkdir -p $(true_wrapper_dir)/$(module_name)/core  # here we store generated files by f90wrap
			mkdir -p $(true_wrapper_dir)/$(module_name)/assim  # here we store generated files by f90wrap
		# COPY WRAPPED NECCESSARY STUFF
		cp  wrappers/_$(wrapping_name).*.so  $(true_wrapper_dir)/$(module_name)/$(wrapping_name)
		cp -r  wrappers/$(wrapping_name)/*  $(true_wrapper_dir)/$(module_name)/$(wrapping_name)

		# Replace f90wrap automatic generation
			#first remove
		cp $(true_wrapper_dir)/$(module_name)/$(wrapping_name)/__init__.py $(true_wrapper_dir)/$(module_name)/$(wrapping_name)/save_init.py
		rm $(true_wrapper_dir)/$(module_name)/$(wrapping_name)/__init__.py

			#then copy wished file
		cp $(lilian_dir)/init_files/wrapping/__init__.py $(true_wrapper_dir)/$(module_name)/$(wrapping_name)  # copy   init file of "wrapped" from sources (AND NOT FROM AUTOMATIC GENERATION)
		cp $(lilian_dir)/init_files/main/__init__.py $(true_wrapper_dir)/$(module_name)												# copy, init file of main module fril sources
		cp -r $(lilian_dir)/core/*  $(true_wrapper_dir)/$(module_name)/core																		# copy main module
		cp -r $(lilian_dir)/assim/*  $(true_wrapper_dir)/$(module_name)/assim																	# copy assim module


		# add necessary files for package generation
		cp $(lilian_dir)/setup.py  $(true_wrapper_dir)

		# change in wrapping dir the  imports of modules  /////
		for file in m_mesh.py m_mpi.py m_model.py m_linear_algebra.py  m_common.py call_model.py  m_adjoint.py m_tap_vars.py; do \
				echo $$file;\
				sed -i -e "s/from wrapping/from dassflow2d.wrapping/g" 	$(true_wrapper_dir)/$(module_name)/$(wrapping_name)/$$file ; \
	  done

# replace  _wrappring by  dassflow2d.wrapping._wrapping   (absolute call)
		for file in m_mesh.py m_mpi.py m_model.py m_linear_algebra.py  m_common.py call_model.py  m_adjoint.py m_tap_vars.py; do \
				echo $$file;\
				sed -i -e "s/import _wrapping/from dassflow2d.wrapping import _wrapping/g" 	$(true_wrapper_dir)/$(module_name)/$(wrapping_name)/$$file ; \
	  done

			# COMPILE package
		pip install --no-compile $(true_wrapper_dir)/.
