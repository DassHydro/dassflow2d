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

#   Reminder of makefile commands
#======================================================================================================================#

# Open console in the directory of this file
# To call make, type "make +commands". Any number of the following commands can be added. They will be executed in order.
# Example : "make clean cleanres cleanmin runexe"

#clean               : removes the exe file in the bin directory
#cleanres, cleanplot : deletes all in the /res and /plot folder in the bin directory
#cleanmin            : deletes all in the /min folder in the bin directory + deletes restart.bin
#cleantap            : deletes all in the /tap folder
#tap_files           : generates files for the /tap folder
#lib                 : generates libraries (depends on user parameters, see below)

#rundirect  : performs a direct mode run
#runmin     : performs an inverse mode run
#runtestadj : performs a test of the adjoint model
#rungrad    : computes the cost gradient

#======================================================================================================================#
#   Model and Compilation Options
#======================================================================================================================#

#Path = 		name of the case directory, default is bin
CASEDIR	= ./bin_A

#Fortran compiler: 	'0' for GNU, '1' for Intel
COMPILO	= 0
#Debug mode: 		'1' to activate. Not activating when computing adjoint model may cause errors.
OPTIM		= 1
#SW model: 		'0' to not compile, '1' to compile SW_MONO, '2' to compile SW_HB, '3' to compile SW_MULTI
MODEL		= 1
#Hydrological model: 	'0' to not compile, '1' to compile
HYDRO		= 0
#Adjoint model: 	'0' to not compile, '1' to compile
ADJOINT	= 1
#Minimization method: 	'1' to compile M1QN3, '2' to compile LBFGSB3
MINMETHOD	= 1
#Solver: 		'1' for MUMPS (recommended), '2' for AGMG, '0' for no recompilation of external libraries
SOLVER		= 0

#Parallel computing: 	'1' to use
MPI		=0
#If MPI=1, number of processors to use
NB_PROC	= 0
DEBUG_WRAPPING = 0

