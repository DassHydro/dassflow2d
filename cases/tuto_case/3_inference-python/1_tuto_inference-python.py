##########################################################
##########################################################
# PERFORM AN INFERENCE RUN USING PYTHON OTPIMIZER (scipy.optimize.minimize, with minmethod= "Newton-CG" )
# Chanel mcdonald?
#
# Use python optimizer !
##########################################################
##########################################################


#===============================================================#
#  Set up python intefare
#===============================================================#

#-------------------#
#  Source librairies
#--------------------#
import dassflow2d as df2d              # df2d : main package
import os                        # os, for shell command execution (mainly for file manipulation)
import shutil
import matplotlib.pyplot as plt
import numpy as np
import h5py
from mpi4py import MPI
#=============================================================================================

# define a class : runmin
# access fortran values for the cost fucntion and its gradient (J, grad J)
# when called, access and return values of J, gradJ in the solver
class runmin():
	def __init__(self, model, init_values):
		#initiliase using the model that is defined using input.txt
		self.model = model
		self.dim_all = df2d.wrapping.m_adjoint.get_dim_all()
		self.jac = np.zeros(self.dim_all, dtype='float')
		self.init_values = init_values
						
	def __call__(self, x):
	
		global nite
		global init_gradj
		
		cost = df2d.wrapping.call_model.func(self = self.model, ctrl_in = x, grad_func = self.jac)  

		#print(nite)
		newj = df2d.wrapping.m_adjoint.get_cost()
		#print("New J " +str(newj))
		


		control_back = self.jac #df2d.m_adjoint.get_array_control_back()
		if nite == 0:
			init_gradj = np.sum(np.square(control_back))**0.5
			#print("self.jac",self.jac)
			#print("init_gradj",init_gradj)
					
		control_back = [i/init_gradj for i in control_back]
		norm_new_gradj = np.sum(np.square(control_back))**0.5
		
		control = df2d.wrapping.m_adjoint.get_array_control()

		plt.plot(range(0,len(control)), control, label = "Ite "+str(nite))
		
		print("#=========================================================#")
		print(" ite " + str(nite) + ", J = " + str(newj) + ", |grad J| = " + str(norm_new_gradj) )
		print(" New control = "+str(control))		
		print("#=========================================================#\n")
		
		nite += 1
		res_cost = (newj, control_back)
		
		return res_cost 
    
#============================================================================================= 


#  Define Parameters 
#----------------------#

dassflow_dir = os.path.abspath(os.path.join(__file__ ,"../../../.."))

case_dir = os.path.join(f"{dassflow_dir}","cases/tuto_case/3_inference-python")
bin_dir = os.path.join(f"{dassflow_dir}","code/bin_A")

print("Dassflow directory is understood as: ", dassflow_dir)
print("Case is copied from: ", case_dir)
print("Running directory is: ", bin_dir)


# delete all files in your simulation directory before starting
os.system(f"rm -r {dassflow_dir}/code/bin_A/*")
# Copy recursively the files provided in DassFlow case repository into your own simulation directory **code/bin_A/**.
os.system(f"cp -r {case_dir}/bin_A/* {dassflow_dir}/code/bin_A")

os.chdir( f"{dassflow_dir}/code/")
os.system(f"make cleanres cleanmin ")

#===============================================================#
#  MAKE INFERENCE
#===============================================================#

# this is necessary for both mpi and non mpi cases in order to initialise some variables
df2d.wrapping.m_mpi.init_mpi()
# get the rank and number of processors
rank = df2d.wrapping.m_mpi.get_proc()
nproc = df2d.wrapping.m_mpi.get_np()
mpi = [rank, nproc]

os.chdir(bin_dir)
df2d.wrapping.read_input( bin_dir + "/input.txt")
my_model = df2d.wrapping.call_model.Model()

my_model.mesh = df2d.wrapping.m_mesh.msh()
my_model.dof = df2d.wrapping.m_model.unk(my_model.mesh)
my_model.dof0 = df2d.wrapping.m_model.unk(my_model.mesh)
df2d.wrapping.call_model.init_solver(my_model)
df2d.wrapping.call_model.init_fortran(my_model)         # UPDATE IN FORTRAN

######### perform inference python
df2d.wrapping.call_model.init_back(my_model)   

import glob
from scipy.optimize import minimize

list_ite = glob.glob('min/*')

#latest = max(list_ite, key=os.path.getmtime).split('/')[1]
#print(f" move file from {os.getcwd()}/min/{latest}  TO {os.getcwd()}/min/lbfgs/{latest}")
#os.rename(f"{os.getcwd()}/min/{latest}", f"{os.getcwd()}/min/lbfgs/{latest}" )

opts = {"disp":True, \
		"gtol" : df2d.wrapping.m_common.get_eps_min(), \
		"maxfun": df2d.wrapping.m_common.get_max_nt_for_adjoint(), \
		"maxiter": 15,#2,  #df2d.m_common.get_restart_min()
		"maxls" : 5, \
        "iprint":5
			}
minmethod = "Newton-CG" #"L-BFGS-B" #"Nelder-Mead"

	# create dass_func instance
#dass_func = df2d.lbfgs.utils.dass_func(model = my_model)

#init_values = my_model.my_friction.manning * 1.5 #(1+np.random())*
#init_values = my_model.get_array_hydrograph() ?? * 1.5 #(1+np.random())
init_values = df2d.wrapping.m_adjoint.get_array_control()

print("Minimization algorithm: " + str(minmethod))
print("Minimization options: " + str(opts))
print("Initial control vector = ", init_values)

plt.plot(range(0,len(init_values)), init_values, label = "First Guess")

runmin_called = runmin(my_model, init_values = init_values)

nite = 0
res_optim1 = minimize(fun=runmin_called,  x0 = init_values, jac = True, method = minmethod, bounds=None, options= opts)

plt.legend()
plt.show()


