#=======================================================#
# Source librairies
#=======================================================#
import dassflow2d as df2d
import numpy as np
import sys
import os



run_type = sys.argv[1] # if runtype = rundirect runmin, classic treatment, if runtype = runminpython



print(run_type)

#=======================================================#
#  run fortran model
#=======================================================#
# run_type= "direct"
if run_type == "min" or run_type == "direct":

	df2d.wrapping.m_mpi.init_mpi()

	# store main path
	dassflow_dir="/home/leo/DISTANT/dassflow2d"
	code_dir =  f"{dassflow_dir}/code"
	bin_dir = f"{code_dir}/bin_A"

	##########
	# Initialise bin
	##########

	os.chdir(code_dir)

	os.system("make cleanres")			   # removes all in bin_dir/res directory
	os.system("make cleanmsh")			   # removes all in bin_dir/msh directory
	os.system("make cleanmin")			   # removes all in bin_dir/min directory

	if os.path.isfile(f"rm {bin_dir}/restart.bin"):
		os.system(f"rm {bin_dir}/restart.bin")   # removes all in bin_dir/msh directory

	os.chdir(bin_dir)

	# ------------ Define Default values------------------

	df2d.wrapping.read_input(f"{bin_dir}/input.txt")

	model = df2d.dassflowmodel(bin_dir =  bin_dir, hdf5_path = f"{bin_dir}/res/simu.hdf5" , run_type = "direct", clean = True, custom_config = None)

	model.init_mesh()

	model.kernel.dof  = df2d.wrapping.m_model.unk(model.kernel.mesh)
	model.kernel.dof0 = model.kernel.dof     

	df2d.wrapping.call_model.init_fortran(model.kernel)
	df2d.wrapping.call_model.run(model.kernel, arg = run_type)
	df2d.wrapping.call_model.clean_model(model.kernel)  

	
	
elif run_type == "lbfgs":

	run_type = "min"
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

			print("#=========================================================#")
			print(" ite " + str(nite) + ", J = " + str(newj) + ", |grad J| = " + str(norm_new_gradj) )
			print(" New control = "+str(control))		
			print("#=========================================================#\n")
			
			nite += 1
			res_cost = (newj, control_back)
			
			return res_cost 
	    
	# this is necessary for both mpi and non mpi cases in order to initialise some variables
	df2d.wrapping.m_mpi.init_mpi()
	# get the rank and number of processors
	rank = df2d.wrapping.m_mpi.get_proc()
	nproc = df2d.wrapping.m_mpi.get_np()
	mpi = [rank, nproc]


	model2 = df2d.dassflowmodel(bin_dir = os.getcwd(), run_type = run_type)
	# then intialise meshing
	model2.init_mesh()
	# build pyvista grid object
	model2.build_grid()
	# initialise dof structure
	model2.init_dof()	
	# initialise remaining structures
	model2.init_fortran()
	
	
	print(model2.kernel.dof.h[:])
	print(model2.kernel.mesh.nc)
	######### perform inference python
	df2d.wrapping.call_model.init_back(model2.kernel)   

	import glob
	from scipy.optimize import minimize

	list_ite = glob.glob('min/*')

	#latest = max(list_ite, key=os.path.getmtime).split('/')[1]
	#print(f" move file from {os.getcwd()}/min/{latest}  TO {os.getcwd()}/min/lbfgs/{latest}")
	#os.rename(f"{os.getcwd()}/min/{latest}", f"{os.getcwd()}/min/lbfgs/{latest}" )

	opts = {"disp":True, \
			"gtol" : df2d.wrapping.m_common.get_eps_min(), \
			"maxfun": df2d.wrapping.m_common.get_max_nt_for_adjoint(), \
			"maxiter": 20,#2,  #df2d.m_common.get_restart_min()
			"maxls" : 20, \
		"iprint":5
				}
	minmethod = "L-BFGS-B" #"L-BFGS-B" #"Nelder-Mead"

		# create dass_func instance
	#dass_func = df2d.lbfgs.utils.dass_func(model = model2)
	#init_values = model2.my_friction.manning * 1.5 #(1+np.random())*
	#init_values = model2.get_array_hydrograph() ?? * 1.5 #(1+np.random())
	init_values = df2d.wrapping.m_adjoint.get_array_control()

	print("Minimization algorithm: " + str(minmethod))
	print("Minimization options: " + str(opts))
	print("Initial control vector = ", init_values)


	runmin_called = runmin(model2.kernel, init_values = init_values)

	nite = 0
	res_optim1 = minimize(fun=runmin_called,  x0 = init_values, jac = True, method = minmethod, bounds=None, options= opts)
	
	
