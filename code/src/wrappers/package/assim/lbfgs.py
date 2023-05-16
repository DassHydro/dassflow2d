import dassflow2d
from . import utils

from scipy.optimize import minimize

import numpy as np
import numpy.linalg as npl

import glob
import os

from mpi4py import MPI

dassflow2d.m_mpi.init_mpi()
rank = dassflow2d.m_mpi.get_proc()
np = dassflow2d.m_mpi.get_np()

class dass_callback_lbfgs(utils.dass_callback):
	def __init__(self,model):
		super().__init__(model)
	def __call__(self,x):
		super().__call__(x)
		
		if rank == 0:
			#move the latest iteration to the directory of lbfgs
			list_ite = glob.glob('min/*')
			latest = max(list_ite, key=os.path.getmtime).split('/')[1]
			print("move file %s to %s"%(latest, 'min/lbfgs'))
			os.rename('min/'+latest, "min/lbfgs/"+latest)
			
			cost = dassflow2d.m_adjoint.get_cost()
			jac_norm = npl.norm(dassflow2d.m_adjoint.get_array_control_back())
			# save the cost function and the gradient norm
			
			with open('min/lbfgs/min_cost','a') as f:
				f.write("%4d %15.7f %15.7f \n" %(self.ite, cost, jac_norm))
			
		


def run(model):

	dassflow2d.wrapping.call_model.init_back(model)   

	# create directory for lbfgs iterations
	try:
		os.makedirs("min/lbfgs")
	except FileExistsError:
		pass
			
	# get the initial guess
	x0 = dassflow2d.m_adjoint.get_array_control()

	# output the first guess
	dassflow2d.m_adjoint.set_ite_min(0)
	dassflow2d.m_adjoint.output_control(model.dof0, model.mesh)

	if rank == 0:
		list_ite = glob.glob('min/*')
		print(list_ite)
		latest = max(list_ite, key=os.path.getmtime).split('/')[1]
		print("move file %s to %s"%(latest, '/min/lbfgs'))
		os.rename(f"{os.getcwd()}/min/{latest}", f"{os.getcwd()}/min/save/{latest}" )
		


	# build the lbfgs options. Cf the documentation
	# "maxcor" : , \ # this parameter is new to lbfgs
	opts = {"disp":True, \
			"gtol" : dassflow2d.m_common.get_eps_min(), \
			"maxfun": dassflow2d.m_common.get_max_nt_for_adjoint(), \
			"maxiter": dassflow2d.m_common.get_restart_min(), \
			"maxls" : 15, \
			}
			
	print("minimization options")
	print(opts)

	# create dass_func instance
	func = dassflow2d.lbfgs.utils.dass_func(model = model)

	callback =  dassflow2d.assim.lbfgs.dass_callback_lbfgs(model)



	# define bounds if you want to use L-BFGS-B
	# the format of bounds
	res_optim = minimize(fun=func, x0 = x0, jac = True, method = 'L-BFGS-B', bounds=None, options= opts, callback = callback)
	
	
	
	
	return res_optim
	
