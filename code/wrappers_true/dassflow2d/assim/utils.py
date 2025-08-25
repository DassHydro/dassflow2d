import numpy as np
import dassflow2d as df2d

class dass_func():
	def __init__(self, model):
		#initiliase using the model that is defined using input.txt
		self.model = model
		self.dim_all = df2d.m_adjoint.get_dim_all()
		self.jac = np.zeros(self.dim_all, dtype='float')
						
	def __call__(self, x ):

		cost = df2d.call_model.func(self = self.model, ctrl_in = x, grad_func = self.jac) 
		return cost, self.jac


class dass_callback():
	def __init__(self, model):
		self.model = model
		self.ite = df2d.m_adjoint.get_ite_min()
		
	def __call__(self,x):
		self.ite += 1
		df2d.m_adjoint.set_ite_min(self.ite)
		# parameter x is not used : 
		# 1. the last call of dass_func() will fill the control vector with new values
		# 2. the control vector is saved in a global variable and inside the subroutine output_control, we fill the dof0 and the mesh with the new values of control vector
		# 3. the instance of type(Model) in this class will thus be updated
		
		# check if x and control is the same ? 
		df2d.m_adjoint.output_control(self.model.dof0, self.model.mesh)
		
		if(self.ite>20):
			return(True)
		else:
			return(x)
