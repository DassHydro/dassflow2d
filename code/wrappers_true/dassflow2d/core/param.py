import dassflow2d
import h5py
import numpy as np


class Friction(object):
        """
        Friction class
        
        contains:
            - self.manning.mesh_value(size=mesh.nc)
            - self.manning.patch_correspondance(size=mesh.nc)
            - self.manning.patch_value(size=nb_patch)
            
#            - self.manning_beta.mesh_value(size=mesh.nc)
#            - self.manning_beta.mesh_correspondance(size=mesh.nc)
#            - self.manning_beta.patch_correspondance(size=nb_patch)
#            - self.manning_beta.patch_value(size=nb_patch)
    
        """            
        def __init__(self):
            self.manning  = dict()
            
        def __force_format_manning__(self):
            """
            force np array of minimal dimention 1 and not 0 
            return : adapted array
            """
            self.manning["mesh_correspondance"]=np.array(self.manning["mesh_correspondance"], ndmin=1,   order = "F", dtype = "float64")
            self.manning["patch_correspondance"]=np.array(self.manning["patch_correspondance"], ndmin=1, order = "F", dtype = "float64")
            self.manning["patch_value"]=np.array(self.manning["patch_value"],                   ndmin=1, order = "F", dtype = "float64")
            
        def __is_FortranPython_same_size__(self):
            """
            check size between fortran and python object
            """
            fkernel_manning = dassflow2d.wrapping.m_model.manning
            py_manning = self.manning["patch_value"]
            
            return(len(fkernel_manning) == len(py_manning))
            
            
            
        def get_manning(self):
            pointeur     = dassflow2d.wrapping.m_model.get_array_land()
            self.manning["mesh_correspondance"] = pointeur.copy()
            
            pointeur     = np.arange(len(dassflow2d.wrapping.m_model.get_array_manning()))+1
            self.manning["patch_correspondance"]    = pointeur.copy()
            
            
            pointeur     = dassflow2d.wrapping.m_model.get_array_manning()
            self.manning["patch_value"]             = pointeur.copy()
            # todo mesh manning value
            
#            self.manning_beta                           = dict()
#            self.manning_beta["mesh_correspondance"]    = dassflow2d.wrapping.m_model.get_array_land()
#            self.manning_beta["patch_correspondance"]   = np.arange(len(dassflow2d.wrapping.m_model.get_array_manning_beta()))+1
#            self.manning_beta["patch_value"]            = dassflow2d.wrapping.m_model.get_array_manning_beta()
#            # todo mesh manning_beta value
        
        def set_manning(self):
            
            self.__force_format_manning__()
            buffer_manning =  self.manning.copy()
            
            if not self.__is_FortranPython_same_size__(): # if python values are not same size of fortran kernel, reallocate fortran
                dassflow2d.wrapping.call_model.reallocate_manning(len(self.manning["patch_value"]))            
                dassflow2d.wrapping.m_model.set_nland(len(self.manning["patch_value"]) )
                self.get_manning() # load new manning value in kernel, new size is uploaded but no values
            
            
            self.manning = buffer_manning.copy()
            # set manning values in fortran kernel
            dassflow2d.wrapping.m_model.manning[...]  = np.asanyarray(self.manning["patch_value"],   order = "F", dtype = "float64")
            #dassflow2d.wrapping.m_model.set_array_manning(self.manning["patch_value"])
            # set land use correspondance in fortran kernel
            dassflow2d.wrapping.m_model.land[...] = np.asanyarray(self.manning["mesh_correspondance"],   order = "F", dtype = "float64")
            #dassflow2d.wrapping.m_model.set_array_land( self.manning["mesh_correspondance"] )
            
            # self.manning["patch_correspondance"]    = np.arange(len(self.manning["patch_value"]))+1
            # self.manning["mesh_value"]     = dassflow2d.wrapping.m_model.get_array_land()

        def save_manning(self, hdf5_path):           
                  
            f = h5py.File(hdf5_path, "a") 
            
            try:    # create input group where will be the meshing data
                        input_group = f.create_group("input")
            except: # if exist, cant be created, then it is loaded
                        input_group = f["input"]
                        
            try:    # create input group where will be the meshing data
                        param_group = input_group.create_group("param")
            except: # if exist, cant be created, then it is loaded
                        param_group = input_group["param"]
                        
            manning_group = param_group.create_group("manning")
            manning_group.create_dataset("patch_value", (len(self.manning["patch_value"])), 'f')
            manning_group.create_dataset("mesh_correspondance", (len(self.manning["mesh_correspondance"])), 'f')
            
            manning_group["patch_value"][:] = self.manning["patch_value"][:]
            manning_group["mesh_correspondance"][:] = self.manning["mesh_correspondance"][:]
            f.close()
            
class Param(dict):
    """
    Parameters data (dictionary)
    
    all key correspond to each parameter accessible.
    
    

    Parameters
    ----------
    input_param: dictionary of all configuration parameters that can be given for a simulation:
        **todo ref sphinx vers input.txt**
    """
    
    def __init__(self):        
        self.friction = Friction()


    def get(self):        
    	""""
    	enable get from fortran kernel for all classes of parameter  availaible
    	"""
    	self.friction.get_manning()       
        

    def set(self):      
    	""""
    	enable set to fortran kernel for all classes of parameter  availaible
    	"""     
    	self.friction.set_manning()
        

    def save(self, hdf5_path):      
    	""""
    	save in hdf5 format all classes of parameter  availaible
    	"""         
    	self.friction.save_manning(hdf5_path=hdf5_path)
