import dassflow2d as df2d
import h5py

class Config(dict):
       """
       configuration data (dictionary)

       all key correspond to a getter and a setter in fortran kernel

       To see all possible key 


       Parameters
       ----------
       input_param: dictionary of all configuration parameters that can be given for a simulation:
              **todo ref sphinx vers input.txt**

       """
    
       def __init__(self) :
              """
              set default values for dictionary
              WARNING: they differ from fortran kernel default value



              Parameters
              ----------
              input_param: dictionary of all configuration parameters that can be given for a simulation:
              **todo ref sphinx vers input.txt**
              """

              self = { "mesh_name":'tofill.geo',
                                                 "ts":0,
                                                 "dta":0,
                                                 "dtw":0,
                                                 "dtp":0,
                                                 "dt":0,
                                                 "temp_scheme":'euler',
                                                 "spatial_scheme":'first_b1',
                                                 "adapt_dt":0,
                                                 "cfl":0.,
                                                 "feedback_inflow":0,
                                                 "coef_feedback":0.,
                                                 "heps":0,
                                                 "friction":0,
                                                 "g":10,
                                                 "w_tecplot":0,
                                                 "w_vtk":0,
                                                 "w_gnuplot":1,
                                                 "w_obs":0,
                                                 "use_obs":0,
                                                 "max_nt_for_adjoint":0,
                                                 "c_manning":0,
                                                 "c_manning_beta":0,
                                                 "c_bathy":0,
                                                 "c_hydrograph":0,
                                                 "c_ratcurve":0,
                                                 "c_rain":0,
                                                 #"c_infil":0,
                                                 "c_ic":0,
                                                 "restart_min":0,
                                                 "eps_min":0}

              
              
       # get fortran kernel values
       def get(self, input_param = None):
              """
              Fill dictionary from fortran kernel values
              """
              if not input_param is None:
                     res = input_param
              else:
                     res = { "mesh_name":'channel.geo',
                                                 "ts":0,
                                                 "dta":0,
                                                 "dtw":0,
                                                 "dtp":0,
                                                 "dt":0,
                                                 "temp_scheme":'euler',
                                                 "spatial_scheme":'first_b1',
                                                 "adapt_dt":0,
                                                 "cfl":0.,
                                                 "feedback_inflow":0,
                                                 "coef_feedback":0.,
                                                 "heps":0,
                                                 "friction":0,
                                                 "g":10,
                                                 "w_tecplot":0,
                                                 "w_vtk":0,
                                                 "w_gnuplot":0,
                                                 "w_obs":0,
                                                 "use_obs":0,
                                                 "max_nt_for_adjoint":0,
                                                 "c_manning":0,
                                                 "c_manning_beta":0,
                                                 "c_bathy":0,
                                                 "c_hydrograph":0,
                                                 "c_ratcurve":0,
                                                 "c_rain":0,
                                                 #"c_infil":0,
                                                 "c_ic":0,
                                                 "restart_min":0,
                                                 "eps_min":0}
              for k in res: 
                     # in module m_common
                     if k == 'mesh_name':
                            res[k] = df2d.wrapping.m_common.get_mesh_name().decode('utf-8')
                     elif k == 'ts':
                            res[k] = df2d.wrapping.m_common.get_ts()
                     elif k == 'dtw':
                            res[k] = df2d.wrapping.m_common.get_dtw()
                     elif k == 'dtp':
                            res[k] = df2d.wrapping.m_common.get_dtp()
                     elif k == 'adapt_dt':
                            res[k] = df2d.wrapping.m_common.get_adapt_dt()
                     elif k == 'dt':
                            res[k] = df2d.wrapping.m_common.get_dt()
                     elif k == 'cfl':
                            res[k] = df2d.wrapping.m_common.get_cfl()
                     elif k == 'temp_scheme':
                            res[k] = df2d.wrapping.m_common.get_temp_scheme().decode('utf-8')
                     elif k == 'spatial_scheme':
                            res[k] = df2d.wrapping.m_common.get_spatial_scheme().decode('utf-8')
                     elif k == 'w_tecplot':
                            res[k] = df2d.wrapping.m_common.get_w_tecplot()
                     elif k == 'w_gnuplot':
                            res[k] = df2d.wrapping.m_common.get_w_gnuplot()
                     elif k == 'w_vtk':
                            res[k] = df2d.wrapping.m_common.get_w_vtk()
                     elif k == 'w_obs':
                            res[k] = df2d.wrapping.m_common.get_w_obs()
                     elif k == 'use_obs':
                            res[k] = df2d.wrapping.m_common.get_use_obs()
                     elif k == 'max_nt_for_adjoint':
                            res[k] = df2d.wrapping.m_common.get_max_nt_for_adjoint()
                     elif k == 'max_nt_for_direct':
                            res[k] = df2d.wrapping.m_common.get_max_nt_for_direct()
                     elif k == 'restart_min':
                            res[k] = df2d.wrapping.m_common.get_restart_min()
                     elif k == 'get_eps_min':
                            res[k] = df2d.wrapping.m_common.get_eps_min()

                     #in module m_model
                     elif k == 'feedback_inflow':
                            res[k] = df2d.wrapping.m_model.get_feedback_inflow()
                     elif k == 'coef_feedback':
                            res[k] = df2d.wrapping.m_model.get_coef_feedback()
                     elif k == 'friction':
                            res[k] = df2d.wrapping.m_model.get_friction()
                     elif k == 'g':
                            res[k] = df2d.wrapping.m_model.get_g()
                     elif k == 'c_manning':
                            res[k] =  df2d.wrapping.m_model.get_c_manning()
                     elif k == 'c_manning_beta':
                            res[k] = df2d.wrapping.m_model.get_c_manning_beta()
                     elif k == 'c_bathy':
                            df2d.wrapping.m_model.get_c_bathy()
                     elif k == 'c_hydrograph':
                            res[k] = df2d.wrapping.m_model.get_c_hydrograph()
                     elif k == 'c_rain':
                            res[k] = df2d.wrapping.m_model.get_c_rain()
                     elif k == 'c_ic':
                            res[k] = df2d.wrapping.m_model.get_c_ic()
                     self[k] = res[k]

              self = res
              return(self)





       def set(self, custom_config = None):
              print(custom_config)
              print(self)
              """
              update values in fortran kernel based on provided input
              ---------------
              param: 
              ---------------
              self : dictionary of class Config

              ---------------
              return: 
              ---------------
              set the values in fortran kernel

              fortran kernel values must be initialise 
              """
              if  custom_config is None:
                     input_param = self
              else:
                     input_param = custom_config

              for k in input_param:
              # in module m_common
                     if k == 'mesh_name':
                            df2d.wrapping.m_common.set_mesh_name(input_param[k])
                     elif k == 'ts':
                            df2d.wrapping.m_common.set_ts(input_param[k])
                            # DTA IN ASSIMILATION FILE (OBSOLETE HERE ?)
                     #       if k == 'dta':
                     #               input_param[k] = df2d.wrapping.m_common.set_dta()
                     elif k == 'dtw':
                            df2d.wrapping.m_common.set_dtw(input_param[k])
                     elif k == 'dtp':
                            df2d.wrapping.m_common.set_dtp(input_param[k])
                     elif k == 'adapt_dt':
                            df2d.wrapping.m_common.set_adapt_dt(input_param[k])
                     elif k == 'dt':
                            df2d.wrapping.m_common.set_dt(input_param[k])
                     elif k == 'cfl':
                            df2d.wrapping.m_common.set_cfl(input_param[k])
                     elif k == 'temp_scheme':
                            df2d.wrapping.m_common.set_temp_scheme(input_param[k])
                     elif k == 'spatial_scheme':
                            df2d.wrapping.m_common.set_spatial_scheme(input_param[k])
                     elif k == 'w_tecplot':
                            df2d.wrapping.m_common.set_w_tecplot(input_param[k])
                     elif k == 'w_gnuplot':
                            df2d.wrapping.m_common.set_w_gnuplot(input_param[k])
                     elif k == 'w_vtk':
                            df2d.wrapping.m_common.set_w_vtk(input_param[k])
                     elif k == 'w_obs':
                            df2d.wrapping.m_common.set_w_obs(input_param[k])
                     elif k == 'use_obs':
                            df2d.wrapping.m_common.set_use_obs(input_param[k])
                     elif k == 'max_nt_for_adjoint':
                            df2d.wrapping.m_common.set_max_nt_for_adjoint(input_param[k])
                     elif k == 'max_nt_for_direct':
                            df2d.wrapping.m_common.set_max_nt_for_direct(input_param[k])
                     elif k == 'restart_min':
                            df2d.wrapping.m_common.set_restart_min(input_param[k])
                     elif k == 'set_eps_min':
                            df2d.wrapping.m_common.set_eps_min(input_param[k])

                     #in m_common linked to infiltration
                     elif k == 'bc_infil':
                            df2d.wrapping.m_common.set_bc_infil(input_param[k])
                     elif k == 'bc_rain':
                            df2d.wrapping.m_common.set_bc_rain(input_param[k])
                     elif k == 'use_Zobs':
                            df2d.wrapping.m_common.set_use_zobs(input_param[k])
                     elif k == 'use_hobs':
                            df2d.wrapping.m_common.set_use_hobs(input_param[k])
                     elif k == 'use_Qobs':
                            df2d.wrapping.m_common.set_use_qobs(input_param[k])
                     elif k == 'use_UVobs':
                            df2d.wrapping.m_common.set_use_uvobs(input_param[k])
                     elif k == 'use_NSE':
                            df2d.wrapping.m_common.set_use_nse(input_param[k])
                            
                     #in m_common linked to XSparams
                     elif k == 'use_xsshp':
                            df2d.wrapping.m_common.set_use_xsshp(input_param[k])
                     elif k == 'xsshp_along_x':
                            df2d.wrapping.m_common.set_xsshp_along_x(input_param[k])
                     elif k == 'xsshp_along_y':
                            df2d.wrapping.m_common.set_xsshp_along_y(input_param[k])
                            
                     elif k == 'use_ptf':
                            df2d.wrapping.m_common.set_use_ptf(input_param[k])

                     #in module m_model
                     elif k == 'feedback_inflow':
                            df2d.wrapping.m_model.set_feedback_inflow(input_param[k])
                     elif k == 'coef_feedback':
                            df2d.wrapping.m_model.set_coef_feedback(input_param[k])
                     elif k == 'friction':
                            df2d.wrapping.m_model.set_friction(input_param[k])
                     elif k == 'g':
                            df2d.wrapping.m_model.set_g(input_param[k])
                     elif k == 'c_manning':
                            df2d.wrapping.m_model.set_c_manning(input_param[k])
                     elif k == 'c_manning_beta':
                            df2d.wrapping.m_model.set_c_manning_beta(input_param[k])
                     elif k == 'c_bathy':
                            df2d.wrapping.m_model.set_c_bathy(input_param[k])
                     elif k == 'c_hydrograph':
                            df2d.wrapping.m_model.set_c_hydrograph(input_param[k])
                     elif k == 'c_rain':
                            df2d.wrapping.m_model.set_c_rain(input_param[k])
                     elif k == 'c_ic':
                            df2d.wrapping.m_model.set_c_ic(input_param[k])    
                     elif k == 'c_shape_s':
                            df2d.wrapping.m_model.set_c_shape_s(input_param[k])
                     elif k == 'c_hmax':
                            df2d.wrapping.m_model.set_c_hmax(input_param[k])
                     elif k == 'c_xcenter':
                            df2d.wrapping.m_model.set_c_xcenter(input_param[k])
                     elif k == 'c_slope_x':
                            df2d.wrapping.m_model.set_c_slope_x(input_param[k])
                     elif k == 'c_slope_x':
                            df2d.wrapping.m_model.set_c_slope_y(input_param[k])
                     elif k == 'regul_bathy':
                            df2d.wrapping.m_model.set_regul_bathy(input_param[k])
                     elif k == 'regul_bathy_grad':
                            df2d.wrapping.m_model.set_regul_bathy_grad(input_param[k])
                     elif k == 'regul_bathy_shape':
                            df2d.wrapping.m_model.set_regul_bathy_shape(input_param[k])

                            #in m_model linked to infiltration
                     elif k == 'c_Ks':
                            df2d.wrapping.m_model.set_c_ks(input_param[k])
                     elif k == 'c_PsiF':
                            df2d.wrapping.m_model.set_c_psif(input_param[k])
                     elif k == 'c_DeltaTheta':
                            df2d.wrapping.m_model.set_c_deltatheta(input_param[k])
                            
                     elif k == 'c_ptf':
                            df2d.wrapping.m_model.set_c_ptf(input_param[k])

              print(f"values {[x for x in input_param.keys()]} set")
              self.get(input_param)
              return()
              
       def save(self, hdf5_path, source_kernel = False ):
              """
              save in hdf5 file specified
              
              saves as an attribute at the root of the hdf5 file.
              
              all key correspond to a getter and a setter in fortran kernel
              
              To see all possible key: ``self.keys()``


              Parameters
              ----------
              
              hdf5_path: str, path to the hdf5 file to save configuration data
              
              source_kernel: boolean, if true update config from fortrna kernel value before save

              """
              if source_kernel:
                     self.get()
              
              f = h5py.File(hdf5_path, "r+")        
              for key,values in self.items() :
                     f.attrs[key]=values
              f.close()
              
       def soure_hdf5(self, hdf5_path ):
              """
              Fill dictionary from  hdf5 filevalues
              Parameters
              ----------
              hdf5_path: str, path to the hdf5 file to load configuration data
              """
              
              self ={ "mesh_name":'channel.geo',
                                                 "ts":0,
                                                 "dta":0,
                                                 "dtw":0,
                                                 "dtp":0,
                                                 "dt":0,
                                                 "temp_scheme":'euler',
                                                 "spatial_scheme":'first_b1',
                                                 "adapt_dt":0,
                                                 "cfl":0.,
                                                 "feedback_inflow":0,
                                                 "coef_feedback":0.,
                                                 "heps":0,
                                                 "friction":0,
                                                 "g":10,
                                                 "w_tecplot":0,
                                                 "w_vtk":0,
                                                 "w_gnuplot":0,
                                                 "w_obs":0,
                                                 "use_obs":0,
                                                 "max_nt_for_adjoint":0,
                                                 "c_manning":0,
                                                 "c_manning_beta":0,
                                                 "c_bathy":0,
                                                 "c_hydrograph":0,
                                                 "c_ratcurve":0,
                                                 "c_rain":0,
                                                 #"c_infil":0,
                                                 "c_ic":0,
                                                 "restart_min":0,
                                                 "eps_min":0}
              f = h5py.File(hdf5_path, "r")
              
              for key in self.keys() :
                     self[key] = f.attrs[key]
              f.close()
