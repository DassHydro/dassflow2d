
import dassflow2d as df2d
# We build the function that enable the  acces to all configuration values defined in fortran kernel:
def get_config(input_param={ "mesh_name":'channel.geo',
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
                                            "eps_min":0,
                                            "bc_infil":0,
                                            "bc_rain":0,
                                            "use_Zobs":0,
                                            "use_Qobs":0,
                                            "use_NSE":0,
                                            "xsshp_along_x":1,
                                            "xsshp_along_y":0,
                                            "use_xsshp":0,
                                            "c_Ks":0,
                                            "c_PsiF":0,
                                            "c_DeltaTheta":0,
                                            "c_shape_s":0,
                                            "c_hmax":0,
                                            "c_xcenter":0,
                                            "c_slope_x":0,
                                            "c_slope_y":0} ):

    res =input_param
#    all_keys = input_param.keys()
    for k in res:
        # in module m_common
       if k == 'mesh_name':
             res[k] = df2d.wrapping.m_common.get_mesh_name().decode('utf-8')
       elif k == 'ts':
             res[k] = df2d.wrapping.m_common.get_ts()
       elif k == 'dtw':
               print("dtw")
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


       #in m_common linked to infiltration
       elif k == 'bc_infil':
               res[k] = df2d.wrapping.m_common.get_bc_infil()
       elif k == 'bc_rain':
               res[k] = df2d.wrapping.m_common.get_bc_rain()
       elif k == 'use_Zobs':
               res[k] = df2d.wrapping.m_common.get_use_zobs()
       elif k == 'use_Qobs':
               res[k] = df2d.wrapping.m_common.get_use_qobs()
       elif k == 'use_NSE':
               res[k] = df2d.wrapping.m_common.get_use_nse()

       #in m_common linked to XSparams
       elif k == 'use_xsshp':
               res[k] = df2d.wrapping.m_common.get_use_xsshp()
       elif k == 'xsshp_along_x':
               res[k] = df2d.wrapping.m_common.get_xsshp_along_x()
       elif k == 'xsshp_along_y':
               res[k] = df2d.wrapping.m_common.get_xsshp_along_y()

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
       elif k == 'c_shape_s':
               res[k] = df2d.wrapping.m_model.get_c_shape_s()
       elif k == 'c_hmax':
               res[k] = df2d.wrapping.m_model.get_c_hmax()
       elif k == 'c_xcenter':
               res[k] = df2d.wrapping.m_model.get_c_xcenter()
       elif k == 'c_slope_x':
               res[k] = df2d.wrapping.m_model.get_c_slope_x()
       elif k == 'c_slope_y':
               res[k] = df2d.wrapping.m_model.get_c_slope_y()

         #in m_model linked to infiltration
       elif k == 'c_Ks':
               res[k] = df2d.wrapping.m_model.get_c_ks()
       elif k == 'c_PsiF':
               res[k] = df2d.wrapping.m_model.get_c_psif()
       elif k == 'c_DeltaTheta':
               res[k] = df2d.wrapping.m_model.get_c_deltatheta()
              
               


    return(res)

# update values in fortran kernel based on provided input
# @param: input_param :: dictionary with keys corresponging to fields to be filled in fortran + values set
# @return: nothing : set the values in fortran kernel
# df2d istance must be on
def set_config(input_param):
#    all_keys = input_param.keys()
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
       elif k == 'use_Qobs':
                df2d.wrapping.m_common.set_use_qobs(input_param[k])
       elif k == 'use_NSE':
                df2d.wrapping.m_common.set_use_nse(input_param[k])
                
              #in m_common linked to XSparams
       elif k == 'use_xsshp':
               df2d.wrapping.m_common.set_use_xsshp(input_param[k])
       elif k == 'xsshp_along_x':
               df2d.wrapping.m_common.set_xsshp_along_x(input_param[k])
       elif k == 'xsshp_along_y':
               df2d.wrapping.m_common.set_xsshp_along_y(input_param[k])

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
               df2d.wrapping.m_model.set_c_rain()
       elif k == 'c_ic':
               df2d.wrapping.m_model.set_c_ic()    
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
                



    print(f"values {[x for x in input_param.keys()]} set")
    return()
