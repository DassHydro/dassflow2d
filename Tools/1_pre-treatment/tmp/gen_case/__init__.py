from __future__ import print_function, absolute_import, division
from . import f90wrap_toplevel
import f90wrap.runtime
import logging

def gen_basic_channel(nx, ny, lx, ly):
    """
    gen_basic_channel(nx, ny, lx, ly)
    
    
    Defined at _gen_channel_case.f90 lines 1-117
    
    Parameters
    ----------
    nx : int
    ny : int
    lx : float
    ly : float
    
    ===================================================================================================================
      SW model parameters
    ===================================================================================================================
    """
    _gen_case.f90wrap_gen_basic_channel(nx=nx, ny=ny, lx=lx, ly=ly)

def gen_bc(in_type, out_type):
    """
    gen_bc(in_type, out_type)
    
    
    Defined at _gen_channel_case.f90 lines 119-132
    
    Parameters
    ----------
    in_type : str
    out_type : str
    
    =================================================================================================
    '
    """
    _gen_case.f90wrap_gen_bc(in_type=in_type, out_type=out_type)

def gen_land_use(manning_alpha, manning_beta):
    """
    gen_land_use(manning_alpha, manning_beta)
    
    
    Defined at _gen_channel_case.f90 lines 134-146
    
    Parameters
    ----------
    manning_alpha : float
    manning_beta : float
    
    =================================================================================================
    '
    """
    _gen_case.f90wrap_gen_land_use(manning_alpha=manning_alpha, \
        manning_beta=manning_beta)

def gen_bc_data(bc_typ, nrow, var1, var2):
    """
    gen_bc_data(bc_typ, nrow, var1, var2)
    
    
    Defined at _gen_channel_case.f90 lines 188-216
    
    Parameters
    ----------
    bc_typ : str
    nrow : int
    var1 : float array
    var2 : float array
    
    =================================================================================================
    '
    """
    _gen_case.f90wrap_gen_bc_data(bc_typ=bc_typ, nrow=nrow, var1=var1, var2=var2)

def gen_obs(nx_obs, ny_obs, xmax_obs, ymax_obs, xmin_obs, ymin_obs, dt_obs):
    """
    gen_obs(nx_obs, ny_obs, xmax_obs, ymax_obs, xmin_obs, ymin_obs, dt_obs)
    
    
    Defined at _gen_channel_case.f90 lines 218-252
    
    Parameters
    ----------
    nx_obs : int
    ny_obs : int
    xmax_obs : float
    ymax_obs : float
    xmin_obs : float
    ymin_obs : float
    dt_obs : float
    
    =================================================================================================
    '
    """
    _gen_case.f90wrap_gen_obs(nx_obs=nx_obs, ny_obs=ny_obs, xmax_obs=xmax_obs, \
        ymax_obs=ymax_obs, xmin_obs=xmin_obs, ymin_obs=ymin_obs, dt_obs=dt_obs)

def h_true_macdo(x, lx, g):
    """
    h_true = h_true_macdo(x, lx, g)
    
    
    Defined at _gen_channel_case.f90 lines 254-257
    
    Parameters
    ----------
    x : float
    lx : float
    g : float
    
    Returns
    -------
    h_true : float
    
    """
    h_true = _gen_case.f90wrap_h_true_macdo(x=x, lx=lx, g=g)
    return h_true

