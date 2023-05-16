# encoding: utf-8
from __future__ import division, print_function, unicode_literals
import os
import matplotlib
if not "DISPLAY" in os.environ:
    matplotlib.use('Agg')
matplotlib.rc('text', usetex=True)
import numpy as np
import matplotlib.pyplot as plt


def plot_iteration_Q(iteration, item, t, control, target=None, prior=None, t_unit=None, tobs=None, river_name=None, outfile=None, external=[]):
  
    # CHECK-UP and type cast
    if not isinstance(t, np.ndarray):
        t = np.array(t)
    if river_name is not None:
        if not isinstance(river_name, str):
            raise ValueError("'river_name' must be a string")
    if outfile is not None:
        if not isinstance(outfile, str):
            raise ValueError("'outfile' must be a string")
    
    # Compute units and multipliers
    if t_unit is None:
        if np.max(t) < 3600:
            t_unit = "s"
        elif np.max(t) < 172800:
            t_unit = "h"
        else:
            t_unit = "days"
    if t_unit == "s":
        t_mult = 1.0
    elif t_unit == "h":
        t_mult = 1.0 / 3600.0
    else:
        t_mult = 1.0 / 86400.0
    # Retrieve control in "real" (physical) space
    current = np.zeros(control.x.size)
    control.get_real_control(current)
    
    # Initialise figure and subplots
    fig, ax1 = plt.subplots(1, 1, sharex=False)
    
    # Plot discharge
    ax1.plot(t*t_mult, target[control.get_item_slice(item)], "r.", label="target")
    if prior is not None:
        ax1.plot(t*t_mult, prior[control.get_item_slice(item)], "k--", label="prior")
    ax1.plot(t*t_mult, current[control.get_item_slice(item)], "b-", label="infered")
    
    for (q_tuple, label) in external:
        t, q = q_tuple
        ax1.plot(t*t_mult, q, "--", label=label)
        
    
    ax1.set_xlabel(r"$t$ $(%s)$" % t_unit)
    ax1.set_ylabel(r"$Q$ $(m^3.s^{-1})$")
    
    if tobs is not None:
        if isinstance(tobs, list):
            for i, tpass in enumerate(tobs):
                for t in tpass:
                    ax1.axvline(x=t*t_mult, color="C{}".format(i), linestyle="-.")
    
    ax1.legend()
    
    qtarget = current[control.get_item_slice(item)]
    qcurrent = target[control.get_item_slice(item)]
    rmse = np.sqrt(np.sum((qcurrent-qtarget)**2) / qtarget.size)
    nrmse = rmse / np.mean(qtarget)
    nse = 1.0 - np.sum((qcurrent-qtarget)**2) / np.sum((qtarget-np.mean(qtarget))**2)
    
    # Set main title
    if iteration < 1:
        iteration_in_title = "final"
    else:
        iteration_in_title = "iteration %i" % iteration
    if river_name is not None:
      fig.suptitle("%s - %s (RMSE=%f, NRMSE=%f, NSE=%f)" % (river_name, iteration_in_title, rmse, nrmse, nse))
    else:
      fig.suptitle("%s (RMSE=%f, NRMSE=%f, NSE=%f)" % (iteration_in_title, rmse, nrmse, nse))
      
    # Adjust plot layout
    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    
    # Show or save figure
    if outfile:
        plt.savefig(outfile)
        plt.close(fig)
    else:
        plt.show()


def plot_iteration_QKhb(iteration, t, xs, control, target=None, prior=None, t_unit=None, xs_unit=None, river_name=None, outfile=None):
  
    # CHECK-UP and type cast
    if not isinstance(t, np.ndarray):
        t = np.array(t)
    if not isinstance(xs, np.ndarray):
        xs = np.array(xs)
    if river_name is not None:
        if not isinstance(river_name, str):
            raise ValueError("'river_name' must be a string")
    if outfile is not None:
        if not isinstance(outfile, str):
            raise ValueError("'outfile' must be a string")
    
    # Compute units and multipliers
    if t_unit is None:
        if np.max(t) < 3600:
            t_unit = "s"
        elif np.max(t) < 172800:
            t_unit = "h"
        else:
            t_unit = "days"
    if t_unit == "s":
        t_mult = 1.0
    elif t_unit == "h":
        t_mult = 1.0 / 3600.0
    else:
        t_mult = 1.0 / 86400.0
    if xs_unit is None:
        if np.max(xs) < 2000:
            xs_unit = "m"
        else:
            xs_unit = "km"
    if xs_unit == "m":
        xs_mult = 1.0
    else:
        xs_mult = 0.001
      
    # Retrieve control in "real" (physical) space
    current = np.zeros(control.x.size)
    control.get_real_control(current)
    
    # Initialise figure and subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=False)
    
    # Plot discharge
    item_index = 0
    if control.items[item_index].id[0:2] == b"BC":
        if target[control.get_item_slice(item_index)].size == t.size:
            ax1.plot(t*t_mult, target[control.get_item_slice(item_index)], "r.", label="target")
            if prior is not None:
                ax1.plot(t*t_mult, prior[control.get_item_slice(item_index)], "k--", label="prior")
            ax1.plot(t*t_mult, current[control.get_item_slice(item_index)], "b-", label="infered")
            ax1.set_xlabel(r"$t$ $(%s)$" % t_unit)
            ax1.set_ylabel(r"$Q$ $(m^3.s^{-1})$")
            ax1.legend()
            item_index +=1
        
    # Plot alpha and beta
    if item_index < len(control.items):
        print(control.items[item_index].id[0:4], control.items[item_index].id[0:4] == "KPAR")
        if control.items[item_index].id[0:4] == b"KPAR":
            if prior is not None:
                ax2.plot(xs*xs_mult, prior[control.get_item_slice(item_index)], "b--")
            ax2.plot(xs*xs_mult, current[control.get_item_slice(item_index)], "b-")
            ax2.set_xlabel(r"$x$ $(%s)$" % xs_unit)
            ax2.set_ylabel(r"$\alpha$ $(m^{1/3}.s^{-1})$")
            ax2.yaxis.label.set_color('b')
            ax2bis = ax2.twinx()
            ax2bis.plot(xs*xs_mult, prior[control.get_item_slice(item_index+1)], "r--")
            ax2bis.plot(xs*xs_mult, current[control.get_item_slice(item_index+1)], "r-")
            ax2bis.set_ylabel(r"$\beta$ $(-)$")
            ax2bis.yaxis.label.set_color('r')
            item_index += 2
    
    # Plot bathy
    if item_index < len(control.items):
        if control.items[item_index].id[0:5] == b"BATHY":
            if prior is not None:
                ax3.plot(xs*xs_mult, prior[control.get_item_slice(item_index)], "k--", label="prior")
            ax3.plot(xs*xs_mult, current[control.get_item_slice(item_index)], "b-", label="current")
            ax3.set_xlabel(r"$x$ $(%s)$" % xs_unit)
            ax3.set_ylabel(r"$b$ $(m)$")
            if prior is not None:
                ax3bis = ax3.twinx()
                ax3bis.plot(xs*xs_mult, current[control.get_item_slice(item_index)] - prior[control.get_item_slice(item_index)], "g--")
                ax3bis.set_ylabel(r"$b-b_0$ $(m)$")
                ax3bis.yaxis.label.set_color('g')
            ax3.legend()
    
    # Set main title
    if iteration < 1:
        iteration_in_title = "final"
    else:
        iteration_in_title = "iteration %i" % iteration
    if river_name is not None:
      fig.suptitle("%s - %s" % (river_name, iteration_in_title))
    else:
      fig.suptitle("%s" % iteration_in_title)
      
    # Adjust plot layout
    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    
    # Show or save figure
    if outfile:
        plt.savefig(outfile)
        plt.close(fig)
    else:
        plt.show()
