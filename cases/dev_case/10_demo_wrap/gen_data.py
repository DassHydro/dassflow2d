#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 11:29:26 2022

@author: livillenave
"""

import numpy as np



g= 10
lx= 1000
q=2
n=0.033


def h_exact(x, g, lx):
    res = (4/g)**(1/3) * ( 1 + 0.5 * np.exp(-16 *(x/lx - 0.5 )**2 ) )
    return(res)
  
def z_exact(x, g, q, lx, ks):    
    for i in range(len(x)-1) :
            z_template[i+1] = z_template[i] +  (q**2/(g*h_exact(x[1], g=g, lx=lx)**3)) * (h_exact(x[i+1], g=g, lx=lx)-h_exact(x[i], g=g, lx=lx)) -                (abs(x[i+1]- x[i])) * q**2/ks**2
    return(z_template)
    
    
x = (np.arange(100)+0.5) *10


h_ex = h_exact(x, g=g, lx=lx)

zb_ex = z_exact(x, g=g,q=q, lx=lx, n = n)
zb_ex = zb_ex - min(zb_ex)
zs_ex = h_ex + zb_ex




test = x[0:2]

z_template = x
z_template[:] = 10


def z_exact(x, g, q, lx):
        
    
    
plt.scatter(x = x,
            y = h_ex)
plt.scatter(x = x,
            y = zb_ex)
plt.scatter(x = x,
            y = zs_ex)





















b=a["spatial_var/bathy"][:,i] 
zs = a["spatial_var/h"][:,i] +  a["spatial_var/bathy"][:,0] 
h = a["spatial_var/h"][:,i]
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.scatter(x = np.arange(len(zs )),
                    y = zs, 
                    linewidths = 1,
                    label = "zs")
ax1.scatter(x = np.arange(len(b )),
                    y = b, 
                    c="red", 
                    linewidths = 1,
                    label = "bathymetry" )
ax2.scatter(x = np.arange(len(h)),
                    y = h, 
                    c="purple", 
                    linewidths = 1,
                    label = "h" ) 
ax2.scatter(x = np.arange(len(h)),
                    y = h_exact(x=x, g=g, lx=lx), 
                    c="green", 
                    linewidths = 1,
                    label = "h_TRUE" )

ax1.set_xlabel('X [m]')
ax1.set_ylabel('elevation (m)', color='g')
ax2.set_ylabel('Water height [m]', color='b')
ax1.legend()
ax2.legend()