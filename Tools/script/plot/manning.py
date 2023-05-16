#!/usr/bin/env python
# encoding: utf-8

#from __future__ import unicode_literals
#from __future__ import print_function
#from __future__ import division
#from StringIO import StringIO
import math
import os, sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from osgeo import ogr,osr # GDAL library
from osgeo.gdalconst import * # GDAL constants
import gdal
from gdalconst import *


import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl


variable_name='manning'
input_file='land_uses'


bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/"

pathOutFile=f'{bin_dir}/plot/min/{variable_name}'
pathTarget=f'{bin_dir}/{input_file}_target.txt'
pathFirstGuess=f'{bin_dir}/{input_file}.txt'

legend=True
legendTime='seconds' # days seconds
plotMin=True
plotIntermediaire=True


if legendTime=='hours' :
   factorTime=1./3600
elif legendTime=='days':
   factorTime=1./(3600*24)
else :
   factorTime=1.

if(variable_name == 'manning'):
    landuse,target, beta_target =np.loadtxt(pathTarget,skiprows=7,unpack=True)
#if(variable_name == 'hydrograph'):
#    landuse,target, beta_target =np.loadtxt(pathTarget,skiprows=8,unpack=True)


#Creation of outputfile 
if (os.path.isdir(pathOutFile)==False):
   os.mkdir(pathOutFile)


#
# First plot : hydrograph.png
#
numberIte=0
iteT,manningT=[],[]
i=0
stop=True

while stop :
   a="{0:03d}".format(i)
   try :
      numberIte=numberIte+1
      landuse,manning=np.loadtxt(f'{bin_dir}/min/{variable_name}.{a}',unpack=True)
      iteT.append(numberIte-1)
      manningT.append(manning)
   except :
      stop=False
   i=i+1


i=i-1

test = np.array(manningT)
for j in range(np.shape(test)[1]):
    a="{0:03d}".format(i+1)
    titlefig = 'Evolution of $K$ with iterations of minimization'     
    plt.plot([0,iteT[-1]],[test[0,j], test[0,j]],label='First guess')
    plt.plot([0,iteT[-1]],[target[j],target[j]],label='Target')
    plt.plot(iteT,test[:,j],label=f'cell {j}')
        
    plt.ylabel('$K$ $(I.S)$')
    plt.xlabel('iterations')
    plt.legend()
    #
    #manningT.append(target)
    #deltay=max(manningT)-min(manningT)
    #ylimT=max(manningT)+deltay*0.05
    #ylimB=min(manningT)-deltay*0.05
    #plt.ylim(ylimB,ylimT)
    
    
    
    c=f'{pathOutFile}/manning_{j}.png'
    plt.savefig(c)
    plt.close()
   # os.system('eog '+str(c) + '&' ) 


