#!/usr/bin/env python
# encoding: utf-8

#from __future__ import unicode_literals
#from __future__ import print_function
#from __future__ import division
import math
import os, sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from osgeo import ogr,osr # GDAL library
from osgeo.gdalconst import * # GDAL constants
import gdal
from gdalconst import *
from function import *

import matplotlib
matplotlib.use('Agg',warn=False)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl


pathPlot='./bin/plot/'
pathOutFile=pathPlot+'./grad/'

legendTime='hours'

titlefig='Hydrogramme gradient'
titlefig=''



pathGrad='./bin/grad/hydrograph001_grad'

if (os.path.isdir(pathPlot)==False):
   os.mkdir(pathPlot)

if (os.path.isdir(pathOutFile)==False):
   os.mkdir(pathOutFile)


factorTime=factorTimeFromLegend(legendTime)
t,q=np.loadtxt(pathGrad,unpack=True)
plt.plot(t*factorTime,q,'b',linewidth=2)

plt.xlabel('$Time$ $( '+str(legendTime)+ ')$')


plt.title(titlefig)


#plt.legend()


nameFile=pathOutFile+'plotGrad.png'
plt.savefig(nameFile)
os.system('eog '+str(nameFile) + '&' ) 
