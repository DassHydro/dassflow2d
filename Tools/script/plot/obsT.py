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

# For each t, plot H in function of numero station

##########
#To change
##########

stations=range(10)

##########
#End to change
##########




pathPlot='./bin/plot/'
pathOutFile=pathPlot+'./obs/'
T,H=[],[]


#Creation of outputfile 
if (os.path.isdir(pathPlot)==False):
   os.mkdir(pathPlot)

if (os.path.isdir(pathOutFile)==False):
   os.mkdir(pathOutFile)


for i in stations:

   a="{0:04d}".format(i+1)
   data=np.loadtxt('./bin/obs/obs_station_'+str(a)+'.plt',skiprows=1)

   t=data[:,0]
   h=data[:,1]
   T.append(t)
   H.append(h)



x=stations
indexTime=range(len(t))
[Ymin,Ymax]=computeYminYmax(H,0.1)


for i in indexTime:
   plt.clf()
   y=[]   
   for j in x:
      h=H[j]
      y.append(h[i])
   time=t[i]
   plt.ylim(Ymin ,Ymax)
   plt.plot(x,y,'b-x')
   plt.xlabel('stations')
   plt.ylabel('$\\bar{H}$ ($m$)')
   titlefig='$\\bar{H}$ observed by stations'
   plt.title(titlefig)
   pathFile=pathOutFile+'fig_'+ str(round(time,2))+'.png'
   plt.savefig(pathFile)

os.system('eog '+str(pathFile) + '&' ) 

