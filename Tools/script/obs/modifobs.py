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


import matplotlib
matplotlib.use('Agg',warn=False)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

numberofobs=21
fmt='%15.8E'
header='# time h_mean u_mean v_mean w_mean'
for i in range(numberofobs):
   #CheminDossierSortie= './station_'+str(i)
   CheminDossierSortie='./bin/obs_shapefile'
   #Creation du fichier de sortie s'il n'existe pas
   if (os.path.isdir(CheminDossierSortie)==False):
      os.mkdir(CheminDossierSortie)


   a="{0:04d}".format(i+1)
   files='./bin/res/obs_station_'+str(a)+'.dat'
   data=np.loadtxt(files,skiprows=1)
   #print data
   #print len(data)
   time=data[:,0]
   #print time
   h=data[:,1]
   u=data[:,2]
   v=data[:,3]
   w=data[:,4]
   
   hnonnul=np.where(h!=0.0)
   
   noise=np.random.normal(0,0.01,len(h))
   h[hnonnul]=h[hnonnul]+0.3+noise[hnonnul]
   w[hnonnul]=w[hnonnul]+0.3+noise[hnonnul]

   data[:,1]=h
   data[:,4]=w
   
   np.savetxt(files,data,header=header,fmt=fmt,comments='')
