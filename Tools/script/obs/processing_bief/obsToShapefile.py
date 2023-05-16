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

numberofobs=22
for i in range(numberofobs):
   #CheminDossierSortie= './station_'+str(i)
   CheminDossierSortie='/home/brisset/Documents/IMT/DassFlow/SWOT_like_data/ouputs/shapefile'
   #Creation du fichier de sortie s'il n'existe pas
   if (os.path.isdir(CheminDossierSortie)==False):
      os.mkdir(CheminDossierSortie)


   a="{0:04d}".format(i+1)
   data=np.loadtxt('/home/brisset/Documents/IMT/DassFlow/SWOT_like_data/ouputs/station_grp_'+str(a)+'.txt',skiprows=1)

   X=data[:,0]
   Y=data[:,1]
   ebathy=data[:,1]


   #	
   # Creation du nouveau shapefile des points
   driver = ogr.GetDriverByName('ESRI Shapefile')
   fshpout=CheminDossierSortie
	
   #Suppression de l'ancien shapefile
   #if os.path.exists(fshpout):
   #   driver.DeleteDataSource(fshpout)
   dataout = driver.CreateDataSource(fshpout)
   proj = osr.SpatialReference()

   proj.ImportFromEPSG(2154) # 4326 = EPSG code for lon/lat WGS84 projection
   layerout = dataout.CreateLayer('obs_'+str(a),proj,geom_type=ogr.wkbPoint)  
   fieldDef1 = ogr.FieldDefn('h_eau', ogr.OFTReal)
   layerout.CreateField(fieldDef1)
   floutDefn = layerout.GetLayerDefn()
   feature_out = ogr.Feature(floutDefn)

   for ipix in range(len(X)):
   #- Create Geometry Point with pixel coordinates
      if float(ebathy[ipix])>0:
         pixel_point = ogr.Geometry(ogr.wkbPoint)
         pixel_point.AddPoint(X[ipix],Y[ipix])
         #- Add the geometry to the feature
         feature_out.SetGeometry(pixel_point)
         #- Set feature attributes
         feature_out.SetField('h_eau', float(ebathy[ipix]))
         layerout.CreateFeature(feature_out)
         #- Delete point geometry
         pixel_point.Destroy()

   feature_out.Destroy()
   dataout.Destroy()
