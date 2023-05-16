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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

from function import *



bin_dir = "/home/livillenave/Documents/distant/dassflow2d-wrap/code/bin_A/"



outputDirectoryQgis=f'{bin_dir}/mesh/'
outputDirectory=outputDirectoryQgis+'stations'


pathobstxt=f'{bin_dir}/obs.txt'



obsfile=open(pathobstxt,'r')
ObsFileLine=obsfile.readlines()
import re 
ObsFileLine = [re.sub("\n", "", x) for x in ObsFileLine]
ObsFileLine[6].split()
j=0
Toread=''
numberLineToRead=0
print( len(ObsFileLine) )

Xobs,Yobs=[],[]



#Recupération des stations
stop=True
i=0
while stop:
   line=ObsFileLine[i]
   print(line)
   if ('stations' in line):# and 'stations_with_grp' not in line) :
      line=cleanString(line)
      line=line.replace('stations','')
      numberStation=int(line)
      compt=numberStation
      while compt>0:
         i=i+1
         line=ObsFileLine[i]
         print(line)
         line_temp=cleanString(line)
         if line_temp!='':
            line=line.split()
            Xobs.append(float(line[0]))
            Yobs.append(float(line[1]))
            compt=compt-1

      stop=False
   i=i+1


#Recupération des stations_with_grp
stop=True
i=0
while stop:
   line=ObsFileLine[i]
   if ('stations_with_grp' in line):
      line=cleanString(line)
      line=line.replace('stations_with_grp','')
      numberofobs=int(line)
      stop=False
   i=i+1





#Creation of output directory Qgis
if os.path.isdir(outputDirectoryQgis)==False:
	os.mkdir(outputDirectoryQgis)

#Creation of output directory shapefile
if os.path.isdir(outputDirectory)==False:
	os.mkdir(outputDirectory)


driver = ogr.GetDriverByName('ESRI Shapefile')



for i in range(numberStation):
   
   a="{0:04d}".format(i+1)

	

   dataout = driver.CreateDataSource(outputDirectory)
   proj = osr.SpatialReference()

   proj.ImportFromEPSG(2154) # 4326 = EPSG code for lon/lat WGS84 projection
   newlayerout = dataout.CreateLayer('obs_'+str(a),proj,geom_type=ogr.wkbPoint)  
   fieldDef1 = ogr.FieldDefn('h_eau', ogr.OFTReal)
   newlayerout.CreateField(fieldDef1)
   floutDefn = newlayerout.GetLayerDefn()
   feature_out = ogr.Feature(floutDefn)


   pixel_point = ogr.Geometry(ogr.wkbPoint)
   pixel_point.AddPoint(Xobs[i],Yobs[i])
   feature_out.SetGeometry(pixel_point)
   feature_out.SetField('h_eau', float(Xobs[i]))
   newlayerout.CreateFeature(feature_out)
   #newlayerout.Destroy()

   feature_out.Destroy()
   dataout.Destroy()

for i in range(numberofobs):


   a="{0:04d}".format(i+1)
   data=np.loadtxt('./bin/station_grp_'+str(a)+'.txt',skiprows=1)

   X=data[:,0]
   Y=data[:,1]
   ebathy=data[:,1]


	

   dataout = driver.CreateDataSource(outputDirectory)
   proj = osr.SpatialReference()

   proj.ImportFromEPSG(2154) # 4326 = EPSG code for lon/lat WGS84 projection
   layerout = dataout.CreateLayer('obs_grp'+str(a),proj,geom_type=ogr.wkbPoint)  
   fieldDef1 = ogr.FieldDefn('h_eau', ogr.OFTReal)
   layerout.CreateField(fieldDef1)
   floutDefn = layerout.GetLayerDefn()
   feature_out = ogr.Feature(floutDefn)

   for ipix in range(len(X)):
   #- Create Geometry Point with pixel coordinates
      if True:
      #if float(ebathy[ipix])>0:
         pixel_point = ogr.Geometry(ogr.wkbPoint)
         pixel_point.AddPoint(X[ipix],Y[ipix])
         #- Add the geometry to the feature
         feature_out.SetGeometry(pixel_point)
         #- Set feature attributes
         feature_out.SetField('h_eau', float(ebathy[ipix]))
         layerout.CreateFeature(feature_out)
         #- Delete point geometry
#         pixel_point.Destroy()

   #feature_out.Destroy()
   dataout.Destroy()
