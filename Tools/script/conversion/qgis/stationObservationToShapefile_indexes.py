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

from function import *


pathInput='./bin/input.txt'
outputDirectoryQgis='./bin/Qgis/'
outputDirectory=outputDirectoryQgis+'stations'


pathobstxt='./bin/obs.txt'


meshFile=readFromInputFile(pathInput,['mesh_name'])
meshFile='./bin/'+str(meshFile[0])


obsfile=open(pathobstxt,'r')
ObsFileLine=obsfile.readlines()

j=0
Toread=''
numberLineToRead=0
print len(ObsFileLine)

Xobs,Yobs=[],[]



#Recupération des stations
stop=True
i=0
while stop:
   line=ObsFileLine[i]
   if ('stations' in line and 'stations_with_grp' not in line) :
      line=cleanString(line)
      line=line.replace('stations','')
      numberStation=int(line)
      compt=numberStation
      while compt>0:
         i=i+1
         line=ObsFileLine[i]
         line_temp=cleanString(line)
         if line_temp!='':
            line=line.split(' ')
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
   layerout = dataout.CreateLayer('obs_'+str(a),proj,geom_type=ogr.wkbPoint)  
   floutDefn = layerout.GetLayerDefn()
   feature_out = ogr.Feature(floutDefn)


   pixel_point = ogr.Geometry(ogr.wkbPoint)
   pixel_point.AddPoint(Xobs[i],Yobs[i])
   feature_out.SetGeometry(pixel_point)
   layerout.CreateFeature(feature_out)
   pixel_point.Destroy()

   feature_out.Destroy()
   dataout.Destroy()

[numberOfNode,numberOfCell,Xnode,Ynode,Connect]=ReadDassFlowMesh(meshFile)

for i in range(numberofobs):

   
  
   a="{0:04d}".format(i+1)
   file=open('./bin/station_grp_'+str(a)+'.txt')
   fileline=file.readlines()
   line=fileline[0]
   line=cleanString(line)
   print 'File: station_grp_'+str(a)+'.txt'
   if line=='indexes':
      print '\t indexes ...'
      data=np.loadtxt('./bin/station_grp_'+str(a)+'.txt',skiprows=1)
      [X,Y]=cellsObservedToPointsObserved(data,Xnode,Ynode,Connect)
   else :
      print '\t points ...'
      data=np.loadtxt('./bin/station_grp_'+str(a)+'.txt',skiprows=1)
      #print data
      X=data[:,0]
      Y=data[:,1]


   if os.path.isfile(outputDirectory+'/obs_grp'+str(a)+'.shp'):
      filepath=outputDirectory+'/obs_grp'+str(a)
      os.remove(filepath+'.dbf')
      os.remove(filepath+'.prj')
      os.remove(filepath+'.shp')
      os.remove(filepath+'.shx')



   dataout = driver.CreateDataSource(outputDirectory)
   proj = osr.SpatialReference()

   proj.ImportFromEPSG(2154) # 4326 = EPSG code for lon/lat WGS84 projection
   layerout = dataout.CreateLayer('obs_grp'+str(a),proj,geom_type=ogr.wkbPoint)  


   floutDefn = layerout.GetLayerDefn()
   feature_out = ogr.Feature(floutDefn)

   for ipix in range(len(X)):
   #- Create Geometry Point with pixel coordinates
      if True:
         pixel_point = ogr.Geometry(ogr.wkbPoint)
         pixel_point.AddPoint(X[ipix],Y[ipix])
         #- Add the geometry to the feature
         feature_out.SetGeometry(pixel_point)
         #- Set feature attributes
         layerout.CreateFeature(feature_out)
         #- Delete point geometry
         pixel_point.Destroy()

   feature_out.Destroy()
   dataout.Destroy()
